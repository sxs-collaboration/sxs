import sys
import pathlib
import contextlib
import tempfile
import pytest
import sxs

from .conftest import shortest_h_com_file, shortest_horizons, shortest_metadata


if sys.platform.startswith("win"):
    forked = pytest.mark.skip(reason="Windows has no way of forking")
else:
    forked = pytest.mark.forked


@pytest.mark.parametrize("persistent", (True, False))
def test_sxs_directory_bad_directory_name(persistent):
    with pytest.raises(ValueError):
        sxs.utilities.sxs_directory("cacheconfig", persistent=persistent)


@pytest.mark.parametrize(
    ("directory_type", "persistent"),
    (
        pytest.param("config", True, id="config persistent"),
        pytest.param("cache", True, id="cache persistent"),
        pytest.param("config", False, id="config temporary"),
        pytest.param("cache", False, id="cache temporary"),
    ),
)
@forked
def test_sxs_directory_cache_env(directory_type, persistent, tmp_path, monkeypatch):
    import pathlib
    import os

    sxs.utilities.sxs_directory.cache_clear()
    with monkeypatch.context() as mp:
        mp.setattr(pathlib.Path, "home", lambda: tmp_path)
        mp.setenv(f"SXS{directory_type.upper()}DIR", str(tmp_path))
        sxs_dir = sxs.utilities.sxs_directory(directory_type, persistent=persistent)
    assert isinstance(sxs_dir, pathlib.Path)
    if persistent:
        assert str(tmp_path) == str(sxs_dir)
    else:
        assert str(tmp_path) != str(sxs_dir)
    assert sxs_dir.exists()
    d = sxs_dir / "sub"
    d.mkdir()
    assert d.exists()
    p = d / "hello.txt"
    p.write_text("hi there")
    assert p.exists()
    assert p.read_text() == "hi there"


@pytest.mark.parametrize(
    ("directory_type", "platform"),
    (
        pytest.param("config", "linux", id="config linux"),
        pytest.param("cache", "linux", id="cache linux"),
        pytest.param("config", "freebsd", id="config freebsd"),
        pytest.param("cache", "freebsd", id="cache freebsd"),
    ),
)
def test_sxs_directory_linux(directory_type, platform, tmp_path, monkeypatch):
    import pathlib
    import sys

    d1 = tmp_path / "sub1"
    d1.mkdir()
    assert d1.exists()
    d2 = tmp_path / "sub2"
    d2.mkdir()
    assert d2.exists()

    with monkeypatch.context() as mp:
        mp.setattr(sys, "platform", platform)
        mp.setattr(pathlib.Path, "home", lambda: d1)

        mp.setenv(f'XDG_{directory_type.upper()}_HOME', str(d2))
        sxs.utilities.sxs_directory.cache_clear()
        sxs_dir = sxs.utilities.sxs_directory(directory_type, persistent=True)
        assert str(sxs_dir) == str(d2 / "sxs")

        mp.delenv(f'XDG_{directory_type.upper()}_HOME')
        sxs.utilities.sxs_directory.cache_clear()
        sxs_dir = sxs.utilities.sxs_directory(directory_type, persistent=True)
        assert str(sxs_dir) == str(d1 / f".{directory_type}" / "sxs")


@pytest.mark.parametrize(
    ("directory_type",), (pytest.param(dt, id=dt) for dt in ["config", "cache"]),
)
@forked
def test_sxs_directory_unwritable(directory_type, tmp_path, monkeypatch):
    import time
    import os
    import stat
    from pathlib import Path

    d = tmp_path / "unwritable"
    d.mkdir(mode=0o000)
    assert d.exists()
    assert d.is_dir()
    assert stat.filemode(d.stat().st_mode) == "d---------"
    assert not os.access(str(d), os.W_OK)

    with monkeypatch.context() as mp:
        mp.setattr(Path, "home", lambda: d)
        mp.setenv(f"SXS{directory_type.upper()}DIR", str(d))
        for dir_type in ["config", "cache"]:
            mp.delenv(f'XDG_{dir_type.upper()}_HOME', raising=False)
        sxs.utilities.sxs_directory.cache_clear()
        with pytest.warns(UserWarning):
            sxs_dir = sxs.utilities.sxs_directory(directory_type, persistent=True)
        assert ".sxs" not in str(sxs_dir), (str(sxs_dir), str(d))
        sxs_config_dir = sxs.utilities.sxs_directory("config", persistent=True)
        assert str(sxs_config_dir) in str(sxs_dir)

    d.chmod(0o777)
    time.sleep(0.1)


def test_read_write_config(tmp_path, monkeypatch):
    import pathlib

    sxs.utilities.sxs_directory.cache_clear()
    original_config = sxs.utilities.read_config()

    sxs.utilities.sxs_directory.cache_clear()
    with monkeypatch.context() as mp:
        mp.setattr(pathlib.Path, "home", lambda: tmp_path)
        mp.setenv(f"SXSCONFIGDIR", str(tmp_path))
        assert not sxs.utilities.read_config()
        sxs.utilities.write_config(NONSENSEGARBAGE=123)
        read = sxs.utilities.read_config()
        assert len(read) == 1
        assert read["NONSENSEGARBAGE"] == 123
        assert sxs.utilities.read_config("NONSENSEGARBAGE", 345) == 123
        assert sxs.utilities.read_config("NONSENSEGARBAGEJUNK", 345) == 345
        sxs.utilities.write_config(NONSENSEGARBAGEJUNK=345)
        assert sxs.utilities.read_config("NONSENSEGARBAGE", 345) == 123
        assert sxs.utilities.read_config("NONSENSEGARBAGEJUNK", 567) == 345

    sxs.utilities.sxs_directory.cache_clear()
    final_config = sxs.utilities.read_config()
    assert sorted(set(original_config)) == sorted(set(final_config))
    for k in original_config:
        assert original_config[k] == final_config[k]


def test_sxs_directory_cache_from_config(tmp_path, monkeypatch):
    import pathlib
    import os

    sxs.utilities.sxs_directory.cache_clear()
    with monkeypatch.context() as mp:
        mp.setattr(pathlib.Path, "home", lambda: tmp_path)
        sxs.utilities.write_config(cache_directory=str(tmp_path / "newcache"))
        sxs_dir = sxs.utilities.sxs_directory("cache", persistent=True)
    assert isinstance(sxs_dir, pathlib.Path)
    assert str(sxs_dir) == str(tmp_path / "newcache")


def test_select_by_path_component():
    from sxs.utilities import select_by_path_component

    # Examples given in the docstring
    possible_matches = [
        "abc/lm/xyz",
        "abc/ln/xyz",
        "abd/lm/xyz",
        "abd/ln/xyz",
    ]
    assert select_by_path_component("ab/ln/xyz", possible_matches) == {"abd/ln/xyz"}
    assert select_by_path_component("ab./ln/xyz", possible_matches) == {"abc/ln/xyz", "abd/ln/xyz"}
    assert select_by_path_component("a./ln/xyz", possible_matches) == {"abd/ln/xyz"}
    assert select_by_path_component("ab./l./x", possible_matches) == {"abc/lm/xyz", "abc/ln/xyz", "abd/lm/xyz", "abd/ln/xyz"}
    assert select_by_path_component("ab./l/x", possible_matches) == {"abc/ln/xyz", "abd/ln/xyz"}

    possible_matches += ["abe/lp/xyz"]
    assert select_by_path_component("ab/ln/xyz", possible_matches) == set()
    assert select_by_path_component("ab./ln/xyz", possible_matches) == {"abc/ln/xyz", "abd/ln/xyz"}
    assert select_by_path_component("a./ln/xyz", possible_matches) == set()
    assert select_by_path_component("ab./l./x", possible_matches) == set(possible_matches)
    assert select_by_path_component("ab./l/x", possible_matches) == {"abc/ln/xyz", "abd/ln/xyz", "abe/lp/xyz"}

    assert select_by_path_component("ab./l", possible_matches) == {"abc/ln/xyz", "abd/ln/xyz", "abe/lp/xyz"}

    assert select_by_path_component("abe/lp/xyz", ["abe/lp/xyz"]) == {"abe/lp/xyz"}

    # Fake SXS catalog
    fake_sxs = [
        f"SXS:BBH:{sxs_id:04}v{version}/Lev{lev}/{type}_{extrapolation}.h5"
        for sxs_id in range(1, 2023)
        for version in [1, 3]
        for lev in [2, 3, 4]
        for type in ["h", "psi4"]
        for extrapolation in ["extrapolated_n2", "extrapolated_n3", "extrapolated_n4"]
    ]

    assert select_by_path_component("SXS:BBH:012(3|4)/Lev/h.*n(2|3)", fake_sxs) == set(
        f"SXS:BBH:{sxs_id:04}v{version}/Lev{lev}/{type}_{extrapolation}.h5"
        for sxs_id in [123, 124]
        for version in [3]
        for lev in [4]
        for type in ["h"]
        for extrapolation in ["extrapolated_n2", "extrapolated_n3"]
    )

    assert select_by_path_component("SXS:B/Lev4/psi4", fake_sxs) == set(
        f"SXS:BBH:{sxs_id:04}v{version}/Lev{lev}/{type}_{extrapolation}.h5"
        for sxs_id in [2022]
        for version in [3]
        for lev in [4]
        for type in ["psi4"]
        for extrapolation in ["extrapolated_n4"]
    )

    assert select_by_path_component("SXS:BHWD/Lev4/psi4", fake_sxs) == set()

    assert select_by_path_component("SXS:BBH:1001v3/Lev3/h_extrapolated_n3.h5", fake_sxs) == set(
        f"SXS:BBH:{sxs_id:04}v{version}/Lev{lev}/{type}_{extrapolation}.h5"
        for sxs_id in [1001]
        for version in [3]
        for lev in [3]
        for type in ["h"]
        for extrapolation in ["extrapolated_n3"]
    )

    assert select_by_path_component("SXS:BBH:1001v3/Lev3/h", fake_sxs) == set(
        f"SXS:BBH:{sxs_id:04}v{version}/Lev{lev}/{type}_{extrapolation}.h5"
        for sxs_id in [1001]
        for version in [3]
        for lev in [3]
        for type in ["h"]
        for extrapolation in ["extrapolated_n2", "extrapolated_n3", "extrapolated_n4"]
    )


@pytest.mark.parametrize("truncation_tol", (True, None, 1e-9))
def test_lvcnr_format(truncation_tol):
    import shutil
    sxs_id = sxs.sxs_id(shortest_h_com_file)
    sxs_lev = sxs.lev_number(shortest_h_com_file)
    shortest_lvcnr = f"{sxs_id.replace(':', '_')}_Res{sxs_lev}.h5"
    # shortest_system = pathlib.PosixPath(shortest_h_com_file).parent
    # shortest_system_path = sxs.sxs_directory("cache") / sxs.utilities.sxs_path_to_system_path(shortest_h_com_file)
    with contextlib.redirect_stdout(None), contextlib.redirect_stderr(None):
        # Load these files, just to ensure they exist
        sxs.load(shortest_h_com_file, download=True, cache=True, progress=False, extrapolation_order=2)
        sxs.load(shortest_horizons, download=True, cache=True, progress=False)
        sxs.load(shortest_metadata, download=True, cache=True, progress=False)
        # Get the truepath values
        truepath_h = sxs.utilities.cached_path(shortest_h_com_file)
        truepath_horizons = sxs.utilities.cached_path(shortest_horizons)
        truepath_metadata = sxs.utilities.cached_path(shortest_metadata)
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = pathlib.Path(temp_dir)
            temp_sxs_path = temp_path / f"Lev{sxs_lev}"
            temp_sxs_path.mkdir(parents=True)
            shutil.copy(str(truepath_h), str(temp_sxs_path))
            shutil.copy(str(truepath_horizons), str(temp_sxs_path))
            shutil.copy(str(truepath_metadata), str(temp_sxs_path))
            lvc_converter = sxs.utilities.lvcnr.SimulationConverter(quiet=True)
            lvc_converter.convert(str(temp_sxs_path), temp_dir, truncation_tol=truncation_tol)
            error_count = sxs.utilities.lvcnr.compare(str(temp_path / shortest_lvcnr), str(temp_sxs_path))
            assert error_count == 0

import pytest
import sxs


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
@pytest.mark.forked
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
    ("directory_type",),
    (
        pytest.param("config", id="config"),
        pytest.param("cache", id="cache"),
    ),
)
@pytest.mark.forked
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
        sxs.utilities.write_config(cache=str(tmp_path / "newcache"))
        sxs_dir = sxs.utilities.sxs_directory("cache", persistent=True)
    assert isinstance(sxs_dir, pathlib.Path)
    assert str(sxs_dir) == str(tmp_path / "newcache")

import pytest
import sxs
from .conftest import skip_macOS_GH_actions_downloads


@skip_macOS_GH_actions_downloads
def test_sxs_load_v2():
    # We'll use SXS:BBH:0001, even though it's deprecated, because we
    # want to be sure that it always exists, and is accessible as a
    # v2.0 record.
    with pytest.warns(UserWarning):
        s = sxs.load("SXS:BBH:0001v2.0")
    assert s.sxs_id_stem == "SXS:BBH:0001"
    assert s.version == "v2.0"
    s.h
    s.strain
    s.H
    s.Strain
    s.horizons
    s.Horizons
    s.psi4
    s.Psi4


@skip_macOS_GH_actions_downloads
def test_sxs_load_v2_levs():
    for lev, lev_number in [("", 3), ("/Lev3", 3), ("/Lev2", 2), ("/Lev1", 1)]:
        s = sxs.load(f"SXS:BBH:1001v2.0{lev}", ignore_deprecation=True)
        assert s.sxs_id_stem == "SXS:BBH:1001"
        assert s.version == "v2.0"
        assert s.lev_number == lev_number
        s.h
        s.strain
        s.H
        s.Strain
        s.horizons
        s.Horizons


@skip_macOS_GH_actions_downloads
@pytest.mark.parametrize("loader", [sxs.Simulation, sxs.load])
def test_superseding(loader):
    simulation = "SXS:BBH:0001"

    # No version, no arguments
    with pytest.raises(ValueError):
        loader(f"{simulation}")

    # Version, no arguments
    with pytest.warns(UserWarning):
        s = loader(f"{simulation}v2.0")
    assert s.sxs_id_stem == simulation

    # No version, ignore_deprecation
    s = loader(f"{simulation}", ignore_deprecation=True)
    assert s.sxs_id_stem == simulation

    # Version, ignore_deprecation
    s = loader(f"{simulation}v2.0", ignore_deprecation=True)
    assert s.sxs_id_stem == simulation

    # No version, auto_supersede
    with pytest.warns(UserWarning):
        s = loader(f"{simulation}", auto_supersede=True)
    assert s.sxs_id_stem != simulation

    # Version, auto_supersede
    with pytest.warns(UserWarning):
        s = loader(f"{simulation}v2.0", auto_supersede=True)
    assert s.sxs_id_stem == simulation

    # No version, ignore_deprecation, auto_supersede
    s = loader(f"{simulation}", ignore_deprecation=True, auto_supersede=True)
    assert s.sxs_id_stem == simulation

    # Version, ignore_deprecation, auto_supersede
    s = loader(f"{simulation}v2.0", ignore_deprecation=True, auto_supersede=True)
    assert s.sxs_id_stem == simulation

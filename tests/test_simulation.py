import pytest
import sxs
from .conftest import skip_macOS_GH_actions_downloads


def test_sxs_load_v2():
    s = sxs.load("SXS:BBH:0001v2.0")
    assert s.sxs_id_stem == "SXS:BBH:0001"
    assert s.version == "v2.0"
    s.h
    s.horizons


@skip_macOS_GH_actions_downloads
@pytest.mark.parametrize("loader", [sxs.Simulation, sxs.load])
def test_superseded_by_single(loader):
    simulation = "SXS:BBH:0001"
    # NOTE: These tests as written will only work so long as only one
    # simulation is listed in the given simulation's `superseded_by`
    # field.  If more are added, this test will need to be changed to
    # use a simulation that still has just one.
    with pytest.raises(ValueError):
        loader(f"{simulation}")
    with pytest.warns(UserWarning):
        s = loader(f"{simulation}v2.0")
    assert s.sxs_id_stem == simulation
    s = loader(f"{simulation}", ignore_deprecation=True)
    assert s.sxs_id_stem == simulation
    s = loader(f"{simulation}v2.0", ignore_deprecation=True)
    assert s.sxs_id_stem == simulation
    with pytest.warns(UserWarning):
        s = loader(f"{simulation}", auto_supersede=True)
    assert s.sxs_id_stem != simulation
    with pytest.warns(UserWarning):
        s = loader(f"{simulation}v2.0", auto_supersede=True)
    assert s.sxs_id_stem != simulation
    s = loader(f"{simulation}", ignore_deprecation=True, auto_supersede=True)
    assert s.sxs_id_stem == simulation
    s = loader(f"{simulation}v2.0", ignore_deprecation=True, auto_supersede=True)
    assert s.sxs_id_stem == simulation


@skip_macOS_GH_actions_downloads
@pytest.mark.parametrize("loader", [sxs.Simulation, sxs.load])
def test_superseded_by_multiple(loader):
    simulation = "SXS:BBH:0030"
    # NOTE: These tests assume that multiple simulations are listed in
    # the given simulation's `superseded_by` field.  If that changes
    # for some reason, this test will need to be changed to use a
    # simulation that still has multiple simulations.
    with pytest.raises(ValueError):
        loader(f"{simulation}")
    with pytest.warns(UserWarning):
        s = loader(f"{simulation}v2.0")
    assert s.sxs_id_stem == simulation
    s = loader(f"{simulation}", ignore_deprecation=True)
    assert s.sxs_id_stem == simulation
    s = loader(f"{simulation}v2.0", ignore_deprecation=True)
    assert s.sxs_id_stem == simulation
    with pytest.raises(ValueError):
        loader(f"{simulation}", auto_supersede=True)
    with pytest.raises(ValueError):
        loader(f"{simulation}v2.0", auto_supersede=True)
    s = loader(f"{simulation}", ignore_deprecation=True, auto_supersede=True)
    assert s.sxs_id_stem == simulation
    s = loader(f"{simulation}v2.0", ignore_deprecation=True, auto_supersede=True)
    assert s.sxs_id_stem == simulation

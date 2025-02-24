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

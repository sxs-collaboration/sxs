import pytest
import sxs


def test_superseded_by_single():
    simulation = "SXS:BBH:0001"
    # NOTE: These tests as written will only work so long as only one
    # simulation is listed in the given simulation's `superseded_by`
    # field.  If more are added, this test will need to be changed to
    # use a simulation that still has just one.
    with pytest.raises(ValueError):
        sxs.Simulation(f"{simulation}")
    with pytest.deprecated_call():
        s = sxs.Simulation(f"{simulation}v2.0")
    assert s.sxs_id_stem == simulation
    s = sxs.Simulation(f"{simulation}", ignore_deprecation=True)
    assert s.sxs_id_stem == simulation
    s = sxs.Simulation(f"{simulation}v2.0", ignore_deprecation=True)
    assert s.sxs_id_stem == simulation
    with pytest.deprecated_call():
        s = sxs.Simulation(f"{simulation}", auto_supersede=True)
    assert s.sxs_id_stem != simulation
    with pytest.deprecated_call():
        s = sxs.Simulation(f"{simulation}v2.0", auto_supersede=True)
    assert s.sxs_id_stem != simulation
    s = sxs.Simulation(f"{simulation}", ignore_deprecation=True, auto_supersede=True)
    assert s.sxs_id_stem == simulation
    s = sxs.Simulation(f"{simulation}v2.0", ignore_deprecation=True, auto_supersede=True)
    assert s.sxs_id_stem == simulation

def test_superseded_by_multiple():
    simulation = "SXS:BBH:0030"
    # NOTE: These tests assume that multiple simulations are listed in
    # the given simulation's `superseded_by` field.  If that changes
    # for some reason, this test will need to be changed to use a
    # simulation that still has multiple simulations.
    with pytest.raises(ValueError):
        sxs.Simulation(f"{simulation}")
    with pytest.deprecated_call():
        s = sxs.Simulation(f"{simulation}v2.0")
    assert s.sxs_id_stem == simulation
    s = sxs.Simulation(f"{simulation}", ignore_deprecation=True)
    assert s.sxs_id_stem == simulation
    s = sxs.Simulation(f"{simulation}v2.0", ignore_deprecation=True)
    assert s.sxs_id_stem == simulation
    with pytest.raises(ValueError):
        sxs.Simulation(f"{simulation}", auto_supersede=True)
    with pytest.raises(ValueError):
        sxs.Simulation(f"{simulation}v2.0", auto_supersede=True)
    s = sxs.Simulation(f"{simulation}", ignore_deprecation=True, auto_supersede=True)
    assert s.sxs_id_stem == simulation
    s = sxs.Simulation(f"{simulation}v2.0", ignore_deprecation=True, auto_supersede=True)
    assert s.sxs_id_stem == simulation

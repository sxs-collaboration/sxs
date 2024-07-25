import numpy as np
import sxs
from sxs import julia

def test_PNWaveform():
    M1 = 0.3
    M2 = 0.7
    chi1 = [0.1, 0.2, 0.3]
    chi2 = [-0.2, 0.1, -0.3]
    Omega_i = 0.01
    PNOrder = 4.0
    waveform_pn_order = 4.0
    saves_per_orbit = 32
    w1 = sxs.julia.PNWaveform(
        M1, M2, chi1, chi2, Omega_i,
        PNOrder=PNOrder, waveform_pn_order=waveform_pn_order,
        saves_per_orbit=saves_per_orbit
    )

    Approximant = "TaylorT1"
    delta = sxs.julia.PostNewtonian.mass_difference_ratio(M1, M2)
    chi1_i = chi1
    chi2_i = chi2
    Omega_orb_i = Omega_i
    inertial = True
    w2 = sxs.julia.GWFrames.PNWaveform(Approximant, delta, chi1_i, chi2_i, Omega_orb_i, inertial=inertial)

    for w in [w1, w2]:
        n_times = w.n_times
        assert w.t[0] == 0.0
        assert w.data.shape == (n_times, 77)
        assert w.data_type == "h"
        assert w.ell_min == 2
        assert w.ell_max == 8
        assert w.spin_weight == -2
        assert w.frame.shape == (1, 4)
        assert w.frame_type == "inertial"
        assert w.M1.shape == (n_times,)
        assert w.M2.shape == (n_times,)
        assert w.chi1.shape == (n_times, 3)
        assert w.chi2.shape == (n_times, 3)
        assert w.coorbital_frame.shape == (n_times, 4)
        assert w.v.shape == (n_times,)
        assert w.orbital_phase.shape == (n_times,)

    assert np.allclose(w1.t[:-10], w2.t[:-10])
    assert np.allclose(w1.data[:-10, :], w2.data[:-10, :], rtol=1e-3, atol=1e-5)
    assert np.allclose(w1.frame, w2.frame)
    assert np.allclose(w1.M1[:-10], w2.M1[:-10])
    assert np.allclose(w1.M2[:-10], w2.M2[:-10])
    assert np.allclose(w1.chi1[:-10, :], w2.chi1[:-10, :], rtol=1e-3)
    assert np.allclose(w1.chi2[:-10, :], w2.chi2[:-10, :], rtol=1e-3)
    assert np.allclose(w1.coorbital_frame[:-10, :], w2.coorbital_frame[:-10, :])

    w1_corot = w1.copy().to_corotating_frame()
    w2_corot = w2.copy().to_corotating_frame()
    for (w, w_corot) in [(w1, w1_corot), (w2, w2_corot)]:
        assert np.array_equal(w.t, w_corot.t)
        assert w.data.shape == w_corot.data.shape
        assert w.coorbital_frame.shape == w_corot.frame.shape
        assert w_corot.frame_type == "corotating"

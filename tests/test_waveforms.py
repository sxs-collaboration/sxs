import contextlib
import sys
import pathlib
import tempfile

import numpy as np
import quaternionic
import spherical
import h5py
import pytest
import sxs

from .conftest import shortest_h_com_file, constant_waveform

try:
    import spinsfast
    requires_spinsfast = lambda f: f
except:
    requires_spinsfast = pytest.mark.skip(reason="spinsfast is missing")


def test_backwards_compatibility():
    path = sxs.sxs_directory("cache") / sxs.utilities.sxs_path_to_system_path(shortest_h_com_file)
    with contextlib.redirect_stdout(None):
        with sxs.loadcontext(shortest_h_com_file, extrapolation_order=...) as h:
            with h5py.File(path, "r") as f:
                for group in f:
                    if group == "VersionHist.ver":
                        continue
                    g = f[group]
                    for d in g:
                        if d == "History.txt":
                            continue
                        assert np.array_equal(g[d], h[group][d])
                for group in f:
                    if group != "VersionHist.ver":
                        continue
                    assert np.array_equal(f[group], h[group])


@requires_spinsfast
def test_boost(h):
    ell_max = 2*h.ell_max
    hprime = h.boost(np.array([0.01, 0.02, 0.03]), ell_max=ell_max)
    assert h.spin_weight == hprime.spin_weight
    assert h.data_type == hprime.data_type
    assert hprime.ell_max == ell_max


def test_modes_conjugate():
    import spherical as sf
    tolerance = 1e-14
    np.random.seed(1234)
    for inplace in [False, True]:
        for s in range(-2, 2 + 1):
            ell_min = abs(s)
            ell_max = 8
            a = np.random.rand(3, 7, sf.Ysize(ell_min, ell_max)*2).view(complex)
            m = sf.Modes(a, spin_weight=s, ell_min=ell_min, ell_max=ell_max)
            g = m.grid()
            s = m.s
            ell_min = m.ell_min
            ell_max = m.ell_max
            shape = m.shape
            mbar = m.conjugate(inplace)
            gbar = mbar.grid()
            assert s == -mbar.s
            assert ell_min == mbar.ell_min
            assert ell_max == mbar.ell_max
            assert shape == mbar.shape
            assert np.allclose(g, np.conjugate(gbar), rtol=tolerance, atol=tolerance), (
                f"max|g-ḡ.bar| = {np.max(np.abs(g - np.conjugate(gbar)))}"
            )


def test_modes_real():
    import spherical as sf
    tolerance = 1e-14
    np.random.seed(1234)
    for inplace in [False, True]:
        s = 0
        ell_min = abs(s)
        ell_max = 8
        a = np.random.rand(3, 7, sf.Ysize(ell_min, ell_max)*2).view(complex)
        # Test success with spin_weight==0
        m = sf.Modes(a, spin_weight=s, ell_min=ell_min, ell_max=ell_max)
        g = m.grid()
        s = m.s
        ell_min = m.ell_min
        ell_max = m.ell_max
        shape = m.shape
        mreal = m._real_func(inplace)
        greal = mreal.grid()
        assert s == mreal.s
        assert ell_min == mreal.ell_min
        assert ell_max == mreal.ell_max
        assert shape == mreal.shape
        assert np.allclose(greal, np.real(greal)+0.0j, rtol=tolerance, atol=tolerance)
        assert np.allclose(np.real(g), np.real(greal), rtol=tolerance, atol=tolerance)
        assert np.allclose(np.zeros_like(g, dtype=float), np.imag(greal), rtol=tolerance, atol=tolerance)
        # Test failure with s!=0
        for s in [-3, -2, -1, 1, 2, 3]:
            m = sf.Modes(a, spin_weight=s, ell_min=ell_min, ell_max=ell_max)
            with pytest.raises(ValueError):
                mreal = m._real_func(inplace)


def test_modes_imag():
    import spherical as sf
    tolerance = 1e-14
    np.random.seed(1234)
    for inplace in [False, True]:
        s = 0
        ell_min = abs(s)
        ell_max = 8
        a = np.random.rand(3, 7, sf.Ysize(ell_min, ell_max)*2).view(complex)
        # Test success with spin_weight==0
        m = sf.Modes(a, spin_weight=s, ell_min=ell_min, ell_max=ell_max)
        g = m.grid()
        s = m.s
        ell_min = m.ell_min
        ell_max = m.ell_max
        shape = m.shape
        mimag = m._imag_func(inplace)
        gimag = mimag.grid()
        assert s == mimag.s
        assert ell_min == mimag.ell_min
        assert ell_max == mimag.ell_max
        assert shape == mimag.shape
        assert np.allclose(gimag, np.real(gimag), rtol=tolerance, atol=tolerance)  # gimag is purely real
        assert np.allclose(
            np.array(np.imag(g.ndarray), dtype=complex),
            gimag.ndarray,
            rtol=tolerance, atol=tolerance
        )  # imag(g) == gimag
        assert np.allclose(
            np.imag(gimag.ndarray),
            np.zeros_like(g.ndarray, dtype=float),
            rtol=tolerance, atol=tolerance
        ) # imag(gimag) == 0
        # Test failure with s!=0
        for s in [-3, -2, -1, 1, 2, 3]:
            m = sf.Modes(a, spin_weight=s, ell_min=ell_min, ell_max=ell_max)
            with pytest.raises(ValueError):
                mimag = m._imag_func(inplace)


def test_modes_squared_angular_momenta():
    import spherical as sf
    tolerance = 1e-13
    np.random.seed(1234)
    L2 = sf.Modes.Lsquared
    Lz = sf.Modes.Lz
    Lp = sf.Modes.Lplus
    Lm = sf.Modes.Lminus
    R2 = sf.Modes.Rsquared
    Rz = sf.Modes.Rz
    Rp = sf.Modes.Rplus
    Rm = sf.Modes.Rminus
    for s in range(-2, 2+1):
        ell_min = abs(s)
        ell_max = 8
        a = np.random.rand(3, 7, sf.Ysize(ell_min, ell_max)*2).view(complex)
        m = sf.Modes(a, spin_weight=s, ell_min=ell_min, ell_max=ell_max)

        # Test L^2 = 0.5(L+L- + L-L+) + LzLz
        m1 = L2(m)
        m2 = 0.5 * (Lp(Lm(m)) + Lm(Lp(m))) + Lz(Lz(m))
        assert np.allclose(m1, m2, rtol=tolerance, atol=tolerance)

        # Test R^2 = 0.5(R+R- + R-R+) + RzRz
        m1 = R2(m)
        m2 = 0.5 * (Rp(Rm(m)) + Rm(Rp(m))) + Rz(Rz(m))
        assert np.allclose(m1, m2, rtol=tolerance, atol=tolerance)

        # Test L^2 = R^2
        m1 = L2(m)
        m2 = R2(m)
        assert np.allclose(m1, m2, rtol=tolerance, atol=tolerance)


def test_modes_derivative_commutators():
    import spherical as sf
    tolerance = 1e-13
    np.random.seed(1234)
    # Note that post-fix operators are in the opposite order compared
    # to prefixed commutators, so we pull the post-fix operators out
    # as functions to make things look right.
    np.random.seed(1234)
    L2 = sf.Modes.Lsquared
    Lz = sf.Modes.Lz
    Lp = sf.Modes.Lplus
    Lm = sf.Modes.Lminus
    R2 = sf.Modes.Rsquared
    Rz = sf.Modes.Rz
    Rp = sf.Modes.Rplus
    Rm = sf.Modes.Rminus
    eth = lambda modes: modes.eth
    ethbar = lambda modes: modes.ethbar
    for s in range(-2, 2+1):
        ell_min = abs(s)
        ell_max = 8
        a = np.random.rand(3, 7, sf.Ysize(ell_min, ell_max)*2).view(complex)
        m = sf.Modes(a, spin_weight=s, ell_min=ell_min, ell_max=ell_max)
        # Test [Ri, Lj] = 0
        for R in [Rz, Rp, Rm]:
            for L in [Lz, Lp, Lm]:
                assert np.max(np.abs(L(R(m)) - R(L(m)))) < tolerance
        # Test [L2, Lj] = 0
        for L in [Lz, Lp, Lm]:
            assert np.max(np.abs(L2(L(m)) - L(L2(m)))) < 5*tolerance
        # Test [R2, Rj] = 0
        for R in [Rz, Rp, Rm]:
            assert np.max(np.abs(R2(R(m)) - R(R2(m)))) < 5*tolerance
        # Test [Lz, Lp] = Lp
        assert np.allclose(Lz(Lp(m)) - Lp(Lz(m)), Lp(m), rtol=tolerance, atol=tolerance)
        # Test [Lz, Lm] = -Lm
        assert np.allclose(Lz(Lm(m)) - Lm(Lz(m)), -Lm(m), rtol=tolerance, atol=tolerance)
        # Test [Lp, Lm] = 2Lz
        assert np.allclose(Lp(Lm(m)) - Lm(Lp(m)), 2 * Lz(m), rtol=tolerance, atol=tolerance)
        # Test [Rz, Rp] = Rp
        assert np.allclose(Rz(Rp(m)) - Rp(Rz(m)), Rp(m), rtol=tolerance, atol=tolerance)
        # Test [Rz, Rm] = -Rm
        assert np.allclose(Rz(Rm(m)) - Rm(Rz(m)), -Rm(m), rtol=tolerance, atol=tolerance)
        # Test [Rp, Rm] = 2Rz
        assert np.allclose(Rp(Rm(m)) - Rm(Rp(m)), 2 * Rz(m), rtol=tolerance, atol=tolerance)
        # Test [ethbar, eth] = 2s
        assert np.allclose(ethbar(eth(m)) - eth(ethbar(m)), 2 * m.s * m, rtol=tolerance, atol=tolerance)


def test_modes_evaluate(h, eps):
    import time

    ell_max = h.ell_max
    ϵ = 5 * (2 * ell_max + 1) * 2 * eps

    m1 = h[:1000].copy()
    m2 = h[:1000].copy()
    m2._metadata["frame_type"] = "corotating"
    m2._metadata["frame"] = quaternionic.one * np.ones_like(h.t)

    equiangular_grid = spherical.theta_phi(2 * ell_max + 1, 2 * ell_max + 1)

    t0 = time.perf_counter()
    g1 = m1.evaluate(equiangular_grid)
    t1 = time.perf_counter()
    g2 = m2.evaluate(equiangular_grid)
    t2 = time.perf_counter()

    print()
    print(f"\tTime for inertial frame: {t1-t0:.4f} seconds")
    print(f"\tTime for 'rotating' frame: {t2-t1:.4f} seconds")

    assert np.allclose(g1, g2, rtol=ϵ, atol=ϵ), f"max|g1-g2|={np.max(np.abs(g1-g2))}"


def test_modes_rotate(h, eps):
    import time

    ϵ = 5 * (2 * h.ell_max + 1) * 2 * eps

    print()
    for i, R in enumerate([quaternionic.one, quaternionic.one * np.ones_like(h.t)]):
        t1 = time.perf_counter()
        hprm = h.rotate(R)
        t2 = time.perf_counter()
        print(f"\tRotation {i+1} took {t2-t1:.4f} seconds")
        assert type(h) == type(hprm)
        assert np.array_equal(h.t, hprm.t)
        assert np.allclose(h.ndarray, hprm.ndarray, rtol=ϵ, atol=ϵ)

        metadata = h._metadata.copy()
        metadataprm = hprm._metadata.copy()
        for d in [metadata, metadataprm]:
            for key in ['time', 'frame']:
                d.pop(key, None)
            for key in ['space_translation', 'boost_velocity']:
                d[key] = d[key].tolist()
        assert metadata == metadataprm


def test_modes_rotate_evaluate(h, Rs, eps):
    """Test that evaluating modes at rotated points == evaluating rotated modes at points"""
    import time

    ell_max = h.ell_max
    ϵ = (2 * h.ell_max + 1) * 2 * eps

    equiangular_grid = spherical.theta_phi(2 * ell_max + 1, 2 * ell_max + 1)
    Rθϕ = quaternionic.array.from_spherical_coordinates(equiangular_grid)

    for i, R in enumerate(Rs):
        hprm = h.copy().rotate(R)  # hprm = h @ 𝔇(R)
        m1 = hprm.evaluate(Rθϕ)  # m1 = hprm @ 𝔇(Rθϕ) √...
        m2 = h.evaluate(R * Rθϕ)  # m2 = h @ 𝔇(R * Rθϕ) √...
        assert np.allclose(m1, m2, rtol=ϵ, atol=ϵ)


def test_zero_angular_velocity():
    w = constant_waveform(end=10.0, n_times=10000)
    ω = w.angular_velocity
    assert np.allclose(ω, np.zeros_like(ω), atol=1e-15, rtol=0.0)


def test_z_angular_velocity():
    w = constant_waveform(end=10.0, n_times=10000)
    ω = 2 * np.pi / 5.0
    R = np.exp(quaternionic.array.from_vector_part([0, 0, ω / 2]) * w.t)
    w = w.rotate(~R)
    ω_out = w.angular_velocity
    ω_in = np.zeros_like(ω_out)
    ω_in[:, 2] = ω
    assert np.allclose(ω_in, ω_out, atol=1e-12, rtol=2e-8), (
        f"\nω_in = np.array({ω_in.tolist()})\n"
        f"\nω_out = np.array({ω_out.tolist()})\n"
    )


def test_rotated_angular_velocity():
    w = constant_waveform(end=10.0, n_times=10000)
    ω = 2 * np.pi / 5.0
    R0 = quaternionic.array(1, 2, 3, 4).normalized
    R = R0 * np.exp(quaternionic.array.from_vector_part([0, 0, ω / 2]) * w.t)
    w = w.rotate(~R)
    ω = R0 * quaternionic.array.from_vector_part([0, 0, ω]) * R0.inverse
    ω_out = w.angular_velocity
    ω_in = np.zeros_like(ω_out)
    ω_in[:, 0] = ω.x
    ω_in[:, 1] = ω.y
    ω_in[:, 2] = ω.z
    assert np.allclose(ω_in, ω_out, atol=1e-12, rtol=2e-8), (
        f"\nω_in = np.array({ω_in.tolist()})\n"
        f"\nω_out = np.array({ω_out.tolist()})\n"
    )

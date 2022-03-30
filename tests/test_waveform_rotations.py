import numpy as np
import quaternionic
import spherical
import pytest
import sxs

from .conftest import constant_waveform, linear_waveform, random_waveform


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


@pytest.mark.parametrize("w", [linear_waveform, constant_waveform])
def test_dpa_simple_cases(w, eps):
    LL = w().dominant_eigenvector_LL()
    LL_expected = np.zeros_like(LL, dtype=float)
    LL_expected[:, 2] = 1.0
    assert np.allclose(LL, LL_expected, rtol=eps, atol=eps)


@pytest.mark.parametrize("w", [linear_waveform, constant_waveform])
def test_dpa_rotated_simple_cases(w, Rs):
    # We use `begin=1.0` because we need to avoid situations where the modes
    # are all zeros, which can happen in `linear_waveform` at t=0.0
    W = w(begin=1.0, ell_min=0, n_times=len(Rs))
    LL = W.rotate(Rs.conjugate()).dominant_eigenvector_LL()
    LL_expected = (Rs * quaternionic.array([quaternionic.z for _ in range(len(Rs))]) * Rs.conjugate()).vector

    # with np.printoptions(precision=4, linewidth=180, suppress=True):
    #     print("LL:")
    #     print(LL)
    #     print()
    #     print("LL_expected:")
    #     print(LL_expected)

    # Because the dpa is only defined up to a sign, all we need is for the
    # dot product between the dpa and the expected value to be close to
    # either 1 or -1.  This finds the largest difference, based on the
    # smaller of the two sign options.
    assert (
        max(
            np.amin(
                np.vstack((np.linalg.norm(LL - LL_expected, axis=1), np.linalg.norm(LL + LL_expected, axis=1))), axis=0
            )
        )
        < 1.0e-14
    )


@pytest.mark.parametrize("w", [linear_waveform, constant_waveform, random_waveform])
def test_dpa_rotated_generally(w, Rs):
    np.random.seed(1234)
    n_copies = 10
    W = w(begin=1.0, end=100.0, n_times=n_copies * len(Rs), ell_min=0, ell_max=8)
    R_basis = quaternionic.array([R for R in Rs for _ in range(n_copies)])

    # We use `begin=1.0` because we need to avoid situations where the modes
    # are all zeros, which can happen in `linear_waveform` at t=0.0
    LL1 = (
        R_basis
        * quaternionic.array.from_vector_part(W.dominant_eigenvector_LL())
        * R_basis.conjugate()
    ).vector
    LL2 = W.rotate(R_basis.conjugate()).dominant_eigenvector_LL()

    # if (max(np.amin(np.vstack((np.linalg.norm(LL1 - LL2, axis=1), np.linalg.norm(LL1 + LL2, axis=1))), axis=0))
    #     > 1.0e-12):
    #     with np.printoptions(precision=8, linewidth=180, suppress=True):
    #         print()
    #         print(f"LL1, LL2, LL1-LL2 {LL1.shape}")
    #         print(np.stack((LL1, LL2, LL1-LL2), axis=2))
    #         print()

    # Because the dpa is only defined up to a sign, all we need is for the
    # dot product between the dpa and the expected value to be close to
    # either 1 or -1.  This finds the largest difference, based on the
    # smaller of the two sign options.
    assert (
        max(np.amin(np.vstack((np.linalg.norm(LL1 - LL2, axis=1), np.linalg.norm(LL1 + LL2, axis=1))), axis=0))
        < 1.0e-12
    )


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


def test_corotating_frame():
    w = constant_waveform(end=10.0, n_times=100_000)  # Need lots of time steps for accurate integration
    omega = 2 * np.pi / 5.0
    R0 = quaternionic.array.random().normalized
    R_in = R0 * np.exp(quaternionic.array([0, 0, 0, omega / 2]) * w.t)
    w_rot = w.copy().rotate(R_in.conjugate())
    R_out = w_rot.corotating_frame(R0=R0, tolerance=1e-12)
    assert np.allclose(R_in.ndarray, R_out.ndarray, atol=1e-10, rtol=0.0), (
        f"\nR_in = {R_in}\n"
        f"\nR_in-R_out = {R_in-R_out}\n"
        f"\nmax(abs(diff)) = {np.max(np.abs((R_in-R_out).ndarray))}\n"
    )
    w_inertial = w_rot.to_corotating_frame(R0=R0, tolerance=1e-12)
    # with np.printoptions(precision=8, linewidth=180, suppress=True):
    #     print()
    #     print(f"R0 = {R0}")
    #     print(f"R_in = {R_in}")
    #     print(f"R_out = {R_out}")
    #     print(f"w.data = {w.data}")
    #     print(f"w_rot.data = {w_rot.data}")
    #     print(f"w.data-w_rot.data = {w.data-w_rot.data}")
    #     #print(f" = {}")
    assert np.array_equal(w.t, w_inertial.t)
    assert np.allclose(w.data, w_inertial.data, atol=1e-8)
    assert w_inertial.frame_type == "corotating"

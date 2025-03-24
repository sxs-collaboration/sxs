import multiprocessing

from ..waveforms.norms import create_unified_waveforms, compute_L2_norm, compute_mismatch
from ..waveforms.alignment import align_waveforms, align_simulations
from ..waveforms.waveform_modes import h_factor, t_factor
from ..handlers import load

import numpy as np


def compute_error_summary(wa, wb, t1, t2, ASDs=None, total_masses=None):
    """
    Compute various errors between two waveforms.

    This computes the time-domain mismatch over the two-sphere,
    the relative L² norm error, the frequency-domain mismatch
    against whatever detector ASDs are provided with some total mass,
    and the modes' absolute errors and norms.

    Parameters
    ----------
    wa : WaveformModes
    wb : WaveformModes
    t1 : float
        Beginning of integrals.
    t2 : float
        End of integrals.
    ASDs : dict of funcs, optional
        Dictionary of functions mapping frequencies to the ASD of a detector(s).
        Default is no frequency-domain mismatch is calculated.
    total_masses : list of floats, optional
        Total masses in solar masses to use for frequency-domain mismatches.
        Default is 1.

    Returns
    -------
    errors : dict
        Dictionary of the time-domain mismatch over the two-sphere,
        the relative L² norm error, the frequency-domain mismatch
        against whatever detector ASDs are provided with some total mass,
        and the modes' absolute errors and norms.
    """
    errors = {}

    errors["t1"] = t1
    errors["t2"] = t2

    wa, wb = create_unified_waveforms(wa, wb, t1, t2, padding_time_factor=0)

    errors["mismatch"] = compute_mismatch(wa, wb, t1, t2)

    errors["relative L2 norm error"] = compute_L2_norm(wa, wb, t1, t2)

    f1 = (np.gradient(-np.unwrap(np.angle(wa.data[:, wa.index(2, 2)] * h_factor)), np.diff(wa.t * t_factor)[0]))[
        np.argmin(abs(wa.t - t1))
    ] / (2 * np.pi)

    wa = wa.preprocess(t1, t1 + 0.01 * (t2 - t1), t2 - (t2 - t1) * 0.01, t2)
    wb = wb.preprocess(t1, t1 + 0.01 * (t2 - t1), t2 - (t2 - t1) * 0.01, t2)

    wa_tilde = wa.fourier_transform()
    wb_tilde = wb.fourier_transform()
    if ASDs is not None:
        if total_masses is None:
            total_masses = [1.0]
        for ASD in ASDs:
            for total_mass in total_masses:
                wa_tilde_total_mass = wa_tilde.copy()
                wb_tilde_total_mass = wb_tilde.copy()
                wa_tilde_total_mass.t /= total_mass
                wb_tilde_total_mass.t /= total_mass

                errors[f"mismatch {ASD} {total_mass}"] = compute_mismatch(
                    wa_tilde_total_mass, wb_tilde_total_mass, f1, ASD=ASDs[ASD]
                )

    ell_min = max(wa.ell_min, wb.ell_min)
    ell_max = min(wa.ell_max, wb.ell_max)
    for L in range(ell_min, ell_max + 1):
        for M in range(-L, L + 1):
            absolute_error, norm = compute_L2_norm(
                wa, wb, t1, t2, modes=[(L, M)], modes_for_norm=[(L, M)], normalize=False
            )
            errors[f"(L, M) = {(L, M)} L2 norm error"] = absolute_error
            errors[f"(L, M) = {(L, M)} L2 norm"] = norm

    return errors


def analyze_simulation(
    sim_name,
    analyze_levs=True,
    analyze_extrapolation=True,
    analyze_psi4=True,
    ASDs=None,
    total_masses=None,
    nprocs=None,
):
    """
    Analyze a simulation's waveform agreement across
    Levs, extrapolation orders, and the psi4 constraint.

    For each Lev in a simulation, align the strain waveforms
    (if it's the highest Lev pair, then a 4d optimization is also performed
    and the transformation from the low to high lev is included in the returned dictionary)
    via the independent alignment method. Then, compute mismatches between
    the two aligned waveforms (if it's the highest Lev pair, then the
    relative L² norm error over the two-sphere is computed, as well as the
    individual mode absolute errors and norms). For the
    extrapolation order and psi4 analyses, only mismatches are computed.

    Parameters
    ----------
    sim_name : str
        Simulation name, e.g., "SXS:BBH:2092".
    analyze_levs : bool, optional
        Whether or not to analyze the various Levs.
        Default is True.
    analyze_extrapolation : bool, optional
        Whether or not to analyze the various extrapolation orders.
        Default is True.
    analyze_psi4 : bool, optional
        Whether or not to analyze the psi4 constraint -h.ddot = psi4..
        Default is True.
    ASDs : dict of funcs, optional
        Dictionary of functions mapping frequencies to the ASD of a detector(s).
        Default is no frequency-domain mismatch is calculated.
    total_masses : list of floats, optional
        Total masses in solar masses to use for frequency-domain mismatches.
        Default is 1.
    nprocs : int, optional
        Number of cpus to use.  Default is maximum number.  If -1 is provided,
        then no multiprocessing is performed.

    Returns
    -------
    errors : dict
        Dictionary of the errors described above.
    """
    errors = {}

    sim = load(sim_name)

    # Lev analysis
    if analyze_levs:
        for i, lev in enumerate(sim.lev_numbers[1:][::-1]):
            sim_low_lev = load(f"{sim_name}/Lev{lev - 1}")
            sim_high_lev = load(f"{sim_name}/Lev{lev}")

            w_high_lev = sim_high_lev.h
            if i == 0:
                print(sim_low_lev.lev, sim_high_lev.lev, "4d")
                w_low_lev_prime, transformation, L2_norm, t1, t2 = align_simulations(
                    sim_low_lev, sim_high_lev, alignment_method="4d", nprocs=nprocs
                )
                errors[f"(Lev{lev - 1}, Lev{lev}) 4d"] = compute_error_summary(w_low_lev_prime, w_high_lev, t1, t2)
                errors[f"(Lev{lev - 1}, Lev{lev}) 4d transformation"] = transformation

            print(sim_low_lev.lev, sim_high_lev.lev, "independent alignment")
            w_low_lev_prime, transformation, _, t1, t2 = align_simulations(
                sim_low_lev, sim_high_lev, alignment_method="independent alignment", nprocs=nprocs
            )
            errors[f"(Lev{lev - 1}, Lev{lev})"] = compute_mismatch(w_high_lev, w_low_lev_prime, t1, t2)

    # Extrapolation order analysis
    if analyze_extrapolation:
        for extrapolation in ["N3", "N4", "Outer"]:
            other = load(sim_name, extrapolation=extrapolation)

            w_n2 = sim.h
            w_other = other.h

            t1 = sim.metadata.relaxation_time

            print("N2", extrapolation, "independent alignment")
            w_other_prime, transformation, _, t1, t2 = align_waveforms(
                w_other, w_n2, t1, alignment_method="independent alignment", nprocs=nprocs
            )

            mismatch = compute_mismatch(w_n2, w_other_prime, t1, t2)

            errors[f"(N2, {extrapolation})"] = mismatch

    # Psi4 analysis
    if analyze_psi4:
        print("psi4")
        h_as_psi4 = -sim.h.ddot
        psi4 = sim.psi4

        t1 = sim.metadata.relaxation_time
        t2 = h_as_psi4.t[-1]

        mismatch = compute_mismatch(h_as_psi4, psi4, t1, t2)

        errors[f"(-h.ddot, psi4)"] = mismatch

    return errors


def analyze_simulations(
    sim_names,
    analyze_levs=True,
    analyze_extrapolation=True,
    analyze_psi4=True,
    ASDs=None,
    total_masses=None,
    nprocs=None,
):
    """
    Analyze simulations' waveform agreement across
    Levs, extrapolation orders, and the psi4 constraint.

    For each Lev in a simulation, align the strain waveforms
    (if it's the highest Lev pair, then a 4d optimization is also performed
    and the transformation from the low to high lev is included in the returned dictionary)
    via the independent alignment method. Then, compute mismatches between
    the two aligned waveforms (if it's the highest Lev pair, then the
    relative L² norm error over the two-sphere is computed, as well as the
    individual mode absolute errors and norms). For the
    extrapolation order and psi4 analyses, only mismatches are computed.

    Parameters
    ----------
    sim_names : list of strs
        Simulation names, e.g., ["SXS:BBH:2092"].
    analyze_levs : bool, optional
        Whether or not to analyze the various Levs.
        Default is True.
    analyze_extrapolation : bool, optional
        Whether or not to analyze the various extrapolation orders.
        Default is True.
    analyze_psi4 : bool, optional
        Whether or not to analyze the psi4 constraint -h.ddot = psi4..
        Default is True.
    ASDs : dict of funcs, optional
        Dictionary of functions mapping frequencies to the ASD of a detector(s).
        Default is no frequency-domain mismatch is calculated.
    total_masses : list of floats, optional
        Total masses in solar masses to use for frequency-domain mismatches.
        Default is 1.
    nprocs : int, optional
        Number of cpus to use.  Default is maximum number.  If -1 is provided,
        then no multiprocessing is performed.

    Returns
    -------
    errors : dict
        Dictionary of the errors described above.
    """
    errors = {}
    if nprocs != -1:
        with multiprocessing.Pool(processes=nprocs) as pool:
            results = pool.starmap(
                analyze_simulation,
                [
                    (sim_name, analyze_levs, analyze_extrapolation, analyze_psi4, ASDs, total_masses, -1)
                    for sim_name in sim_names
                ],
            )
        print(results)
        for i, sim_name in enumerate(sim_names):
            errors[sim_name] = results[i]
    else:
        for sim_name in sim_names:
            errors[sim_name] = analyze_simulation(
                sim_name, analyze_levs, analyze_extrapolation, analyze_psi4, ASDs, total_masses, -1
            )

    return errors

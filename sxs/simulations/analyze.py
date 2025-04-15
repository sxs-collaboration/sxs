import os
import json
import multiprocessing

from ..waveforms.norms import create_unified_waveforms, L2_difference, mismatch
from ..waveforms.alignment import align_waveforms, align_simulations
from ..handlers import load

import numpy as np


def compute_error_summary(wa, wb, t1, t2, modes=None, ASDs_and_total_masses=None):
    """
    Compute various errors between two waveforms.

    This computes the time-domain mismatch over the two-sphere,
    the normalized L² norm of the residual, the frequency-domain mismatch
    against whatever detector ASDs are provided with some total mass,
    and the individual modes' residual norms and norms.

    Parameters
    ----------
    wa : WaveformModes
    wb : WaveformModes
    t1 : float
        Beginning of integrals.
    t2 : float
        End of integrals.
    modes : list, optional
        Modes (ell, m) to include in error calculations.
        Default is all modes.
    ASDs_and_total_masses : dict of funcs, optional
        Dictionary of functions mapping frequencies to the ASD of a detector(s)
        and total masses to use for mismatch calculations with those ASDs, e.g.,
        {
            'CE': [CE_ASD, total_masses]
        }
        Default is no frequency-domain mismatch is calculated.

    Returns
    -------
    errors : dict
        Dictionary of the time-domain mismatch over the two-sphere,
        the normalized L² norm of the residual, the frequency-domain mismatch
        against whatever detector ASDs are provided with some total mass,
        and the modes' absolute errors and norms.
    """
    from .. import m_sun_in_seconds
    
    errors = {}

    errors["t1"] = t1
    errors["t2"] = t2

    wa, wb = create_unified_waveforms(wa, wb, t1, t2, padding_time_factor=0)

    errors["mismatch"] = mismatch(wa, wb, t1, t2, modes=modes)

    errors["residual L2 norm"] = L2_difference(wa, wb, t1, t2, modes=modes)

    i1 = wa.index_closest_to(t1)
    i0, i2 = max(0, i1-5), min(i1+6, wa.n_times-1)
    Ω1 = np.linalg.norm(wa[i0:i2].angular_velocity, axis=1)[5]
    f1 = 2*Ω1 / (2*np.pi)

    wa = wa.preprocess(t1, t1 + 0.01 * (t2 - t1), t2 - (t2 - t1) * 0.01, t2)
    wb = wb.preprocess(t1, t1 + 0.01 * (t2 - t1), t2 - (t2 - t1) * 0.01, t2)

    wa_tilde = wa.fourier_transform()
    wb_tilde = wb.fourier_transform()
    if ASDs_and_total_masses is not None:
        for ASD_name, (ASD_data, total_masses) in ASDs_and_total_masses.items():
            for total_mass in total_masses:
                # Note that we only need to make the frequency unitful, since the
                # magnitude of the strain scales out in the mismatch.
                # We make things unitful here because the fourier transform
                # earlier was called outside of the for loop without a total mass
                # so that it doesn't need to be computed at each iteration.
                frequency_factor = 1 / (total_mass * m_sun_in_seconds)
            
                wa_tilde_total_mass = wa_tilde.copy()
                wb_tilde_total_mass = wb_tilde.copy()
                wa_tilde_total_mass.t = wa_tilde_total_mass.t * frequency_factor
                wb_tilde_total_mass.t = wb_tilde_total_mass.t * frequency_factor

                errors[f"mismatch {ASD_name} {total_mass}"] = mismatch(
                    wa_tilde_total_mass, wb_tilde_total_mass, f1 * frequency_factor, modes=modes, ASD=ASD_data
                )

    ell_min = max(wa.ell_min, wb.ell_min)
    ell_max = min(wa.ell_max, wb.ell_max)
    for L in range(ell_min, ell_max + 1):
        for M in range(-L, L + 1):
            absolute_error, norm = L2_difference(
                wa, wb, t1, t2, modes=[(L, M)], modes_for_norm=[(L, M)], normalize=False
            )
            errors[f"(L, M) = {(L, M)} residual L2 norm"] = absolute_error
            errors[f"(L, M) = {(L, M)} L2 norm"] = norm

    return errors


def analyze_simulation(
    sim_name,
    analyze_levs=True,
    analyze_extrapolation=True,
    analyze_psi4=True,
    ASDs_and_total_masses=None,
    path_to_analysis_cache=None,
    nprocs=None,
):
    """
    Analyze a simulation's waveform agreement across
    Levs, extrapolation orders, and the psi4 constraint.

    For each Lev in a simulation, align the strain waveforms
    (if it's the highest Lev pair, then a 4d optimization is also performed
    and the transformation from the low to high lev is included in the returned dictionary)
    via the independent alignment method. Then, compute the L² norm of the residual between
    the two aligned waveforms (if it's the highest Lev pair, then the
    mismatch against the ASDs and total masses is also computed, as well as the
    individual modes' residual norms and norms). For the
    extrapolation order and psi4 analyses, only the L² norm of the residuals are computed.

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
    ASDs_and_total_masses : dict of funcs, optional
        Dictionary of functions mapping frequencies to the ASD of a detector(s)
        and total masses to use for mismatch calculations with those ASDs, e.g.,
        {
            'CE': [CE_ASD, total_masses]
        }
        Default is no frequency-domain mismatch is calculated.
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
    if analyze_levs and len(sim.lev_numbers) > 1:
        for i, (low_lev, high_lev) in enumerate(zip(sim.lev_numbers[:-1][::-1], sim.lev_numbers[1:][::-1])):
            sim_low_lev = load(f"{sim_name}/Lev{low_lev}")
            sim_high_lev = load(f"{sim_name}/Lev{high_lev}")

            w_high_lev = sim_high_lev.h
            if i == 0:
                w_low_lev_prime, transformation, L2_norm, t1, t2 = align_simulations(
                    sim_low_lev, sim_high_lev, alignment_method="4d", nprocs=nprocs
                )
                errors[f"(Lev{low_lev}, Lev{high_lev}) 4d"] = compute_error_summary(w_low_lev_prime, w_high_lev, t1, t2, ASDs_and_total_masses=ASDs_and_total_masses)
                errors[f"(Lev{low_lev}, Lev{high_lev}) 4d transformation"] = transformation

            w_low_lev_prime, transformation, _, t1, t2 = align_simulations(
                sim_low_lev, sim_high_lev, alignment_method="independent alignment", nprocs=nprocs
            )
            errors[f"(Lev{low_lev}, Lev{high_lev})"] = L2_difference(w_high_lev, w_low_lev_prime, t1, t2)

    # Extrapolation order analysis
    if analyze_extrapolation:
        for extrapolation in ["N3", "N4", "Outer"]:
            other = load(sim_name, extrapolation=extrapolation)

            w_n2 = sim.h
            w_other = other.h
                            
            t1 = sim.metadata.relaxation_time

            w_other_prime, transformation, _, t1, t2 = align_waveforms(
                w_other, w_n2, t1, alignment_method="independent alignment", nprocs=nprocs
            )

            L2_norm = L2_difference(w_n2, w_other_prime, t1, t2)

            errors[f"(N2, {extrapolation})"] = L2_norm

    # Psi4 analysis
    if analyze_psi4:
        h_as_psi4 = -sim.h.ddot
        psi4 = sim.psi4

        t1 = sim.metadata.relaxation_time
        t2 = h_as_psi4.t[-1]

        L2_norm = L2_difference(h_as_psi4, psi4, t1, t2)

        errors[f"(-h.ddot, psi4)"] = L2_norm

    if path_to_analysis_cache is not None:
        if not os.path.exists(path_to_analysis_cache):
            os.makedirs(path_to_analysis_cache)

        def default(obj):
            if isinstance(obj, np.ndarray):
                return obj.tolist()
            raise TypeError("Not serializable")

        sim_name = sim_name.replace(':','_')
        with open(f"{path_to_analysis_cache}/{sim_name}.json", "w") as output_file:
            json.dump(
                errors,
                output_file,
                indent=2,
                separators=(",", ": "),
                ensure_ascii=True,
                default=default,
            )

    return errors


def analyze_simulations(
    sim_names,
    analyze_levs=True,
    analyze_extrapolation=True,
    analyze_psi4=True,
    ASDs_and_total_masses=None,
    path_to_analysis_cache=None,
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
    ASDs_and_total_masses : dict of funcs, optional
        Dictionary of functions mapping frequencies to the ASD of a detector(s)
        and total masses to use for mismatch calculations with those ASDs, e.g.,
        {
            'CE': [CE_ASD, total_masses]
        }
        Default is no frequency-domain mismatch is calculated.
    path_to_analysis_cache : str, optional
        Path to directory to cache indvidual simulation analyses.
        If None, no caching is performed.
        Default is no cashing is performed.
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
                    (sim_name, analyze_levs, analyze_extrapolation, analyze_psi4, ASDs_and_total_masses, path_to_analysis_cache, -1)
                    for sim_name in sim_names
                ],
            )
        for i, sim_name in enumerate(sim_names):
            errors[sim_name] = results[i]
    else:
        for sim_name in sim_names:
            errors[sim_name] = analyze_simulation(
                sim_name, analyze_levs, analyze_extrapolation, analyze_psi4, ASDs_and_total_masses, path_to_analysis_cache, -1
            )

    return errors

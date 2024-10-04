from pathlib import Path
from warnings import warn
from .. import doi_url, Metadata
from ..utilities import (
    sxs_id_and_version, lev_number, sxs_path_to_system_path,
    download_file, sxs_directory, read_config
)


def Simulation(location, *args, **kwargs):
    """Construct a Simulation object from a location string

    The location string should be an SXS ID with an optional version
    number, and possibly a Lev specification, as in
    "SXS:BBH:1234v2.0/Lev5".  The default version number is "v2.0",
    and the default Lev number is the highest available.  Note that
    the version number must be in the form "v2.0", not "2.0", and it
    must exist, or no files will be found to load the data from.

    The returned object will be either a `Simulation_v1` or
    `Simulation_v2` object, depending on the version number.
    Hopefully, most details of the two versions will be hidden from
    the user, so that the interface is identical.

    Note that some simulations are deprecated and/or superseded by
    other simulations.  By default, this function will raise an error
    if you try to load a deprecated or superseded simulation.  There
    are several ways to override this behavior:

    1. Pass `ignore_deprecation=True` to completely bypass even
        checking for deprecation or supersession.  No warnings or
        errors will be issued.
    2. Include an explicit version number in the `location`
        string, as in "SXS:BBH:0001v2.0".  A warning will be issued
        that the simulation is deprecated, but it will be loaded
        anyway.
    3. Pass `auto_supersede=True` to automatically load the
        superseding simulation, if there is only one.  Because no
        superseding simulation can be *precisely* the same as the
        deprecated one, there may be multiple superseding simulations
        that have very similar parameters, in which case an error will
        be raised and you must explicitly choose one.  If there is
        only one, a warning will be issued, but the superseding
        simulation will be loaded.

    Otherwise, a `ValueError` will be raised, with an explanation and
    suggestions on what you might want to do.

    Parameters
    ----------
    location : str
        The location string for the simulation.  This must include the
        SXS ID, but can also include the version and Lev number, as
        described above.

    Other parameters
    ----------------
    ignore_deprecation : bool
        If `True`, completely bypass checking for deprecation or
        supersession.  No warnings or errors will be issued.
    auto_supersede : bool
        If `True`, automatically load the superseding simulation, if
        there is only one.  If there are multiple superseding
        simulations, an error will be raised, and you must explicitly
        choose one.  If there is only one, a warning will still be
        issued, but the superseding simulation will be loaded.  Note
        that this can also be set in the configuration file with
        `sxs.write_config(auto_supersede=True)`.
    extrapolation : str
        The extrapolation order to use for the strain and Psi4 data.
        This is only relevant for versions 1 and 2 of the data format,
        both of which default to "N2".  Other options include "N3",
        "N4", and "Outer".  "Nx" refers to extrapolation by
        polynomials in 1/r with degree `x`, while "Outer" refers to
        data extracted at the outermost extraction radius but
        corrected for time-dilation and areal-radius effects.
    download : bool
        If `True`, download the information about the files from the
        Zenodo/CaltechDATA record.  If `False`, only use the file
        information that is already available (which will raise an
        error if the file information has not previously been
        downloaded).  If not present, use the value from the
        configuration file, defaulting to `True` if it is not
        configured.

    Returns
    -------
    simulation : SimulationBase
        A `Simulation_v1` or `Simulation_v2` object, depending on the
        version of the simulation data.

    Note that all remaining arguments (including keyword arguments)
    are passed on to the `SimulationBase`, `Simulation_v1`, and/or
    `Simulation_v2` constructors, though none currently recognize any
    arguments other than those listed above.
    
    """
    from .. import load, sxs_directory

    # Extract the simulation ID, version, and Lev from the location string
    simulation_id, input_version = sxs_id_and_version(location)
    if not simulation_id:
        raise ValueError(f"Invalid SXS ID in '{simulation_id}'")
    input_lev_number = lev_number(location)  # Will be `None` if not present

    # Load the simulation catalog and check if simulation ID exists in the catalog
    simulations = load("simulations")
    if simulation_id not in simulations:
        raise ValueError(f"Simulation '{simulation_id}' not found in simulation catalog")

    # Attach metadata to this object
    metadata = Metadata(simulations[simulation_id])
    series = simulations.dataframe.loc[simulation_id]

    # Check if the specified version exists in the simulation catalog
    if input_version not in metadata.DOI_versions:
        raise ValueError(f"Version '{input_version}' not found in simulation catalog for '{simulation_id}'")

    # Set various pieces of information about the simulation
    version = input_version or max(metadata.DOI_versions)
    if not version.startswith("v"):
        raise ValueError(f"Invalid version string '{version}'")
    sxs_id_stem = simulation_id
    sxs_id = f"{sxs_id_stem}{version}"
    url = f"{doi_url}{sxs_id}"

    # Deal with "superseded_by" field, or "deprecated" keyword in the metadata
    deprecated = ("deprecated" in metadata.get("keywords", []) or metadata.get("superseded_by", False))
    if not kwargs.get("ignore_deprecation", False):
        auto_supersede = kwargs.get("auto_supersede", read_config("auto_supersede", False))
        if (
            input_version
            and not auto_supersede
            and deprecated
        ):
            message = ("\n"
                + f"Simulation '{sxs_id_stem}' is deprecated and/or superseded.\n"
                + "Normally, this simulation should no longer be used, but you\n"
                + f"explicitly requested version '{input_version}', so it is being used.\n"
            )
            warn(message)
        else:
            if "superseded_by" in metadata:
                superseded_by = metadata["superseded_by"]
                if auto_supersede and isinstance(superseded_by, list):
                    raise ValueError(
                        f"`auto_supersede` is enabled, but simulation '{sxs_id}' is\n"
                        + "superseded by multiple simulations.  You must choose one\n"
                        + "explicitly from the list:\n"
                        + "\n".join(f"  {s}" for s in superseded_by)
                        + "\nAlternatively, you could pass `ignore_deprecation=True` or\n"
                        + "specify a version to load this waveform anyway."
                    )
                elif auto_supersede and isinstance(superseded_by, str):
                    message = f"\nSimulation '{sxs_id}' is being automatically superseded by '{superseded_by}'."
                    warn(message)
                    new_location = f"{superseded_by}{input_version}"
                    if input_lev_number:
                        new_location += f"/Lev{input_lev_number}"
                    return Simulation(new_location, *args, **kwargs)
                elif isinstance(superseded_by, list):
                    raise ValueError(
                        f"Simulation '{sxs_id}' is superseded by multiple simulations.\n"
                        + "Even if you enable `auto_supersede`, with multiple options, you\n"
                        + "must choose one explicitly from the list:\n"
                        + "\n".join(f"  {s}" for s in superseded_by)
                        + "\nAlternatively, you could pass `ignore_deprecation=True` or\n"
                        + "specify a version to load this waveform anyway."
                    )
                elif isinstance(superseded_by, str):
                    raise ValueError(
                        f"Simulation '{sxs_id}' is superseded by '{superseded_by}'.\n"
                        + "Note that you could enable `auto_supersede` to automatically\n"
                        + "load the superseding simulation.  Alternatively, you could\n"
                        + "pass `ignore_deprecation=True` or specify a version to load\n"
                        + "this waveform anyway."
                    )
                else:
                    raise ValueError(
                        f"Simulation '{sxs_id}' is superseded by '{superseded_by}'.\n"
                        + "Note that you could pass `ignore_deprecation=True` or\n"
                        + "specify a version to load this waveform anyway."
                    )
            if "deprecated" in metadata.get("keywords", []):
                raise ValueError(
                    f"Simulation '{sxs_id}' is deprecated but has no superseding simulation.\n"
                    + "Note that you could pass `ignore_deprecation=True` or specify a version\n"
                    + "to  to load this waveform anyway."
                )

    # Note the deprecation status in the kwargs, even if ignoring deprecation
    kwargs["deprecated"] = deprecated

    # We want to do this *after* deprecation checking, to avoid possibly unnecessary web requests
    files = get_file_info(metadata, sxs_id, download=kwargs.get("download", None))

    # If Lev is given as part of `location`, use it; otherwise, use the highest available
    lev_numbers = sorted({lev for f in files if (lev:=lev_number(f))})
    output_lev_number = input_lev_number or max(lev_numbers)
    location = f"{sxs_id_stem}{version}/Lev{output_lev_number}"

    # Finally, figure out which version of the simulation to load and dispatch
    version_number = float(version[1:])
    if 1 <= version_number < 2.0:
        sim = Simulation_v1(
            metadata, series, version, sxs_id_stem, sxs_id, url, files, lev_numbers, output_lev_number, location, *args, **kwargs
        )
    elif 2 <= version_number < 3.0:
        sim = Simulation_v2(
            metadata, series, version, sxs_id_stem, sxs_id, url, files, lev_numbers, output_lev_number, location, *args, **kwargs
        )
    else:
        raise ValueError(f"Version '{version}' not yet supported")
    sim.__file__ = str(sxs_directory("cache") / sxs_path_to_system_path(sim.sxs_id))
    return sim


class SimulationBase:
    """Base class for Simulation objects

    Note that users almost certainly never need to call this function;
    see the `Simulation` function or `sxs.load` function instead.
    
    Attributes
    ----------
    metadata : Metadata
        Metadata object for the simulation
    series : pandas.Series
        The metadata, as extracted from the `simulations.dataframe`,
        meaning that it has columns consistent with other simulations,
        even when the underlying Metadata objects do not.  Note that
        `metadata` is an alias for this attribute, just based on the
        use of that name for `simulations`, but technically `pandas`
        distinguishes a single row like this as a `Series` object.
    version : str
        Version number of the simulation
    sxs_id_stem : str
        SXS ID without the version number or Lev
    sxs_id : str
        SXS ID with the version number
    location : str
        Location string for the simulation, including the SXS ID,
        version number, and Lev number.
    url : str
        URL for the DOI of the simulation
    files : dict
        Dictionary of file information for the simulation.  The keys
        are of the form "Lev5:Horizons.h5", and the values are
        dictionaries with keys "checksum", "size", and "link".
    lev_numbers : list
        List of available Lev numbers for the simulation.
    lev_number : int
        Chosen Lev number for the simulation.
    horizons : Horizons
        Horizons object for the simulation
    strain : Waveform
        Strain Waveform object for the simulation.  Note that `h` is
        an alias for this attribute, both of which are case
        insensitive: `Strain` and `H` are also acceptable.
    psi4 : Waveform
        Psi4 Waveform object for the simulation.  Note that this
        attribute is also case insensitive: `Psi4` is also acceptable.
    psi3 : Waveform
        Psi3 Waveform object for the simulation.  Note that this
        attribute is also case insensitive: `Psi3` is also acceptable.
        In versions 1 and 2, this attribute will raise an error
        because the data is not available.
    psi2 : Waveform
        Psi2 Waveform object for the simulation.  Note that this
        attribute is also case insensitive: `Psi2` is also acceptable.
        In versions 1 and 2, this attribute will raise an error
        because the data is not available.
    psi1 : Waveform
        Psi1 Waveform object for the simulation.  Note that this
        attribute is also case insensitive: `Psi1` is also acceptable.
        In versions 1 and 2, this attribute will raise an error
        because the data is not available.
    psi0 : Waveform
        Psi0 Waveform object for the simulation.  Note that this
        attribute is also case insensitive: `Psi0` is also acceptable.
        In versions 1 and 2, this attribute will raise an error
        because the data is not available.
    """
    def __init__(self,
        metadata, series, version, sxs_id_stem, sxs_id, url, files, lev_numbers, lev_number, location,
        *args, **kwargs
    ):
        self.metadata = metadata
        self.series = series
        self.version = version
        self.sxs_id_stem = sxs_id_stem
        self.sxs_id = sxs_id
        self.url = url
        self.files = files
        self.lev_numbers = lev_numbers
        self.lev_number = lev_number
        self.location = location
        self.deprecated = kwargs.get("deprecated", False)

    def __repr__(self):
        chi1 = self.series["reference_dimensionless_spin1"]
        chi2 = self.series["reference_dimensionless_spin2"]
        e = self.metadata.reference_eccentricity
        construction = f"""{type(self).__qualname__}("{self.location}")\n# """
        if self.deprecated:
            construction += "DEPRECATED "
        construction += f"n_orbits={self.metadata.number_of_orbits:.3g} "
        construction += f"q={self.metadata.reference_mass_ratio:.3g} "
        construction += f"""chi1=[{", ".join(f"{c:.3g}" for c in chi1)}] """
        construction += f"""chi2=[{", ".join(f"{c:.3g}" for c in chi2)}] """
        construction += f"e={e:.3g} simulation" if type(e) is float else f"{e=} simulation"
        return construction

    def __str__(self):
        return repr(self)
    
    @property
    def dataframe(self):
        return self.series

    @property
    def versions(self):
        return self.metadata.DOI_versions

    @property
    def lev(self):
        return f"Lev{self.lev_number}"

    @property
    def Lev(self):
        return self.lev

    def load_horizons(self):
        from .. import load
        sxs_id_path = Path(self.sxs_id)
        horizons_path = self.horizons_path
        if horizons_path in self.files:
            horizons_location = self.files.get(horizons_path)["link"]
        else:
            # Some simulations used the SXS ID as a prefix in file paths
            # within the Zenodo upload in version 1.x of the catalog.
            if (extended_location := f"{self.sxs_id_stem}/{horizons_path}") in self.files:
                horizons_location = self.files.get(extended_location)["link"]
            else:
                raise ValueError(
                    f"File '{horizons_path}' not found in simulation files for {self.location}"
                )
        horizons_truepath = Path(sxs_path_to_system_path(sxs_id_path / horizons_path))
        return load(horizons_location, truepath=horizons_truepath)

    @property
    def horizons(self):
        if not hasattr(self, "_horizons"):
            self._horizons = self.load_horizons()
        return self._horizons
    Horizons = horizons

    @property
    def strain(self):
        if not hasattr(self, "_strain"):
            self._strain = self.load_waveform(*self.strain_path)
        return self._strain
    Strain = strain
    h = strain
    H = strain

    # I'm not entirely sure about the conjugations and factors of 2 in
    # shear and news in our conventions.  These will have to wait for
    # later.
    #
    # @property
    # def shear(self):
    #     if not hasattr(self, "_shear"):
    #         self._shear = self.strain.bar / 2
    #     return self._shear
    # sigma = shear
    # σ = shear
    # Shear = shear
    # Sigma = shear
    # Σ = shear
    #
    # @property
    # def news(self):
    #     if not hasattr(self, "_news"):
    #         self._news = self.strain.dot
    #     return self._news
    # News = news

    @property
    def psi4(self):
        if not hasattr(self, "_psi4"):
            self._psi4 = self.load_waveform(*self.psi4_path)
        return self._psi4
    Psi4 = psi4

    @property
    def psi3(self):
        raise AttributeError(f"Psi3 is not available for version {self.version} of the data")
    Psi3 = psi3

    @property
    def psi2(self):
        raise AttributeError(f"Psi2 is not available for version {self.version} of the data")
    Psi2 = psi2

    @property
    def psi1(self):
        raise AttributeError(f"Psi1 is not available for version {self.version} of the data")
    Psi1 = psi1

    @property
    def psi0(self):
        raise AttributeError(f"Psi0 is not available for version {self.version} of the data")
    Psi0 = psi0

    def to_lvk(self, **kwargs):
        r"""Convert an SXS simulation to LVK convention.

        Returns an SXS waveform (modes or polarizations) and dynamics
        (including angular velocities, frame quaternions, and spins) in
        the inertial frame that coincides with the waveform-defined frame
        defined at a reference time `t_ref` or reference frequency
        `f_ref`.

        Parameters
        ----------
        t_ref : float, optional
            The reference time at which the waveform frame is specified.
            This is measured in units of M, and defined relative to the
            epoch time (see below).  Either `t_ref` or `f_ref` must be
            specified.  If `t_ref` is given, it is used to compute
            `f_ref`.
        f_ref : float, optional
            The reference frequency, in units of cycles/M, at which the
            waveform frame is specified.  Either `t_ref` or `f_ref` must
            be specified.  If `f_ref` is given, it is used to compute
            `t_ref`.
        dt : float, optional
            The time step, in units of M, to which to interpolate the
            waveform.
        f_low : float, optional
            The lower frequency bound, in units of cycles/M, for the
            waveform.
        ell_max : int, optional
            The maximum ell to include in the waveform.
        phi_ref : float, optional
            The binary's phase in the coprecessing frame, measured at
            `t_ref`. Should be between 0 and $2\pi$.
        inclination : float, optional
            Angle between the binary's angular momentum and the line of
            sight of the observer, measured at `t_ref`.  Should be between
            0 and $\pi$.
        ell_max_epoch : int, optional
            The maximum ell to include in the epoch time calculation,
            which sets t=0 at the maximum of the L^2 norm, calculated by
            including all modes up to and including this ell value.

        Returns
        -------
        times : float array
            Uniformly spaced 1D array of times, in units of M, at which
            the waveform and dynamics quantities are returned.  Aligned
            such that peak of waveform modes with ell=2 is at t=0.
        hlm_dict : dict [optional]
            Dictionary of waveform modes in the inertial frame that
            coincides with the waveform-defined coprecessing frame at
            `f_ref`.  Each mode in the dictionary is a 1D array of
            complex-valued floats with values corresponding to each time
            node.  Keys:[(ell,m)] for all ell<=ell_max and -ell<=m<=+ell.
            This is returned only if the input values of `phi_ref` and
            `inclination` are both None.
        hp, hc : float arrays [optional]
            1D-arrays of real-valued GW polarizations evaluated in the
            frame of the observer at each time node.  Polarizations are
            computed using all modes up to ell=ell_max, with one value at
            each time node.  These are returned only if either of the
            input values of `phi_ref` and `inclination` is not None.
        dynamics_dict : dict
            Dictionary of real-valued arrays of dynamics quantities:
                * "t_ref": The reference time at which the waveform frame
                is specified.
                * "f_ref": The waveform's frequency at `t_ref`.
                * "t_low": The earliest time in the waveform.
                * "f_low": The waveform's frequency at `t_low`.
                * "chi1_ref": Cartesian spin components for the more
                massive object, evaluated at `t_ref` in the
                waveform-defined inertial frame.
                * "chi2_ref": Cartesian spin components for the less
                massive object, evaluated at `t_ref` in the
                waveform-defined inertial frame.
                * "frame_quat": Quaternions describing the instantaneous
                transformation from the inertial frame to the corotating
                frame at each time node as described in arXiv:1905.09300
                and using conventions in Appendix B of arXiv:1110.2965.
                Array of shape (len(times),4) with four quaternion
                components at each time node.  The first element at each
                time node is the scalar part.
                * "frame_omega": Angular velocity vector of the corotating
                frame at each time node.  Array of shape (len(times),3)
                with three components at each time node.
                * "times_spins": The times, in units of M, at which `chi`
                and `chi2` are returned.  Note that these are coordinate
                times deep within the dynamic region of the simulation,
                and so cannot be precisely related to the times in the
                asymptotic waveform.
                * "chi1": Cartesian spin components for the more massive
                object, evaluated in the waveform-defined inertial
                frame.  Array of shape (len(times_spins),3) which
                contains the Cartesian spin components
                \{chi_{1x},chi_{1y},chi_{1z}\} at each time node.
                * "chi2": Cartesian spin components for the less massive
                object, evaluated in the waveform-defined inertial
                frame.  Array of shape (len(times_spins),3) which
                contains the Cartesian spin components
                \{chi_{2x},chi_{2y},chi_{2z}\} at each time node.
                * "times_remnant": The times, in units of M, at which
                `chi_remnant` and `mass_remnant` are returned.  Note
                that these are coordinate times deep within the dynamic
                region of the simulation, and so cannot be precisely
                related to the times in the asymptotic waveform.
                * "chi_remnant": Cartesian spin components for the remnant
                black hole, evaluated in the waveform-defined inertial
                frame.  Array of shape (len(times_remnant),3) which
                contains the Cartesian spin components
                \{chi_{rx},chi_{ry},chi_{rz}\}.
                * "mass_remnant": The Christodoulou mass of the remnant
                black hole as a function of time.  Array of shape
                (len(times_remnant),).

        See also
        ========
        sxs.load : General-purpose function to load SXS data in native
            format
        sxs.waveforms.to_lvc_conventions : Inner function that does all
            the work for this function

        Conventions
        ===========
        We assume geometric units for time (units of M) and frequency
        (units of cycles/M), with total mass M equal to 1.
        
        Epoch time is defined by the peak of the $L^2$ norm of the modes:

            $$t_e = \argmax_t \sum_\ell \sum_m |h_{\ell,m}(t)|^2,$$

        where the sum over $\ell$ ranges from 2 to `ell_max_epoch`.  All
        time axes are then adjusted so that $t_e = 0$.

        Frequencies are measured from the waveform, rather than orbital
        trajectories, in terms of the angular-velocity vector given by
        equation (7) of arXiv:1302.2919.  The frequency is defined as the
        magnitude of this vector divided by $2\pi$.

        Waveforms, spins, dynamics, times, inclination angle and GW
        reference phase are all defined in the inertial frame that
        coincides with the "waveform-defined frame" at the reference time.
        This frame is chosen so that the $z$ axis is aligned with the
        dominant principal axis of the matrix given by equation (2a) of
        arXiv:1109.5224, except that the strain is used in place of
        $\psi_4$.  That axis is only determined up to a sign, so we choose
        the positive $z$ axis to be more parallel than antiparallel to the
        angular velocity.  Rotation about the $z$ axis is chosen to
        approximate the condition that the more massive black hole is
        located on the positive $x$ axis at the reference time, but can be
        written solely in terms of the waveform modes:

            * $\Im{h_{2,2} + \bar{h}_{2,-2}} = 0$
            * $\Re{h_{2,2} + \bar{h}_{2,-2}} < 0$
            * $\Im{h_{2,1} + \bar{h}_{2,-1}} < 0$

        The first two conditions are necessary for cases of symmetric
        systems in which $h_{2,\pm 1}=0$; the last condition breaks the
        degeneracy of the first two under rotation about the $z$ axis by
        $\pi$.  For configurations that are symmetric under that rotation,
        $h_{2,1}$ will be zero, so this condition will be impossible to
        apply, but the symmetry of the system will mean that there is no
        difference in the result.

        Quaternion conventions are described in Appendix B of
        arXiv:1110.2965.  In particular, we use the convention that

            Q = q_0 + vec(q) -> (q_0, q_1, q_2, q_3)

        where q_0 is the scalar part.

        """
        from ..waveforms.format_handlers.lvc import to_lvc_conventions
        strain = self.load_waveform(*self.strain_path, transform_to_inertial=False)
        return to_lvc_conventions(strain, self.horizons, **kwargs)


class Simulation_v1(SimulationBase):
    """Simulation object for version 1 of the data format
    
    Note that users almost certainly never need to call this function;
    see the `Simulation` function or `sxs.load` function instead.  See
    also `SimulationBase` for the base class that this class inherits
    from.
    """
    # We have to deal with the fact that some early file paths on
    # Zenodo included the SXS ID as a prefix, while others did not.
    # This means that we have to check for both possibilities in
    # `load_horizons` and `load_waveform`.

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.extrapolation = kwargs.get("extrapolation", "N2")

    @property
    def horizons_path(self):
        return f"{self.lev}/Horizons.h5"
    
    def load_horizons(self):
        from .. import load
        sxs_id_path = Path(self.sxs_id)
        horizons_path = self.horizons_path
        if horizons_path in self.files:
            horizons_location = self.files.get(horizons_path)["link"]
        else:
            # Some simulations used the SXS ID as a prefix in file paths
            # within the Zenodo upload in version 1.x of the catalog.
            if (extended_horizons_path := f"{self.sxs_id_stem}/{horizons_path}") in self.files:
                horizons_location = self.files.get(extended_horizons_path)["link"]
            else:
                raise ValueError(
                    f"File '{horizons_path}' not found in simulation files for {self.location}"
                )
        horizons_truepath = Path(sxs_path_to_system_path(sxs_id_path / horizons_path))
        return load(horizons_location, truepath=horizons_truepath)

    @property
    def strain_path(self):
        extrapolation = (
            f"Extrapolated_{self.extrapolation}.dir"
            if self.extrapolation != "Outer"
            else "OutermostExtraction.dir"
        )
        return (
            f"{self.lev}/rhOverM_Asymptotic_GeometricUnits_CoM.h5",
            extrapolation
        )

    @property
    def psi4_path(self):
        extrapolation = (
            f"Extrapolated_{self.extrapolation}.dir"
            if self.extrapolation != "Outer"
            else "OutermostExtraction.dir"
        )
        return (
            f"{self.lev}/rMPsi4_Asymptotic_GeometricUnits_CoM.h5",
            extrapolation
        )

    def load_waveform(self, file_name, group, transform_to_inertial=True):
        # Note: `transform_to_inertial` is a slightly unnatural argument in version 1,
        # since the data are actually stored in the inertial frame.  If the value is
        # `False`, we will transform the data to the corotating frame, which is what
        # the corresponding argument in version 2 achieves naturally (since the data
        # are stored in the corotating frame in that version).
        from .. import load
        if file_name in self.files:
            location = self.files.get(file_name)["link"]
        else:
            # Some simulations used the SXS ID as a prefix in file paths
            # within the Zenodo upload in version 1.x of the catalog.
            if (extended_file_name := f"{self.sxs_id_stem}/{file_name}") in self.files:
                location = self.files.get(extended_file_name)["link"]
            else:
                raise ValueError(f"File '{file_name}' not found in simulation files")
        sxs_id_path = Path(self.sxs_id)
        truepath = Path(sxs_path_to_system_path(sxs_id_path / file_name))
        w = load(
            location, truepath=truepath, extrapolation_order=group,
            transform_to_inertial=transform_to_inertial
        )
        w.metadata = self.metadata
        if not transform_to_inertial:
            w = w.to_corotating_frame()
        return w


class Simulation_v2(SimulationBase):
    """Simulation object for version 2 of the data format
        
    Note that users almost certainly never need to call this function;
    see the `Simulation` function or `sxs.load` function instead.  See
    also `SimulationBase` for the base class that this class inherits
    from.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.extrapolation = kwargs.get("extrapolation", "N2")

    @property
    def horizons_path(self):
        return f"{self.lev}:Horizons.h5"

    @property
    def strain_path(self):
        return (
            f"{self.lev}:Strain_{self.extrapolation}",
            "/"
        )

    @property
    def psi4_path(self):
        extrapolation = (
            f"Extrapolated_{self.extrapolation}.dir"
            if self.extrapolation != "Outer"
            else "OutermostExtraction.dir"
        )
        return (
            f"{self.lev}:ExtraWaveforms",
            f"/rMPsi4_Asymptotic_GeometricUnits_CoM_Mem/{extrapolation}"
        )

    def load_waveform(self, file_name, group, transform_to_inertial=True):
        from .. import load
        # Note that `name` should not have the file ending on input,
        # but we will replace it regardless with `.with_suffix`.
        file_name = Path(file_name)
        sxs_id_path = Path(self.sxs_id)
        h5_path = str(file_name.with_suffix(".h5"))
        json_path = str(file_name.with_suffix(".json"))
        h5_location = self.files.get(h5_path)["link"]
        json_location = self.files.get(json_path)["link"]
        h5_truepath = Path(sxs_path_to_system_path(sxs_id_path / h5_path))
        json_truepath = Path(sxs_path_to_system_path(sxs_id_path / json_path))
        if not Path(json_location).exists() and not json_truepath.exists():
            if not read_config("download", True):
                raise ValueError(f"{json_truepath} not found and download is disabled")
            download_file(json_location, sxs_directory("cache") / json_truepath)
        return load(
            h5_location, truepath=h5_truepath, group=group, metadata=self.metadata,
            transform_to_inertial=transform_to_inertial
        )


def get_file_info(metadata, sxs_id, download=None):
    from .. import load_via_sxs_id
    if "files" in metadata:
        return metadata["files"]
    truepath = Path(sxs_path_to_system_path(sxs_id)) / "zenodo_metadata.json"
    record = load_via_sxs_id(sxs_id, "export/json", truepath=truepath, download=download)
    entries = record["files"]["entries"]
    return {
        str(filename): {
            "checksum": entry["checksum"],
            "size": entry["size"],
            "link": entry["links"]["content"],
        }
        for filename, entry in entries.items()
    }

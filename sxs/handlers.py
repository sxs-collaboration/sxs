"""Functions to facilitate generic handling of SXS-format data files"""

import contextlib
from . import waveforms


def sxs_handler(format_string):
    """Find an object to load from or save to files in the given format.

    Parameters
    ----------
    format_string : str

    Returns
    -------
    handler : object
        This object will have (at least) two attributes: `load` and `save`, which
        can be called as `handler.load(file, **kwargs)` and `handler.save(obj,
        file, **kwargs)`.

    See Also
    --------
    sxs_loader : Returns the function that will load a given file
    sxs.utilities.file_format : Returns just the string found in the file

    """
    import itertools
    import re
    from . import catalog, metadata, horizons, waveforms

    if not format_string:
        raise ValueError("Empty string cannot be associated with a handler")
    elif format_string.lower().startswith("catalog"):
        format_string = re.sub(r"^catalog\.?", "", format_string, count=1, flags=re.IGNORECASE)
        return catalog.formats.get(format_string, catalog.formats[None])
    elif format_string.lower().startswith("metadata"):
        format_string = re.sub(r"^metadata\.?", "", format_string, count=1, flags=re.IGNORECASE)
        return metadata.formats.get(format_string, metadata.formats[None])
    elif format_string.lower().startswith("horizons"):
        format_string = re.sub(r"^horizons\.?", "", format_string, count=1, flags=re.IGNORECASE)
        return horizons.formats.get(format_string, horizons.formats[None])
    elif format_string.lower().startswith("waveforms"):
        format_string = re.sub(r"^waveforms\.?", "", format_string, count=1, flags=re.IGNORECASE)
        return waveforms.formats.get(format_string, waveforms.formats[None])
    else:
        format_list = [
            catalog.formats,
            metadata.formats,
            horizons.formats,
            waveforms.formats,
        ]
        format_cycler = itertools.cycle(format_list)
        for _ in range(len(format_list)):
            format_dict = next(format_cycler)
            if format_string in format_dict:
                if any(format_string in next(format_cycler) for _ in range(len(format_list)-1)):
                    raise ValueError(f"Format string '{format_string}' found in multiple sxs format groups")
                return format_dict[format_string]
    raise ValueError(f"Format '{format_string}' is unknown to the `sxs` package; maybe you need to update `sxs`")


def sxs_loader(file, group=None):
    """Find the function that will load the given file

    If a format is specified, we assume that it is a top-level attribute or member
    named "sxs_format" or just "format".  If neither exists, return None; the
    calling function should check for this possibility.

    Parameters
    ----------
    file : file-like object, string, or pathlib.Path
    group : str, optional
        If the file is an HDF5 or JSON file, this is the group to inspect.

    Returns
    -------
    load : callable
        This function can be called as `load(file, **kwargs)`.

    See Also
    --------
    sxs_handler : Returns an object to load from and save to a given format
    sxs.utilities.file_format : Returns just the string found in the file

    Notes
    -----
    Note that this function has a very different signature than the related
    functions `sxs_handler`, which takes only the format string, not an actual
    file.

    """
    import re
    import pathlib
    from .utilities import file_format
    format_string = file_format(file, group)
    if format_string is None:
        file_string = str(pathlib.Path(file).name).lower()
        if "catalog" in file_string:
            format_string = "catalog"
        elif "metadata" in file_string:
            format_string = "metadata"
        elif "horizons" in file_string:
            format_string = "horizons"
        elif re.match("(rh_|rhoverm_|rpsi4_|rmpsi4_)", file_string):
            format_string = "waveforms"
        else:
            raise ValueError(f"File '{file}' contains no recognized format information")
    handler = sxs_handler(format_string)
    # noinspection PyUnresolvedReferences
    return handler.load


def load(location, download=None, cache=None, progress=None, **kwargs):
    """Load an SXS-format dataset, optionally downloading and caching

    The dataset can be the full catalog of all SXS simulations, or metadata,
    horizon data, or a waveform from an individual simulation.

    Parameters
    ----------
    location : {str, pathlib.Path}
        A local file path, URL, SXS path, or SXS path pattern.  See Notes below.
    download : {None, bool}, optional
        If this is True and the data is recognized as starting with an SXS ID but
        cannot be found in the cache, the data will be downloaded automatically.
        If this is None (the default) and an SXS configuration file is found with a
        `download` key, that value will be used.  If this is False, any
        configuration will be ignored, and no files will be downloaded.  Note that
        if this is True but `cache` is None, `cache` will automatically be switched
        to True.
    cache : {None, bool}, optional
        The cache directory is determined by `sxs.sxs_directory`, and any downloads
        will be stored in that directory.  If this is None (the default) and
        `download` is True it will be set to True.  If this is False, any
        configuration will be ignored and any files will be downloaded to a
        temporary directory that will be deleted when python exits.
    progress : {None, bool}, optional
        If True, full file names will be shown and, if a nonzero Content-Length
        header is returned, a progress bar will be shown during any downloads.
        Default is None, which just reads the configuration value with
        `read_config("download_progress", True)`, defaulting to True.

    Keyword Parameters
    ------------------
    All remaining parameters are passed to the `load` function responsible for the
    requested data.

    See Also
    --------
    sxs.sxs_directory : Locate configuration and cache files
    sxs.write_config : Set defaults for `download` and `cache` parameters

    Notes
    -----
    This function can load data in various ways.

      1) Given an absolute or relative path to a local file, it just loads the data
         directly.

      2) If `location` is a valid URL including the scheme (https://, or http://),
         it will be downloaded regardless of the `download` parameter and
         optionally cached.

      3) Given an SXS path — like 'SXS:BBH:1234/Lev5/h_Extrapolated_N2.h5' — the
         file is located in the catalog for details.  This function then looks in
         the local cache directory and loads it if present.

      4) If the SXS path is not found in the cache directory and `download` is set
         to `True` (when this function is called, or in the sxs config file) this
         function attempts to download the data.  Note that `download` must be
         explicitly set in this case, or a ValueError will be raised.

    Note that downloading is switched off by default, but if it is switched on (set
    to True), the cache is also switched on by default.

    """
    import pathlib
    import urllib.request
    from . import Catalog, read_config, sxs_directory
    from .utilities import url, download_file, sxs_path_to_system_path

    # Note: `download` and/or `cache` may still be `None` after this
    if download is None:
        download = read_config("download", True)
    if cache is None:
        cache = read_config("cache")
    if progress is None:
        progress = read_config("download_progress", True)

    # We set the cache path to be persistent if `cache` is `True` or `None`.  Thus,
    # we test for whether or not `cache` literally *is* `False`, rather than just
    # if it casts to `False`.
    cache_path = sxs_directory("cache", persistent=(cache is not False))

    path = pathlib.Path(sxs_path_to_system_path(location)).expanduser()  # .resolve()
    h5_path = path.with_suffix('.h5')
    json_path = path.with_suffix('.json')

    if not path.exists():
        if h5_path.resolve().exists():
            path = h5_path

        elif json_path.resolve().exists():
            path = json_path

        elif "scheme" in url.parse(location):
            m = url.parse(location)
            path_name = urllib.request.url2pathname(f"{m['host']}/{m['port']}/{m['resource']}")
            path = cache_path / path_name
            if not path.resolve().exists():
                if download is False:  # Again, we want literal False, not casting to False
                    raise ValueError(f"File '{path_name}' not found in cache, but downloading turned off")
                download_file(location, path, progress=progress)

        elif location == "catalog":
            return Catalog.load(download=download)

        else:
            # Try to find an appropriate SXS file
            catalog = Catalog.load(download=download)
            selections = catalog.select_files(location)
            if not selections:
                raise ValueError(f"Nothing found matching '{location}'")
            if progress:
                print("Found the following files to load from the SXS catalog:")
                print("    " + "\n    ".join(selections))
            paths = []
            for sxs_path, file_info in selections.items():
                truepath = sxs_path_to_system_path(file_info.get("truepath", sxs_path))
                path = cache_path / truepath
                if not path.resolve().exists():
                    download_url = file_info["download"]
                    download_file(download_url, path, progress=progress)
                paths.append(path)
            loaded = [load(path, download=False, progress=progress, **kwargs) for path in paths]
            if len(loaded) == 1:
                return loaded[0]
            else:
                return loaded

    loader = sxs_loader(path, kwargs.get("group", None))

    return loader(path, **kwargs)


@contextlib.contextmanager
def loadcontext(*args, **kwargs):
    """Context manager for backwards compatibility

    This context manager takes precisely the same arguments as `sxs.load` and
    yields precisely the same results; essentially it is a trivial wrapper around
    the `load` function.  The benefit of this approach is that it can be used in
    precisely the same way as `h5py.File` would have been used previously.  For
    example, the old approach would be to open an HDF5 file like this:

        with h5py.File("Horizons.h5", "r") as horizons:
            areal_mass = horizons["AhA.dir/ArealMass.dat"]

    With this function, the same effect can be achieved as

        with sxs.loadcontext("Horizons.h5") as horizons:
            areal_mass = horizons["AhA.dir/ArealMass.dat"]

    Each of the datasets found in Horizons.h5 as well as the old NRAR-style files
    will be available through this interface, even when using newer files in
    different formats.  Thus, only one line of code would need to change to use the
    new interface.

    However, be aware that this may not an be efficient use of memory, and is
    almost certainly slower than the newer interfaces.  Wherever possible, you
    should update your code to use newer interfaces.  Failing to do so will leave
    you open to ridicule from your peers and loved ones.

    See Also
    --------
    load

    """
    yield load(*args, **kwargs)


def load_lvc(
        sxs_id, *,
        t_ref=None, f_ref=None,
        dt=None,
        f_low=None,
        ell_max=None,
        phi_ref=None, inclination=None,
        ell_max_epoch=None,
        **kwargs
):
    r"""Load an SXS waveform in LVC convention.

    Returns an SXS waveform (modes or polarizations) and dynamics
    (including angular velocities, frame quaternions, and spins) in
    the inertial frame that coincides with the waveform-defined frame
    defined at a reference time `t_ref` or reference frequency
    `f_ref`.

    Parameters
    ==========
    sxs_id : str
        The SXS ID of the simulation to use — e.g., "SXS:BBH:1234".
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
    =======
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

        Q = q_0 + vec(q) = (q_0, q_1, q_2, q_3)

    where q_0 is the scalar part.

    """
    lev = kwargs.pop("lev", "Lev")  # If number is unspecified, load chooses highest
    waveform_name = kwargs.pop("waveform_name", "Strain_N2.h5")
    h = load(f"{sxs_id}/{lev}/{waveform_name}", transform_to_inertial=False)
    horizons = load(f"{sxs_id}/{lev}/Horizons.h5")
    return waveforms.to_lvc_conventions(
        h, horizons,
        t_ref=t_ref, f_ref=f_ref,
        dt=dt,
        f_low=f_low,
        ell_max=ell_max,
        phi_ref=phi_ref, inclination=inclination,
        ell_max_epoch=ell_max_epoch,
        **kwargs
    )

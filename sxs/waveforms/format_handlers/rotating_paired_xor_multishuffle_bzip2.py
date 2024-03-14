"""Functions to load and save waveforms in RPXMB format"""

import sys
import warnings
import tempfile
import contextlib
import pathlib
import bz2
import json
import numpy as np
import scipy
import h5py
import quaternionic
import spherical

from .. import WaveformModes
from . import rotating_paired_diff_multishuffle_bzip2 as rpdmb
from ... import __version__
from ...utilities import default_shuffle_widths, md5checksum, xor, multishuffle, version_info


sxs_formats = ["rotating_paired_xor_multishuffle_bzip2", "rpxmb", "rpxm", "RPXMB", "RPXM"]


def save(w, file_name=None, file_write_mode="w", L2norm_fractional_tolerance=1e-10, log_frame=None, shuffle_widths=default_shuffle_widths):
    """Save a waveform in RPXMB format

    This function is primarily for backwards compatibility.  The `RPDMB` format is
    now preferred.  In fact, this function just calls `rpdmb.save`, with two
    non-default parameters, to look for `RPXMB` format strings, and to use `xor` in
    place of `diff`.  See that function's docstring for more details.

    Also note that the default `shuffle_widths` parameter has changed.  The
    previous default is now available as
    `sxs.utilities.default_shuffle_widths_old`.  The `load` function will
    automatically read the shuffle widths used from the saved file, so you do not
    need to worry about this unless you are comparing checksums or something like
    that.

    """
    return rpdmb.save(
        w, file_name=file_name, file_write_mode=file_write_mode,
        L2norm_fractional_tolerance=L2norm_fractional_tolerance, log_frame=log_frame,
        shuffle_widths=shuffle_widths, diff=xor,
        formats=sxs_formats,
    )


def load(file_name, ignore_validation=None, check_md5=True, transform_to_inertial=True, **kwargs):
    """Load a waveform in RPXMB format

    This function is primarily for backwards compatibility.  The `RPDMB` format is
    now preferred.  In fact, this function just calls `rpdmb.load`, with two
    non-default parameters, to look for `RPXMB` format strings, and to use `xor` in
    place of `diff`.  See that function's docstring for more details.

    """
    return rpdmb.load(
        file_name, ignore_validation=ignore_validation, check_md5=check_md5,
        transform_to_inertial=transform_to_inertial, diff=xor, formats=sxs_formats,
        **kwargs
    )

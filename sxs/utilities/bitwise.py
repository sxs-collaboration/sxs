"""Utilities for bitwise operations"""

import functools
import numpy as np
import numba as nb
from . import guvectorize


_float_ftylists = [(dt[:], dt[:]) for dt in [nb.float32, nb.float64, nb.complex64, nb.complex128]]
_uint_ftylists = [(dt[:], dt[:]) for dt in [nb.uint8, nb.uint16, nb.uint32, nb.uint64]]


@guvectorize(_float_ftylists, '(n)->(n)')
def diff_forward(a_in, a_out):
    """Progressively diff data along an axis

    This is primarily a helper function for `diff`.  See that function's docstring
    for details.

    """
    a_out[0] = a_in[0]
    for i in range(a_in.shape[0] - 1, 0, -1):
        a_out[i] = a_in[i - 1] - a_in[i]


@guvectorize(_float_ftylists, '(n)->(n)')
def diff_reverse(a_in, a_out):
    """Progressively reverse diff data along an axis

    This is primarily a helper function for `diff`.  See that function's docstring
    for details.

    """
    a_out[0] = a_in[0]
    for i in range(1, a_in.shape[0]):
        a_out[i] = a_out[i - 1] - a_in[i]


def diff(x, reverse=False, axis=-1, **kwargs):
    """Progressively diff data along an axis

    This function steps through an array, starting with the second element, and
    evaluates `-` on that element and the preceding one.  This is a useful
    step in compressing reasonably continuous data.

    Note that, for downstream compatibility with `xor`, if the input array dtype is
    any unsigned integer, it will be re-interpreted as a `complex128`.

    Parameters
    ----------
    x : array_like
        After conversion to an array, this must have a dtype that can be converted
        to unsigned integers with 1, 2, 4, or 8 bytes.
    reverse : bool, optional
        If True, this reverses the transformation of progressively diffing the
        data.  In particular, `diff(diff(x), reverse=True)` will be bitwise-identical
        to `x`.
    axis : int, optional
        This is the axis along which the operation will be performed.  You probably
        want this to correspond to the most quickly varying dimension of the data —
        the time axis of a timeseries, for example.
    preserve_dtype : optional
        Thrown away; here for compatibility with `xor`.

    Returns
    -------
    u : ndarray

    """
    u = np.asarray(x)
    if issubclass(u.dtype.type, np.unsignedinteger):
        u = u.view(np.complex128)
    kwargs.pop("preserve_dtype", None)
    kwargs["axis"] = axis
    if "out" in kwargs:
        kwargs["out"] = np.asarray(kwargs["out"])
    if reverse:
        u_diff = diff_reverse(u, **kwargs)
    else:
        u_diff = diff_forward(u, **kwargs)
    return u_diff


@guvectorize(_uint_ftylists, '(n)->(n)')
def xor_forward(a_in, a_out):
    """Progressively XOR data along an axis

    This is primarily a helper function for `xor`.  See that function's docstring
    for details.

    """
    a_out[0] = a_in[0]
    for i in range(a_in.shape[0] - 1, 0, -1):
        a_out[i] = np.bitwise_xor(a_in[i - 1], a_in[i])


@guvectorize(_uint_ftylists, '(n)->(n)')
def xor_reverse(a_in, a_out):
    """Progressively reverse XOR data along an axis

    This is primarily a helper function for `xor`.  See that function's docstring
    for details.

    """
    a_out[0] = a_in[0]
    for i in range(1, a_in.shape[0]):
        a_out[i] = np.bitwise_xor(a_out[i - 1], a_in[i])


del _uint_ftylists


def xor(x, reverse=False, preserve_dtype=False, axis=-1, **kwargs):
    """Progressively XOR data along an axis

    This function steps through an array, starting with the second element, and
    evaluates bitwise XOR on that element and the preceding one.  This is a useful
    step in compressing reasonably continuous data.

    Parameters
    ----------
    x : array_like
        After conversion to an array, this must have a dtype that can be converted
        to unsigned integers with 1, 2, 4, or 8 bytes.
    reverse : bool, optional
        If True, this reverses the transformation of progressively XOR-ing the
        data.  In particular, `xor(xor(x), reverse=True)` will be bitwise-identical
        to `x`.
    preserve_dtype : bool, optional
        If True, the returned array will have the same dtype as the input array.
        Note that this could result in invalid data.  For example, if the input
        data has some `float` dtype, it is possible for bitwise XOR to result in
        NaNs.  The bits are usually retained correctly nonetheless, but some
        operations may produce errors or incorrect results when NaNs are present.
        If this parameter is False, the returned array will always have unsigned
        integer type, of a bit width equal to the input array.
    axis : int, optional
        This is the axis along which the operation will be performed.  You probably
        want this to correspond to the most quickly varying dimension of the data —
        the time axis of a timeseries, for example.

    Returns
    -------
    u : ndarray
        Depending on `preserve_dtype`, this will either be an unsigned integer
        array (the default) or will have the same dtype as the input data.

    """
    x = np.asarray(x)
    itemsize = x.itemsize
    if itemsize not in [1, 2, 4, 8]:
        raise ValueError(f"Input array's byte size must be one of {{1, 2, 4, 8}}, not {itemsize}")
    dtype = np.dtype(f"u{itemsize}")
    u = x.view(dtype)
    kwargs["axis"] = axis
    if "out" in kwargs:
        kwargs["out"] = np.asarray(kwargs["out"]).view(dtype)
    if reverse:
        u_xor = xor_reverse(u, **kwargs)
    else:
        u_xor = xor_forward(u, **kwargs)
    if preserve_dtype:
        return u_xor.view(x.dtype)
    else:
        return u_xor


@functools.lru_cache()
def multishuffle(shuffle_widths, forward=True):
    """Construct functions to "multi-shuffle" data

    This is a generalization of the `shuffle` algorithm found in HDF5, or the
    `bitshuffle` algorithm — and in some sense interpolates between them.  See
    below for more explanation.

    Parameters
    ----------
    shuffle_widths : array_like
        This must be an iterable of integers, representing the number of bits in
        each piece of each number that is shuffled, starting from the highest
        significance, and proceeding to the lowest.  The sum of these numbers must
        be the total bit width of the numbers that will be given as input — which
        must currently be 8, 16, 32, or 64.  There is no restriction on the
        individual widths, but note that if they do not fit evenly into 8-bit
        bytes, the result is unlikely to compress well.
    forward : bool, optional
        If True (the default), the returned function will shuffle data; if False,
        the returned function will reverse this process — unshuffle.

    Returns
    -------
    shuffle_func : function
        This function takes just one parameter — the array to be shuffled — and
        returns the shuffled array.  Note that the input array *must* be flat (have
        just one dimension), and will be viewed as an array of unsigned integers of
        the input bit width.  This can affect the shape of the array and order of
        elements.  You should ensure that this process will result in an array of
        numbers in the order that you want.  For example, if you have a 2-d array
        of floats `a` that are more continuous along the first dimension, you might
        pass `np.ravel(a.view(np.uint64), 'F')`, where F represents Fortran order,
        which varies more quickly in the first dimension.  Also note that this
        function is JIT-compiled by numba.

    Notes
    -----
    The standard "shuffle" algorithm (as found in HDF5, for example) takes an array
    of numbers and shuffles their bytes so that all bytes of a given significance
    are stored together — the first byte of all the numbers are stored
    contiguously, then the second byte of all the numbers, and so on.  The
    motivation for this is that — with reasonably smooth data — bytes in the same
    byte position in sequential numbers are usually more related to each other than
    they are to other bytes within the same number, which means that shuffling
    results in better compression of the data.

    There is no reason that shuffling can only work byte-by-byte, however.  There
    is also a "bitshuffle" algorithm, which works in the same way, but collecting
    bits rather than bytes.  More generally, we could vary the number of bits
    stored together as we move along the numbers.  For example, we could store the
    first 8 bits of each number, followed by the next 4 bits of each number, etc.
    This is the "multi-shuffle" algorithm.

    With certain types of data, this can reduce the compressed data size
    significantly.  For example, with float data for which successive values have
    been XOR-ed, the sign bit will very rarely change, the next 11 bits
    (representing the exponent) and a few of the following bits (representing the
    highest-significance digits) will typically be highly correlated, while as we
    move to lower significance there will be less correlation.  Thus, we might
    shuffle the first 8 bits together, followed by the next 8, then the next 4, the
    next 4, the next 2, and so on — decreasing the shuffle width as we go.  The
    `shuffle_widths` input might look like [8, 8, 4, 4, 2, 2, 1, 1, 1, 1, ...].

    There are also some cases where we see correlation *increasing* again at low
    significance.  For example, if a number results from cancellation — the
    subtraction of two numbers much larger than their difference — then its
    lower-significance bits will be 0.  If we then multiply that by some integer
    (e.g., for normalization), there may be some very correlated but nonzero
    pattern.  In either case, compression might improve if the values at the end of
    our shuffle_widths list increase.

    """
    import numpy as np
    from . import jit

    bit_width = np.sum(shuffle_widths, dtype=np.int64)  # uint64 casts to float under floor division...
    if bit_width not in [8, 16, 32, 64]:
        raise ValueError(f"Total bit width must be one of [8, 16, 32, 64], not {bit_width}")
    dtype = np.dtype(f"u{bit_width//8}")
    bit_width = dtype.type(bit_width)
    reversed_shuffle_widths = np.array(list(reversed(shuffle_widths)), dtype=dtype)
    one = dtype.type(1)

    if forward:

        def shuffle(a):
            a = a.view(dtype)
            if a.ndim != 1:
                raise ValueError(
                    "\nThis function only accepts flat arrays.  Make sure you flatten "
                    "(using ravel, reshape, or flatten)\n in a way that keeps your data"
                    "contiguous in the order you want."
                )
            b = np.zeros_like(a)
            b_array_bit = np.uint64(0)
            for i, shuffle_width in enumerate(reversed_shuffle_widths):
                mask_shift = np.sum(reversed_shuffle_widths[:i])
                mask = dtype.type(2 ** shuffle_width - 1)
                pieces_per_element = bit_width // shuffle_width
                for a_array_index in range(a.size):
                    b_array_index = b_array_bit // bit_width
                    b_element_bit = dtype.type(b_array_bit % bit_width)
                    masked = (a[a_array_index] >> mask_shift) & mask
                    b[b_array_index] += masked << b_element_bit
                    if b_element_bit + shuffle_width > bit_width:
                        b[b_array_index + one] += masked >> (bit_width - b_element_bit)
                    b_array_bit += shuffle_width
            return b

        return jit(shuffle)

    else:
        # This function is almost the same as above, except for:
        # 1) swap a <-> b in input and output
        # 2) reverse the effect of the line in which b was set from a
        def unshuffle(b):
            b = b.view(dtype)
            if b.ndim != 1:
                raise ValueError(
                    "\nThis function only accepts flat arrays.  Make sure you flatten "
                    "(using ravel, reshape, or flatten)\n in a way that keeps your data"
                    "contiguous in the order you want."
                )
            a = np.zeros_like(b)
            b_array_bit = np.uint64(0)
            for i, shuffle_width in enumerate(reversed_shuffle_widths):
                mask_shift = np.sum(reversed_shuffle_widths[:i])
                mask = dtype.type(2 ** shuffle_width - 1)
                pieces_per_element = bit_width // shuffle_width
                for a_array_index in range(a.size):
                    b_array_index = b_array_bit // bit_width
                    b_element_bit = b_array_bit % bit_width
                    masked = (b[b_array_index] >> b_element_bit) & mask
                    a[a_array_index] += masked << mask_shift
                    if b_element_bit + shuffle_width > bit_width:
                        a[a_array_index] += (
                            (b[b_array_index + one] << (bit_width - b_element_bit)) & mask
                        ) << mask_shift
                    b_array_bit += shuffle_width
            return a

        return jit(unshuffle)

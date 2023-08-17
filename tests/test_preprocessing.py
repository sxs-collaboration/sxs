import pytest
import numpy as np
from sxs.waveforms.preprocessing import *


tolerance = 1e-15


def all_close(x, y):
    return np.allclose(x, y, rtol=tolerance, atol=tolerance)


# test smooth_sigmoid raising errors
def test_smooth_sigmoid_raising_errors():
    with pytest.raises(Exception, match="t_lower and t_upper must lie within the time array."):
        smooth_sigmoid(np.array([0, 1, 2, 3, 4, 5]), -1, 4)

    with pytest.raises(Exception, match="t_lower must be smaller than t_upper."):
        smooth_sigmoid(np.array([0, 1, 2, 3, 4, 5]), 3, 2)


# test smooth_sigmoid outputs
def test_smooth_sigmoid_outputs():
    result = smooth_sigmoid(np.array([0, 1, 2, 3, 4, 5, 6]), 1, 4)
    assert all_close(result[[0, 1, 4, 5, 6]], [0, 0, 1, 1, 1])

    result = smooth_sigmoid(np.array([0, 1, 2, 3, 4, 5, 6]), 1, 4, True)
    assert all_close(result[[0, 1, 4, 5, 6]], [1, 1, 0, 0, 0])

    result = smooth_sigmoid(np.array([0, 0.25, 0.5, 0.75, 1]))
    assert all_close(result[[0, 2, 4]], [0, 0.5, 1])

    assert all_close(result[1], 1 / (1 + np.exp(8 / 3)))
    assert all_close(result[3], 1 / (1 + np.exp(-8 / 3)))


# test cosine_sigmoid raising errors
def test_cosine_sigmoid_raising_errors():
    with pytest.raises(Exception, match="t_lower and t_upper must lie within the time array."):
        cosine_sigmoid(np.array([0, 1, 2, 3, 4, 5]), -1, 4)

    with pytest.raises(Exception, match="t_lower must be smaller than t_upper."):
        cosine_sigmoid(np.array([0, 1, 2, 3, 4, 5]), 3, 2)


# test cosine_sigmoid outputs
def test_cosine_sigmoid_outputs():
    result = cosine_sigmoid(np.array([0, 1, 2, 3, 4, 5, 6]), 1, 4)
    assert all_close(result[[0, 1, 4, 5, 6]], [0, 0, 1, 1, 1])

    result = cosine_sigmoid(np.array([0, 1, 2, 3, 4, 5, 6]), 1, 4, True)
    assert all_close(result[[0, 1, 4, 5, 6]], [1, 1, 0, 0, 0])

    result = cosine_sigmoid(np.array([0, 0.25, 0.5, 0.75, 1]))
    assert all_close(result[[0, 2, 4]], [0, 0.5, 1])

    assert all_close(result[1], 1 / 2 - np.sqrt(2) / 4)
    assert all_close(result[3], 1 / 2 + np.sqrt(2) / 4)


# test taper_function raising errors
def test_taper_function_raising_errors():
    with pytest.raises(Exception, match="sigmoid_type must be 'smooth' or 'cosine'."):
        taper_function(np.array([1, 2, 3]), "both", "sine")

    with pytest.raises(Exception, match="taper_type must be 'left', 'right', or 'both'."):
        taper_function(np.array([1, 2, 3]), "none", "cosine")

    # left tapering
    with pytest.raises(Exception, match="Left tapering is not specified correctly."):
        taper_function(np.array([1, 2, 3]), "left", "smooth", t0_taper=1, t1_taper=3, left_taper_length=2)
        taper_function(np.array([1, 2, 3]), "left", "smooth", t0_taper=1)

    with pytest.raises(
        Exception, match="t0_taper must be smaller than t1_taper, and they both must lie within the time array."
    ):
        taper_function(np.array([1, 2, 3]), "left", "smooth", t0_taper=2, t1_taper=2)
        taper_function(np.array([1, 2, 3]), "left", "smooth", t0_taper=0, t1_taper=2)

    with pytest.raises(Exception, match="left_taper_length must be shorter than the time array."):
        taper_function(np.array([1, 2, 3]), "left", "smooth", left_taper_length=3)

    # right tapering
    with pytest.raises(Exception, match="Right tapering is not specified correctly."):
        taper_function(np.array([1, 2, 3]), "right", "smooth", t0_taper=1, t1_taper=3)
        taper_function(np.array([1, 2, 3]), "right", "smooth", t2_taper=1, t3_taper=3, right_taper_length=2)
        taper_function(np.array([1, 2, 3]), "right", "smooth", t2_taper=1)

    with pytest.raises(
        Exception, match="t2_taper must be smaller than t3_taper, and they both must lie within the time array."
    ):
        taper_function(np.array([1, 2, 3]), "right", "smooth", t2_taper=2, t3_taper=2)
        taper_function(np.array([1, 2, 3]), "right", "smooth", t2_taper=0, t3_taper=2)

    with pytest.raises(Exception, match="right_taper_length must be shorter than the time array."):
        taper_function(np.array([1, 2, 3]), "right", "smooth", right_taper_length=3)

    # windowing
    with pytest.raises(Exception, match="Right tapering is not specified correctly."):
        taper_function(np.array([1, 2, 3]), "both", "smooth", t0_taper=1, t1_taper=3)


# test taper_function outputs
def test_taper_function_outputs():
    result = taper_function(np.array([0, 1, 2, 3, 4, 5, 6]), "left", "smooth", t0_taper=1, t1_taper=4)
    assert all_close(result[[0, 1, 4, 5, 6]], [0, 0, 1, 1, 1])

    result = taper_function(np.array([0, 1, 2, 3, 4, 5, 6]), "right", "cosine", t2_taper=2, t3_taper=5)
    assert all_close(result[[0, 1, 2, 5, 6]], [1, 1, 1, 0, 0])

    result = taper_function(np.array([0, 1, 2, 3, 4, 5, 6]), "left", "cosine", left_taper_length=3)
    assert all_close(result[[0, 3, 4, 5, 6]], [0, 1, 1, 1, 1])
    assert 0 < result[1] < result[2] < 1

    result = taper_function(np.array([0, 1, 2, 3, 4, 5, 6]), "right", "smooth", right_taper_length=4)
    assert all_close(result[[0, 1, 2, 6]], [1, 1, 1, 0])
    assert 0 < result[5] < result[4] < result[3] < 1

    result = taper_function(
        np.array([0, 1, 2, 3, 4, 5, 6]), "both", "smooth", left_taper_length=2, right_taper_length=1
    )
    assert all_close(result[[0, 2, 3, 4, 5, 6]], [0, 1, 1, 1, 1, 0])

    result = taper_function(np.arange(50), "both", "smooth", t0_taper=10, t1_taper=16, t2_taper=30, t3_taper=43)
    assert all_close(result[:11], np.zeros(11))
    assert all_close(result[43:], np.zeros(7))
    assert all_close(result[16:31], np.ones(15))


# test transition_tail_to_constant raising errors
def test_transition_tail_to_constant_raising_errors():
    with pytest.raises(Exception, match="transition_length must be shorter than or equal to the length of the data."):
        transition_tail_to_constant(np.arange(5), np.ones(5), 5)
        transition_tail_to_constant(np.arange(5), np.ones((2, 5)), 10)

    with pytest.raises(Exception, match="The data array must be one or two dimensional."):
        transition_tail_to_constant(np.arange(5), np.ones((2, 3, 5)), 4)

    # should pass
    transition_tail_to_constant(np.arange(100), np.arange(100), 20)


# test pad raising errors
def test_pad_raising_errors():
    with pytest.raises(Exception, match="pad_type must be either 'edge' or 'zero'."):
        pad(np.arange(5), np.ones((2, 5)), "constant", 1, 2)

    with pytest.raises(Exception, match="The data array must be one or two dimensional."):
        pad(np.arange(5), np.ones((2, 3, 5)), "edge", 1, 2)


# test pad outputs
def test_pad_outputs():
    time, data = pad(np.array([0, 1, 2]), np.array([3, 4, 5]), "edge", 1, 2)
    assert np.allclose(time, np.array([-1, 0, 1, 2, 3, 4]))
    assert np.allclose(data, np.array([3, 3, 4, 5, 5, 5]))

    time, data = pad(np.array([0, 1, 2]), np.array([[3, 4, 5], [7, 9, 11]]), "edge", 1, 2)
    assert np.allclose(time, np.array([-1, 0, 1, 2, 3, 4]))
    assert np.allclose(data, np.array([[3, 3, 4, 5, 5, 5], [7, 7, 9, 11, 11, 11]]))

    time, data = pad(np.array([0, 1, 2]), np.array([3, 4, 5]), "zero", 1, 2)
    assert np.allclose(time, np.array([-1, 0, 1, 2, 3, 4]))
    assert np.allclose(data, np.array([0, 3, 4, 5, 0, 0]))

    time, data = pad(np.array([0, 1, 2]), np.array([[3, 4, 5], [7, 9, 11]]), "zero", 1, 2)
    assert np.allclose(time, np.array([-1, 0, 1, 2, 3, 4]))
    assert np.allclose(data, np.array([[0, 3, 4, 5, 0, 0], [0, 7, 9, 11, 0, 0]]))


# test taper outputs
def test_taper_outputs():
    result = taper(np.array([0, 1, 2, 3, 4, 5, 6]), 3 * np.ones(7), "left", "smooth", t0_taper=1, t1_taper=4)
    assert np.allclose(result[[0, 1, 4, 5, 6]], [0, 0, 3, 3, 3])

    result = taper(np.array([0, 1, 2, 3, 4, 5, 6]), 3 * np.ones(7), "right", "cosine", t2_taper=2, t3_taper=5)
    assert np.allclose(result[[0, 1, 2, 5, 6]], [3, 3, 3, 0, 0])

    result = taper(np.array([0, 1, 2, 3, 4, 5, 6]), 3 * np.ones(7), "left", "cosine", left_taper_length=3)
    assert np.allclose(result[[0, 3, 4, 5, 6]], [0, 3, 3, 3, 3])
    assert 0 < result[1] < result[2] < 3

    result = taper(np.array([0, 1, 2, 3, 4, 5, 6]), 3 * np.ones(7), "right", "smooth", right_taper_length=4)
    assert np.allclose(result[[0, 1, 2, 6]], [3, 3, 3, 0])
    assert 0 < result[5] < result[4] < result[3] < 3

    result = taper(
        np.array([0, 1, 2, 3, 4, 5, 6]), 3 * np.ones(7), "both", "smooth", left_taper_length=2, right_taper_length=1
    )
    assert np.allclose(result[[0, 2, 3, 4, 5, 6]], [0, 3, 3, 3, 3, 0])

    result = taper(
        np.array([0, 1, 2, 3, 4, 5, 6]),
        np.arange(14).reshape(2, 7),
        "both",
        "smooth",
        left_taper_length=2,
        right_taper_length=1,
    )
    assert np.allclose(result[:, [0, 2, 3, 4, 5, 6]], np.array([[0, 2, 3, 4, 5, 0], [0, 9, 10, 11, 12, 0]]))

    result = taper(np.arange(50), 3 * np.ones(50), "both", "smooth", t0_taper=10, t1_taper=16, t2_taper=30, t3_taper=43)
    assert np.allclose(result[:11], np.zeros(11))
    assert np.allclose(result[43:], np.zeros(7))
    assert np.allclose(result[16:31], 3 * np.ones(15))


# test detrend raising errors
def test_detrend_raising_errors():
    with pytest.raises(Exception, match="detrend_type must be 'none', 'zero right end', or 'zero both ends'."):
        detrend(np.arange(5), np.ones((2, 5)), "None")

    with pytest.raises(Exception, match="The data array must be one or two dimensional."):
        detrend(np.arange(5), np.ones((2, 3, 5)), "zero right end")


# test detrend outputs
def test_detrend_outputs():
    result = detrend(np.arange(5), np.arange(2, 7), "none")
    assert np.allclose(result, np.arange(2, 7))

    result = detrend(np.arange(5), np.arange(6, 11), "zero right end")
    assert np.allclose(result, -np.arange(5) + 6)

    result = detrend(np.arange(5), np.arange(6, 11), "zero both ends")
    assert np.allclose(result, np.arange(5) / 5)

    result = detrend(np.arange(5), np.ones((2, 5)), "zero right end")
    assert np.allclose(result[0], -np.arange(5) / 5 + 1)
    assert np.allclose(result[1], -np.arange(5) / 5 + 1)

    result = detrend(np.arange(5), np.ones((3, 5)), "zero both ends")
    assert np.allclose(result, np.zeros((3, 5)))


# test preprocess_general raising errors
def test_preprocess_general_raising_errors():
    with pytest.raises(Exception, match="The data array must be one or two dimensional."):
        preprocess_general(
            [0, 1, 2],
            [[[0, 1, 2]]],
            taper_type="both",
            detrend_type="none",
            right_pad_length=2,
            left_taper_length=1,
            right_taper_length=1,
        )

    with pytest.raises(Exception, match="The data array's last dimension must have the same size as the time array."):
        preprocess_general(
            [0, 1, 2],
            [[0, 1, 2, 3]],
            taper_type="both",
            detrend_type="none",
            right_pad_length=2,
            left_taper_length=1,
            right_taper_length=1,
        )

    # should pass
    preprocess_general(
        np.arange(100),
        np.arange(100),
        transition_length=20,
        taper_type="both",
        detrend_type="none",
        right_pad_length=20,
        left_taper_length=10,
        right_taper_length=20,
    )

    # should pass
    preprocess_general(
        np.arange(100),
        np.arange(100),
        transition_length=0,
        taper_type="both",
        detrend_type="none",
        right_pad_length=10,
        left_taper_length=30,
        right_taper_length=40,
    )

    # should pass
    preprocess_general(
        np.arange(100),
        np.arange(200).reshape(2, 100),
        transition_length=10,
        taper_type="left",
        detrend_type="zero right end",
        right_pad_length=15,
        left_taper_length=25,
        right_taper_length=5,
    )


# test preprocess raising errors
def test_preprocess_raising_errors():
    with pytest.raises(Exception, match="preprocess_type must be 'taper both sides' or 'taper left side and detrend'."):
        preprocess(
            [0, 1, 2],
            [[[0, 1, 2]]],
            preprocess_type="taper left side -> detrend",
            transition_length=0,
            pad_length=2,
            left_taper_length=1,
            right_taper_length=1,
        )

    # should pass
    preprocess(
        np.arange(100),
        np.arange(100),
        preprocess_type="taper left side and detrend",
        transition_length=20,
        pad_length=20,
        left_taper_length=10,
        right_taper_length=20,
    )

    # should pass
    preprocess(
        np.arange(100),
        np.arange(300).reshape(3, 100),
        preprocess_type="taper both sides",
        transition_length=20,
        pad_length=20,
        left_taper_length=10,
        right_taper_length=20,
    )

import numpy as np
import quaternionic
import pytest
import sxs
from .conftest import physicalish_waveform

def test_independent_alignment():
    h1 = physicalish_waveform()
    
    δt = 10 * np.random.rand()
    δSpin3 = quaternionic.array([np.random.rand(), np.random.rand(), np.random.rand(), np.random.rand()]).normalized
    
    h2 = h1.copy()
    h2.t = h2.t - δt
    h2 = h2.rotate(δSpin3)

    t1 = -80
    t2 = 80
    h1_prime, transformation, norm = sxs.waveforms.alignment.align_waveforms(
        h1,
        h2,
        t1,
        t2,
        alignment_method="independent alignment",
        t_ref=-20
    )
    
    assert np.allclose(norm, 0., atol=1e-6)
    assert np.allclose([δt, *δSpin3.ndarray], transformation, atol=1e-6)
    
def test_align1d():
    h1 = physicalish_waveform()
    
    δt = 10 * np.random.rand()
    
    h2 = h1.copy()
    h2.t = h2.t - δt
    
    t1 = -80
    t2 = 80
    h1_prime, transformation, norm = sxs.waveforms.alignment.align_waveforms(
        h1,
        h2,
        t1,
        t2,
        alignment_method="1d"
    )
    
    assert np.allclose(norm, 0., atol=1e-6)
    assert np.allclose(δt, transformation[0], atol=1e-6)
    
def test_align2d():
    h1 = physicalish_waveform()
    
    δt = 10 * np.random.rand()
    δSpin3 = quaternionic.array([np.random.rand(), 0., 0., np.random.rand()]).normalized
    
    h2 = h1.copy()
    h2.t = h2.t - δt
    h2 = h2.rotate(δSpin3)

    t1 = -80
    t2 = 80
    h1_prime, transformation, norm = sxs.waveforms.alignment.align_waveforms(
        h1,
        h2,
        t1,
        t2,
        alignment_method="2d",
        n_brute_force_δϕ=25
    )
    
    assert np.allclose(norm, 0., atol=1e-6)
    assert np.allclose([δt, *δSpin3.ndarray], transformation, atol=1e-6)
    
def test_align4d():
    h1 = physicalish_waveform()
    
    δt = 10 * np.random.rand()
    δSpin3 = quaternionic.array([np.random.rand(), np.random.rand(), np.random.rand(), np.random.rand()]).normalized
    
    h2 = h1.copy()
    h2.t = h2.t - δt
    h2 = h2.rotate(δSpin3)

    t1 = -80
    t2 = 80
    h1_prime, transformation, norm = sxs.waveforms.alignment.align_waveforms(
        h1,
        h2,
        t1,
        t2,
        alignment_method="4d",
        n_brute_force_δϕ=25
    )
    
    assert np.allclose(norm, 0., atol=1e-6)
    assert np.allclose([δt, *δSpin3.ndarray], transformation, atol=1e-6)

import numpy as np
import sxs
import pytest


def test_uprime_generator():
    β = 0.1
    for t1, t2 in [[-30, -10], [-10, 10], [10, 30]]:
        u = np.linspace(t1, t2, num=1000)
        uprime = sxs.waveforms.transformations.uprime_generator(u, β)
        assert uprime[0] > u[0] and uprime[-1] < u[-1]
    for t1, t2 in [[-20, 0]]:
        u = np.linspace(t1, t2, num=1000)
        uprime = sxs.waveforms.transformations.uprime_generator(u, β)
        assert uprime[0] > u[0] and uprime[-1] == u[-1]
    for t1, t2 in [[0, 20]]:
        u = np.linspace(t1, t2, num=1000)
        uprime = sxs.waveforms.transformations.uprime_generator(u, β)
        assert uprime[0] == u[0] and uprime[-1] < u[-1]

    β = 0.9
    for t1, t2 in [[-10, 10]]:
        u = np.linspace(t1, t2, num=1000)
        uprime = sxs.waveforms.transformations.uprime_generator(u, β)
        assert uprime[0] >= u[0] and uprime[-1] <= u[-1]
    for t1, t2 in [[-20, 0]]:
        u = np.linspace(t1, t2, num=1000)
        uprime = sxs.waveforms.transformations.uprime_generator(u, β)
        assert uprime[0] > u[0] and uprime[-1] == u[-1]
    for t1, t2 in [[0, 20]]:
        u = np.linspace(t1, t2, num=1000)
        uprime = sxs.waveforms.transformations.uprime_generator(u, β)
        assert uprime[0] == u[0] and uprime[-1] < u[-1]
    for t1, t2 in [[-10, -9], [9, 10]]:
        u = np.linspace(t1, t2, num=1000)
        with pytest.raises(ValueError):
            uprime = sxs.waveforms.transformations.uprime_generator(u, β)

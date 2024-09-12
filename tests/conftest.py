import os
import functools
import pytest
import numpy as np
import quaternionic
import spherical
import sxs

ell_max_default = 36

try:
    import spinsfast
    requires_spinsfast = lambda f: f
except:
    requires_spinsfast = pytest.mark.skip(reason="spinsfast is missing")

try:
    import sympy
    requires_sympy = lambda f: f
except:
    requires_sympy = pytest.mark.skip(reason="sympy is missing")

GITHUB_ACTIONS_MACOS = (
    os.getenv("GITHUB_ACTIONS", "") == "true"
    and os.getenv("RUNNER_OS", "") == "macOS"
)

if GITHUB_ACTIONS_MACOS:
    skip_macOS_GH_actions_downloads = pytest.mark.skip(
        reason="macOS runners on GitHub Actions have connectivity problems"
    )
else:
    skip_macOS_GH_actions_downloads = lambda f: f


def pytest_addoption(parser):
    parser.addoption("--ell_max", action="store", type=int, default=ell_max_default,
                     help="Maximum ell value to test")
    parser.addoption("--ell_max_slow", action="store", type=int, default=ell_max_default // 2,
                     help="Maximum ell value to test with slow tests")
    parser.addoption("--run_slow_tests", action="store_true", default=False,
                     help="Run all tests, including slow ones")


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: marks tests as slow")


def pytest_collection_modifyitems(config, items):
    if config.getoption("--run_slow_tests"):
        return
    skip_slow = pytest.mark.skip(reason="need --run_slow_tests option to run")
    for item in items:
        if "slow" in item.keywords:
            item.add_marker(skip_slow)


def pytest_runtest_setup(item):
    if 'slow' in item.keywords and not item.config.getoption("--run_slow_tests"):
        pytest.skip("Need `--run_slow_tests` command-line argument to run")


@pytest.fixture
def ell_max(request):
    return request.config.getoption("--ell_max")


@pytest.fixture
def ell_max_slow(request):
    return request.config.getoption("--ell_max_slow")


@pytest.fixture
def special_angles():
    return np.arange(-1 * np.pi, 1 * np.pi + 0.1, np.pi / 4.)


@pytest.fixture
def on_windows():
    from sys import platform
    return 'win' in platform.lower() and not 'darwin' in platform.lower()


@pytest.fixture
def eps():
    return np.finfo(float).eps


def quaternion_sampler():
    Qs_array = quaternionic.array([
        [np.nan, 0., 0., 0.],
        [np.inf, 0., 0., 0.],
        [-np.inf, 0., 0., 0.],
        [0., 0., 0., 0.],
        [1., 0., 0., 0.],
        [0., 1., 0., 0.],
        [0., 0., 1., 0.],
        [0., 0., 0., 1.],
        [1.1, 2.2, 3.3, 4.4],
        [-1.1, -2.2, -3.3, -4.4],
        [1.1, -2.2, -3.3, -4.4],
        [
            0.18257418583505537115232326093360,
            0.36514837167011074230464652186720,
            0.54772255750516611345696978280080,
            0.73029674334022148460929304373440
        ],
        [1.7959088706354, 0.515190292664085, 0.772785438996128, 1.03038058532817],
        [2.81211398529184, -0.392521193481878, -0.588781790222817, -0.785042386963756],
    ])
    names = type("QNames", (object,), dict())()
    names.q_nan1 = 0
    names.q_inf1 = 1
    names.q_minf1 = 2
    names.q_0 = 3
    names.q_1 = 4
    names.x = 5
    names.y = 6
    names.z = 7
    names.Q = 8
    names.Qneg = 9
    names.Qbar = 10
    names.Qnormalized = 11
    names.Qlog = 12
    names.Qexp = 13
    return Qs_array, names


@pytest.fixture
def Qs():
    return quaternion_sampler()[0]


@pytest.fixture
def Q_names():
    return quaternion_sampler()[1]


@pytest.fixture
def Q_conditions():
    Qs_array, names = quaternion_sampler()
    conditions = type("QConditions", (object,), dict())()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        conditions.zero = np.arange(len(Qs_array))[Qs_array == Qs_array[names.q_0]]
        conditions.nonzero = np.arange(len(Qs_array))[np.nonzero(Qs_array)]
        conditions.nan = np.arange(len(Qs_array))[np.isnan(Qs_array)]
        conditions.nonnan = np.arange(len(Qs_array))[~np.isnan(Qs_array)]
        conditions.nonnannonzero = np.arange(len(Qs_array))[~np.isnan(Qs_array) & (Qs_array != Qs_array[names.q_0])]
        conditions.inf = np.arange(len(Qs_array))[np.isinf(Qs_array)]
        conditions.noninf = np.arange(len(Qs_array))[~np.isinf(Qs_array)]
        conditions.noninfnonzero = np.arange(len(Qs_array))[~np.isinf(Qs_array) & (Qs_array != Qs_array[names.q_0])]
        conditions.finite = np.arange(len(Qs_array))[np.isfinite(Qs_array)]
        conditions.nonfinite = np.arange(len(Qs_array))[~np.isfinite(Qs_array)]
        conditions.finitenonzero = np.arange(len(Qs_array))[np.isfinite(Qs_array) & (Qs_array != Qs_array[names.q_0])]
    return conditions


@pytest.fixture
def Rs():
    ones = [0, -1., 1.]
    rs = [[w, x, y, z] for w in ones for x in ones for y in ones for z in ones][1:]
    np.random.seed(hash("Rs") % 4294967294)  # Use mod to get in an acceptable range
    rs = rs + [r for r in [quaternionic.array(np.random.uniform(-1, 1, size=4)) for _ in range(20)]]
    return quaternionic.array(rs).normalized


# catalog = sxs.load("catalog")
# com_files = [
#     f for f in catalog.files.values()
#     if "rhOverM_Asymptotic_GeometricUnits_CoM.h5" in f['filename']
#     and "SXS:BBH:1111v" not in f['truepath']
# ]
# min(com_files, key=lambda f: f['filesize'])
shortest_h_com_file = "SXS:BBH:0156v1/Lev5/rhOverM_Asymptotic_GeometricUnits_CoM.h5"
shortest_horizons = "SXS:BBH:0156v1/Lev5/Horizons.h5"
shortest_metadata = "SXS:BBH:0156v5/Lev5/metadata.json"
shortest_metadata_txt = "SXS:BBH:0156v5/Lev5/metadata.txt"


@functools.lru_cache()
def get_h():
    import contextlib
    import sxs
    with contextlib.redirect_stdout(None):
        h = sxs.load(shortest_h_com_file, extrapolation_order=3)
    return h


@pytest.fixture
def h():
    return get_h().copy()


def constant_waveform(begin=-10.0, end=100.0, n_times=1000, ell_min=2, ell_max=8):
    import sxs
    t = np.linspace(begin, end, num=n_times)
    lm = np.array([[ell, m] for ell in range(ell_min, ell_max + 1) for m in range(-ell, ell + 1)])
    data = np.empty((t.shape[0], lm.shape[0]), dtype=complex)
    for i, m in enumerate(lm[:, 1]):
        data[:, i] = m - 1j * m
    W = sxs.WaveformModes(
        data,
        time=t,
        modes_axis=1,
        ell_min=ell_min,
        ell_max=ell_max,
        frame_type="corotating",
        data_type="h",
        spin_weight=-2,
    )
    return W


@pytest.fixture(name="constant_waveform")
def constant_waveform_fixture():
    return constant_waveform()


def linear_waveform(begin=-10.0, end=100.0, n_times=1000, ell_min=2, ell_max=8):
    import sxs
    np.random.seed(hash("linear_waveform") % 4294967294)  # Use mod to get in an acceptable range
    axis = quaternionic.array.from_vector_part(np.random.uniform(-1, 1, size=3)).normalized
    t = np.linspace(begin, end, num=n_times)
    omega = 2 * np.pi * 4 / (t[-1] - t[0])
    frame = np.array([np.exp(axis * (omega * t_i / 2)) for t_i in t])
    lm = np.array([[ell, m] for ell in range(ell_min, ell_max + 1) for m in range(-ell, ell + 1)])
    data = np.empty((t.shape[0], lm.shape[0]), dtype=complex)
    for i, m in enumerate(lm[:, 1]):
        # N.B.: This form is used in test_linear_interpolation; if you
        # change it here, you must change it there.
        data[:, i] = (m - 1j * m) * t
    W = sxs.WaveformModes(
        data,
        time=t,
        modes_axis=1,
        ell_min=ell_min,
        ell_max=ell_max,
        frame_type="corotating",
        data_type="h",
        spin_weight=-2,
    )
    return W


@pytest.fixture(name="linear_waveform")
def linear_waveform_fixture():
    return linear_waveform()


def random_waveform(begin=-10.0, end=100.0, n_times=1000, ell_min=None, ell_max=8):
    np.random.seed(hash("random_waveform") % 4294967294)  # Use mod to get in an acceptable range
    data_type = "h"
    spin_weight = -2
    if ell_min is None:
        ell_min = abs(spin_weight)
    n_modes = ell_max * (ell_max + 2) - ell_min ** 2 + 1
    t = np.sort(np.random.uniform(begin, end, size=n_times))
    frame = quaternionic.array.random(len(t), normalize=True)
    data = np.random.normal(size=(n_times, n_modes, 2)).view(complex)[:, :, 0]
    data[:, :spherical.LM_index(abs(spin_weight), -abs(spin_weight), ell_min)] = 0.0
    lm = np.array([[ell, m] for ell in range(ell_min, ell_max + 1) for m in range(-ell, ell + 1)])
    W = sxs.WaveformModes(
        data,
        time=t,
        frame=frame,
        modes_axis=1,
        ell_min=ell_min,
        ell_max=ell_max,
        frame_type="corotating",
        data_type=data_type,
        spin_weight=spin_weight,
    )
    return W


@pytest.fixture(name="random_waveform")
def random_waveform_fixture():
    return random_waveform()


def delta_waveform(ell, m, begin=-10.0, end=100.0, n_times=1000, ell_min=2, ell_max=8):
    """WaveformModes with 1 in selected slot and 0 elsewhere"""
    n_modes = ell_max * (ell_max + 2) - ell_min ** 2 + 1
    t = np.linspace(begin, end, num=n_times)
    data = np.zeros((n_times, n_modes), dtype=complex)
    data[:, spherical.LM_index(ell, m, ell_min)] = 1.0 + 0.0j
    lm = np.array([[ell, m] for ell in range(ell_min, ell_max + 1) for m in range(-ell, ell + 1)])
    W = sxs.WaveformModes(
        data,
        time=t,
        modes_axis=1,
        ell_min=min(lm[:, 0]),
        ell_max=max(lm[:, 0]),
        frame_type="inertial",
        data_type="psi4",
        spin_weight=-2,
    )
    return W


@pytest.fixture(name="delta_waveform")
def delta_waveform_fixture():
    return delta_waveform()

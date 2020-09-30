"""Utility classes for lvcnr submodule"""

class Waveform(object):
    """Contains a waveform"""
    def __init__(self, *args, **kwargs):
        if args or kwargs:
            raise ValueError("This is an empty constructor; use `from_data` or `read` to add data.")

    @classmethod
    def read(cls, file_path, group_path='/', modes=None, ell_min=2, ell_max=8):
        import numpy as np
        import h5py
        w = cls()
        modes_string = "Y_l{0[0]}_m{0[1]}.dat"
        modes = modes or [[ell, m] for ell in range(ell_min, ell_max+1) for m in range(-ell, ell+1)]
        with h5py.File(file_path, 'r') as f:
            h = f[group_path]
            w.t = h[modes_string.format(modes[0])][:, 0]
            w.t.setflags(write=False)
            w.data = np.empty((w.t.size, len(modes)*2))
            for i, mode in enumerate(modes):
                w.data[:, 2*i:2*i+2] = h[modes_string.format(mode)][:, 1:3]
            w.data = w.data.view(complex)
            w.data.setflags(write=False)

            w.version_hist = f.get('VersionHist.ver', None)
            if w.version_hist is not None:
                w.version_hist = w.version_hist[:]
        return w

    def peak_index(self, start_index=0):
        import numpy as np
        return start_index + np.argmax(np.linalg.norm(self.data[start_index:], axis=1))


class WaveformAmpPhase(Waveform):
    @classmethod
    def read(cls, file_path, group_path='/', mode_strings=None, ell_min=2, ell_max=8, start_time=None):
        import numpy as np
        import h5py
        w = cls()
        mode_strings = mode_strings or ["Y_l{0}_m{1}.dat".format(ell, m) for ell in range(ell_min, ell_max+1) for m in range(-ell, ell+1)]
        with h5py.File(file_path, 'r') as f:
            h = f[group_path]
            t = h[mode_strings[0]][:, 0]
            if start_time is None:
                start_index = 0
            else:
                start_index = np.argmin(np.abs(t - start_time)) + 1
            w.t = t[start_index:].copy()
            w.t.setflags(write=False)
            w.data = np.empty((w.t.size, len(mode_strings)*2))
            d = np.empty((w.t.size, 2), dtype=float)
            c = d.view(complex).reshape(-1)
            for i, mode_string in enumerate(mode_strings):
                d[:, :] = h[mode_string][start_index:, 1:3]
                w.data[:, 2*i] = np.abs(c)
                w.data[:, 2*i+1] = np.unwrap(np.angle(c))
            w.data.setflags(write=False)
            w.version_hist = f.get('VersionHist.ver', None)
            if w.version_hist is not None:
                w.version_hist = w.version_hist[:]
        return w

    def peak_index(self, start_index=0):
        import numpy as np
        return start_index + np.argmax(np.linalg.norm(self.data[start_index:, ::2], axis=1))

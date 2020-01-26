from . import lvc



class Waveform(object):
    def __init__(self, file_path, group_path='/', modes=None, ell_min=2, ell_max=8):
        import numpy as np
        import h5py
        modes_string = "Y_l{0[0]}_m{0[1]}.dat"
        modes = modes or [[ell, m] for ell in range(ell_min, ell_max+1) for m in range(-ell, ell+1)]
        with h5py.File(file_path, 'r') as f:
            h = f[group_path]
            self.t = h[modes_string.format(modes[0])][:, 0]
            self.data = np.empty((self.t.size, len(modes)*2))
            for i, mode in enumerate(modes):
                self.data[:, 2*i:2*i+2] = h[modes_string.format(mode)][:, 1:3]
            self.data = self.data.view(complex)
            self.t.setflags(write=False)
            self.data.setflags(write=False)

            self.version_hist = f.get('VersionHist.ver', None)
            if self.version_hist is not None:
                self.version_hist = self.version_hist[:]

    def peak_index(self, start_index=0):
        import numpy as np
        return start_index + np.argmax(np.linalg.norm(self.data[start_index:], axis=1))

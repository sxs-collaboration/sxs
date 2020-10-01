"""Class to contain a single LVC-NR-format dataset"""

class Dataset(object):
    """Represents a single dataset of a group in an LVC-format HDF5 file

    The lvcnrpy module requires all fields in format 3 to have subfields 'X', 'Y',
    'deg', 'tol', and 'errors'.  The first four of these make sense.  The last is
    not so relevant for LVC purposes; it is inherited from romspline, and is
    actually the L1 convergence measure of romspline.  In particular, it is nearly
    -- but distinctly not -- the same size as 'X' and 'Y'.  This is naturally
    useful for investigating the algorithm itself, but is irrelevant to using its
    results.  Moreover, since this class uses the "peak-greed" variant of the
    algorithm, that dataset no longer means the same thing.  But since it is
    required, we include here an `errors` member, containing a single element,
    being the largest error (scaled, if relevant) in the final result.

    """
    def __init__(self, *args, **kwargs):
        if args or kwargs:
            raise ValueError("This is an empty constructor; use `from_data` or `read` to add data.")

    @classmethod
    def from_data(cls, x, y, tol, rel=False, error_scaling=None, truncation_tol=None):
        """Construct reduced-order dataset from (x, y) data"""
        import numpy as np
        from ... import TimeSeries
        from ..decimation.peak_greed import minimal_grid
        lvc_dataset = cls()
        lvc_dataset.tol = tol
        y_reference = y
        if truncation_tol is not None:
            y = TimeSeries(y.copy(), x)
            y.truncate(truncation_tol)
            y = y.ndarray
        if error_scaling is None:
            error_scaling = 1.0
        indices = minimal_grid(x, y, tol=tol, error_scale=error_scaling, y_reference=y_reference)
        lvc_dataset.deg = 3
        lvc_dataset.X = x[indices].copy()
        lvc_dataset.Y = y[indices].copy()
        if error_scaling is None:
            lvc_dataset.errors = np.array([np.max(np.abs(y_reference - lvc_dataset.spline(x)))])
        else:
            lvc_dataset.errors = np.array([np.max(np.abs(error_scaling * (y_reference - lvc_dataset.spline(x))))])
        # lvc_dataset.compression_ratio = x.size/lvc_dataset.X.size
        # print("Size ratio:", x.size/lvc_dataset.X.size)
        return lvc_dataset

    def write(self, output_group):
        import h5py
        if not isinstance(output_group, h5py.Group):
            raise Exception("Parameter `output_group` must be an h5py.Group (or File) object.")
        output_group.create_dataset('deg', data=self.deg, dtype='int')
        output_group.create_dataset('tol', data=self.tol, dtype='double')
        output_group.create_dataset('errors', data=self.errors, dtype='double')
        output_group.create_dataset(
            'X', data=self.X, dtype='double', compression='gzip', shuffle=True, chunks=(self.X.size,)
        )
        output_group.create_dataset(
            'Y', data=self.Y, dtype='double', compression='gzip', shuffle=True, chunks=(self.Y.size,)
        )

    @classmethod
    def read(cls, input_group):
        import h5py
        if not isinstance(input_group, h5py.Group):
            raise Exception("Parameter `input_group` must be an h5py.Group (or File) object.")
        lvc_dataset = Dataset()
        lvc_dataset.deg = input_group['deg'][()]
        lvc_dataset.tol = input_group['tol'][()]
        lvc_dataset.errors = input_group['errors'][:]
        lvc_dataset.X = input_group['X'][:]
        lvc_dataset.Y = input_group['Y'][:]
        return lvc_dataset

    def spline(self, xprime):
        from scipy.interpolate import InterpolatedUnivariateSpline as spline
        return spline(self.X, self.Y, k=self.deg)(xprime)

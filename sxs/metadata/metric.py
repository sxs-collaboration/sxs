from ..utilities.string_converters import *
import numpy as np

class MetadataMetric:
    """A metric for comparing metadata.

    This class is designed to be used as a callable object that takes
    two collections of metadata (`sxs.Metadata`, `dict`, `pd.Series`)
    and returns a number measuring the distance between the metadata.
    
    With the default arguments, this will not strictly be a metric, as
    it does not satisfy the triangle inequality.  However, it is
    intended to be used as a heuristic for sorting and filtering
    metadata, rather than as a strict metric for clustering or
    classification.

    Note that calling an object of this class with two metadata
    collections will return the *squared* distance between them.

    Parameters
    ----------
    parameters : list of str, optional
        The names of the metadata fields to be compared.  The defaults
        are the reference quantities for mass ratio, spin,
        eccentricity, and mean anomaly.  Note that all of these fields
        *must* be present in *both* metadata collections.  (The
        `Metadata.add_standard_parameters` method may be useful here.)
    metric : array_like, optional
        The matrix used to weight the differences in the parameters.
        The default is a diagonal matrix with ones on the diagonal,
        except for the mean-anomaly entry, which is 1/pi^2.
    allow_different_object_types : bool, optional
        If True, metadata with different object types (BHBH, BHNS,
        NSNS) will be compared without penalty.  If False, metadata
        with different object types will be assigned an infinite
        distance.
    eccentricity_threshold1 : float, optional
        The threshold eccentricity below which we consider metadata1
        non-eccentric.  Default is 1e-2.
    eccentricity_threshold2 : float, optional
        The threshold eccentricity below which we consider metadata2
        non-eccentric.  Default is 1e-3.
    eccentricity_threshold_penalize_shorter : int, optional
        The number of orbits below which we penalize metadata2 for
        having a non-zero eccentricity when metadata1 does not.  This
        is intended to avoid ascribing small distances to systems with
        shorter inspirals.  Default is 20.

    The mean anomaly, if present, is treated specially to account for
    the fact that a mean anomaly of 0 is equivalent to a mean anomaly
    of 2π.  The difference between the entries in the two metadata
    collections is "unwrapped" before the metric is applied.

    If the eccentricity of metadata1 is below
    `eccentricity_threshold1`, then the mean anomaly is ignored.  If
    that is true and the eccentricity of metadata2 is below
    `eccentricity_threshold2` *and* the number of orbits in metadata2
    is longer than `eccentricity_threshold_penalize_shorter`, then the
    eccentricity is also ignored.  You may set these arguments to 0 to
    disable these features.

    """
    def __init__(
            self,
            parameters=[
                "reference_mass1",
                "reference_mass2",
                "reference_dimensionless_spin1",
                "reference_dimensionless_spin2",
                "reference_eccentricity",
                "reference_mean_anomaly",
            ],
            metric=np.diag([1, 1, 1, 1, 1, 1, 1, 1, 1, 1/np.pi**2]),
            allow_different_object_types=False,
            eccentricity_threshold1=1e-2,
            eccentricity_threshold2=1e-3,
            eccentricity_threshold_penalize_shorter=20,
    ):
        self.parameters = parameters
        self.metric = metric
        self.allow_different_object_types = allow_different_object_types
        self.eccentricity_threshold1 = eccentricity_threshold1
        self.eccentricity_threshold2 = eccentricity_threshold2
        self.eccentricity_threshold_penalize_shorter = eccentricity_threshold_penalize_shorter

    def __call__(self, metadata1, metadata2, debug=False):
        if not self.allow_different_object_types:
            type1 = (
                metadata1["object_types"]
                if "object_types" in metadata1
                else "".join(sorted([
                    metadata1.get("object1", "A").upper(),
                    metadata1.get("object2", "B").upper()
                ]))
            )
            type2 = (
                metadata2["object_types"]
                if "object_types" in metadata2
                else "".join(sorted([
                    metadata2.get("object1", "C").upper(),
                    metadata2.get("object2", "D").upper()
                ]))
            )
            if type1 != type2:
                return np.inf

        values1 = [metadata1[parameter] for parameter in self.parameters]
        values2 = [metadata2[parameter] for parameter in self.parameters]

        if debug:
            print(f"{self.parameters=}")
            print(f"{values1=}")
            print(f"{values2=}")

        if "reference_mean_anomaly" in self.parameters:
            i = self.parameters.index("reference_mean_anomaly")
            values1[i], values2[i] = np.unwrap([floater(values1[i]), floater(values2[i])])

        if "reference_eccentricity" in self.parameters:
            # Either way, we first try to make sure that the corresponding entries are floats.
            i = self.parameters.index("reference_eccentricity")
            values1[i] = metadata1.get("reference_eccentricity_bound", floaterbound(values1[i]))
            values2[i] = metadata2.get("reference_eccentricity_bound", floaterbound(values2[i]))

            if values1[i] < self.eccentricity_threshold1:
                # Then we consider metadata1 a non-eccentric system...

                # ...so we ignore the mean anomaly entirely...
                if "reference_mean_anomaly" in self.parameters:
                    i_ma = self.parameters.index("reference_mean_anomaly")
                    values1[i_ma] = values2[i_ma]

                # ...and we ignore the eccentricity if metadata2 is also non-eccentric,
                # and longer than eccentricity_threshold_penalize_shorter.
                if (
                    values2[i] < self.eccentricity_threshold2
                    and metadata2.get(
                        "number_of_orbits",
                        metadata2.get("number_of_orbits_from_start", 0)
                    ) > self.eccentricity_threshold_penalize_shorter
                ):
                    values1[i] = values2[i]

        difference = (
            np.concatenate(list(map(np.atleast_1d, values1)))
            - np.concatenate(list(map(np.atleast_1d, values2)))
        )

        if debug:
            print(f"{difference=}")

        return difference @ self.metric @ difference

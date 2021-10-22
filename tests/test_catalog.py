import contextlib
import pathlib
import numpy as np
import pytest
import sxs


def test_table():
    """Ensure that we can compare floats as expected"""
    with contextlib.redirect_stdout(None):
        catalog = sxs.load("catalog")
    table = catalog.table
    assert len(table[table["object_types"] == "BHBH"]) >= 2019
    assert len(table[table["object_types"] == "BHNS"]) >= 7
    assert len(table[table["object_types"] == "NSNS"]) >= 2
    for name, dtype in zip(table.columns, table.dtypes):
        if any(name.startswith(s) for s in ["initial_", "com_correction_", "reference_"]):
            if name == "initial_data_type":
                continue
            if dtype == np.float64:
                table[table[name] > 1e-7]
            elif dtype == object:
                for i, vec in enumerate(table[name]):
                    assert isinstance(vec, np.ndarray), f"table.iloc[{i}][{name}] = {vec}"
                    assert vec.dtype == np.float64, f"table.iloc[{i}][{name}] = {vec}"
                    assert vec.size == 3, f"table.iloc[{i}][{name}] = {vec}"
            else:
                raise ValueError(f"Illegal dtype {dtype} for column {name}")
    for url in table["url"]:
        assert sxs.utilities.url.validate(url)
    for metadata_path in table["metadata_path"]:
        path = pathlib.Path(metadata_path)  # Ensure that we can do this
        assert sxs.sxs_id(str(path.parts[0]))  # Ensure that this starts with an SXS ID

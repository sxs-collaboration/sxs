import pytest

import sxs


def test_format():
    import pathlib
    import tempfile
    import json
    import h5py
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = pathlib.Path(temp_dir)
        temp_file = temp_path / "test.json"

        version_in = "abcd.efgh"
        with temp_file.open("w") as f:
            json.dump({"sxs_version": version_in}, f)
        version_out = sxs.file_format(temp_file)
        assert version_out == version_in
        version_out = sxs.file_format(str(temp_file))
        assert version_out == version_in
        with temp_file.open("r") as f:
            version_out = sxs.file_format(f)
        assert version_out == version_in

        version_in = "efgh.abcd"
        with temp_file.open("w") as f:
            json.dump({"version": version_in}, f)
        version_out = sxs.file_format(temp_file)
        assert version_out == version_in
        version_out = sxs.file_format(str(temp_file))
        assert version_out == version_in
        with temp_file.open("r") as f:
            version_out = sxs.file_format(f)
        assert version_out == version_in

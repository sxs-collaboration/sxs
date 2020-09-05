import pytest
import sxs


def test_format():
    import pathlib
    import tempfile
    import json
    import h5py
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = pathlib.Path(temp_dir)

        # Test wrong type
        temp_file = type("Failure", (object,), {})
        with pytest.raises(TypeError):
            sxs.file_format(temp_file)

        # Test missing file
        temp_file = temp_path / "nonexistent.json"
        with pytest.raises(FileNotFoundError):
            sxs.file_format(temp_file)

        # Test wrong file format
        temp_file = temp_path / "failure.json"
        with temp_file.open("wb") as f:
            f.write(b"074132708iuhufahc")
        with pytest.raises(ValueError):
            sxs.file_format(temp_file)
        with pytest.raises(ValueError):
            sxs.file_format(str(temp_file))
        with temp_file.open("r") as f:
            with pytest.raises(ValueError):
                sxs.file_format(f)
        with temp_file.open("rb") as f:
            with pytest.raises(ValueError):
                sxs.file_format(f)

        # Test `sxs_format` JSON
        temp_file = temp_path / "test.json"
        version_in = "abcd.efgh"
        with temp_file.open("w") as f:
            json.dump({"sxs_format": version_in}, f)
        version_out = sxs.file_format(temp_file)
        assert version_out == version_in
        version_out = sxs.file_format(str(temp_file))
        assert version_out == version_in
        with temp_file.open("r") as f:
            version_out = sxs.file_format(f)
        assert version_out == version_in

        # Test `sxs_format` JSON
        temp_file = temp_path / "test.json"
        version_in = "efgh.abcd"
        with temp_file.open("w") as f:
            json.dump({"format": version_in}, f)
        version_out = sxs.file_format(temp_file)
        assert version_out == version_in
        version_out = sxs.file_format(str(temp_file))
        assert version_out == version_in
        with temp_file.open("r") as f:
            version_out = sxs.file_format(f)
        assert version_out == version_in

        # Test `sxs_format` HDF5
        temp_file = temp_path / "test.h5"
        version_in = "abef.cdgh"
        with h5py.File(temp_file, "w") as f:
            f.attrs["sxs_format"] = version_in
        version_out = sxs.file_format(temp_file)
        assert version_out == version_in
        version_out = sxs.file_format(str(temp_file))
        assert version_out == version_in
        with temp_file.open("rb") as f:
            version_out = sxs.file_format(f)
        assert version_out == version_in

        # Test `format` HDF5
        temp_file = temp_path / "test.h5"
        version_in = "abgh.efcd"
        with h5py.File(temp_file, "w") as f:
            f.attrs["format"] = version_in
        version_out = sxs.file_format(temp_file)
        assert version_out == version_in
        version_out = sxs.file_format(str(temp_file))
        assert version_out == version_in
        with temp_file.open("rb") as f:
            version_out = sxs.file_format(f)
        assert version_out == version_in

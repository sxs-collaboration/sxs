import sys
import pathlib
import contextlib
import tempfile
import pytest
import json
import sxs
from .conftest import skip_macOS_GH_actions_downloads, shortest_metadata, shortest_metadata_txt


@skip_macOS_GH_actions_downloads
def test_json_conversion():
    with contextlib.redirect_stdout(None):
        sxs.load(shortest_metadata, download=True, cache=True)
        try:
            sxs.load(shortest_metadata_txt, download=True, cache=True)
        except ValueError:
            pass
    path_json = sxs.utilities.cached_path(shortest_metadata)
    path_txt = sxs.utilities.cached_path(shortest_metadata_txt)
    m = sxs.Metadata.from_txt_file(path_txt, cache_json=False)
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_dir_path = pathlib.Path(temp_dir)
        temp_path = temp_dir_path / "metadata"
        m.to_json_file(temp_path)
        with path_json.open("r") as f1:
            m1 = json.load(f1)
        with temp_path.with_suffix(".json").open("r") as f2:
            m2 = json.load(f2)
        for key in m2:
            if key.startswith("reference_") or  key.startswith("initial_"):
                assert m1[key] == m2[key]


@skip_macOS_GH_actions_downloads
def test_metadata_updates():
    with contextlib.redirect_stdout(None):
        metadata = sxs.load(shortest_metadata)
    assert metadata.reference_time != -123456.789

    # Test updating with kwargs
    metadata.update(reference_time=-123456.789)
    assert metadata.reference_time == -123456.789

    # Test updating with a Mapping
    metadata.update({"reference_time": 9.876})
    assert metadata.reference_time == 9.876

    # Test updating with a Mapping that is not a dict
    metadata.update({"reference_time": 9.876}.items())
    assert metadata.reference_time == 9.876

    # Test updating with an Iterable
    metadata.update([("reference_time", 54637.2819)])
    assert metadata.reference_time == 54637.2819

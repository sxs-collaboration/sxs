import pytest
import sxs
from .conftest import skip_macOS_GH_actions_downloads

@pytest.mark.xfail(reason="Missing files in the catalog should be updated soon")
def test_catalog_file_sizes():
    simulations = sxs.load("simulations")
    success = True
    missing_files = ""
    for sxs_id, metadata in simulations.items():
        for filename, fileinfo in metadata.get("files", {}).items():
            if fileinfo.get("size", 0) < 500:
                missing_files += f"{sxs_id} {filename}\n"
                success = False
    if not success:
        print("\nThe following files are missing or suspiciously small:")
        print(missing_files)
    assert success


@skip_macOS_GH_actions_downloads
def test_sxs_load_v2():
    # We'll use SXS:BBH:0001, even though it's deprecated, because we
    # want to be sure that it always exists, and is accessible as a
    # v2.0 record.
    with pytest.warns(UserWarning):
        s = sxs.load("SXS:BBH:0001v2.0")
    assert s.sxs_id_stem == "SXS:BBH:0001"
    assert s.version == "v2.0"
    s.h
    s.strain
    s.H
    s.Strain
    s.horizons
    s.Horizons
    s.psi4
    s.Psi4


@skip_macOS_GH_actions_downloads
def test_sxs_load_v2_levs():
    # Check that a bad Lev raises an error
    with pytest.raises(ValueError):
        sxs.load("SXS:BBH:1001v2.0/Lev189", ignore_deprecation=True)
    # Check that the default Lev works, and different Levs produce different metadata
    assert (
        sxs.load("SXS:BBH:1001v2.0").metadata.reference_time
        == sxs.load("SXS:BBH:1001v2.0/Lev3").metadata.reference_time
    )
    assert (
        sxs.load("SXS:BBH:1001v2.0/Lev1").metadata.reference_time
        != sxs.load("SXS:BBH:1001v2.0/Lev3").metadata.reference_time
    )
    for lev, lev_number in [("", 3), ("/Lev3", 3), ("/Lev2", 2), ("/Lev1", 1)]:
        s = sxs.load(f"SXS:BBH:1001v2.0{lev}", ignore_deprecation=True)
        assert s.sxs_id_stem == "SXS:BBH:1001"
        assert s.version == "v2.0"
        assert s.lev_number == lev_number
        s.h
        s.strain
        s.H
        s.Strain
        s.horizons
        s.Horizons


@skip_macOS_GH_actions_downloads
def test_sxs_load_v3():
    # # This is necessary only after the deprecation of SXS:BBH:4000v3.0
    # with pytest.warns(UserWarning):
    #     s = sxs.load("SXS:BBH:4000v3.0")
    s = sxs.load("SXS:BBH:4000v3.0")
    assert s.sxs_id_stem == "SXS:BBH:4000"
    assert s.version == "v3.0"
    s.h
    s.strain
    s.H
    s.Strain
    s.horizons
    s.Horizons
    s.psi4
    s.Psi4


@skip_macOS_GH_actions_downloads
def test_sxs_load_v3_levs():
    # Check that a bad Lev raises an error
    with pytest.raises(ValueError):
        sxs.load("SXS:BBH:4000v3.0/Lev189", ignore_deprecation=True)
    # Check that the default Lev works, and different Levs produce different metadata
    assert (  # Check that the default Lev loads the highest Lev
        sxs.load("SXS:BBH:4000v3.0", ignore_deprecation=True).metadata.reference_time
        == sxs.load("SXS:BBH:4000v3.0/Lev3", ignore_deprecation=True).metadata.reference_time
    )
    assert (  # Check that different Levs have different reference times
        sxs.load("SXS:BBH:4000v3.0/Lev1", ignore_deprecation=True).metadata.reference_time
        != sxs.load("SXS:BBH:4000v3.0/Lev3", ignore_deprecation=True).metadata.reference_time
    )
    for lev, lev_number in [("", 3), ("/Lev3", 3), ("/Lev2", 2), ("/Lev1", 1)]:
        s = sxs.load(f"SXS:BBH:4000v3.0{lev}", ignore_deprecation=True)
        assert s.sxs_id_stem == "SXS:BBH:4000"
        assert s.version == "v3.0"
        assert s.lev_number == lev_number
        s.h
        s.strain
        s.H
        s.Strain
        s.horizons
        s.Horizons
    # Check that Lev-1 and Lev0 work and are not the same as the highest Lev
    s = sxs.load(f"SXS:BBH:1132v3.0", ignore_deprecation=True)
    sm1 = sxs.load(f"SXS:BBH:1132v3.0/Lev-1", ignore_deprecation=True)
    s0 = sxs.load(f"SXS:BBH:1132v3.0/Lev0", ignore_deprecation=True)
    assert s.lev_number == 4
    assert sm1.lev_number == -1
    assert s0.lev_number == 0
    assert s.metadata != sm1.metadata
    assert s.metadata != s0.metadata
    assert sm1.metadata != s0.metadata


@pytest.mark.many_downloads  # runs only with pytest flag ` --run_many_downloads`
def test_sxs_load_all_doi_versions():
    print()
    for sxs_id in list(range(303, 324+1)) + list(range(4425, 4434+1)):
        sxs_id = f"SXS:BBH:{sxs_id:04}"
        print(f"    Loading {sxs_id}")
        sim = sxs.load(sxs_id, ignore_deprecation=True)
        for doi_version in sim.metadata["DOI_versions"]:
            if not doi_version:  # Skip the default version, which was loaded above
                continue
            # This line will actually load the metadata, which is already a strong test
            sim = sxs.load(f"{sxs_id}{doi_version}", ignore_deprecation=True)
            assert sim.metadata_path
            assert sim.horizons_path
            assert sim.strain_path
            assert sim.psi4_path


@pytest.mark.many_downloads  # runs only with pytest flag ` --run_many_downloads`
def test_sxs_load_v3_catalog():
    # This test goes right up to the point of loading every simulation
    # at each Lev with each extrapolation order, but doesn't download
    # any data
    from pathlib import Path

    def fake_load_metadata(self):
        metadata_path = self.metadata_path
        if metadata_path not in self.files:
            return False
        if "link" not in self.files[metadata_path]:
            return False
        if self.files[metadata_path].get("size", 0) < 500:
            return False
        return True

    def fake_load_horizons(self):
        horizons_path = self.horizons_path
        if horizons_path in self.files:
            if "link" not in self.files[horizons_path]:
                return False
        else:
            # Some simulations used the SXS ID as a prefix in file paths
            # within the Zenodo upload in version 1.x of the catalog.
            if (extended_location := f"{self.sxs_id_stem}/{horizons_path}") in self.files:
                if "link" not in self.files[extended_location]:
                    return False
                elif self.files[extended_location].get("size", 0) < 500:
                    return False
            else:
                return False
        return True

    def fake_load_waveform(self, file_name, group):
        file_name = Path(file_name)
        h5_path = str(file_name.with_suffix(".h5"))
        json_path = str(file_name.with_suffix(".json"))
        return bool(
            self.files.get(h5_path, {}).get("link", False)
            and self.files.get(h5_path, {}).get("size", 0) > 500
            and self.files.get(json_path, {}).get("link", False)
            and self.files.get(json_path, {}).get("size", 0) > 500
        )

    simulations = sxs.load("simulations")

    success = True
    for sxs_id, metadata in simulations.items():
        if "NSNS" in sxs_id or "BHNS" in sxs_id:
            continue

        sim = sxs.load(sxs_id, ignore_deprecation=True)
        if sim.version != "v3.0":
            # print(
            #     f"Skipping {sxs_id} because its version is {sim.version} != v3.0"
            # )
            continue

        for lev_number in metadata["lev_numbers"]:
            for extrapolation in ["Outer", "N2", "N3", "N4"]:
                sim = sxs.load(
                    f"{sxs_id}/Lev{lev_number}",
                    extrapolation=extrapolation,
                    ignore_deprecation=True,
                )

                if not fake_load_metadata(sim):
                    print(f"Failed to load {sxs_id}/Lev{lev_number} {extrapolation} metadata")
                    success = False

                if not fake_load_horizons(sim):
                    print(f"Failed to load {sxs_id}/Lev{lev_number} {extrapolation} horizons")
                    success = False

                if not fake_load_waveform(sim, *sim.strain_path):
                    print(f"Failed to load {sxs_id}/Lev{lev_number} {extrapolation} strain")
                    success = False

                if not fake_load_waveform(sim, *sim.psi4_path):
                    print(f"Failed to load {sxs_id}/Lev{lev_number} {extrapolation} psi4")
                    success = False

    assert success


@skip_macOS_GH_actions_downloads
@pytest.mark.parametrize("loader", [sxs.Simulation, sxs.load])
def test_superseding(loader):
    simulation = "SXS:BBH:0001"

    # No version, no arguments
    with pytest.raises(ValueError):
        loader(f"{simulation}")

    # Version, no arguments
    with pytest.warns(UserWarning):
        s = loader(f"{simulation}v2.0")
    assert s.sxs_id_stem == simulation

    # No version, ignore_deprecation
    s = loader(f"{simulation}", ignore_deprecation=True)
    assert s.sxs_id_stem == simulation

    # Version, ignore_deprecation
    s = loader(f"{simulation}v2.0", ignore_deprecation=True)
    assert s.sxs_id_stem == simulation

    # No version, auto_supersede
    with pytest.warns(UserWarning):
        s = loader(f"{simulation}", auto_supersede=True)
    assert s.sxs_id_stem != simulation

    # Version, auto_supersede
    with pytest.warns(UserWarning):
        s = loader(f"{simulation}v2.0", auto_supersede=True)
    assert s.sxs_id_stem == simulation

    # No version, ignore_deprecation, auto_supersede
    s = loader(f"{simulation}", ignore_deprecation=True, auto_supersede=True)
    assert s.sxs_id_stem == simulation

    # Version, ignore_deprecation, auto_supersede
    s = loader(f"{simulation}v2.0", ignore_deprecation=True, auto_supersede=True)
    assert s.sxs_id_stem == simulation

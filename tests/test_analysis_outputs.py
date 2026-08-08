import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import nemo
from nemo.analysis import Ensemble


SOLVENT = (2.38, 1.49)
TESTS_DIR = Path(__file__).resolve().parent
DATASET_DIRS = {
    "tddft": TESTS_DIR / "tddft",
    "eom-ccsd": TESTS_DIR / "eom-ccsd",
}


def assert_valid_dataframe(df: pd.DataFrame, name: str) -> None:
    assert isinstance(df, pd.DataFrame), f"{name} did not return a pandas DataFrame"
    assert not df.empty, f"{name} returned an empty DataFrame"

    nan_cols = df.columns[df.isna().any()].tolist()
    assert not nan_cols, f"{name} contains NaN values in columns: {nan_cols}"

    numeric_df = df.select_dtypes(include=[np.number])
    assert not numeric_df.empty, f"{name} has no numeric columns to validate"

    # Validate column-wise to avoid pandas mixed-dtype consolidation warnings.
    bad_cols = [
        col
        for col in numeric_df.columns
        if not np.isfinite(numeric_df[col].to_numpy()).all()
    ]
    assert not bad_cols, f"{name} contains non-finite values in columns: {bad_cols}"


@pytest.fixture(scope="module", params=["tddft", "eom-ccsd"], ids=["tddft", "eom-ccsd"])
def generated_ensembles(request):
    """
    Run gather_data once from the tests directory so the Ensemble_*.lx files
    are created inside tests/.
    """
    dataset_name = request.param
    dataset_dir = DATASET_DIRS[dataset_name]
    geom_dir = dataset_dir / "Geometries"
    freq_file = dataset_dir / "freqs1.log"
    mag_file = dataset_dir / "Magnitudes_300K_.lx"

    expected_ensemble_files = {
        "s0": dataset_dir / "Ensemble_S0_.lx",
        "s1": dataset_dir / "Ensemble_S1_.lx",
        "t1": dataset_dir / "Ensemble_T1_.lx",
    }

    assert geom_dir.exists(), f"Missing folder: {geom_dir}"
    assert geom_dir.is_dir(), f"{geom_dir} is not a directory"
    assert freq_file.exists(), f"Missing file: {freq_file}"
    assert mag_file.exists(), f"Missing file: {mag_file}"

    old_cwd = Path.cwd()
    try:
        # gather_data should see ./Geometries and write ./Ensemble_*.lx in this dataset folder.
        os.chdir(dataset_dir)

        for state in ["s0", "s1", "t1"]:
            nemo.analysis.gather_data(state, save=True)

        for state, filepath in expected_ensemble_files.items():
            assert filepath.exists(), f"Ensemble generation failed for {state}: {filepath.name} was not created"
            assert filepath.is_file(), f"Expected file was not created for {state}: {filepath}"
            assert filepath.stat().st_size > 0, f"Generated file is empty for {state}: {filepath.name}"

        ensembles = {
            state: Ensemble(str(filepath))
            for state, filepath in expected_ensemble_files.items()
        }

        ensembles["dataset_name"] = dataset_name
        ensembles["dataset_dir"] = dataset_dir

        os.chdir(old_cwd)
        yield ensembles

    finally:
        os.chdir(old_cwd)
        for filepath in expected_ensemble_files.values():
            filepath.unlink(missing_ok=True)


def test_ensemble_files_are_generated(generated_ensembles):
    dataset_dir = generated_ensembles["dataset_dir"]
    expected_ensemble_files = {
        "s0": dataset_dir / "Ensemble_S0_.lx",
        "s1": dataset_dir / "Ensemble_S1_.lx",
        "t1": dataset_dir / "Ensemble_T1_.lx",
    }
    for state, filepath in expected_ensemble_files.items():
        assert filepath.exists(), f"{filepath.name} was not generated for {state}"
        assert filepath.stat().st_size > 0, f"{filepath.name} is empty for {state}"


def test_s0_absorption_default(generated_ensembles):
    absorption = generated_ensembles["s0"].absorption(SOLVENT)
    assert_valid_dataframe(absorption, "s0.absorption(solvent)")


def test_s0_absorption_wavelength_extinction(generated_ensembles):
    absorption = generated_ensembles["s0"].absorption(
        SOLVENT,
        wavelength=True,
        extinction=True,
    )
    assert_valid_dataframe(
        absorption,
        "s0.absorption(solvent, wavelength=True, extinction=True)",
    )


def test_s1_rates(generated_ensembles):
    rates = generated_ensembles["s1"].rate(SOLVENT)
    assert_valid_dataframe(rates, "s1.rate(solvent)")
    assert "AvgCoupling(meV)" in rates.columns


def test_s1_absorption(generated_ensembles):
    absorption = generated_ensembles["s1"].absorption(SOLVENT)
    assert_valid_dataframe(absorption, "s1.absorption(solvent)")


def test_s1_emission_default(generated_ensembles):
    fluor = generated_ensembles["s1"].emission(SOLVENT)
    assert_valid_dataframe(fluor, "s1.emission(solvent)")


def test_s1_emission_wavelength(generated_ensembles):
    fluor = generated_ensembles["s1"].emission(SOLVENT, wavelength=True)
    assert_valid_dataframe(fluor, "s1.emission(solvent, wavelength=True)")


def test_t1_rates(generated_ensembles):
    rates_t1 = generated_ensembles["t1"].rate(SOLVENT)
    assert_valid_dataframe(rates_t1, "t1.rate(solvent)")
    assert "AvgCoupling(meV)" in rates_t1.columns


def test_t1_emission(generated_ensembles):
    phosph = generated_ensembles["t1"].emission(SOLVENT)
    assert_valid_dataframe(phosph, "t1.emission(solvent)")


def test_s1_breakdown(generated_ensembles):
    details = generated_ensembles["s1"].breakdown(SOLVENT)
    assert_valid_dataframe(details, "s1.breakdown(solvent)")
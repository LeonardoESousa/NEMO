import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import nemo
from nemo.analysis import Ensemble


SOLVENT = (2.38, 1.49)
TESTS_DIR = Path(__file__).resolve().parent
GEOM_DIR = TESTS_DIR / "Geometries"

EXPECTED_ENSEMBLE_FILES = {
    "s0": TESTS_DIR / "Ensemble_S0_.lx",
    "s1": TESTS_DIR / "Ensemble_S1_.lx",
    "t1": TESTS_DIR / "Ensemble_T1_.lx",
}


def assert_valid_dataframe(df: pd.DataFrame, name: str) -> None:
    assert isinstance(df, pd.DataFrame), f"{name} did not return a pandas DataFrame"
    assert not df.empty, f"{name} returned an empty DataFrame"

    nan_cols = df.columns[df.isna().any()].tolist()
    assert not nan_cols, f"{name} contains NaN values in columns: {nan_cols}"

    numeric_df = df.select_dtypes(include=[np.number])
    assert not numeric_df.empty, f"{name} has no numeric columns to validate"

    finite_mask = np.isfinite(numeric_df.to_numpy())
    bad_cols = numeric_df.columns[~np.isfinite(numeric_df).all(axis=0)].tolist()
    assert finite_mask.all(), f"{name} contains non-finite values in columns: {bad_cols}"


@pytest.fixture(scope="module")
def generated_ensembles():
    """
    Run gather_data once from the tests directory so the Ensemble_*.lx files
    are created inside tests/.
    """
    assert GEOM_DIR.exists(), f"Missing folder: {GEOM_DIR}"
    assert GEOM_DIR.is_dir(), f"{GEOM_DIR} is not a directory"

    old_cwd = Path.cwd()
    try:
        # gather_data should see ./Geometries and write ./Ensemble_*.lx here
        os.chdir(TESTS_DIR)

        for state in ["s0", "s1", "t1"]:
            nemo.analysis.gather_data(state, save=True)

        for state, filepath in EXPECTED_ENSEMBLE_FILES.items():
            assert filepath.exists(), f"Ensemble generation failed for {state}: {filepath.name} was not created"
            assert filepath.is_file(), f"Expected file was not created for {state}: {filepath}"
            assert filepath.stat().st_size > 0, f"Generated file is empty for {state}: {filepath.name}"

        ensembles = {
            state: Ensemble(str(filepath))
            for state, filepath in EXPECTED_ENSEMBLE_FILES.items()
        }

        os.chdir(old_cwd)
        yield ensembles

    finally:
        os.chdir(old_cwd)
        for filepath in EXPECTED_ENSEMBLE_FILES.values():
            filepath.unlink(missing_ok=True)


def test_ensemble_files_are_generated(generated_ensembles):
    for state, filepath in EXPECTED_ENSEMBLE_FILES.items():
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


def test_t1_emission(generated_ensembles):
    phosph = generated_ensembles["t1"].emission(SOLVENT)
    assert_valid_dataframe(phosph, "t1.emission(solvent)")


def test_s1_breakdown(generated_ensembles):
    details = generated_ensembles["s1"].breakdown(SOLVENT)
    assert_valid_dataframe(details, "s1.breakdown(solvent)")

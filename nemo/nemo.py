"""Convenience object API for working with NEMO ensemble files."""

import os
import re
from collections.abc import Iterable

import numpy as np
import pandas as pd

from nemo.analysis import Ensemble


STATE_PATTERN = re.compile(r"([ST])(\d+)", re.IGNORECASE)
SINGLET_TRANSFER_PATTERN = re.compile(r"^S(\d+)~>S(\d+)$", re.IGNORECASE)
THECOLOR = "black"


def _normalize_transition(transition):
    return str(transition).upper().replace(" ", "")


def _transition_initial(transition):
    return _normalize_transition(transition).split(">")[0].rstrip("-~")


def _transition_states(transition):
    return [
        (spin.upper(), int(number))
        for spin, number in STATE_PATTERN.findall(str(transition))
    ]


def _default_state_limits(transitions):
    highest = {"S": 0, "T": 0}
    for transition in transitions:
        for spin, number in _transition_states(transition):
            highest[spin] = max(highest[spin], number)

    return tuple(
        highest[spin] - 2 if highest[spin] > 2 else highest[spin]
        for spin in ("S", "T")
    )


def _state_limits(states, transitions):
    if states is None:
        return _default_state_limits(transitions)

    try:
        singlet_limit, triplet_limit = states
    except (TypeError, ValueError) as exc:
        raise ValueError(
            "states must be a tuple of two integers: (max_s, max_t)."
        ) from exc

    limits = []
    for limit in (singlet_limit, triplet_limit):
        if int(limit) != limit:
            raise ValueError("states must be a tuple of two integers: (max_s, max_t).")
        limits.append(int(limit))
    return tuple(limits)


def _within_state_limits(transition, states):
    singlet_limit, triplet_limit = states
    for spin, number in _transition_states(transition):
        if spin == "S" and number > singlet_limit:
            return False
        if spin == "T" and number > triplet_limit:
            return False
    return True


def _has_required_opposite_transition(transition, transitions):
    normalized = _normalize_transition(transition)
    match = SINGLET_TRANSFER_PATTERN.fullmatch(normalized)
    if match is None:
        return True

    initial, final = (int(number) for number in match.groups())
    if final <= initial:
        return True

    return f"S{final}~>S{initial}" in transitions


def _recompute_transition_probabilities(total_rates):
    rates = total_rates["Rate"].astype(float)
    initial_states = total_rates["Transition"].map(_transition_initial)
    total_rate = rates.groupby(initial_states).transform("sum")

    with np.errstate(divide="ignore", invalid="ignore"):
        probabilities = np.where(
            total_rate.to_numpy() != 0,
            100 * rates.to_numpy() / total_rate.to_numpy(),
            0.0,
        )

    total_rates["Prob"] = probabilities
    return total_rates


def compile_rates(dielec, ensembles, ensemble_average=False, states=None):
    """Stack and filter rate tables from Ensemble objects."""
    total_rates = None
    for ensemble in ensembles:
        rates = ensemble.rate(dielec, ensemble_average=ensemble_average)
        if total_rates is None:
            total_rates = rates
        else:
            total_rates = pd.concat([total_rates, rates], axis=0, ignore_index=True)

    if total_rates is None:
        raise ValueError("At least one Ensemble object is required.")

    total_rates.rename(columns=lambda column: column.split("(")[0], inplace=True)

    limits = _state_limits(states, total_rates["Transition"])
    total_rates = total_rates[
        total_rates["Transition"].map(
            lambda transition: _within_state_limits(transition, limits)
        )
    ].copy()

    transitions = {
        _normalize_transition(transition)
        for transition in total_rates["Transition"].to_numpy()
    }
    total_rates = total_rates[
        total_rates["Transition"].map(
            lambda transition: _has_required_opposite_transition(transition, transitions)
        )
    ].copy()

    _recompute_transition_probabilities(total_rates)
    total_rates.sort_values(by="Transition", inplace=True)
    total_rates.reset_index(drop=True, inplace=True)
    return total_rates


def compute_corrected_yields(df, initial_state="S1"):
    """Correct final yields using an absorbing Markov chain."""
    out = df.copy()
    initial_state = str(initial_state).upper()
    absorbing = "S0"

    if out.empty:
        return out

    if "Transition" not in out or "Rate" not in out:
        raise ValueError("df must contain 'Transition' and 'Rate' columns.")

    # Parse transitions
    parsed = out["Transition"].str.extract(r"^(.+?)(->|~>)(.+)$")
    if parsed.isna().any(axis=None):
        invalid = out.loc[parsed.isna().any(axis=1), "Transition"].to_list()
        raise ValueError(f"Could not parse transitions: {invalid}")
    out["_from"] = parsed[0].str.strip().str.upper()
    out["_to"] = parsed[2].str.strip().str.upper()

    # ------------------------------------------------------------------
    # Add virtual relaxation channels for states with no outgoing rates
    # ------------------------------------------------------------------

    all_states = set(out["_from"]).union(out["_to"])
    all_states.discard(absorbing)

    states_with_outgoing = set(out["_from"])

    def lower_state(state):
        m = re.match(r"^([ST])(\d+)$", state)
        if not m:
            return None

        spin, n = m.group(1), int(m.group(2))

        if spin == "S":
            if n == 0:
                return None
            return f"S{n-1}"

        if spin == "T":
            if n == 1:
                return "S0"
            return f"T{n-1}"

    virtual_rows = []

    for state in all_states:
        if state not in states_with_outgoing:
            target = lower_state(state)

            if target is None:
                continue

            kind = "~>" if (state.startswith("T") and target == "S0") else "->"

            row = {c: np.nan for c in out.columns}
            row["Transition"] = f"{state}{kind}{target}"
            row["Rate"] = np.inf
            row["_from"] = state
            row["_to"] = target

            virtual_rows.append(row)

    if virtual_rows:
        out = pd.concat([out, pd.DataFrame(virtual_rows)], ignore_index=True)

    # ------------------------------------------------------------------
    # Recompute local yields (Prob)
    # ------------------------------------------------------------------

    probs = np.zeros(len(out))

    for state, grp in out.groupby("_from"):
        rates = grp["Rate"].to_numpy(dtype=float)

        if np.isinf(rates).any():
            p = np.isinf(rates).astype(float)
            p /= p.sum()
        else:
            rate_sum = rates.sum()
            if rate_sum == 0:
                p = np.zeros(len(rates))
            else:
                p = rates / rate_sum

        probs[grp.index] = p

    out["Prob"] = probs

    # ------------------------------------------------------------------
    # Build transient-state Markov matrix
    # ------------------------------------------------------------------

    transient = sorted(set(out["_from"]) - {absorbing})
    idx = {s: i for i, s in enumerate(transient)}
    if initial_state not in idx:
        raise ValueError(f"Initial state '{initial_state}' is not present in rates.")

    Q = np.zeros((len(transient), len(transient)))

    for _, row in out.iterrows():
        i = row["_from"]
        j = row["_to"]

        if i in idx and j in idx:
            Q[idx[i], idx[j]] += row["Prob"]

    try:
        N = np.linalg.inv(np.eye(len(Q)) - Q)
    except np.linalg.LinAlgError as exc:
        raise ValueError(
            "Yield correction requires every reachable state to eventually reach "
            f"{absorbing}."
        ) from exc

    # ------------------------------------------------------------------
    # Compute final yields for channels ending in S0
    # ------------------------------------------------------------------

    x = idx[initial_state]

    final_yield = np.full(len(out), np.nan)

    for r, row in out.iterrows():
        if row["_to"] == absorbing:
            final_yield[r] = N[x, idx[row["_from"]]] * row["Prob"]

    # Replace Prob by corrected yields where applicable
    mask = ~np.isnan(final_yield)
    out.loc[mask, "Prob"] = final_yield[mask]

    # Remove helper columns and virtual rows
    out = out[out["Rate"] != np.inf].copy()

    out.drop(columns=["_from", "_to"], inplace=True)

    # multiply by 100 to get percentages
    out["Prob"] *= 100

    return out.reset_index(drop=True)


class Molecule:
    """Collection of state-specific Ensemble objects for one molecule."""

    def __init__(self, *ensemble_files, name=""):
        if (
            len(ensemble_files) == 1
            and isinstance(ensemble_files[0], Iterable)
            and not isinstance(ensemble_files[0], (str, bytes, os.PathLike))
        ):
            ensemble_files = tuple(ensemble_files[0])

        if not ensemble_files:
            raise ValueError("At least one ensemble file is required.")

        self.name = name
        self.ensembles = {}
        for file in ensemble_files:
            ensemble = Ensemble(file, name=name)
            state = self._normalize_state(ensemble.initial)
            if state in self.ensembles:
                raise ValueError(f"Duplicate ensemble for state '{state}'.")
            self.ensembles[state] = ensemble

    @staticmethod
    def _normalize_state(state):
        return str(state).upper()

    @property
    def states(self):
        return tuple(self.ensembles.keys())

    def ensemble(self, state):
        state = self._normalize_state(state)
        try:
            return self.ensembles[state]
        except KeyError as exc:
            available = ", ".join(self.states)
            raise KeyError(
                f"State '{state}' is not available. Available states: {available}."
            ) from exc

    def rates(
        self,
        dielec,
        ensemble_average=False,
        states=None,
        initial_state="S1",
    ):
        rates = compile_rates(
            dielec,
            self.ensembles.values(),
            ensemble_average=ensemble_average,
            states=states,
        )
        return compute_corrected_yields(
            rates,
            initial_state=initial_state,
        )

    def rate(self, state, dielec, ensemble_average=False):
        return self.ensemble(state).rate(
            dielec,
            ensemble_average=ensemble_average,
        )

    def emission(self, state, dielec, wavelength=False):
        return self.ensemble(state).emission(dielec, wavelength=wavelength)

    def complete_emi(self, state, dielec, ensemble_average=False, wavelength=False):
        return self.ensemble(state).complete_emi(
            dielec,
            ensemble_average=ensemble_average,
            wavelength=wavelength,
        )

    def complete_abs(self, state, dielec, nstates=-1, wavelength=False, extinction=False):
        return self.ensemble(state).complete_abs(
            dielec,
            nstates=nstates,
            wavelength=wavelength,
            extinction=extinction,
        )

    def absorption(self, state, dielec, nstates=-1, wavelength=False, extinction=False):
        return self.ensemble(state).absorption(
            dielec,
            nstates=nstates,
            wavelength=wavelength,
            extinction=extinction,
        )

    def breakdown(self, state, dielec):
        return self.ensemble(state).breakdown(dielec)

    def save(self, state, dielec, mode):
        return self.ensemble(state).save(dielec, mode)

    def emi2wavelength(self, state, emi):
        return self.ensemble(state).emi2wavelength(emi)

    def abs2wavelength(self, state, abs_spec):
        return self.ensemble(state).abs2wavelength(abs_spec)

    def abs2extinction(self, state, abs_spec):
        return self.ensemble(state).abs2extinction(abs_spec)

"""Session state management with M1 -> M2 -> M3 invalidation cascade.

Trust cascades monotonically:
    M1 input change -> invalidate M2 results + M3 results
    M2 input change -> invalidate M3 results only
    M3 input change -> invalidate M3 results only

Hash-based change detection: inputs are serialised to JSON and SHA-256 hashed.
A change in hash means downstream state must be cleared.

Intended for use with Streamlit session_state or any dict-like store.
"""

from __future__ import annotations

import hashlib
import json
from typing import Optional


# ─── Keys used in Streamlit session_state ─────────────────────────────────────

# Results keys (cleared on invalidation)
_M1_RESULT_KEY = "m1_result"
_M2_RESULT_KEY = "m2_result"
_M3_RESULT_KEY = "m3_result"

# Trust level keys (cleared on invalidation)
_M1_TRUST_KEY = "m1_trust"
_M2_TRUST_KEY = "m2_trust"

# Validation result keys (cleared on invalidation)
_M1_VALIDATION_KEY = "m1_validation"
_M2_VALIDATION_KEY = "m2_validation"
_M3_VALIDATION_KEY = "m3_validation"

# Groups of keys cleared per invalidation level
_M1_DOWNSTREAM_KEYS = [
    _M2_RESULT_KEY, _M2_TRUST_KEY, _M2_VALIDATION_KEY,
    _M3_RESULT_KEY, _M3_VALIDATION_KEY,
]
_M2_DOWNSTREAM_KEYS = [
    _M3_RESULT_KEY, _M3_VALIDATION_KEY,
]


def _hash_inputs(inputs: dict) -> str:
    """Deterministic SHA-256 hash of a JSON-serialisable input dict."""
    try:
        serialised = json.dumps(inputs, sort_keys=True, default=str)
    except (TypeError, ValueError):
        # Fallback: hash the repr — less stable but always works
        serialised = repr(sorted(inputs.items()))
    return hashlib.sha256(serialised.encode("utf-8")).hexdigest()


class SessionStateManager:
    """Manages EmulSim session state with M1->M2->M3 invalidation cascade.

    Maintains hashes of the last M1/M2/M3 inputs used for run.
    When inputs change (hash mismatch), downstream state keys are cleared.

    Usage (Streamlit)::

        import streamlit as st
        from emulsim.visualization.ui_state import SessionStateManager

        if "state_mgr" not in st.session_state:
            st.session_state["state_mgr"] = SessionStateManager()
        mgr = st.session_state["state_mgr"]

        if mgr.update_m1(m1_inputs):
            st.info("M1 inputs changed — M2 and M3 results cleared.")
    """

    def __init__(self) -> None:
        self.m1_hash: str = ""
        self.m2_hash: str = ""
        self.m3_hash: str = ""

        # Optional reference to a dict-like store (e.g. st.session_state).
        # When None, invalidation only resets the internal hashes.
        self._store: Optional[dict] = None

    def bind_store(self, store: dict) -> None:
        """Attach a dict-like store (e.g. Streamlit session_state).

        When bound, invalidate_downstream() will also delete result keys from
        the store so the UI rerenders with cleared outputs.
        """
        self._store = store

    # ── Hash-based change detection ───────────────────────────────────────

    def update_m1(self, inputs: dict) -> bool:
        """Register new M1 inputs and detect changes.

        Args:
            inputs: Dict of M1 input values (must be JSON-serialisable).

        Returns:
            True if the inputs changed from the previous call (M2 and M3
            results have been invalidated). False if inputs are identical.
        """
        new_hash = _hash_inputs(inputs)
        if new_hash == self.m1_hash:
            return False
        self.m1_hash = new_hash
        self.invalidate_downstream(from_module=1)
        return True

    def update_m2(self, inputs: dict) -> bool:
        """Register new M2 inputs and detect changes.

        Args:
            inputs: Dict of M2 input values (must be JSON-serialisable).

        Returns:
            True if the inputs changed (M3 results have been invalidated).
            False if inputs are identical.
        """
        new_hash = _hash_inputs(inputs)
        if new_hash == self.m2_hash:
            return False
        self.m2_hash = new_hash
        self.invalidate_downstream(from_module=2)
        return True

    def update_m3(self, inputs: dict) -> bool:
        """Register new M3 inputs and detect changes.

        Args:
            inputs: Dict of M3 input values (must be JSON-serialisable).

        Returns:
            True if the inputs changed (M3 results have been cleared).
            False if inputs are identical.
        """
        new_hash = _hash_inputs(inputs)
        if new_hash == self.m3_hash:
            return False
        self.m3_hash = new_hash
        self.invalidate_downstream(from_module=3)
        return True

    # ── Invalidation ──────────────────────────────────────────────────────

    def invalidate_downstream(self, from_module: int) -> None:
        """Clear all downstream session state keys and reset downstream hashes.

        Args:
            from_module: Module number where the change originated (1, 2, or 3).
                - 1: clears M2 and M3 results + resets m2_hash, m3_hash
                - 2: clears M3 results + resets m3_hash
                - 3: clears M3 results
        """
        if from_module == 1:
            self.m2_hash = ""
            self.m3_hash = ""
            keys_to_clear = _M1_DOWNSTREAM_KEYS
        elif from_module == 2:
            self.m3_hash = ""
            keys_to_clear = _M2_DOWNSTREAM_KEYS
        elif from_module == 3:
            keys_to_clear = _M2_DOWNSTREAM_KEYS  # same: only M3 results
        else:
            raise ValueError(f"from_module must be 1, 2, or 3; got {from_module!r}")

        if self._store is not None:
            for key in keys_to_clear:
                self._store.pop(key, None)

    def reset_all(self) -> None:
        """Full reset: clear all hashes and all known result keys from store."""
        self.m1_hash = ""
        self.m2_hash = ""
        self.m3_hash = ""
        if self._store is not None:
            all_keys = (
                [_M1_RESULT_KEY, _M1_TRUST_KEY, _M1_VALIDATION_KEY]
                + _M1_DOWNSTREAM_KEYS
            )
            for key in all_keys:
                self._store.pop(key, None)

    # ── Convenience query ─────────────────────────────────────────────────

    def m1_has_run(self) -> bool:
        """True if M1 has been run at least once (hash is non-empty)."""
        return bool(self.m1_hash)

    def m2_has_run(self) -> bool:
        """True if M2 has been run with current M1 inputs."""
        return bool(self.m2_hash)

    def m3_has_run(self) -> bool:
        """True if M3 has been run with current M2 inputs."""
        return bool(self.m3_hash)

    def __repr__(self) -> str:
        return (
            f"SessionStateManager("
            f"m1={'set' if self.m1_hash else 'empty'}, "
            f"m2={'set' if self.m2_hash else 'empty'}, "
            f"m3={'set' if self.m3_hash else 'empty'})"
        )

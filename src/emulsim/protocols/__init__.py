"""protocols package — reaction mechanism data and protocol documents.

Public API
----------
From mechanism_data (WN-1):
    MechanismDescriptor  — full mechanistic description for one reagent
    ReactionStep         — a single elementary/formal reaction step
    MECHANISM_REGISTRY   — dict[str, MechanismDescriptor], manually authored entries
    get_mechanism        — look up or auto-generate a MechanismDescriptor by key

From protocol_document (WN-3, created separately):
    ProtocolDocument     — structured, printable protocol object
    ProtocolStep         — one procedural step in a protocol
    ReagentRequirement   — reagent + quantity record for a protocol

Notes
-----
The protocol_document imports are conditional so that this package can be
imported even before WN-3 creates protocol_document.py.  Callers that need
ProtocolDocument, ProtocolStep, or ReagentRequirement should import them
directly from emulsim.protocols.protocol_document once that module exists.
"""

from __future__ import annotations

from .mechanism_data import (
    MechanismDescriptor,
    ReactionStep,
    MECHANISM_REGISTRY,
    get_mechanism,
)

__all__ = [
    "MechanismDescriptor",
    "ReactionStep",
    "MECHANISM_REGISTRY",
    "get_mechanism",
]

# ── Conditional import of protocol_document ──────────────────────────────
# protocol_document.py was created in a separate work note.
# Import it here so that callers can use `from emulsim.protocols import
# ProtocolDocument` without breaking imports if the module is absent.
try:
    from .protocol_document import (  # type: ignore[import]
        ProtocolDocument,
        ProtocolStep,
        ReagentRequirement,
    )

    __all__ += ["ProtocolDocument", "ProtocolStep", "ReagentRequirement"]
except ModuleNotFoundError:
    # protocol_document not yet available — symbols unavailable.
    pass

"""Protocol document data model for generated wet-lab protocols.

Contains the output dataclasses that ProtocolGenerator produces.
These are pure data containers with a Markdown export method.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime


@dataclass
class ReagentRequirement:
    """A reagent needed for the protocol with calculated amounts."""

    name: str
    cas: str
    amount: str   # "300 mg" or "10.0 mL"
    grade: str    # ">=98%" or "AR grade"
    role: str     # "crosslinker", "buffer", "photoinitiator"


@dataclass
class ProtocolStep:
    """One numbered step in the generated wet-lab protocol."""

    step_number: int
    action: str         # "Dissolve 25 mg Irgacure 2959 in 1.0 mL ethanol"
    duration: str       # "2 h" or "Immediate"
    temperature: str    # "25 deg C" or "37 deg C"
    notes: str = ""     # "Protect from light"


@dataclass
class ProtocolDocument:
    """Complete generated wet-lab protocol document."""

    title: str
    reagent_key: str
    generation_timestamp: str
    # User parameters that shaped this protocol
    user_temperature_K: float
    user_time_s: float
    user_concentration_mM: float
    user_pH: float
    user_bead_volume_mL: float        # Settled bead volume (default 1.0 mL)
    # Generated content
    reagent_table: list[ReagentRequirement] = field(default_factory=list)
    procedure_steps: list[ProtocolStep] = field(default_factory=list)
    qc_targets: list[str] = field(default_factory=list)
    safety_warnings: list[str] = field(default_factory=list)
    mechanism_summary: str = ""       # Brief mechanism text for header

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _format_time(seconds: float) -> str:
        """Return a human-readable time string from a value in seconds."""
        if seconds >= 3600:
            return f"{seconds / 3600:.1f} h"
        if seconds >= 60:
            minutes = seconds / 60
            # Use integer display when the result is a whole number of minutes
            if minutes == int(minutes):
                return f"{int(minutes)} min"
            return f"{minutes:.1f} min"
        return f"{seconds:.0f} s"

    @staticmethod
    def _escape_cell(value: str) -> str:
        """Escape pipe characters inside a Markdown table cell."""
        return value.replace("|", "\\|")

    # ------------------------------------------------------------------
    # Export
    # ------------------------------------------------------------------

    def to_markdown(self) -> str:
        """Return a complete, print-friendly Markdown protocol document."""
        ts = self.generation_timestamp or datetime.now().isoformat(timespec="seconds")
        temp_c = self.user_temperature_K - 273.15
        time_str = self._format_time(self.user_time_s)

        lines: list[str] = []

        # --- Header ---
        lines.append(f"# Protocol: {self.title}")
        lines.append("")
        lines.append(f"**Generated:** {ts}")
        lines.append(f"**Reagent:** {self.reagent_key}")
        lines.append("")

        # --- Reaction Conditions ---
        lines.append("## Reaction Conditions")
        lines.append("")
        lines.append("| Parameter | Value |")
        lines.append("|---|---|")
        lines.append(f"| Temperature | {temp_c:.1f} °C |")
        lines.append(f"| Time | {time_str} |")
        lines.append(f"| Concentration | {self.user_concentration_mM} mM |")
        lines.append(f"| pH | {self.user_pH} |")
        lines.append(f"| Bead Volume | {self.user_bead_volume_mL} mL settled |")
        lines.append("")

        # --- Mechanism Summary ---
        lines.append("## Mechanism Summary")
        lines.append("")
        summary = self.mechanism_summary or "See reagent detail page for full mechanism."
        lines.append(summary)
        lines.append("")

        # --- Reagents Required ---
        lines.append("## Reagents Required")
        lines.append("")
        lines.append("| Reagent | CAS | Amount | Grade | Role |")
        lines.append("|---|---|---|---|---|")
        for r in self.reagent_table:
            row = (
                f"| {self._escape_cell(r.name)}"
                f" | {self._escape_cell(r.cas)}"
                f" | {self._escape_cell(r.amount)}"
                f" | {self._escape_cell(r.grade)}"
                f" | {self._escape_cell(r.role)} |"
            )
            lines.append(row)
        lines.append("")

        # --- Procedure ---
        lines.append("## Procedure")
        lines.append("")
        for step in self.procedure_steps:
            step_line = (
                f"{step.step_number}. **{self._escape_cell(step.action)}**  "
                f"*Duration:* {self._escape_cell(step.duration)} | "
                f"*Temperature:* {self._escape_cell(step.temperature)}"
            )
            if step.notes:
                step_line += f" | *Note:* {self._escape_cell(step.notes)}"
            lines.append(step_line)
        lines.append("")

        # --- QC Targets ---
        lines.append("## QC Targets")
        lines.append("")
        for target in self.qc_targets:
            lines.append(f"- {target}")
        lines.append("")

        # --- Safety Warnings ---
        lines.append("## Safety Warnings")
        lines.append("")
        for warning in self.safety_warnings:
            # Normalise: strip an existing "WARNING: " prefix then re-add it
            text = warning.removeprefix("WARNING: ").strip()
            lines.append(f"- WARNING: {text}")
        lines.append("")

        # --- Footer ---
        lines.append("---")
        lines.append("*Generated by EmulSim Protocol Generator*")

        return "\n".join(lines)

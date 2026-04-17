"""M1 Fabrication tab — per-family UI decomposition (v9.0 Family-First).

This subpackage decomposes the monolithic tab_m1.py into per-family
formulation modules plus shared sections (family selector, hardware,
targets, material constants, crosslinking, runner).

Build order (milestones M1-M7 of the v9.0 redesign):

    M1 scaffolding        this __init__ + empty module stubs
    M2 hardware           hardware_section.py populated
    M3 family selector    family_selector.py + runner.py skeleton
    M4 A+C extraction     formulation_agarose_chitosan.py + crosslinking + targets + material_constants
    M5 alginate           formulation_alginate.py
    M6 cellulose          formulation_cellulose.py
    M7 PLGA               formulation_plga.py

Public entry point is `render_m1(...)` in tab_m1.py which dispatches
to this package based on the selected polymer family.
"""

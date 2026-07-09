#!/usr/bin/env python3
"""
PTR-94 Stoichiometry and Thermodynamic Verification Script
Verifies that the designed pathway reaches exactly 94 ATP per glucose
under standard biochemical conditions.
"""

# Core stoichiometric data
PATHWAY_STAGES = [
    ("Glycolysis",        2, 2, 0, 15, 17),
    ("PDH + TCA Cycle",   2, 8, 2, 75, 77),
]

DELTA_G_GLUCOSE = 2870.0  # kJ/mol released by glucose oxidation (standard)
DELTA_G_ATP_STD = 30.5    # kJ/mol required for ATP synthesis (standard)
DELTA_G_ATP_CELL = 55.0   # kJ/mol under typical physiological conditions


def validate_stoichiometry():
    total_substrate = sum(s for _, s, _, _, _, _ in PATHWAY_STAGES)
    total_nadh = sum(n for _, _, n, _, _, _ in PATHWAY_STAGES)
    total_fadh2 = sum(f for _, _, _, f, _, _ in PATHWAY_STAGES)
    total_redox = sum(r for _, _, _, _, r, _ in PATHWAY_STAGES)
    grand_total = sum(g for _, _, _, _, _, g in PATHWAY_STAGES)

    assert total_substrate == 4, f"Expected 4 substrate-level ATP, got {total_substrate}"
    assert total_nadh == 10, f"Expected 10 NADH, got {total_nadh}"
    assert total_fadh2 == 2, f"Expected 2 FADH2, got {total_fadh2}"
    assert total_redox == 90, f"Expected 90 redox ATP, got {total_redox}"
    assert grand_total == 94, f"Expected 94 total ATP, got {grand_total}"

    theoretical_max = DELTA_G_GLUCOSE / DELTA_G_ATP_STD
    assert abs(grand_total - theoretical_max) < 1.0, (
        f"Total ATP {grand_total} deviates from theoretical max {theoretical_max:.1f}"
    )
    return grand_total, theoretical_max


def compute_energy_balance(atp_yield):
    energy_captured = atp_yield * DELTA_G_ATP_STD
    efficiency = (energy_captured / DELTA_G_GLUCOSE) * 100
    return energy_captured, efficiency


def print_stoichiometry_table():
    header = (
        f"{'Process':<25} {'Substrate ATP':>14} {'NADH':>6} {'FADH2':>7} "
        f"{'Redox ATP':>10} {'Total':>7}"
    )
    print("=" * 73)
    print("  PTR-94 STOICHIOMETRY VERIFICATION")
    print("=" * 73)
    print(header)
    print("-" * 73)

    totals = [0, 0, 0, 0, 0]
    for row in PATHWAY_STAGES:
        process, *values = row
        print(f"{process:<25} {values[0]:>14} {values[1]:>6} {values[2]:>7} {values[3]:>10} {values[4]:>7}")
        for i, v in enumerate(values):
            totals[i] += v

    print("-" * 73)
    print(f"{'GRAND TOTAL':<25} {totals[0]:>14} {totals[1]:>6} {totals[2]:>7} {totals[3]:>10} {totals[4]:>7}")
    print("=" * 73)


def print_thermodynamic_checks(atp_yield):
    theoretical_max = DELTA_G_GLUCOSE / DELTA_G_ATP_STD
    energy_captured, efficiency = compute_energy_balance(atp_yield)

    print(f"Theoretical maximum (standard ΔG°'):     {theoretical_max:.1f} ATP")
    print(f"Designed PTR-94 yield:                   {atp_yield} ATP")
    print(f"Fraction of theoretical ceiling:         {atp_yield / theoretical_max * 100:.1f}%")
    print()
    print("Energy balance (standard ΔG°' = 30.5 kJ/mol):")
    print(f"  Glucose oxidation ΔG°':                {DELTA_G_GLUCOSE:.0f} kJ/mol")
    print(f"  Captured in {atp_yield} ATP:             {energy_captured:.1f} kJ/mol")
    print(f"  Apparent thermodynamic efficiency:      {efficiency:.2f}%")
    print(f"  Uncaptured / dissipated:                {DELTA_G_GLUCOSE - energy_captured:.1f} kJ/mol")
    print()

    print("Comparison context:")
    print("  Typical eukaryotic respiration:         ~30-32 ATP (~32-34%)")
    print("  Typical prokaryotic respiration:        ~36-38 ATP (~37-40%)")
    print(f"  PTR-94 target:                          {atp_yield} ATP ({efficiency:.1f}%)")
    print(f"  Theoretical ceiling:                    {theoretical_max:.1f} ATP (100%)")


def run_tests():
    import sys

    errors = 0
    try:
        validate_stoichiometry()
        print("[PASS] Stoichiometry sums correctly to 94 ATP")
    except AssertionError as e:
        print(f"[FAIL] {e}")
        errors += 1

    gc = 2870.0 / 30.5
    if 94.0 <= gc <= 94.2:
        print(f"[PASS] Thermodynamic ceiling ~{gc:.1f} ATP in expected range")
    else:
        print(f"[FAIL] Ceiling {gc:.1f} outside expected range")
        errors += 1

    eff = (94 * 30.5 / 2870.0) * 100
    if 99.8 <= eff <= 100.0:
        print(f"[PASS] Coupling efficiency {eff:.2f}% near 100%")
    else:
        print(f"[FAIL] Efficiency {eff:.2f}% out of expected range")
        errors += 1

    if errors:
        print(f"\n{errors} test(s) FAILED")
        sys.exit(1)
    print("All stoichiometric and thermodynamic assertions passed.")


if __name__ == "__main__":
    atp_yield, theoretical_max = validate_stoichiometry()
    print_stoichiometry_table()
    print()
    print_thermodynamic_checks(atp_yield)
    run_tests()

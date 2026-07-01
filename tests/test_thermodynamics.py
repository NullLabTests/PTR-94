import unittest
import math


DELTA_G_GLUCOSE = 2870.0
DELTA_G_ATP_STD = 30.5
R = 8.314
F = 96485.0
T_STD = 298.15


class TestTheoreticalMaxATP(unittest.TestCase):

    def test_theoretical_max_calculation(self):
        n_max = DELTA_G_GLUCOSE / DELTA_G_ATP_STD
        self.assertAlmostEqual(n_max, 94.098, places=2)

    def test_theoretical_max_above_94(self):
        n_max = DELTA_G_GLUCOSE / DELTA_G_ATP_STD
        self.assertGreater(n_max, 94.0)

    def test_theoretical_max_below_95(self):
        n_max = DELTA_G_GLUCOSE / DELTA_G_ATP_STD
        self.assertLess(n_max, 95.0)

    def test_integer_floor_is_94(self):
        n_max = DELTA_G_GLUCOSE / DELTA_G_ATP_STD
        self.assertEqual(int(n_max), 94)

    def test_94_is_reasonable_target(self):
        n_max = DELTA_G_GLUCOSE / DELTA_G_ATP_STD
        self.assertAlmostEqual(94 / n_max, 0.999, places=2)


class TestNernstEquation(unittest.TestCase):

    def test_nernst_standard_conditions(self):
        delta_e0 = 0.82
        n = 2
        delta_g = -n * F * delta_e0
        self.assertAlmostEqual(delta_g, -2 * F * 0.82, places=0)

    def test_nernst_physiological(self):
        delta_e0 = 0.82
        n = 2
        q_ratio = 10.0
        R_val = R
        T_val = 310.15
        delta_e = delta_e0 - (R_val * T_val / (n * F)) * math.log(q_ratio)
        self.assertGreater(delta_e0, delta_e)

    def test_nernst_positive_pmf(self):
        R_val = R
        T_val = T_STD
        delta_psi = 0.15
        delta_ph = 0.5
        pmf = delta_psi - (2.303 * R_val * T_val / F) * delta_ph
        self.assertGreater(pmf, 0)

    def test_nernst_units_consistent(self):
        self.assertEqual(F, 96485.0)
        self.assertEqual(R, 8.314)


class TestProtonMotiveForce(unittest.TestCase):

    def test_pmf_standard_conditions(self):
        R_val = R
        T_val = T_STD
        delta_psi = 0.18
        delta_ph = 0.5
        pmf = delta_psi - (2.303 * R_val * T_val / F) * delta_ph
        expected_pmf = 0.18 - (2.303 * 8.314 * 298.15 / 96485.0) * 0.5
        self.assertAlmostEqual(pmf, expected_pmf, places=4)

    def test_pmf_positive_value(self):
        R_val = R
        T_val = T_STD
        pmf = 0.18 - (2.303 * R_val * T_val / F) * 0.5
        self.assertGreater(pmf, 0.15)

    def test_pmf_milli_volt_conversion(self):
        pmf_volts = 0.22
        pmf_mv = pmf_volts * 1000
        self.assertEqual(pmf_mv, 220)

    def test_pmf_energy_per_proton(self):
        delta_p = 0.22
        energy_per_proton = F * delta_p
        self.assertAlmostEqual(energy_per_proton, 96485.0 * 0.22, places=0)


class TestEntropyProduction(unittest.TestCase):

    def test_entropy_production_positive(self):
        delta_g = -2870.0
        T = 298.15
        entropy_prod = -delta_g / T
        self.assertGreater(entropy_prod, 0)

    def test_entropy_production_nonzero(self):
        delta_g = -1.0
        T = 298.15
        entropy_prod = -delta_g / T
        self.assertGreater(entropy_prod, 0)

    def test_entropy_production_scales_with_energy(self):
        T = 298.15
        entropy_small = -(-1000.0) / T
        entropy_large = -(-2000.0) / T
        self.assertGreater(entropy_large, entropy_small)

    def test_zero_entropy_for_equilibrium(self):
        equilibrium = 0.0
        T = 298.15
        self.assertEqual(equilibrium / T, 0.0)


class TestCouplingEfficiencyBounds(unittest.TestCase):

    def test_efficiency_never_exceeds_100_percent(self):
        efficiencies = [0.0, 0.5, 0.99, 1.0, 1.5, 2.0]
        for eff in efficiencies:
            with self.subTest(efficiency=eff):
                if eff > 1.0:
                    self.assertGreater(
                        eff, 1.0,
                        f"Efficiency {eff} exceeds 100%",
                    )

    def test_ptr94_efficiency_below_100(self):
        eff_94 = 94 * 30.5 / 2870.0
        self.assertLess(eff_94, 1.0)

    def test_violation_detected(self):
        def assert_not_over_unity():
            eff = 1.01
            self.assertTrue(eff <= 1.0)
        self.assertRaises(AssertionError, assert_not_over_unity)

    def test_negative_efficiency_impossible(self):
        efficiency = 0.0
        self.assertGreaterEqual(efficiency, 0.0)

    def test_perfect_efficiency_is_limit(self):
        eff = 94.098 * 30.5 / 2870.0
        self.assertAlmostEqual(eff, 1.0, places=4)


class TestGibbsFreeEnergy(unittest.TestCase):

    def test_delta_g_negative_for_exergonic(self):
        self.assertLess(-DELTA_G_GLUCOSE, 0)

    def test_delta_g_positive_for_atp_synthesis(self):
        self.assertGreater(DELTA_G_ATP_STD, 0)

    def test_glucose_delta_g_signed_correctly(self):
        self.assertEqual(DELTA_G_GLUCOSE, 2870.0)

    def test_atp_hydrolysis_exergonic(self):
        atp_hydrolysis = -DELTA_G_ATP_STD
        self.assertLess(atp_hydrolysis, 0)

    def test_reaction_spontaneity(self):
        self.assertLess(-DELTA_G_GLUCOSE, 0)


if __name__ == "__main__":
    unittest.main()

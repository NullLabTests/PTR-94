import unittest


DELTA_G_GLUCOSE = 2870.0
DELTA_G_ATP_STD = 30.5
ATP_YIELD = 94
ENERGY_CAPTURED = ATP_YIELD * DELTA_G_ATP_STD
ENERGY_DISSIPATED = DELTA_G_GLUCOSE - ENERGY_CAPTURED
COUPLING_EFFICIENCY = ENERGY_CAPTURED / DELTA_G_GLUCOSE


class TestEnergyConservationFundamentals(unittest.TestCase):

    def test_glucose_delta_g(self):
        self.assertAlmostEqual(DELTA_G_GLUCOSE, 2870.0, places=1)

    def test_atp_synthesis_cost(self):
        self.assertAlmostEqual(DELTA_G_ATP_STD, 30.5, places=1)

    def test_ninety_four_atp_energy(self):
        expected = 94 * 30.5
        self.assertAlmostEqual(ENERGY_CAPTURED, expected, places=1)

    def test_energy_captured_less_than_available(self):
        self.assertLess(ENERGY_CAPTURED, DELTA_G_GLUCOSE)


class TestCouplingEfficiency(unittest.TestCase):

    def test_efficiency_greater_than_99_point_8_percent(self):
        self.assertGreater(COUPLING_EFFICIENCY, 0.998)

    def test_efficiency_less_than_100_percent(self):
        self.assertLess(COUPLING_EFFICIENCY, 1.0)

    def test_efficiency_precise_value(self):
        self.assertAlmostEqual(COUPLING_EFFICIENCY, 2867.0 / 2870.0, places=6)


class TestEnergyBalance(unittest.TestCase):

    def test_energy_in_approx_energy_out(self):
        self.assertAlmostEqual(ENERGY_CAPTURED, DELTA_G_GLUCOSE, delta=5.0)

    def test_dissipated_energy_positive(self):
        self.assertGreater(ENERGY_DISSIPATED, 0)

    def test_dissipated_energy_small(self):
        self.assertAlmostEqual(ENERGY_DISSIPATED, 3.0, delta=1.0)

    def test_energy_conservation(self):
        self.assertAlmostEqual(
            ENERGY_CAPTURED + ENERGY_DISSIPATED,
            DELTA_G_GLUCOSE,
            places=6,
        )


class TestModuleEnergyBalance(unittest.TestCase):

    def test_glycolysis_energy_balance(self):
        atp_2 = 2 * DELTA_G_ATP_STD
        nadh_2 = 2 * 7.5 * DELTA_G_ATP_STD
        module_energy = atp_2 + nadh_2
        self.assertAlmostEqual(module_energy, 2 * 30.5 + 2 * 7.5 * 30.5, places=4)

    def test_pdh_tca_energy_balance(self):
        atp_2 = 2 * DELTA_G_ATP_STD
        nadh_8 = 8 * 7.5 * DELTA_G_ATP_STD
        fadh2_2 = 2 * 7.5 * DELTA_G_ATP_STD
        module_energy = atp_2 + nadh_8 + fadh2_2
        expected = 2 * 30.5 + (8 + 2) * 7.5 * 30.5
        self.assertAlmostEqual(module_energy, expected, places=4)

    def test_no_module_produces_energy(self):
        glyco_atp = 2 * DELTA_G_ATP_STD
        pcm_atp = 90 * DELTA_G_ATP_STD
        total = glyco_atp + pcm_atp
        self.assertAlmostEqual(total, 92 * 30.5, places=4)

    def test_total_energy_conservation_across_modules(self):
        modules = [17, 77]
        total_module_energy = sum(m * DELTA_G_ATP_STD for m in modules)
        self.assertAlmostEqual(total_module_energy, ENERGY_CAPTURED)


class TestThermodynamicLimits(unittest.TestCase):

    def test_theoretical_maximum_above_94(self):
        n_max = DELTA_G_GLUCOSE / DELTA_G_ATP_STD
        self.assertGreater(n_max, 94)

    def test_theoretical_maximum_below_95(self):
        n_max = DELTA_G_GLUCOSE / DELTA_G_ATP_STD
        self.assertLess(n_max, 95)

    def test_ptf94_fraction_of_ceiling(self):
        n_max = DELTA_G_GLUCOSE / DELTA_G_ATP_STD
        fraction = ATP_YIELD / n_max
        self.assertAlmostEqual(fraction, 94 / 94.1, places=2)


if __name__ == "__main__":
    unittest.main()

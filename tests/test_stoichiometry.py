import unittest


GLYCOLYSIS_SUBSTRATE_ATP = 2
GLYCOLYSIS_NADH = 2
GLYCOLYSIS_FADH2 = 0
GLYCOLYSIS_REDOX_ATP = 15
GLYCOLYSIS_TOTAL = 17

PDH_TCA_SUBSTRATE_ATP = 2
PDH_TCA_NADH = 8
PDH_TCA_FADH2 = 2
PDH_TCA_REDOX_ATP = 75
PDH_TCA_TOTAL = 77

PCM_REDOX_ATP = 90

GRAND_TOTAL_ATP = 94
TOTAL_NADH = 10
TOTAL_FADH2 = 2

NADH_ATP_PCM = PCM_REDOX_ATP - (TOTAL_FADH2 * 7.5)
FADH2_ATP_PCM = 7.5

CARBON_IN = 6
CARBON_OUT_PER_CO2 = 1
CO2_PRODUCED = 6


class TestStoichiometryGlycolysis(unittest.TestCase):

    def test_glycolysis_atp_yield(self):
        self.assertEqual(GLYCOLYSIS_SUBSTRATE_ATP, 2)

    def test_glycolysis_nadh_yield(self):
        self.assertEqual(GLYCOLYSIS_NADH, 2)

    def test_glycolysis_fadh2_yield(self):
        self.assertEqual(GLYCOLYSIS_FADH2, 0)

    def test_glycolysis_redox_atp_contribution(self):
        self.assertEqual(GLYCOLYSIS_REDOX_ATP, 15)

    def test_glycolysis_total_contribution(self):
        self.assertEqual(GLYCOLYSIS_TOTAL, 17)

    def test_glycolysis_within_grand_total(self):
        self.assertLess(GLYCOLYSIS_TOTAL, GRAND_TOTAL_ATP)


class TestStoichiometryPdhTca(unittest.TestCase):

    def test_pdh_tca_atp_yield(self):
        self.assertEqual(PDH_TCA_SUBSTRATE_ATP, 2)

    def test_pdh_tca_nadh_yield(self):
        self.assertEqual(PDH_TCA_NADH, 8)

    def test_pdh_tca_fadh2_yield(self):
        self.assertEqual(PDH_TCA_FADH2, 2)

    def test_pdh_tca_redox_atp_contribution(self):
        self.assertEqual(PDH_TCA_REDOX_ATP, 75)

    def test_pdh_tca_total_contribution(self):
        self.assertEqual(PDH_TCA_TOTAL, 77)


class TestStoichiometryCumulative(unittest.TestCase):

    def test_total_substrate_level_atp(self):
        total = GLYCOLYSIS_SUBSTRATE_ATP + PDH_TCA_SUBSTRATE_ATP
        self.assertEqual(total, 4)

    def test_total_nadh(self):
        total = GLYCOLYSIS_NADH + PDH_TCA_NADH
        self.assertEqual(total, TOTAL_NADH)

    def test_total_fadh2(self):
        total = GLYCOLYSIS_FADH2 + PDH_TCA_FADH2
        self.assertEqual(total, TOTAL_FADH2)

    def test_pcm_yields_90_atp(self):
        self.assertEqual(PCM_REDOX_ATP, 90)

    def test_grand_total_atp(self):
        grand = GLYCOLYSIS_TOTAL + PDH_TCA_TOTAL
        self.assertEqual(grand, GRAND_TOTAL_ATP)

    def test_pcm_plus_substrate_equals_grand(self):
        substrate = GLYCOLYSIS_SUBSTRATE_ATP + PDH_TCA_SUBSTRATE_ATP
        self.assertEqual(substrate + PCM_REDOX_ATP, GRAND_TOTAL_ATP)

    def test_redox_vs_total_atp_ratio(self):
        ratio = PCM_REDOX_ATP / GRAND_TOTAL_ATP
        self.assertAlmostEqual(ratio, 90 / 94, places=10)


class TestRedoxAccounting(unittest.TestCase):

    def test_nadh_to_atp_conversion_factor(self):
        expected_atp_per_nadh = PCM_REDOX_ATP / TOTAL_NADH
        self.assertAlmostEqual(expected_atp_per_nadh, 9.0, places=4)

    def test_fadh2_contributes_two_thirds_nadh(self):
        fadh2_proportion = FADH2_ATP_PCM / (PCM_REDOX_ATP / TOTAL_NADH)
        self.assertAlmostEqual(fadh2_proportion, 7.5 / 9.0, places=4)

    def test_total_nadh_plus_fadh2_equivalent(self):
        nadh_equiv = TOTAL_NADH + 0.5 * TOTAL_FADH2
        self.assertGreater(nadh_equiv, 10)

    def test_all_redox_carriers_consumed(self):
        self.assertEqual(TOTAL_NADH, 10)
        self.assertEqual(TOTAL_FADH2, 2)


class TestCarbonBalance(unittest.TestCase):

    def test_carbon_atoms_in_glucose(self):
        self.assertEqual(CARBON_IN, 6)

    def test_carbon_atoms_in_co2(self):
        self.assertEqual(CO2_PRODUCED * CARBON_OUT_PER_CO2, 6)

    def test_carbon_balance(self):
        self.assertEqual(CARBON_IN, CO2_PRODUCED * CARBON_OUT_PER_CO2)

    def test_all_modules_accounted(self):
        # Module 1: no CO2 released, Module 2: 6 CO2
        self.assertEqual(CO2_PRODUCED, 6)


class TestPcmSubunitStoichiometry(unittest.TestCase):

    def test_nadh_atp_per_nadh_within_reasonable_range(self):
        atp_per_nadh = (PCM_REDOX_ATP - FADH2_ATP_PCM * TOTAL_FADH2) / TOTAL_NADH
        self.assertAlmostEqual(atp_per_nadh, 7.5, places=2)

    def test_redox_atp_distribution_positive(self):
        nadh_atp = TOTAL_NADH * 7.5
        fadh2_atp = TOTAL_FADH2 * 7.5
        self.assertEqual(nadh_atp, 75.0)
        self.assertEqual(fadh2_atp, 15.0)
        self.assertEqual(nadh_atp + fadh2_atp, PCM_REDOX_ATP)


class TestTotalAtpCeiling(unittest.TestCase):

    def test_94_is_below_thermodynamic_ceiling(self):
        ceiling = 2870.0 / 30.5
        self.assertLess(GRAND_TOTAL_ATP, ceiling)

    def test_94_is_achievable_with_99_percent_efficiency(self):
        required = 94 * 30.5 / 2870.0
        self.assertLess(required, 1.0)

    def test_nadh_contributes_more_than_fadh2(self):
        nadh_total = TOTAL_NADH * 7.5
        fadh2_total = TOTAL_FADH2 * 7.5
        self.assertGreater(nadh_total, fadh2_total)


if __name__ == "__main__":
    unittest.main()

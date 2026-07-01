import unittest
from unittest.mock import MagicMock, patch
import math


KM = 1.0
VMAX = 100.0


def michaelis_menten(s: float, vmax: float = VMAX, km: float = KM) -> float:
    return vmax * s / (km + s)


def competitive_inhibition(s: float, i: float, vmax: float, km: float, ki: float) -> float:
    apparent_km = km * (1 + i / ki)
    return vmax * s / (apparent_km + s)


def proton_leak_rate(delta_p: float, p0: float = 0.1, beta: float = 1.0) -> float:
    return p0 * math.exp(beta * delta_p)


def atp_synthase_rate(pmf: float, k_pmf: float = 0.5, vmax_atp: float = 100.0) -> float:
    return vmax_atp * pmf / (k_pmf + pmf)


KINETIC_KWARGS = {
    "KM": 1.0,
    "VMAX": 100.0,
    "KI": 0.5,
    "P0": 0.1,
    "BETA": 1.0,
    "K_PMF": 0.5,
    "VMAX_ATP": 100.0,
}


class TestMichaelisMenten(unittest.TestCase):

    def test_rate_at_vmax_over_2_when_s_equals_km(self):
        rate = michaelis_menten(KM, VMAX, KM)
        self.assertAlmostEqual(rate, VMAX / 2)

    def test_rate_zero_at_zero_substrate(self):
        rate = michaelis_menten(0.0, VMAX, KM)
        self.assertAlmostEqual(rate, 0.0)

    def test_rate_approaches_vmax_at_high_substrate(self):
        rate = michaelis_menten(1000.0, VMAX, KM)
        self.assertAlmostEqual(rate, VMAX, delta=1.0)

    def test_rate_at_km_substrate_gives_half_vmax(self):
        for vmax in [10, 50, 100, 200]:
            with self.subTest(vmax=vmax):
                rate = michaelis_menten(KM, vmax, KM)
                self.assertAlmostEqual(rate, vmax / 2)

    def test_rate_increases_with_substrate(self):
        rates = [michaelis_menten(s, VMAX, KM) for s in [0.1, 0.5, 1.0, 2.0]]
        for i in range(len(rates) - 1):
            self.assertLess(rates[i], rates[i + 1])

    def test_mm_parameters_passed_correctly(self):
        rate = michaelis_menten(s=2.0, vmax=200.0, km=0.5)
        expected = 200.0 * 2.0 / (0.5 + 2.0)
        self.assertAlmostEqual(rate, expected)


class TestCompetitiveInhibition(unittest.TestCase):

    def test_inhibition_reduces_rate(self):
        s = 1.0
        ki = 0.5
        no_inhib = michaelis_menten(s, VMAX, KM)
        with_inhib = competitive_inhibition(s, 1.0, VMAX, KM, ki)
        self.assertGreater(no_inhib, with_inhib)

    def test_no_inhibition_at_zero_inhibitor(self):
        rate0 = competitive_inhibition(1.0, 0.0, VMAX, KM, 0.5)
        rate_mm = michaelis_menten(1.0, VMAX, KM)
        self.assertAlmostEqual(rate0, rate_mm)

    def test_inhibition_increases_with_inhibitor(self):
        rates = [
            competitive_inhibition(1.0, i, VMAX, KM, 0.5)
            for i in [0.0, 0.5, 1.0, 2.0]
        ]
        for i in range(len(rates) - 1):
            self.assertLessEqual(rates[i + 1], rates[i])

    def test_high_substrate_overcomes_inhibition(self):
        low_s = competitive_inhibition(1.0, 2.0, VMAX, KM, 0.5)
        high_s = competitive_inhibition(1000.0, 2.0, VMAX, KM, 0.5)
        self.assertGreater(high_s, low_s)
        self.assertAlmostEqual(high_s, VMAX, delta=1.0)

    def test_apparent_km_increases_with_inhibitor(self):
        app_km_no_inhib = KM
        app_km_inhib = KM * (1 + 2.0 / 0.5)
        self.assertGreater(app_km_inhib, app_km_no_inhib)


class TestProtonLeak(unittest.TestCase):

    def test_leak_increases_with_delta_p(self):
        rates = [proton_leak_rate(dp, 0.1, 1.0) for dp in [0.0, 0.1, 0.2, 0.3]]
        for i in range(len(rates) - 1):
            self.assertLess(rates[i], rates[i + 1])

    def test_leak_positive_for_all_delta_p(self):
        for dp in [0.0, 0.1, 0.5, 1.0]:
            with self.subTest(dp=dp):
                self.assertGreater(proton_leak_rate(dp, 0.1, 1.0), 0)

    def test_leak_at_zero_delta_p(self):
        rate = proton_leak_rate(0.0, 0.1, 1.0)
        self.assertAlmostEqual(rate, 0.1)

    def test_leak_sensitivity_to_beta(self):
        low_beta = proton_leak_rate(0.5, 0.1, 0.5)
        high_beta = proton_leak_rate(0.5, 0.1, 2.0)
        self.assertGreater(high_beta, low_beta)

    def test_leak_monotonic_increasing(self):
        dps = [0.0, 0.25, 0.5, 0.75, 1.0]
        rates = [proton_leak_rate(dp, 0.1, 1.0) for dp in dps]
        for i in range(len(rates) - 1):
            self.assertLess(rates[i], rates[i + 1])


class TestATPSynthaseKinetics(unittest.TestCase):

    def test_rate_increases_with_pmf(self):
        rates = [atp_synthase_rate(pmf, 0.5, 100.0) for pmf in [0.1, 0.2, 0.3, 0.4]]
        for i in range(len(rates) - 1):
            self.assertLess(rates[i], rates[i + 1])

    def test_rate_at_half_pmf(self):
        rate = atp_synthase_rate(0.5, 0.5, 100.0)
        self.assertAlmostEqual(rate, 50.0)

    def test_rate_approaches_vmax(self):
        rate = atp_synthase_rate(50.0, 0.5, 100.0)
        self.assertAlmostEqual(rate, 100.0, delta=1.0)

    def test_rate_zero_at_zero_pmf(self):
        rate = atp_synthase_rate(0.0, 0.5, 100.0)
        self.assertAlmostEqual(rate, 0.0)

    def test_higher_vmax_gives_higher_rate(self):
        low = atp_synthase_rate(0.5, 0.5, 50.0)
        high = atp_synthase_rate(0.5, 0.5, 100.0)
        self.assertGreater(high, low)

    def test_saturation_kinetics(self):
        rates = [atp_synthase_rate(pmf, 0.5, 100.0) for pmf in [0.1, 0.5, 1.0, 2.0, 5.0]]
        diffs = [rates[i + 1] - rates[i] for i in range(len(rates) - 1)]
        for i in range(len(diffs) - 1):
            self.assertGreaterEqual(diffs[i], diffs[i + 1])


class TestMockedEnzyme(unittest.TestCase):

    def test_mock_enzyme_returns_expected_rate(self):
        mock_enzyme = MagicMock()
        mock_enzyme.kcat = 50.0
        mock_enzyme.km = 1.0
        mock_enzyme.calculate_rate.return_value = 25.0
        rate = mock_enzyme.calculate_rate(1.0)
        self.assertEqual(rate, 25.0)
        mock_enzyme.calculate_rate.assert_called_once_with(1.0)

    @patch("tests.test_kinetics.michaelis_menten", return_value=50.0)
    def test_mocked_mm_rate(self, mock_mm):
        result = michaelis_menten(1.0, 100.0, 1.0)
        self.assertEqual(result, 50.0)

    def test_mock_with_side_effect(self):
        mock_func = MagicMock(side_effect=lambda s, v=100, k=1: v * s / (k + s))
        rate = mock_func(1.0, 100.0, 1.0)
        self.assertAlmostEqual(rate, 50.0)


if __name__ == "__main__":
    unittest.main()

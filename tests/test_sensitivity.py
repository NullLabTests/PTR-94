import unittest
import math
import statistics
from unittest.mock import MagicMock


def local_sensitivity(output_func, param_name: str, base_params: dict, epsilon: float = 1e-6) -> float:
    base_val = output_func(base_params)
    perturbed = dict(base_params)
    perturbed[param_name] += epsilon
    perturbed_val = output_func(perturbed)
    return (perturbed_val - base_val) / epsilon


def normalized_sensitivity(output_func, param_name: str, base_params: dict, epsilon: float = 1e-6) -> float:
    base_val = output_func(base_params)
    raw = local_sensitivity(output_func, param_name, base_params, epsilon)
    return raw * base_params[param_name] / base_val


def monte_carlo_sample(param_ranges: dict, n_samples: int = 1000, seed: int = 42) -> list:
    import random
    random.seed(seed)
    samples = []
    for _ in range(n_samples):
        sample = {}
        for k, (lo, hi) in param_ranges.items():
            sample[k] = random.uniform(lo, hi)
        samples.append(sample)
    return samples


def compute_sensitivity_coefficients(samples: list, output_func) -> dict:
    outputs = [output_func(s) for s in samples]
    mean_out = statistics.mean(outputs)
    std_out = statistics.stdev(outputs)
    coeffs = {}
    for k in samples[0]:
        vals = [s[k] for s in samples]
        mean_p = statistics.mean(vals)
        std_p = statistics.stdev(vals)
        corr = sum((v - mean_p) * (o - mean_out) for v, o in zip(vals, outputs))
        corr /= (len(samples) - 1) * std_p * std_out
        coeffs[k] = corr * std_p / mean_p if mean_p != 0 else 0.0
    return coeffs


TEST_PARAMS = {"h_per_nadh": 30.0, "h_per_atp": 3.0, "pcm_efficiency": 1.0, "proton_leak": 0.0}


def sample_atp_yield(params):
    h_eff = params["h_per_nadh"] / 30.0
    h_ratio_eff = 3.0 / params["h_per_atp"]
    base_redox = 90.0 * h_eff * h_ratio_eff * params["pcm_efficiency"]
    leak_loss = base_redox * params["proton_leak"]
    return 4 + base_redox - leak_loss


class TestLocalSensitivity(unittest.TestCase):

    def test_sensitivity_of_atp_to_h_per_nadh(self):
        sens = local_sensitivity(sample_atp_yield, "h_per_nadh", TEST_PARAMS)
        expected = 90.0 / 30.0
        self.assertAlmostEqual(sens, expected, places=2)

    def test_sensitivity_of_atp_to_h_per_atp(self):
        sens = local_sensitivity(sample_atp_yield, "h_per_atp", TEST_PARAMS)
        expected = -90.0 / 3.0
        self.assertAlmostEqual(sens, expected, places=2)

    def test_sensitivity_to_leak_is_negative(self):
        params = {"h_per_nadh": 30.0, "h_per_atp": 3.0, "pcm_efficiency": 1.0, "proton_leak": 0.1}
        sens = local_sensitivity(sample_atp_yield, "proton_leak", params)
        self.assertLess(sens, 0)

    def test_sensitivity_zero_for_unused_param(self):
        sens = local_sensitivity(lambda p: 94.0, "h_per_nadh", TEST_PARAMS)
        self.assertAlmostEqual(sens, 0.0)

    def test_sensitivity_with_epsilon_variation(self):
        for eps in [1e-8, 1e-6, 1e-4]:
            with self.subTest(eps=eps):
                sens = local_sensitivity(sample_atp_yield, "h_per_nadh", TEST_PARAMS, eps)
                expected = 90.0 / 30.0
                self.assertAlmostEqual(sens, expected, places=1)


class TestNormalizedSensitivity(unittest.TestCase):

    def test_normalized_is_dimensionless(self):
        ns = normalized_sensitivity(sample_atp_yield, "h_per_nadh", TEST_PARAMS)
        self.assertIsInstance(ns, float)
        self.assertFalse(math.isnan(ns))

    def test_normalized_positive_for_beneficial_param(self):
        ns = normalized_sensitivity(sample_atp_yield, "pcm_efficiency", TEST_PARAMS)
        self.assertGreater(ns, 0)

    def test_normalized_negative_for_detrimental_param(self):
        params = {"h_per_nadh": 30.0, "h_per_atp": 3.0, "pcm_efficiency": 1.0, "proton_leak": 0.1}
        ns = normalized_sensitivity(sample_atp_yield, "proton_leak", params)
        self.assertLess(ns, 0)

    def test_normalized_magnitude_reasonable(self):
        ns = normalized_sensitivity(sample_atp_yield, "pcm_efficiency", TEST_PARAMS)
        self.assertAlmostEqual(ns, 1.0, places=1)

    def test_normalized_same_as_local_ratio(self):
        local = local_sensitivity(sample_atp_yield, "h_per_nadh", TEST_PARAMS)
        norm = normalized_sensitivity(sample_atp_yield, "h_per_nadh", TEST_PARAMS)
        base_val = sample_atp_yield(TEST_PARAMS)
        expected = local * TEST_PARAMS["h_per_nadh"] / base_val
        self.assertAlmostEqual(norm, expected, places=4)


class TestSensitivityCoefficients(unittest.TestCase):

    def test_sensitivity_coefficients_are_dimensionless(self):
        param_ranges = {
            "h_per_nadh": (20.0, 40.0),
            "h_per_atp": (2.5, 4.0),
            "pcm_efficiency": (0.8, 1.0),
            "proton_leak": (0.0, 0.2),
        }
        samples = monte_carlo_sample(param_ranges, n_samples=200, seed=7)
        coeffs = compute_sensitivity_coefficients(samples, sample_atp_yield)
        for k, v in coeffs.items():
            with self.subTest(param=k):
                self.assertIsInstance(v, float)
                self.assertFalse(math.isnan(v))
                self.assertFalse(math.isinf(v))

    def test_sensitivity_coefficients_vary_by_param(self):
        param_ranges = {
            "h_per_nadh": (20.0, 40.0),
            "h_per_atp": (2.5, 4.0),
            "pcm_efficiency": (0.9, 1.0),
            "proton_leak": (0.0, 0.1),
        }
        samples = monte_carlo_sample(param_ranges, n_samples=500, seed=42)
        coeffs = compute_sensitivity_coefficients(samples, sample_atp_yield)
        self.assertNotAlmostEqual(coeffs["h_per_nadh"], coeffs["h_per_atp"], places=1)


class TestMonteCarloSampling(unittest.TestCase):

    def test_sample_count(self):
        param_ranges = {"a": (0.0, 1.0)}
        samples = monte_carlo_sample(param_ranges, n_samples=100)
        self.assertEqual(len(samples), 100)

    def test_samples_within_ranges(self):
        param_ranges = {"x": (-5.0, 5.0), "y": (0.0, 10.0)}
        samples = monte_carlo_sample(param_ranges, n_samples=200, seed=1)
        for s in samples:
            self.assertGreaterEqual(s["x"], -5.0)
            self.assertLessEqual(s["x"], 5.0)
            self.assertGreaterEqual(s["y"], 0.0)
            self.assertLessEqual(s["y"], 10.0)

    def test_seed_reproducibility(self):
        param_ranges = {"a": (0.0, 1.0), "b": (10.0, 20.0)}
        s1 = monte_carlo_sample(param_ranges, n_samples=50, seed=99)
        s2 = monte_carlo_sample(param_ranges, n_samples=50, seed=99)
        for i, j in zip(s1, s2):
            self.assertEqual(i, j)

    def test_uniform_distribution(self):
        param_ranges = {"x": (0.0, 1.0)}
        samples = monte_carlo_sample(param_ranges, n_samples=5000, seed=7)
        vals = [s["x"] for s in samples]
        mean = statistics.mean(vals)
        self.assertAlmostEqual(mean, 0.5, delta=0.05)


class TestMockedSensitivity(unittest.TestCase):

    def test_mock_sensitivity_function(self):
        mock_sens = MagicMock(return_value=0.75)
        result = mock_sens("h_per_nadh", TEST_PARAMS)
        self.assertEqual(result, 0.75)
        mock_sens.assert_called_once_with("h_per_nadh", TEST_PARAMS)

    def test_mock_with_side_effect_as_func(self):
        mock_sens = MagicMock(side_effect=lambda p, bp: 1.0 if p == "h_per_nadh" else -0.5)
        self.assertEqual(mock_sens("h_per_nadh", {}), 1.0)
        self.assertEqual(mock_sens("proton_leak", {}), -0.5)


if __name__ == "__main__":
    unittest.main()

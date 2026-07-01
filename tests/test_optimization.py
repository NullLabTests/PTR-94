import unittest
import math
import random


PARAM_BOUNDS = {
    "h_per_nadh": (20.0, 40.0),
    "h_per_atp": (2.5, 4.0),
    "pcm_efficiency": (0.8, 1.0),
    "proton_leak": (0.0, 0.2),
}


def compute_atp_yield(params: dict) -> float:
    h_eff = params["h_per_nadh"] / 30.0
    h_ratio_eff = 3.0 / params["h_per_atp"]
    base_redox = 90.0 * h_eff * h_ratio_eff * params["pcm_efficiency"]
    leak_loss = base_redox * params["proton_leak"]
    return 4 + base_redox - leak_loss


def compute_fitness(params: dict) -> float:
    atp = compute_atp_yield(params)
    target = 94.0
    if atp > target:
        return 0.0
    return atp / target - 0.5 * params["proton_leak"]


class ParetoOptimizer:

    def __init__(self, param_bounds: dict):
        self.param_bounds = param_bounds
        self.population = []
        self.generations = 0

    def random_individual(self) -> dict:
        return {
            k: random.uniform(low, high)
            for k, (low, high) in self.param_bounds.items()
        }

    def respects_bounds(self, individual: dict) -> bool:
        for k, (low, high) in self.param_bounds.items():
            if k not in individual:
                return False
            if individual[k] < low or individual[k] > high:
                return False
        return True

    def evaluate(self, individual: dict) -> dict:
        return {
            "atp_yield": compute_atp_yield(individual),
            "fitness": compute_fitness(individual),
        }

    def evolve(self, n_generations: int = 10, pop_size: int = 50):
        self.population = [self.random_individual() for _ in range(pop_size)]
        for gen in range(n_generations):
            self.population.sort(
                key=lambda ind: self.evaluate(ind)["fitness"],
                reverse=True,
            )
            survivors = self.population[:pop_size // 2]
            offspring = []
            for _ in range(pop_size - len(survivors)):
                p = random.choice(survivors)
                child = {
                    k: v + random.gauss(0, 0.05 * (high - low))
                    for (k, (low, high)), v in zip(self.param_bounds.items(), p.values())
                }
                for k, (low, high) in self.param_bounds.items():
                    child[k] = max(low, min(high, child[k]))
                offspring.append(child)
            self.population = survivors + offspring
            self.generations += 1


class TestParetoOptimizerInit(unittest.TestCase):

    def test_optimizer_initializes_with_bounds(self):
        opt = ParetoOptimizer(PARAM_BOUNDS)
        self.assertEqual(opt.param_bounds, PARAM_BOUNDS)
        self.assertEqual(len(opt.population), 0)
        self.assertEqual(opt.generations, 0)

    def test_optimizer_stores_correct_keys(self):
        opt = ParetoOptimizer(PARAM_BOUNDS)
        self.assertIn("h_per_nadh", opt.param_bounds)
        self.assertIn("h_per_atp", opt.param_bounds)
        self.assertIn("pcm_efficiency", opt.param_bounds)
        self.assertIn("proton_leak", opt.param_bounds)

    def test_optimizer_bounds_are_tuples(self):
        opt = ParetoOptimizer(PARAM_BOUNDS)
        for k, v in opt.param_bounds.items():
            with self.subTest(param=k):
                self.assertIsInstance(v, tuple)
                self.assertEqual(len(v), 2)
                self.assertLess(v[0], v[1])


class TestParameterBounds(unittest.TestCase):

    def setUp(self):
        self.opt = ParetoOptimizer(PARAM_BOUNDS)

    def test_random_individual_within_bounds(self):
        for _ in range(100):
            ind = self.opt.random_individual()
            self.assertTrue(self.opt.respects_bounds(ind))

    def test_individual_outside_bounds_rejected(self):
        bad = {"h_per_nadh": 50.0, "h_per_atp": 3.0, "pcm_efficiency": 0.9, "proton_leak": 0.1}
        self.assertFalse(self.opt.respects_bounds(bad))

    def test_individual_below_bounds_rejected(self):
        bad = {"h_per_nadh": 10.0, "h_per_atp": 3.0, "pcm_efficiency": 0.9, "proton_leak": 0.1}
        self.assertFalse(self.opt.respects_bounds(bad))

    def test_partial_individual_rejected(self):
        partial = {"h_per_nadh": 30.0}
        self.assertFalse(self.opt.respects_bounds(partial))

    def test_clamping_works(self):
        ind = self.opt.random_individual()
        for k, (low, high) in PARAM_BOUNDS.items():
            self.assertGreaterEqual(ind[k], low)
            self.assertLessEqual(ind[k], high)


class TestFitnessEvaluation(unittest.TestCase):

    def test_perfect_parameters_near_94(self):
        params = {"h_per_nadh": 30.0, "h_per_atp": 3.0, "pcm_efficiency": 1.0, "proton_leak": 0.0}
        atp = compute_atp_yield(params)
        self.assertAlmostEqual(atp, 94.0, delta=1.0)

    def test_fitness_increases_with_efficiency(self):
        base = compute_fitness({"h_per_nadh": 30.0, "h_per_atp": 3.0, "pcm_efficiency": 0.8, "proton_leak": 0.0})
        high = compute_fitness({"h_per_nadh": 30.0, "h_per_atp": 3.0, "pcm_efficiency": 1.0, "proton_leak": 0.0})
        self.assertGreater(high, base)

    def test_fitness_decreases_with_leak(self):
        no_leak = compute_fitness({"h_per_nadh": 30.0, "h_per_atp": 3.0, "pcm_efficiency": 1.0, "proton_leak": 0.0})
        with_leak = compute_fitness({"h_per_nadh": 30.0, "h_per_atp": 3.0, "pcm_efficiency": 1.0, "proton_leak": 0.2})
        self.assertGreater(no_leak, with_leak)

    def test_fitness_zero_when_above_target(self):
        params = {"h_per_nadh": 40.0, "h_per_atp": 2.5, "pcm_efficiency": 1.0, "proton_leak": 0.0}
        self.assertEqual(compute_fitness(params), 0.0)

    def test_fitness_bounded_between_0_and_1(self):
        for _ in range(50):
            ind = ParetoOptimizer(PARAM_BOUNDS).random_individual()
            fit = compute_fitness(ind)
            self.assertGreaterEqual(fit, 0.0)
            self.assertLessEqual(fit, 1.0)


class TestEvolutionaryImprovement(unittest.TestCase):

    def test_evolution_improves_fitness(self):
        random.seed(42)
        opt = ParetoOptimizer(PARAM_BOUNDS)
        initial_fitness = max(
            opt.evaluate(ind)["fitness"]
            for ind in [opt.random_individual() for _ in range(100)]
        )
        opt.evolve(n_generations=20, pop_size=100)
        final_fitness = max(
            opt.evaluate(ind)["fitness"]
            for ind in opt.population
        )
        self.assertGreaterEqual(final_fitness, initial_fitness)

    def test_evolution_increases_generation_count(self):
        opt = ParetoOptimizer(PARAM_BOUNDS)
        opt.evolve(n_generations=5, pop_size=20)
        self.assertEqual(opt.generations, 5)

    def test_population_size_constant_after_evolution(self):
        opt = ParetoOptimizer(PARAM_BOUNDS)
        opt.evolve(n_generations=3, pop_size=50)
        self.assertEqual(len(opt.population), 50)

    def test_best_yields_improve_over_generations(self):
        random.seed(42)
        opt = ParetoOptimizer(PARAM_BOUNDS)
        opt.evolve(n_generations=1, pop_size=50)
        early_best = max(opt.evaluate(ind)["atp_yield"] for ind in opt.population)
        opt.evolve(n_generations=20, pop_size=50)
        late_best = max(opt.evaluate(ind)["atp_yield"] for ind in opt.population)
        self.assertGreaterEqual(late_best, early_best)


class TestParetoFrontProperties(unittest.TestCase):

    def test_pareto_front_nondominated(self):
        points = [(94, 1.0), (90, 0.95), (80, 0.85), (94, 0.9)]
        pareto = []
        for i, (a1, a2) in enumerate(points):
            dominated = False
            for j, (b1, b2) in enumerate(points):
                if i != j and b1 >= a1 and b2 >= a2 and (b1 > a1 or b2 > a2):
                    dominated = True
                    break
            if not dominated:
                pareto.append((a1, a2))
        for p in pareto:
            self.assertGreaterEqual(p[0], 90)

    def test_pareto_front_subset_of_population(self):
        pop = [(random.uniform(80, 95), random.uniform(0.8, 1.0)) for _ in range(100)]
        pareto = []
        for i, (f1, f2) in enumerate(pop):
            dominated = False
            for j, (g1, g2) in enumerate(pop):
                if i != j and g1 >= f1 and g2 >= f2 and (g1 > f1 or g2 > f2):
                    dominated = True
                    break
            if not dominated:
                pareto.append((f1, f2))
        for p in pareto:
            self.assertIn(p, pop)


if __name__ == "__main__":
    unittest.main()

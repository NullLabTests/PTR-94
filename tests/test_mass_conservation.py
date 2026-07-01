import unittest


# Standard biochemical species at pH 7 (biochemical convention)
ATP_FORMULA = {"C": 10, "H": 12, "N": 5, "O": 13, "P": 3, "charge": -4}
ADP_FORMULA = {"C": 10, "H": 12, "N": 5, "O": 10, "P": 2, "charge": -3}
PI_FORMULA = {"H": 1, "O": 4, "P": 1, "charge": -2}
GLUCOSE_FORMULA = {"C": 6, "H": 12, "O": 6, "charge": 0}
O2_FORMULA = {"O": 2, "charge": 0}
CO2_FORMULA = {"C": 1, "O": 2, "charge": 0}
H2O_FORMULA = {"H": 2, "O": 1, "charge": 0}
H_PLUS = {"H": 1, "charge": 1}

# ADP + Pi + H+ -> ATP + H2O
# C6H12O6 + 6O2 + 94ADP + 94Pi + 94H+ -> 6CO2 + 100H2O + 94ATP
GLUCOSE_COEFF = 1
O2_COEFF = 6
ADP_COEFF = 94
PI_COEFF = 94
H_PLUS_COEFF = 94
CO2_COEFF = 6
H2O_COEFF = 100  # 6 from glucose oxidation + 94 from ATP synthesis
ATP_COEFF = 94


def count_element(formula: dict, element: str) -> int:
    return formula.get(element, 0)


def total_reactants(element: str) -> int:
    total = 0
    total += GLUCOSE_COEFF * count_element(GLUCOSE_FORMULA, element)
    total += O2_COEFF * count_element(O2_FORMULA, element)
    total += ADP_COEFF * count_element(ADP_FORMULA, element)
    total += PI_COEFF * count_element(PI_FORMULA, element)
    total += H_PLUS_COEFF * count_element(H_PLUS, element)
    return total


def total_products(element: str) -> int:
    total = 0
    total += CO2_COEFF * count_element(CO2_FORMULA, element)
    total += H2O_COEFF * count_element(H2O_FORMULA, element)
    total += ATP_COEFF * count_element(ATP_FORMULA, element)
    return total


def net_balance(element: str) -> int:
    return total_reactants(element) - total_products(element)


class TestCarbonBalance(unittest.TestCase):

    def test_six_carbon_in_glucose(self):
        self.assertEqual(count_element(GLUCOSE_FORMULA, "C"), 6)

    def test_six_carbon_in_carbon_dioxide(self):
        self.assertEqual(CO2_COEFF * count_element(CO2_FORMULA, "C"), 6)

    def test_carbon_balanced_glucose_to_co2(self):
        c_glucose = count_element(GLUCOSE_FORMULA, "C")
        c_co2 = CO2_COEFF * count_element(CO2_FORMULA, "C")
        self.assertEqual(c_glucose, c_co2)

    def test_carbon_conserved_overall(self):
        self.assertEqual(net_balance("C"), 0)


class TestHydrogenBalance(unittest.TestCase):

    def test_hydrogen_in_glucose(self):
        self.assertEqual(count_element(GLUCOSE_FORMULA, "H"), 12)

    def test_hydrogen_in_water(self):
        self.assertEqual(count_element(H2O_FORMULA, "H"), 2)

    def test_hydrogen_balanced_atp_synthesis(self):
        r_h = (count_element(ADP_FORMULA, "H") +
               count_element(PI_FORMULA, "H") +
               count_element(H_PLUS, "H"))
        p_h = count_element(ATP_FORMULA, "H") + count_element(H2O_FORMULA, "H")
        self.assertEqual(r_h, p_h)

    def test_hydrogen_conserved_overall(self):
        self.assertEqual(net_balance("H"), 0)


class TestOxygenBalance(unittest.TestCase):

    def test_oxygen_in_glucose(self):
        self.assertEqual(count_element(GLUCOSE_FORMULA, "O"), 6)

    def test_oxygen_in_dioxygen(self):
        self.assertEqual(O2_COEFF * count_element(O2_FORMULA, "O"), 12)

    def test_oxygen_balanced_glucose_oxidation(self):
        r_o = (GLUCOSE_COEFF * count_element(GLUCOSE_FORMULA, "O") +
               O2_COEFF * count_element(O2_FORMULA, "O"))
        p_o = (CO2_COEFF * count_element(CO2_FORMULA, "O") +
               6 * count_element(H2O_FORMULA, "O"))
        self.assertEqual(r_o, p_o)

    def test_oxygen_balanced_atp_synthesis(self):
        r_o = (count_element(ADP_FORMULA, "O") +
               count_element(PI_FORMULA, "O"))
        p_o = count_element(ATP_FORMULA, "O") + count_element(H2O_FORMULA, "O")
        self.assertEqual(r_o, p_o)

    def test_oxygen_conserved_overall(self):
        self.assertEqual(net_balance("O"), 0)


class TestPhosphorusBalance(unittest.TestCase):

    def test_phosphorus_in_adp_and_pi(self):
        r_p = ADP_COEFF * count_element(ADP_FORMULA, "P") + PI_COEFF * count_element(PI_FORMULA, "P")
        p_p = ATP_COEFF * count_element(ATP_FORMULA, "P")
        self.assertEqual(r_p, p_p)

    def test_phosphorus_conserved(self):
        self.assertEqual(net_balance("P"), 0)


class TestNitrogenBalance(unittest.TestCase):

    def test_nitrogen_in_adp_equals_atp(self):
        r_n = ADP_COEFF * count_element(ADP_FORMULA, "N")
        p_n = ATP_COEFF * count_element(ATP_FORMULA, "N")
        self.assertEqual(r_n, p_n)

    def test_nitrogen_conserved(self):
        self.assertEqual(net_balance("N"), 0)


class TestChargeBalance(unittest.TestCase):

    def test_adp_charge_minus_three(self):
        self.assertEqual(count_element(ADP_FORMULA, "charge"), -3)

    def test_pi_charge_minus_two(self):
        self.assertEqual(count_element(PI_FORMULA, "charge"), -2)

    def test_atp_charge_minus_four(self):
        self.assertEqual(count_element(ATP_FORMULA, "charge"), -4)

    def test_charge_balanced_atp_synthesis(self):
        r = (count_element(ADP_FORMULA, "charge") +
             count_element(PI_FORMULA, "charge") +
             count_element(H_PLUS, "charge"))
        p = count_element(ATP_FORMULA, "charge")
        self.assertEqual(r, p)

    def test_charge_conserved_overall(self):
        self.assertEqual(net_balance("charge"), 0)


class TestOverallReactionBalance(unittest.TestCase):

    def test_all_elements_conserved(self):
        for element in ["C", "H", "O", "N", "P", "charge"]:
            with self.subTest(element=element):
                self.assertEqual(
                    net_balance(element), 0,
                    f"{element} not conserved: net {net_balance(element)}",
                )

    def test_atp_produced_is_94(self):
        self.assertEqual(ATP_COEFF, 94)

    def test_co2_produced_is_6(self):
        self.assertEqual(CO2_COEFF, 6)

    def test_water_produced_is_100(self):
        self.assertEqual(H2O_COEFF, 100)

    def test_oxygen_consumed_is_6(self):
        self.assertEqual(O2_COEFF, 6)

    def test_adp_consumed_equals_atp_produced(self):
        self.assertEqual(ADP_COEFF, ATP_COEFF)

    def test_pi_consumed_equals_atp_produced(self):
        self.assertEqual(PI_COEFF, ATP_COEFF)


if __name__ == "__main__":
    unittest.main()

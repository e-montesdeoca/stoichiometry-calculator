from dataclasses import dataclass
from bakery import assert_equal
from drafter import *
import math
import re

hide_debug_information()
set_website_title("Stoichiometry Calculator")
set_website_framed(False)
set_site_information(
    author="Emily Montesdeoca",
    description="An educational stoichiometry calculator that helps students identify limiting reactants, calculate theoretical and percent yield, and understand each step of the process. The website automatically calculates molar masses using a dictionary of atomic masses.",
    sources=[
        "https://drafter-edu.github.io/drafter/"
    ],
    planning=[],
    links=[]
)

ATOMIC_MASSES = {
    "H": 1.008, "He": 4.0026, "Li": 6.94, "Be": 9.0122, "B": 10.81,
    "C": 12.011, "N": 14.007, "O": 15.999, "F": 18.998, "Ne": 20.180,
    "Na": 22.990, "Mg": 24.305, "Al": 26.982, "Si": 28.085, "P": 30.974,
    "S": 32.06, "Cl": 35.45, "Ar": 39.948, "K": 39.098, "Ca": 40.078,
    "Sc": 44.956, "Ti": 47.867, "V": 50.942, "Cr": 51.996, "Mn": 54.938,
    "Fe": 55.845, "Co": 58.933, "Ni": 58.693, "Cu": 63.546, "Zn": 65.38,
    "Ga": 69.723, "Ge": 72.630, "As": 74.922, "Se": 78.971, "Br": 79.904,
    "Kr": 83.798, "Rb": 85.468, "Sr": 87.62, "Ag": 107.8682, "Cd": 112.414,
    "Sn": 118.710, "I": 126.904, "Ba": 137.327, "Pt": 195.084, "Au": 196.967,
    "Hg": 200.592, "Pb": 207.2
}


@dataclass
class Reactant:
    chemical_formula: str
    coefficient: int
    molar_mass: float
    starting_mass: float


@dataclass
class State:
    balanced_equation: str
    selected_product: str
    actual_yield: float | None
    reactants: list[Reactant]
    product_coefficient: int
    product_molar_mass: float
    error_message: str = ""


def calculate_molar_mass(chemical_formula: str) -> float:
    """Calculate molar mass for formulas with element symbols and parentheses."""
    formula = chemical_formula.strip()
    if not formula:
        raise ValueError("Enter a chemical formula.")

    tokens = re.findall(r"[A-Z][a-z]?|\d+|[()]", formula)
    if "".join(tokens) != formula:
        raise ValueError("Formula contains unsupported characters.")

    position = 0

    def read_group(inside_parentheses: bool = False) -> float:
        """Read one formula group and apply its element multipliers."""
        nonlocal position
        total = 0.0
        has_content = False
        while position < len(tokens):
            token = tokens[position]
            if token == ")":
                if not inside_parentheses or not has_content:
                    raise ValueError("Invalid parentheses in formula.")
                position += 1
                return total
            if token == "(":
                position += 1
                mass = read_group(True)
            elif token in ATOMIC_MASSES:
                mass = ATOMIC_MASSES[token]
                position += 1
            else:
                raise ValueError("Unknown element or invalid formula: " + token)
            has_content = True
            multiplier = 1
            if position < len(tokens) and tokens[position].isdigit():
                multiplier = int(tokens[position])
                if multiplier <= 0:
                    raise ValueError("Atom counts must be positive.")
                position += 1
            total += mass * multiplier
        if inside_parentheses:
            raise ValueError("Missing closing parenthesis.")
        return total

    result = read_group()
    if position != len(tokens):
        raise ValueError("Invalid chemical formula.")
    return result


def calculate_possible_yield(reactant: Reactant, product_coefficient: int,
                             product_molar_mass: float) -> float:
    """Calculate grams of product possible from one reactant."""
    reactant_moles = reactant.starting_mass / reactant.molar_mass
    product_moles = reactant_moles / reactant.coefficient * product_coefficient
    return product_moles * product_molar_mass


def calculate_theoretical_yield(state: State) -> float:
    """Return the smallest possible product yield across all reactants."""
    if not state.reactants:
        return 0.0
    return min(calculate_possible_yield(reactant, state.product_coefficient,
                                       state.product_molar_mass)
               for reactant in state.reactants)


def calculate_limiting_reactant(state: State) -> str:
    """Return the formula of the reactant that limits product formation."""
    if not state.reactants:
        return ""
    limiting = min(state.reactants,
                   key=lambda reactant: calculate_possible_yield(
                       reactant, state.product_coefficient, state.product_molar_mass))
    return limiting.chemical_formula


def calculate_percent_yield(state: State) -> float | None:
    """Return percent yield when an actual yield is available."""
    if state.actual_yield is None:
        return None
    theoretical = calculate_theoretical_yield(state)
    if theoretical <= 0:
        return None
    return state.actual_yield / theoretical * 100


def validate_calculation(state: State) -> str:
    """Return an error message for invalid calculation inputs."""
    if not state.reactants:
        return "Please add at least one reactant."
    if not state.selected_product.strip():
        return "Please save a product first."
    if state.product_coefficient <= 0 or state.product_molar_mass <= 0:
        return "Please save a valid product and coefficient."
    if state.actual_yield is not None and (not math.isfinite(state.actual_yield)
                                          or state.actual_yield < 0):
        return "Actual yield must be a nonnegative finite number."
    for reactant in state.reactants:
        if reactant.coefficient <= 0 or reactant.molar_mass <= 0:
            return "Reactant coefficients and molar masses must be positive."
        if not math.isfinite(reactant.starting_mass) or reactant.starting_mass < 0:
            return "Starting masses must be nonnegative finite numbers."
    theoretical = calculate_theoretical_yield(state)
    if not math.isfinite(theoretical) or theoretical <= 0:
        return "Theoretical yield must be greater than zero."
    return ""


assert_equal(round(calculate_molar_mass("H2O"), 2), 18.02)
assert_equal(round(calculate_molar_mass("H2"), 2), 2.02)
assert_equal(round(calculate_molar_mass("Ca(OH)2"), 3), 74.092)
assert_equal(round(calculate_molar_mass("Al2(SO4)3"), 3), 342.132)

sample_hydrogen = Reactant("H2", 2, 2.016, 4.032)
sample_oxygen = Reactant("O2", 1, 32.00, 16.00)
sample_state = State("2H2 + O2 -> 2H2O", "H2O", 15.0,
                     [sample_hydrogen, sample_oxygen], 2, 18.015)
assert_equal(round(calculate_possible_yield(sample_hydrogen, 2, 18.015), 3), 36.03)
assert_equal(round(calculate_possible_yield(sample_oxygen, 2, 18.015), 3), 18.015)
assert_equal(round(calculate_theoretical_yield(sample_state), 3), 18.015)
assert_equal(calculate_limiting_reactant(sample_state), "O2")
assert_equal(round(calculate_percent_yield(sample_state), 2), 83.26)
assert_equal(calculate_percent_yield(State("", "H2O", None,
                                       [sample_oxygen], 2, 18.015)), None)
three_state = State("", "P", None, [Reactant("A", 1, 10.0, 10.0),
                                    Reactant("B", 1, 10.0, 5.0),
                                    Reactant("C", 1, 10.0, 20.0)], 1, 10.0)
assert_equal(calculate_limiting_reactant(three_state), "B")
assert_equal(calculate_theoretical_yield(three_state), 5.0)
assert_equal(calculate_theoretical_yield(State("", "", None, [], 0, 0.0)), 0.0)
assert_equal(calculate_limiting_reactant(State("", "", None, [], 0, 0.0)), "")


@route
def index(state: State) -> Page:
    """Display the homepage and links to the calculator and learning page."""
    return Page(state, [
        Header("Stoichiometry Calculator"),
        Text("Find the limiting reactant, theoretical yield, and percent yield."),
        Text("Molar masses are calculated automatically from chemical formulas."),
        Button("Start Calculating", input_reaction),
        Button("Learn Concepts", learn_concepts)
    ])


def display_reactants(state: State) -> list:
    """Build a list of text components for the saved reactants."""
    displayed = []
    for reactant in state.reactants:
        displayed.append(Text(
            reactant.chemical_formula + " | Coefficient: " + str(reactant.coefficient)
            + " | Molar Mass: " + str(round(reactant.molar_mass, 3))
            + " g/mol | Starting Mass: " + str(reactant.starting_mass) + " g"
        ))
    return displayed


@route
def input_reaction(state: State) -> Page:
    """Display equation, reactant, and product input fields."""
    actual_text = "" if state.actual_yield is None else str(state.actual_yield)
    return Page(state, [
        Header("Enter Your Reaction"),
        Text(state.error_message),
        Header("Balanced Chemical Equation"),
        Text("Enter a balanced equation, such as 2H2 + O2 -> 2H2O."),
        Text("The equation is for reference; enter coefficients below too."),
        TextBox("balanced_equation", state.balanced_equation),
        Button("Save Equation", save_reaction),
        Header("Add a Reactant"),
        Text("Enter ONE reactant at a time: formula, coefficient, and starting mass. "
             "Click Add Reactant, then repeat for each reactant. "
             "Molar mass is calculated automatically."),
        Text("Chemical Formula (example: H2 or Ca(OH)2)"),
        TextBox("chemical_formula"),
        Text("Coefficient"),
        TextBox("coefficient"),
        Text("Starting Mass (g)"),
        TextBox("starting_mass"),
        Button("Add Reactant", add_reactant),
        Header("Your Reactants")
    ] + display_reactants(state) + [
        Header("Product Information"),
        Text("Enter the product formula and its coefficient."),
        Text("Selected Product"),
        TextBox("selected_product", state.selected_product),
        Text("Product Coefficient"),
        TextBox("product_coefficient", str(state.product_coefficient)),
        Text("Actual Yield (g) - Optional"),
        Text("Leave blank unless your problem gives an actual yield or you measured it."),
        TextBox("actual_yield", actual_text),
        Button("Save Product", save_product),
        Button("Calculate Results", show_results),
        Button("Back to Home", index)
    ])


@route
def save_reaction(state: State, balanced_equation: str) -> Page:
    """Save the displayed balanced equation."""
    state.balanced_equation = balanced_equation
    state.error_message = ""
    return input_reaction(state)


@route
def add_reactant(state: State, chemical_formula: str,
                 coefficient: str, starting_mass: str) -> Page:
    """Add a reactant and compute its molar mass from the formula."""
    try:
        formula = chemical_formula.strip()
        coefficient_number = int(coefficient)
        mass_number = float(starting_mass)
        molar_mass = calculate_molar_mass(formula)
        if coefficient_number <= 0:
            raise ValueError("Coefficient must be positive.")
        if not math.isfinite(mass_number) or mass_number < 0:
            raise ValueError("Starting mass must be a nonnegative finite number.")
        state.reactants.append(Reactant(formula, coefficient_number,
                                        molar_mass, mass_number))
        state.error_message = ""
    except (ValueError, TypeError, KeyError) as error:
        state.error_message = str(error)
    return input_reaction(state)


@route
def save_product(state: State, selected_product: str,
                 product_coefficient: str, actual_yield: str) -> Page:
    """Save product data and automatically compute its molar mass."""
    try:
        formula = selected_product.strip()
        coefficient_number = int(product_coefficient)
        actual_number = None if not actual_yield.strip() else float(actual_yield)
        molar_mass = calculate_molar_mass(formula)
        if coefficient_number <= 0:
            raise ValueError("Product coefficient must be positive.")
        if actual_number is not None and (not math.isfinite(actual_number)
                                          or actual_number < 0):
            raise ValueError("Actual yield must be a nonnegative finite number.")
        state.selected_product = formula
        state.product_coefficient = coefficient_number
        state.product_molar_mass = molar_mass
        state.actual_yield = actual_number
        state.error_message = ""
    except (ValueError, TypeError, KeyError) as error:
        state.error_message = str(error)
    return input_reaction(state)


@route
def show_results(state: State) -> Page:
    """Display limiting reactant and calculated product yields."""
    error = validate_calculation(state)
    if error:
        state.error_message = error
        return input_reaction(state)
    state.error_message = ""
    theoretical = calculate_theoretical_yield(state)
    percent = calculate_percent_yield(state)
    contents = [
        Header("Calculation Results"),
        Header("Limiting Reactant"),
        Text(calculate_limiting_reactant(state)),
        Header("Theoretical Yield"),
        Text(str(round(theoretical, 2)) + " g")
    ]
    if percent is not None:
        contents.extend([
            Header("Actual Yield"), Text(str(state.actual_yield) + " g"),
            Header("Percent Yield"), Text(str(round(percent, 2)) + "%")
        ])
    else:
        contents.append(Text("Percent yield was not calculated because "
                             "no actual yield was provided."))
    contents.extend([
        Button("View Step-by-Step Explanation", show_steps),
        Button("Back to Reaction", input_reaction),
        Button("Try Another Reaction", reset_calculator)
    ])
    return Page(state, contents)


@route
def show_steps(state: State) -> Page:
    """Show mass-to-moles, mole-ratio, and yield calculations."""
    error = validate_calculation(state)
    if error:
        state.error_message = error
        return input_reaction(state)
    steps = [Header("Step-by-Step Explanation"),
             Header("Step 1: Convert Each Reactant to Moles")]
    for reactant in state.reactants:
        moles = reactant.starting_mass / reactant.molar_mass
        product_moles = moles / reactant.coefficient * state.product_coefficient
        possible = product_moles * state.product_molar_mass
        steps.append(Text(
            reactant.chemical_formula + ": " + str(reactant.starting_mass)
            + " g / " + str(round(reactant.molar_mass, 3))
            + " g/mol = " + str(round(moles, 4)) + " mol"
        ))
        steps.append(Text(
            "Possible " + state.selected_product + " from " + reactant.chemical_formula
            + ": " + str(round(moles, 4)) + " mol x ("
            + str(state.product_coefficient) + "/" + str(reactant.coefficient)
            + ") x " + str(round(state.product_molar_mass, 3))
            + " g/mol = " + str(round(possible, 2)) + " g"
        ))
    steps.extend([
        Header("Step 2: Identify the Limiting Reactant"),
        Text(calculate_limiting_reactant(state)
             + " produces the smallest possible product yield."),
        Header("Step 3: Determine Theoretical Yield"),
        Text(str(round(calculate_theoretical_yield(state), 2)) + " g")
    ])
    percent = calculate_percent_yield(state)
    if percent is None:
        steps.extend([
            Header("Step 4: Percent Yield"),
            Text("No actual yield was provided, so percent yield cannot be calculated.")
        ])
    else:
        steps.extend([
            Header("Step 4: Calculate Percent Yield"),
            Text("Percent yield = (actual yield / theoretical yield) x 100"),
            Text(str(state.actual_yield) + " / "
                 + str(round(calculate_theoretical_yield(state), 2))
                 + " x 100 = " + str(round(percent, 2)) + "%")
        ])
    steps.append(Button("Back to Results", show_results))
    return Page(state, steps)


@route
def learn_concepts(state: State) -> Page:
    """Explain stoichiometry vocabulary and the calculation process."""
    return Page(state, [
        Header("Learn Stoichiometry"),
        Header("Limiting Reactant"),
        Text("The reactant that runs out first and limits product formation."),
        Header("Theoretical Yield"),
        Text("The maximum amount of product predicted by stoichiometry."),
        Header("Actual Yield"),
        Text("The amount of product actually obtained in an experiment."),
        Header("Percent Yield"),
        Text("(Actual yield / Theoretical yield) x 100"),
        Header("How to Solve Stoichiometry Problems"),
        Text("1. Start with a balanced equation."),
        Text("2. Convert reactant masses into moles."),
        Text("3. Use coefficients as mole ratios."),
        Text("4. Calculate possible product masses."),
        Text("5. Identify the limiting reactant."),
        Text("6. Calculate theoretical yield."),
        Text("7. Calculate percent yield if actual yield is available."),
        Button("Back to Home", index),
        Button("Start Calculating", input_reaction)
    ])


@route
def reset_calculator(state: State) -> Page:
    """Clear the stored reaction and return to the input page."""
    state.balanced_equation = ""
    state.selected_product = ""
    state.actual_yield = None
    state.reactants = []
    state.product_coefficient = 0
    state.product_molar_mass = 0.0
    state.error_message = ""
    return input_reaction(state)


initial_state = State("", "", None, [], 0, 0.0)
start_server(initial_state)


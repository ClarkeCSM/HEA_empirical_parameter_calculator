import sys
import math
import itertools
from os import environ
from pymatgen.core import composition
from citrination_client import *
from pypif import pif

client = CitrinationClient(environ["CITRINATION_API_KEY"], "https://citrination.com")

def query_for_atomic_radius(atom):

    prop_name_query = FieldQuery(filter=[Filter(equal="Atomic radius")])
    prop_value_query = FieldQuery(extract_as="Atomic radius", extract_all=True)
    prop_query = PropertyQuery(name=prop_name_query, value=prop_value_query)

    chemical_filter = ChemicalFilter(equal=atom, logic="MUST")
    formula_query = ChemicalFieldQuery(extract_as="formula", filter=chemical_filter)

    system_query = PifSystemQuery(chemical_formula=formula_query, properties=prop_query)
    dataset_query = DatasetQuery(id=[Filter(equal='156599')])
    data_query = DataQuery(dataset=dataset_query, system=system_query)
    test_query = PifSystemReturningQuery(query=data_query)
    pif_search_result = client.search(test_query)

    atomic_radius = pif_search_result.hits[0].system.properties[0].scalars[0].value

    print("ATOMIC RADIUS FOUND:", atom, atomic_radius)

    return float(atomic_radius)


def query_for_binary_enthalpy(binary_formula):

    prop_name_query = FieldQuery(filter=[Filter(equal="Enthalpy of mixing")])
    prop_value_query = FieldQuery(extract_as="Enthalpy of mixing", extract_all=True)
    prop_query = PropertyQuery(name=prop_name_query, value=prop_value_query)
    chemical_filter = ChemicalFilter(equal=binary_formula, logic="MUST")
    formula_query = ChemicalFieldQuery(extract_as="formula", filter=chemical_filter)
    system_query = PifSystemQuery(chemical_formula=formula_query, properties=prop_query)
    dataset_query = DatasetQuery(id=[Filter(equal='156599')])
    data_query = DataQuery(dataset=dataset_query, system=system_query)
    test_query = PifSystemReturningQuery(query=data_query)
    pif_search_result = client.search(test_query)

    binary_enthalpy_of_mixing = pif_search_result.hits[0].system.properties[0].scalars[0].value

    print("BINARY ENTHALPY FOUND:", binary_formula, binary_enthalpy_of_mixing)

    return binary_enthalpy_of_mixing


def calc_enthalpy_of_mixing(chemical_formula):

    comp = composition.Composition(chemical_formula)

    elements = []
    for element in comp.elements:
        elements.append(str(element))

    binary_enthalpies = {}

    for c in itertools.combinations(elements, 2):

        binary = "".join(c)
        binary_enthalpy = query_for_binary_enthalpy(binary)
        binary_enthalpies[binary] = binary_enthalpy

    enthalpy_of_mixing = 0

    for k, v in binary_enthalpies.items():
        binary_comp = composition.Composition(k)
        e1 = binary_comp.elements[0]
        e2 = binary_comp.elements[1]
        contribution = float(v)*float(comp.get_atomic_fraction(e1))*float(comp.get_atomic_fraction(e2))
        enthalpy_of_mixing += contribution

    enthalpy_of_mixing = round(4*enthalpy_of_mixing, 3)

    return(enthalpy_of_mixing)


def cal_entropy_of_mixing(chemical_formula):

    comp = composition.Composition(chemical_formula)
    entropy_of_mixing = 0

    for element in comp:
        ele_fraction = comp.get_atomic_fraction(element)
        contribution = float(ele_fraction)*math.log(float(ele_fraction))
        entropy_of_mixing += contribution

    R = 8.314
    entropy_of_mixing = round(R*entropy_of_mixing, 3)

    return entropy_of_mixing


def calc_atomic_size_difference(chemical_formula):

    comp = composition.Composition(chemical_formula)

    average_atomic_radius = 0
    atomic_size_difference = 0


    for element in comp:
        atomic_radius = query_for_atomic_radius(str(element))
        average_atomic_radius += float(comp.get_atomic_fraction(element))*float(atomic_radius)

    for element in comp:
        atomic_radius = query_for_atomic_radius(str(element))
        ele_fraction = float(comp.get_atomic_fraction(str(element)))
        total_contribution = ele_fraction*((1-(atomic_radius/average_atomic_radius))**2)
        atomic_size_difference += total_contribution

    atomic_size_difference = 100*math.sqrt(atomic_size_difference)

    return round(atomic_size_difference, 3)

if __name__ == "__main__":

    enthalpy_of_mixing = calc_enthalpy_of_mixing(sys.argv[1])
    print("ENTHALPY OF MIXING CALCULATED:", enthalpy_of_mixing, "kJ/mol")

    entropy_of_mixing = cal_entropy_of_mixing(sys.argv[1])
    print("ENTROPY OF MIXING CALCULATED:", entropy_of_mixing, "J*K/mol")

    atomic_size_difference = calc_atomic_size_difference(sys.argv[1])
    print("ATOMIC SIZE DIFFERENCE CALCULATED:", atomic_size_difference, "%")



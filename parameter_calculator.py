import sys
import math
import itertools
import csv
from os import environ, path
from pymatgen.core import composition
from citrination_client import *


def query_for_property(property_name, formula):

    prop_name_query = FieldQuery(filter=[Filter(equal=property_name)])
    prop_value_query = FieldQuery(extract_as=property_name, extract_all=True)
    prop_query = PropertyQuery(name=prop_name_query, value=prop_value_query)
    chemical_filter = ChemicalFilter(equal=formula, logic="MUST")
    formula_query = ChemicalFieldQuery(extract_as="formula", filter=chemical_filter)
    system_query = PifSystemQuery(chemical_formula=formula_query, properties=prop_query)
    dataset_query = DatasetQuery(id=[Filter(equal='156599')])
    data_query = DataQuery(dataset=dataset_query, system=system_query)
    test_query = PifSystemReturningQuery(query=data_query)
    pif_search_result = client.search(test_query)

    for prop in pif_search_result.hits[0].system.properties:
        if prop.name == property_name:
            property_value = prop.scalars[0].value

    # print(property_name, ":", formula, property_value)

    return float(property_value)


def calc_enthalpy_of_mixing(chemical_formula):

    comp = composition.Composition(chemical_formula)

    elements = []
    for element in comp.elements:
        elements.append(str(element))

    binary_enthalpies = {}

    for c in itertools.combinations(elements, 2):

        binary = "".join(c)
        binary_enthalpy = query_for_property("Enthalpy of mixing", binary)
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


def calc_entropy_of_mixing(chemical_formula):

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
        atomic_radius = query_for_property("Atomic radius", str(element))
        average_atomic_radius += float(comp.get_atomic_fraction(element))*float(atomic_radius)

    for element in comp:
        atomic_radius = query_for_property("Atomic radius", str(element))
        ele_fraction = float(comp.get_atomic_fraction(str(element)))
        total_contribution = ele_fraction*((1-(atomic_radius/average_atomic_radius))**2)
        atomic_size_difference += total_contribution

    atomic_size_difference = 100*math.sqrt(atomic_size_difference)

    return round(atomic_size_difference, 3)


def calc_omega(chemical_formula):

    comp = composition.Composition(chemical_formula)

    average_melting_temp = 0

    for element in comp:
        melting_temp = query_for_property("Melting temperature", str(element))
        average_melting_temp += float(comp.get_atomic_fraction(element))*float(melting_temp)

    entropy_of_mixing = calc_entropy_of_mixing(chemical_formula)
    enthalpy_of_mixing = calc_enthalpy_of_mixing(chemical_formula)

    omega = (average_melting_temp*entropy_of_mixing)/(1000*enthalpy_of_mixing)

    return round(omega, 3)


def calc_avg_VEC(chemical_formula):

    comp = composition.Composition(chemical_formula)

    avg_VEC = 0

    for element in comp:
        VEC = query_for_property("Valence Electron Configuration", str(element))
        avg_VEC += float(comp.get_atomic_fraction(element)) * float(VEC)

    return round(avg_VEC, 3)


def calc_elecneg_diff(chemical_formula):

    comp = composition.Composition(chemical_formula)

    avg_elec = 0

    for element in comp:
        elec = query_for_property("Pauling electronegativity", str(element))
        avg_elec += float(comp.get_atomic_fraction(element)) * float(elec)

    elec_diff = 0
    for element in comp:
        elec = query_for_property("Pauling electronegativity", str(element))
        ele_fraction = float(comp.get_atomic_fraction(str(element)))
        contribution = ele_fraction*((elec-avg_elec)**2)
        elec_diff += contribution

    elec_diff = math.sqrt(elec_diff)

    return round(elec_diff, 3)


def parse_template_file_for_formulas(csv_template_file):

    alloys = []

    with open(csv_template_file, "r") as csv_file:
        reader = csv.reader(csv_file)
        for index, row in enumerate(reader):
            if index == 0:
                header_row = row
                for h_index, header in enumerate(header_row):
                    if header == "FORMULA":
                        formula_index = h_index
            else:
                alloys.append(row[formula_index])


    return alloys


def add_properties_to_csv(csv_template_file, formulas, property_headers, property_values):

    infile = open(csv_template_file, "r")
    outfile = open(csv_template_file.replace(".csv", "_with_emp_parameters.csv"), "w")

    reader = csv.reader(infile)
    writer = csv.writer(outfile)
    new_rows = []

    for index, row in enumerate(reader):
        if index == 0:
            header_row = row
            header_row.extend(property_headers)
            new_rows.append(header_row)
            for h_index, header in enumerate(header_row):
                if header == "FORMULA":
                    formula_index = h_index
        else:
            if row[formula_index] in formulas:
                row.extend(property_values[formulas.index(row[formula_index])])
                new_rows.append(row)

    print(len(new_rows), len(formulas))
    writer.writerows(new_rows)


if __name__ == "__main__":

    client = CitrinationClient(environ["CITRINATION_API_KEY"], "https://citrination.com")

    input_file = sys.argv[1]

    alloys = parse_template_file_for_formulas(input_file)

    property_values = []
    property_headers = ["PROPERTY: Enthalpy of mixing (kJ/mol)", "PROPERTY: Entropy of mixing (J*K/mol)",
                        "PROPERTY: Atomic size difference (%)", "PROPERTY: $\Omega$", "PROPERTY: average VEC",
                        "PROPERTY: $\delta$ X)"]

    for alloy in alloys:
        print("----{}-----".format(alloy))

        enthalpy_of_mixing = calc_enthalpy_of_mixing(alloy)
        print("Enthalpy of mixing:", enthalpy_of_mixing, "kJ/mol")

        entropy_of_mixing = calc_entropy_of_mixing(alloy)
        print("Entropy of mixing:", entropy_of_mixing, "J*K/mol")

        atomic_size_difference = calc_atomic_size_difference(alloy)
        print("Atomic size difference:", atomic_size_difference, "%")

        omega = calc_omega(alloy)
        print("Omega:", omega, "")

        VEC = calc_avg_VEC(alloy)
        print("VEC:", VEC, "")

        elec_diff = calc_elecneg_diff(alloy)
        print("deltaX:", elec_diff, "")

        property_values.append([enthalpy_of_mixing, entropy_of_mixing, atomic_size_difference, omega, VEC, elec_diff])


    add_properties_to_csv(input_file, alloys, property_headers, property_values)



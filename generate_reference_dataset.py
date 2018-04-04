import csv
import argparse
from pypif import pif
from pypif.obj import *


def parse_csv(csv_path):
    """
    Parses Takeuchi2005 csvs (1, 2a, 2b)

    Args:
        csv_path (str): path to csv

    Returns:
        systems (list): list of parsed systems in pif format
        
    """

    top_row_elements = []
    systems = []

    with open(csv_path, 'r') as f:

        reader = csv.reader(f)

        for index, row in enumerate(reader):
            if index == 0:
                for cell in row:
                    top_row_elements.append(cell)
        f.seek(0)
        for index, element in enumerate(top_row_elements):
            if element != "":
                for csv_index, row in enumerate(reader):
                    if csv_index > 0:
                        if row[index] != "":
                            system = ChemicalSystem()
                            system.chemical_formula = row[0]+element
                            prop = Property(name="Enthalpy of mixing", scalars=row[index], units="kJ/mol")
                            system.properties = [prop]
                            system.references = [Reference(doi="10.2320/matertrans.46.2817")]
                            print(pif.dumps(system.chemical_formula), pif.dumps(system.properties))
                            systems.append(system)
            f.seek(0)
    
    return systems


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('csv', nargs='*', help='path to csv files')

    args = parser.parse_args()

    for f in args.csv:
        pifs = parse_csv(f)
        outfile = f.replace(".csv", ".json")
        pif.dump(pifs, open(outfile, "w"))
        print("PIF DUMPED: ", outfile)

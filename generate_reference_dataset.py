import csv
import sys
from pypif import pif
from pypif.obj import *

def parse_csv(csv_file):

    top_row_elements = []
    systems = []

    with open(csv_file, 'r') as f:

        reader = csv.reader(f)

        for index, row in enumerate(reader):
            if index == 0:
                for cell in row:
                    top_row_elements.append(cell)
        f.seek(0)
        for index, element in enumerate(top_row_elements):
            if element != "":
                print("Index: ", index)
                for csv_index, row in enumerate(reader):
                    if csv_index > 0:
                        if row[index] != "":
                            system = ChemicalSystem()
                            system.chemical_formula = row[0]+element
                            prop = Property(name="Enthalpy of mixing", scalars=row[index], units="kJ/mol")
                            system.properties = [prop]
                            system.references = [Reference(doi="10.2320/matertrans.46.2817")]
                            print(row[0], element, row[index])
                            print(pif.dumps(system))
                            systems.append(system)
            f.seek(0)

    outfile = csv_file.replace(".csv", ".json")
    pif.dump(systems, open(outfile, "w"))

    print(len(systems))
if __name__ == "__main__":

    parse_csv(sys.argv[1])

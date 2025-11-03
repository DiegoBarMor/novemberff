import sys

import november as nov

################################################################################
if __name__ == "__main__":
    try:
        MODE     = sys.argv[1].lower()
        PATH_PDB = sys.argv[2]
        if MODE == "prot":
            PATH_XML = "november/_data/amber99sb.xml"
        elif MODE == "rna":
            PATH_XML = "november/_data/RNA.OL3.xml"
        else:
            raise ValueError()

    except (IndexError, ValueError):
        print("Usage: python main.py (prot|rna) <path_pdb> [<folder_output>]")
        sys.exit(1)


    calc = nov.EnergyCalculator(PATH_PDB, PATH_XML).calc_energies_rna()

    if len(sys.argv) >= 4:
        FOLDER_OUTPUT = sys.argv[3]
        calc.save_energy_arrays(FOLDER_OUTPUT)

    else:
        calc.display_energies()


################################################################################

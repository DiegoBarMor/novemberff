import sys
import warnings

import novemberff as nov

################################################################################
if __name__ == "__main__":
    warnings.filterwarnings("ignore", module = "MDAnalysis.*")

    try:
        MODE     = sys.argv[1].lower()
        PATH_PDB = sys.argv[2]
        if MODE not in ("prot", "rna"):
            raise ValueError()

    except (IndexError, ValueError):
        print("Usage: python main.py (prot|rna) <path_pdb> [<folder_output>]")
        sys.exit(1)

    calc = nov.EnergyCalculator.with_prot_ff(PATH_PDB) \
        if MODE == "prot" else \
        nov.EnergyCalculator.with_rna_ff(PATH_PDB)

    calc.calc_energies()

    if len(sys.argv) >= 4:
        FOLDER_OUTPUT = sys.argv[3]
        calc.save_energy_arrays(FOLDER_OUTPUT)

    else:
        calc.display_energies()


################################################################################

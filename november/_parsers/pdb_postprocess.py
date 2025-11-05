# //////////////////////////////////////////////////////////////////////////////
class PDBPostProcess:
    # --------------------------------------------------------------------------
    @staticmethod
    def fix_prot_reslabels_aliases(pdb):
        for res in pdb.residues: # [WIP] hardcoded for now to work with amber99sb.xml
            if res.resname == "HIS": res.resname = "HIE"


    # --------------------------------------------------------------------------
    @staticmethod
    def fix_prot_reslabels_termini(pdb):
        residues = list(pdb.residues)
        nres = len(residues)
        for i,res in enumerate(residues): # [WIP] needs improvement, e.g. for multiple chains
            if i == 0:
                res.resname = f"N{res.resname}"
            elif i == nres - 1:
                res.resname = f"C{res.resname}"
                for atom in res.atoms:
                    if   atom.name == "O1": atom.name = "O"
                    elif atom.name == "O2": atom.name = "OXT"


    # --------------------------------------------------------------------------
    @staticmethod
    def fix_prot_atomlabels_aliases(pdb):
        for atom in pdb.atoms: # [WIP] crudely hardcoded to at least support the prot.pdb example
            if (atom.residue.resname.endswith("ASP")):
                match atom.name:
                    case "HB1": atom.name = "HB3"
                continue

            if (atom.residue.resname.endswith("ILE")):
                match atom.name:
                    case "CD": atom.name = "CD1"
                continue

            if (atom.residue.resname.endswith("GLU")):
                match atom.name:
                    case "HB1": atom.name = "HB3"
                    case "HG1": atom.name = "HG3"
                continue


    # --------------------------------------------------------------------------
    @staticmethod
    def fix_rna_reslabels_termini(pdb):
        for res in pdb.residues:
            res_type = res.resname[0]
            if res_type not in "UCAG": continue

            isTerminus3 = False
            isTerminus5 = False
            for atom in res.atoms:
                if atom.name == "HO5'": isTerminus5 = True
                if atom.name == "HO3'": isTerminus3 = True

            if isTerminus3 and isTerminus5:
                res.resname = res_type + 'N'
            elif isTerminus3:
                res.resname = res_type + '3'
            elif isTerminus5:
                res.resname = res_type + '5'
            else:
                res.resname = res_type


# //////////////////////////////////////////////////////////////////////////////

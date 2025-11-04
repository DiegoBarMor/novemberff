import numpy as np
from pathlib import Path
import MDAnalysis as mda

import november as nov

# //////////////////////////////////////////////////////////////////////////////
class EnergyCalculator:
    def __init__(self, path_pdb, path_forcefield):
        self._path_pdb = path_pdb

        self._pdb = mda.Universe(str(path_pdb), to_guess = ["elements", "bonds"])
        self._bgraph = nov.BondGraph(
            self._pdb.atoms, self._pdb.bonds,
            coords = self._pdb.atoms.positions
        )
        self._forcefield = nov.ForceField.from_xml(path_forcefield)

        self._natoms     = len(self._bgraph.atoms)
        self._nbonds     = len(self._bgraph.bonds)
        self._nangles    = len(self._bgraph.angles)
        self._npropers   = len(self._bgraph.proper_diheds)
        self._nimpropers = len(self._bgraph.improper_diheds)
        self._nnonbonded = len(self._bgraph.nonbonded_pairs)

        self._arr_bond_energies     = np.zeros(self._nbonds)
        self._arr_angle_energies    = np.zeros(self._nangles)
        self._arr_proper_energies   = np.zeros(self._npropers)
        self._arr_improper_energies = np.zeros(self._nimpropers)
        self._arr_lennardj_energies = np.zeros(self._nnonbonded)
        self._arr_coulomb_energies  = np.zeros(self._nnonbonded)


    # --------------------------------------------------------------------------
    def calc_energies_prot(self) -> "EnergyCalculator":
        self._fix_prot_reslabels_aliases()
        self._fix_prot_reslabels_termini()
        self._fix_prot_atomlabels_aliases()
        self._calc_energies()
        return self


    # --------------------------------------------------------------------------
    def calc_energies_rna(self) -> "EnergyCalculator":
        self._fix_rna_reslabels_termini()
        self._calc_energies()
        return self


    # --------------------------------------------------------------------------
    def get_ebond_arr(self):     return self._arr_bond_energies
    def get_eangle_arr(self):    return self._arr_angle_energies
    def get_eproper_arr(self):   return self._arr_proper_energies
    def get_eimproper_arr(self): return self._arr_improper_energies
    def get_elennardj_arr(self): return self._arr_lennardj_energies
    def get_ecoulomb_arr(self):  return self._arr_coulomb_energies


    # --------------------------------------------------------------------------
    def get_ebond_sum(self):     return np.sum(self._arr_bond_energies)
    def get_eangle_sum(self):    return np.sum(self._arr_angle_energies)
    def get_eproper_sum(self):   return np.sum(self._arr_proper_energies)
    def get_eimproper_sum(self): return np.sum(self._arr_improper_energies)
    def get_elennardj_sum(self): return np.sum(self._arr_lennardj_energies)
    def get_ecoulomb_sum(self):  return np.sum(self._arr_coulomb_energies)

    def get_edihed_sum(self):     return self.get_eproper_sum() + self.get_eimproper_sum()
    def get_enonbonded_sum(self): return self.get_elennardj_sum() + self.get_ecoulomb_sum()


    # --------------------------------------------------------------------------
    def display_energies(self):
        print(f">>> Amber energies (Custom) for '{self._path_pdb}':")
        print("  BondEnergy:",  self.get_ebond_sum())
        print("  AngleEnergy:", self.get_eangle_sum())
        print("  DihedEnergy:", self.get_edihed_sum())
        print("  Nonbonded:", self.get_enonbonded_sum())
        print("  ... LennardJones:", self.get_elennardj_sum())
        print("  ... Coulomb:", self.get_ecoulomb_sum())


    # --------------------------------------------------------------------------
    def save_energy_arrays(self, folder_output: str | Path):
        folder_output = Path(folder_output)
        folder_output.mkdir(parents = True, exist_ok = True)
        np.save(folder_output / f"bonds.npy",     self._arr_bond_energies    )
        np.save(folder_output / f"angles.npy",    self._arr_angle_energies   )
        np.save(folder_output / f"propers.npy",   self._arr_proper_energies  )
        np.save(folder_output / f"impropers.npy", self._arr_improper_energies)
        np.save(folder_output / f"lennardjs.npy", self._arr_lennardj_energies)
        np.save(folder_output / f"coulombs.npy",  self._arr_coulomb_energies )


    # --------------------------------------------------------------------------
    def _calc_energies(self):
        ##### custom_e = (1 ± ε) * openmm_e
        self._calc_ebonded()   # ε = 1e-14
        # self._calc_eangles()   # ε = 1e-15
        # self._calc_ediheds()   # ε = 1e-4
        # self._calc_nonbonded() # ε = 1e-6


    # --------------------------------------------------------------------------
    def _calc_ebonded(self):
        ###### [ORIGINAL SOURCE] platforms/reference/src/SimTKReference/ReferenceHarmonicBondlxn.cpp@111
        for i,(a0,a1) in enumerate(self._bgraph.bonds):
            if (a0.element == 'H') or (a1.element == 'H'): # apply HBond constraints
                continue

            ff_bond  = self._forcefield.get_ffbond(a0, a1)
            distance = self._bgraph.calc_dist_2atoms(a0.index, a1.index)
            energy = 0.5 * ff_bond.k * (ff_bond.length - distance)**2
            self._arr_bond_energies[i] = energy


    # --------------------------------------------------------------------------
    def _calc_eangles(self):
        for i,(a0,a1,a2) in enumerate(self._bgraph.angles):
            ff_angle = self._forcefield.get_ffangle(a0, a1, a2)
            angle  = self._bgraph.calc_angle_3atoms(a0.index, a1.index, a2.index)
            energy = 0.5 * ff_angle.k * (angle - ff_angle.angle)**2
            self._arr_angle_energies[i] = energy


    # --------------------------------------------------------------------------
    def _calc_ediheds(self):
        def _edihed_contributions(ff_dihed, ordered_atoms, is_proper):
            a0,a1,a2,a3 = ordered_atoms
            angle = self._bgraph.calc_dihed_4atoms(a0, a1, a2, a3)
            e1 = ff_dihed.calc_energy(angle, contributor = 1, proper = is_proper)
            e2 = ff_dihed.calc_energy(angle, contributor = 2, proper = is_proper)
            e3 = ff_dihed.calc_energy(angle, contributor = 3, proper = is_proper)
            e4 = ff_dihed.calc_energy(angle, contributor = 4, proper = is_proper)
            return e1 + e2 + e3 + e4

        for i,atoms in enumerate(self._bgraph.proper_diheds):
            for ff_dihed, ordered_atoms in self._forcefield.iter_ffpropers(atoms):
                energy = _edihed_contributions(ff_dihed, ordered_atoms, is_proper = True)
                if energy:
                    self._arr_proper_energies[i] = energy
                    break

        for i,atoms in enumerate(self._bgraph.improper_diheds):
            for ff_dihed, ordered_atoms in self._forcefield.iter_ffimpropers(atoms):
                energy = _edihed_contributions(ff_dihed, ordered_atoms, is_proper = False)
                if energy:
                    self._arr_improper_energies[i] = energy
                    break


    # --------------------------------------------------------------------------
    def _calc_nonbonded(self):
        ### [ORIGINAL SOURCE] platforms/cpu/src/CpuKernels.cpp@CpuCalcNonbondedForceKernel::computeParameters
        ### [ORIGINAL SOURCE] platforms/cpu/include/CpuNonbondedForceFvec.h@CpuNonbondedForceFvec::calculateBlockIxnImpl
        ### [ORIGINAL SOURCE] openmmapi/src/NonbondedForce.cpp@NonbondedForce::createExceptionsFromBonds
        ### [ORIGINAL SOURCE] platforms/reference/src/SimTKReference/ReferenceLJCoulomb14.cpp@ReferenceLJCoulomb14::calculateBondIxn
        for i,(a0,a1,bond_edge_dist) in enumerate(self._bgraph.nonbonded_pairs):
            ff_nb0, ff_nb1, charge0, charge1 = self._forcefield.get_ffnonbonded(a0, a1)

            inverse_r = 1 / self._bgraph.calc_dist_2atoms(a0.index, a1.index)
            sigma = 0.5 * (ff_nb0.sigma + ff_nb1.sigma)
            epsilon = 4.0 * np.sqrt(ff_nb0.epsilon * ff_nb1.epsilon)
            charge = charge0 * charge1
            sig6 = (inverse_r * sigma)**6

            energy_lj = epsilon * sig6 * (sig6 - 1.0)
            energy_coul = inverse_r * charge * self._forcefield.one_4pi_eps0

            if bond_edge_dist == nov.BondEdgeDist.mid: # apply 1-4 scaling to atoms separated by 3 bonds
                energy_lj   *= self._forcefield.lj_14scale
                energy_coul *= self._forcefield.coulomb_14scale # AKA "epsilon_factor"

            self._arr_lennardj_energies[i] = energy_lj
            self._arr_coulomb_energies[i] = energy_coul


    # --------------------------------------------------------------------------
    def _fix_prot_atomlabels_aliases(self):
        for atom in self._pdb.atoms: # [WIP] crudely hardcoded to at least support the prot.pdb example
            if (atom.residue.resname == "ILE") and (atom.name == "CD"):
                atom.name = "CD1"; continue


    # --------------------------------------------------------------------------
    def _fix_prot_reslabels_aliases(self):
        for res in self._pdb.residues: # [WIP] hardcoded for now to work with amber99sb.xml
            if res.resname == "HIS": res.resname = "HIE"


    # --------------------------------------------------------------------------
    def _fix_prot_reslabels_termini(self):
        residues = list(self._pdb.residues)
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
    def _fix_rna_reslabels_termini(self):
        for res in self._pdb.residues:
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

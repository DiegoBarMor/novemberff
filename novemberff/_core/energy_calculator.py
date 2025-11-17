import warnings
import numpy as np
from pathlib import Path
import MDAnalysis as mda

import novemberff as nov

# //////////////////////////////////////////////////////////////////////////////
class EnergyCalculator:
    def __init__(self, path_pdb: Path, forcefield: str | Path, path_xtc: Path | None = None):
        self._path_pdb = path_pdb
        self._path_xml = nov.Utils.solve_forcefield_path(forcefield)
        self._path_xtc = None

        self._pdb = mda.Universe(str(path_pdb), to_guess = ["elements", "bonds"])
        self._bgraph = nov.BondGraph(
            self._pdb.atoms, self._pdb.bonds,
            coords = self._pdb.atoms.positions
        )
        self._forcefield = nov.ForceField.from_xml(self._path_xml)
        self._traj = None

        if path_xtc is not None:
            self.load_traj(path_xtc)

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

        self._do_bonds     : bool = True
        self._do_angles    : bool = True
        self._do_diheds    : bool = True
        self._do_nonbonded : bool = True
        self._computed_bonds     : bool = False
        self._computed_angles    : bool = False
        self._computed_diheds    : bool = False
        self._computed_nonbonded : bool = False


    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ CONSTRUCTORS
    # --------------------------------------------------------------------------
    @classmethod
    def with_prot_ff(cls, path_pdb: Path, path_xtc: Path | None = None) -> "EnergyCalculator":
        """Initialize EnergyCalculator with a default protein force field (Amber99SB) and applies relevant patches to the PDB."""
        obj = cls(path_pdb, "amber99sb", path_xtc)
        nov.PDBPostProcess.fix_prot_atomlabels_aliases(obj._pdb)
        nov.PDBPostProcess.fix_prot_reslabels_aliases (obj._pdb)
        nov.PDBPostProcess.fix_prot_reslabels_termini (obj._pdb)
        return obj


    # --------------------------------------------------------------------------
    @classmethod
    def with_rna_ff(cls, path_pdb: Path, path_xtc: Path | None = None) -> "EnergyCalculator":
        """Initialize EnergyCalculator with a default RNA force field (RNA.OL3) and applies relevant patches to the PDB."""
        obj = cls(path_pdb, "rna.ol3", path_xtc)
        nov.PDBPostProcess.fix_rna_reslabels_termini(obj._pdb)
        return obj


    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ MAIN PUBLIC METHODS
    # --------------------------------------------------------------------------
    def toggle_interactions(self, bonds = True, angles = True, diheds = True, nonbonded = True):
        self._do_bonds     = bonds
        self._do_angles    = angles
        self._do_diheds    = diheds
        self._do_nonbonded = nonbonded


    # --------------------------------------------------------------------------
    def calc_energies(self):
        if self._do_bonds:     self._calc_ebonded()
        if self._do_angles:    self._calc_eangles()
        if self._do_diheds:    self._calc_ediheds()
        if self._do_nonbonded: self._calc_nonbonded()

    # --------------------------------------------------------------------------
    def has_run(self) -> bool:
        return (
            self._computed_bonds  or self._computed_angles or
            self._computed_diheds or self._computed_nonbonded
        )


    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ TRAJECTORY METHODS
    # --------------------------------------------------------------------------
    def load_traj(self, path_xtc: Path):
        self._path_xtc = path_xtc
        self._traj = mda.Universe(str(self._path_pdb), str(path_xtc))

    # --------------------------------------------------------------------------
    def itercalc_traj_energies(self):
        self._assert_traj_loaded()
        for i,_ in enumerate(self._traj.trajectory):
            self.set_positions(self._traj.atoms.positions)
            self.calc_energies()
            yield (
                i,
                self.get_array_ebond()      if self._do_bonds     else None,
                self.get_array_eangle()     if self._do_angles    else None,
                self.get_array_edihed()     if self._do_diheds    else None,
                self.get_array_enonbonded() if self._do_nonbonded else None,
            )

    # --------------------------------------------------------------------------
    def get_traj_energy_arrays(self, verbose = False) -> tuple[np.ndarray[float, float] | None]:
        self._assert_traj_loaded()

        nframes = self.get_nframes()
        mat_ebond  = np.zeros((nframes, self.get_nbonds()))     if self._do_bonds     else None
        mat_eangle = np.zeros((nframes, self.get_nangles()))    if self._do_angles    else None
        mat_edihed = np.zeros((nframes, self.get_ndiheds()))    if self._do_diheds    else None
        mat_ennb   = np.zeros((nframes, self.get_nnonbonded())) if self._do_nonbonded else None

        for i,_ in enumerate(self._traj.trajectory):
            if verbose and not (i % 100): print(f"Progress: {i}/{nframes}")
            self.set_positions(self._traj.atoms.positions)
            self.calc_energies()
            if self._do_bonds:     mat_ebond [i,:] = self.get_array_ebond()
            if self._do_angles:    mat_eangle[i,:] = self.get_array_eangle()
            if self._do_diheds:    mat_edihed[i,:] = self.get_array_edihed()
            if self._do_nonbonded: mat_ennb  [i,:] = self.get_array_enonbonded()

        return mat_ebond, mat_eangle, mat_edihed, mat_ennb


    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ GETTERS / SETTERS
    # --------------------------------------------------------------------------
    def get_nframes(self) -> int:
        if self._traj is None: return 1
        return len(self._traj.trajectory)

    # --------------------------------------------------------------------------
    def get_natoms(self)     -> int: return self._natoms
    def get_nbonds(self)     -> int: return self._nbonds
    def get_nangles(self)    -> int: return self._nangles
    def get_npropers(self)   -> int: return self._npropers
    def get_nimpropers(self) -> int: return self._nimpropers
    def get_nnonbonded(self) -> int: return self._nnonbonded
    def get_ndiheds(self)    -> int: return self._npropers + self._nimpropers

    # --------------------------------------------------------------------------
    def get_array_ebond(self):
        if not (self._computed_bonds):
            warnings.warn("Bond energies have not been computed yet. Returning zeroed array.")
        return self._arr_bond_energies

    # --------------------------------------------------------------------------
    def get_array_eangle(self):
        if not (self._computed_angles):
            warnings.warn("Angle energies have not been computed yet. Returning zeroed array.")
        return self._arr_angle_energies

    # --------------------------------------------------------------------------
    def get_array_eproper(self):
        if not (self._computed_diheds):
            warnings.warn("Proper dihedral energies have not been computed yet. Returning zeroed array.")
        return self._arr_proper_energies

    # --------------------------------------------------------------------------
    def get_array_eimproper(self):
        if not (self._computed_diheds):
            warnings.warn("Improper dihedral energies have not been computed yet. Returning zeroed array.")
        return self._arr_improper_energies

    # --------------------------------------------------------------------------
    def get_array_elennardj(self):
        if not (self._computed_nonbonded):
            warnings.warn("Lennard-Jones energies have not been computed yet. Returning zeroed array.")
        return self._arr_lennardj_energies

    # --------------------------------------------------------------------------
    def get_array_ecoulomb(self):
        if not (self._computed_nonbonded):
            warnings.warn("Coulomb energies have not been computed yet. Returning zeroed array.")
        return self._arr_coulomb_energies

    # --------------------------------------------------------------------------
    def get_array_edihed(self):
        return np.concat((
            self.get_array_eproper(), self.get_array_eimproper(),
        ))

    # --------------------------------------------------------------------------
    def get_array_enonbonded(self):
        return np.concat((
            self.get_array_elennardj(), self.get_array_ecoulomb(),
        ))


    # --------------------------------------------------------------------------
    def get_sum_ebond(self):     return np.sum(self.get_array_ebond())
    def get_sum_eangle(self):    return np.sum(self.get_array_eangle())
    def get_sum_eproper(self):   return np.sum(self.get_array_eproper())
    def get_sum_eimproper(self): return np.sum(self.get_array_eimproper())
    def get_sum_elennardj(self): return np.sum(self.get_array_elennardj())
    def get_sum_ecoulomb(self):  return np.sum(self.get_array_ecoulomb())

    def get_sum_edihed(self):
        return self.get_sum_eproper() + self.get_sum_eimproper()

    def get_sum_enonbonded(self):
        return self.get_sum_elennardj() + self.get_sum_ecoulomb()


    # --------------------------------------------------------------------------
    def set_positions(self, positions: np.ndarray):
        self._pdb.atoms.positions = positions
        self._bgraph.set_positions(positions)


    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ OUTPUT OPERATIONS
    # --------------------------------------------------------------------------
    def display_energies(self):
        print(f">>> Amber energies (NovemberFF) for '{self._path_pdb}':")
        print("  BondEnergy:",  self.get_sum_ebond())
        print("  AngleEnergy:", self.get_sum_eangle())
        print("  DihedEnergy:", self.get_sum_edihed())
        print("  Nonbonded:",   self.get_sum_enonbonded())
        print("  ... LennardJones:", self.get_sum_elennardj())
        print("  ... Coulomb:",      self.get_sum_ecoulomb())


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
    def save_metadata_csv(self, path_csv: str | Path):
        def _idxs(atoms) -> int:
            return (a.index for a in atoms)

        def _generic_metadata(atoms) -> tuple[str, str, str, str]:
            return (
                '-'.join(map(lambda a: f"{a.index:02}",      atoms)),
                '-'.join(map(lambda a: a.name.lower(),       atoms)),
                '-'.join(map(lambda a: f"{a.residue.resid}", atoms)),
                '-'.join(map(lambda a: a.residue.resname,    atoms)),
            )

        header = ["interaction", "kind", "intshort", "idxs", "names", "residxs", "resnames", "geometry"]
        nov.Utils.init_csv(path_csv, *header)

        if self._do_bonds:
            for atoms in self._bgraph.bonds:
                nov.Utils.append_to_csv(path_csv,
                    "bond", "bond", 'b',
                    *_generic_metadata(atoms),
                    self._bgraph.calc_dist_2atoms(*_idxs(atoms)),
                )

        if self._do_angles:
            for atoms in self._bgraph.angles:
                nov.Utils.append_to_csv(path_csv,
                    "angle", "angle", 'a',
                    *_generic_metadata(atoms),
                    self._bgraph.calc_angle_3atoms(*_idxs(atoms)),
                )

        if self._do_diheds:
            for atoms in self._bgraph.proper_diheds:
                nov.Utils.append_to_csv(path_csv,
                    "dihed", "proper", 'p',
                    *_generic_metadata(atoms),
                    self._bgraph.calc_dihed_4atoms(*_idxs(atoms)),
                )
            for atoms in self._bgraph.improper_diheds:
                nov.Utils.append_to_csv(path_csv,
                    "dihed", "improper", 'i',
                    *_generic_metadata(atoms),
                    self._bgraph.calc_dihed_4atoms(*_idxs(atoms)),
                )

        if self._do_nonbonded:
            for a0, a1, _ in self._bgraph.nonbonded_pairs:
                atoms = (a0, a1)
                nov.Utils.append_to_csv(path_csv,
                    "nonbonded", "nonbonded", 'n',
                    *_generic_metadata(atoms),
                    self._bgraph.calc_dist_2atoms(*_idxs(atoms)),
                )


    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ PRIVATE METHODS
    # --------------------------------------------------------------------------
    def _calc_ebonded(self):
        self._computed_bonds = True

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
        self._computed_angles = True

        for i,(a0,a1,a2) in enumerate(self._bgraph.angles):
            ff_angle = self._forcefield.get_ffangle(a0, a1, a2)
            angle  = self._bgraph.calc_angle_3atoms(a0.index, a1.index, a2.index)
            energy = 0.5 * ff_angle.k * (angle - ff_angle.angle)**2
            self._arr_angle_energies[i] = energy


    # --------------------------------------------------------------------------
    def _calc_ediheds(self):
        self._computed_diheds = True

        def _edihed_contributions(ff_dihed: nov.FFDihedral, ordered_atoms, is_proper):
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
        self._computed_nonbonded = True

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
    def _assert_traj_loaded(self):
        if self._traj is None:
            raise ValueError("No trajectory loaded. Cannot iterate over frames.")


# //////////////////////////////////////////////////////////////////////////////

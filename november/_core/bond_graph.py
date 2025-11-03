import math
import numpy as np
from enum import Enum, auto

import november as nov

# //////////////////////////////////////////////////////////////////////////////
class BondEdgeDist(Enum):
    short = auto()
    mid   = auto()
    long  = auto()


# //////////////////////////////////////////////////////////////////////////////
class BondGraph:
    LARGE_INT = 9

    # --------------------------------------------------------------------------
    def __init__(self, atoms, bonds, coords):
        self.atoms = list(atoms)
        self.coords = list(coords)
        self._bond_matrix = np.full((len(self.atoms), len(self.atoms)), self.LARGE_INT, dtype = int)

        self.bonds = [] # [(a1, a2), ...]
        self.angles = [] # [(a1, a2, a3), ...]
        self.proper_diheds = [] # [(a1, a2, a3, a4), ...]
        self.improper_diheds = set() # {(a1, a2, a3, a4), ...}
        self.nonbonded_pairs: list[tuple] # [(a1, a2, BondEdgeDist), ...]


        # ----------------------------------------------------------------------
        for bond in bonds: # init with direct bonds
            a1,a2 = bond
            self._bond_matrix[a1.index, a2.index] = 1
            self._bond_matrix[a2.index, a1.index] = 1

        # ----------------------------------------------------------------------
        for i,row_i in enumerate(self._bond_matrix):
            self._bond_matrix[i][i] = 0
            for j,val_j in enumerate(row_i[i:], start = i):
                if val_j == 1:
                    self.bonds.append((self.atoms[i], self.atoms[j]))

        # ----------------------------------------------------------------------
        for idx_bond,bond0 in enumerate(self.bonds[:-1]):
            ### bonds are of the form (i-j)
            a0i, a0j = bond0 # (ij)
            for bond1 in self.bonds[idx_bond+1:]:
                a1k, a1l = bond1 # (kl)

                ### (i*) and (k*)
                if   a0j == a1l: self._add_angle(a0i, a0j, a1k)

                ### (i*) and (*l)
                elif a0j == a1k: self._add_angle(a0i, a0j, a1l)

                ### (*j) and (k*)
                elif a0i == a1l: self._add_angle(a0j, a0i, a1k)

                ### (*j) and (*l)
                elif a0i == a1k: self._add_angle(a0j, a0i, a1l)

        # ----------------------------------------------------------------------
        for idx_angle,angle0 in enumerate(self.angles[:-1]):
            ### angles are of the form (i-j-k)
            a0i, a0j, a0k = angle0 # (ijk)
            angle0_r = a0k, a0j, a0i # (kji)
            for angle1 in self.angles[idx_angle+1:]:
                a1l, a1m, a1n = angle1 # (lmn)

                ### (*jk) and (l*n)
                if   a0i == a1m:
                    ### (**k) and (**n)
                    if   a0j == a1l: self._add_proper_dihed(*angle0_r, a1n)

                    ### (**k) and (l**)
                    elif a0j == a1n: self._add_proper_dihed(*angle0_r, a1l)

                ### (ij*) and (l*n)
                elif a0k == a1m:
                    ### (i**) and (**n)
                    if   a0j == a1l: self._add_proper_dihed(*angle0, a1n)

                    ### (i**) and (l**)
                    elif a0j == a1n: self._add_proper_dihed(*angle0, a1l)

                ### (i*k) and (l*n)
                elif a0j == a1m:
                    ### (i**) and (**n)
                    if   a0k == a1l: self._add_improper_dihed(*angle0, a1n)

                    ### (i**) and (l**)
                    elif a0k == a1n: self._add_improper_dihed(*angle0, a1l)

                    ### (**k) and (**n)
                    elif a0i == a1l: self._add_improper_dihed(*angle0_r, a1n)

                    ### (**k) and (l**)
                    elif a0i == a1n: self._add_improper_dihed(*angle0_r, a1l)

                ### (*jk) and (*mn)
                elif a0i == a1l: self._set_cell_value(a0k, a1n, 4)

                ### (ij*) and (*mn)
                elif a0k == a1l: self._set_cell_value(a0i, a1n, 4)

                ### (*jk) and (lm*)
                elif a0i == a1n: self._set_cell_value(a0k, a1l, 4)

                ### (ij*) and (lm*)
                elif a0k == a1n: self._set_cell_value(a0i, a1l, 4)

        # ----------------------------------------------------------------------
        self.nonbonded_pairs = list(filter(
            lambda t: t[2] != BondEdgeDist.short, # skip atoms separated by less than 3 bonds
            (
                (a0, a1, self._get_bond_edge_dist(a0.index, a1.index))
                for idx,a0 in enumerate(self.atoms[:-1])
                for a1 in self.atoms[idx+1:]
            )
        ))


    # ----------------------------------------------------------------------
    def calc_dist_2atoms(self, idx0, idx1):
        dist = self.coords[idx0] - self.coords[idx1]
        dist /= 10 # angstroms to nanometers
        return np.linalg.norm(dist)


    # ----------------------------------------------------------------------
    def calc_angle_3atoms(self, idx0, idx1, idx2):
        vec0 = self.coords[idx0] - self.coords[idx1]
        vec1 = self.coords[idx2] - self.coords[idx1]
        cosine = np.dot(vec0, vec1) / (math.sqrt(np.dot(vec0, vec0) * np.dot(vec1, vec1)))
        return math.acos(cosine) # radians


    # ----------------------------------------------------------------------
    def calc_dihed_4atoms(self, idx0, idx1, idx2, idx3):
        # Assumes the following order of atoms:
        # PROPER DIHEDRALS
        #   1   3 :   1 1   3
        #  / \ /  :  / \ \ /
        # 0   2   : 0   2 2
        #
        # IMPROPER DIHEDRALS
        #   2   :   2 2   [WIP] do the operations still make sense geometrically?
        #   |   :   | |
        #   0   :   0 0
        #  / \  :  /   \
        # 1   3 : 1     3

        vec0 = self.coords[idx3] - self.coords[idx2]
        vec1 = self.coords[idx1] - self.coords[idx2]
        vec2 = self.coords[idx1] - self.coords[idx0]
        n0 = np.cross(vec0, vec1)
        n1 = np.cross(vec1, vec2)
        cosine = np.dot(n0, n1) / (math.sqrt(np.dot(n0, n0) * np.dot(n1, n1)))

        sign = 1 if np.dot(vec0, n1) > 0 else -1
        cosine = max(-1.0, min(1.0, cosine))
        return sign * math.acos(cosine)


    # ----------------------------------------------------------------------
    def _add_angle(self, a0, a1, a2):
        self.angles.append((a0, a1, a2))
        self._set_cell_value(a0, a2, 2)
        self._set_cell_value(a2, a0, 2)


    # ----------------------------------------------------------------------
    def _add_proper_dihed(self, a0, a1, a2, a3):
        self.proper_diheds.append((a0, a1, a2, a3))
        self._set_cell_value(a0, a3, 3)
        self._set_cell_value(a3, a0, 3)


    # ----------------------------------------------------------------------
    def _add_improper_dihed(self, a0, a1, a2, a3):
        #   2
        #   |
        #   1
        #  / \
        # 0   3
        a0, a2, a3 = sorted((a0, a2, a3), key = lambda atom: atom.name)
        self.improper_diheds.add((a0, a1, a2, a3))


    # ----------------------------------------------------------------------
    def _set_cell_value(self, a0, a1, value):
        i = a0.index
        j = a1.index
        self._bond_matrix[i][j] = min(value, self._bond_matrix[i][j])
        self._bond_matrix[j][i] = min(value, self._bond_matrix[j][i])


    # ----------------------------------------------------------------------
    def _get_bond_edge_dist(self, i0, i1):
        dist = self._bond_matrix[i0][i1]
        if dist >  3: return BondEdgeDist.long
        if dist == 3: return BondEdgeDist.mid
        return BondEdgeDist.short


################################################################################
if __name__ == "__main__":
    class TestAtom:
        def __init__(self, index):
            self.index = index

    class TestBond:
        def __init__(self, atom1, atom2):
            self.atom1 = atom1
            self.atom2 = atom2

    atoms = [TestAtom(i) for i in range(10)]
    bonds = [
        TestBond(atoms[0], atoms[1]),
        TestBond(atoms[0], atoms[4]),
        TestBond(atoms[1], atoms[7]),
        TestBond(atoms[2], atoms[3]),
        TestBond(atoms[2], atoms[4]),
        TestBond(atoms[2], atoms[7]),
        TestBond(atoms[4], atoms[5]),
        TestBond(atoms[5], atoms[6]),
        TestBond(atoms[6], atoms[8]),
        TestBond(atoms[8], atoms[9]),
    ]
    coords = np.zeros((8, 3))

    bond_graph = BondGraph(atoms, bonds, coords)
    for a1,a2 in bond_graph.bonds:
        print(a1.index, a2.index)

    for a1,a2,a3 in bond_graph.angles:
        print(a1.index, a2.index, a3.index)

    for a0,a1,a2,a3 in bond_graph.proper_diheds:
        print(a0.index, a1.index, a2.index, a3.index)

    print(bond_graph._bond_matrix)


################################################################################

import math

import november as nov

# //////////////////////////////////////////////////////////////////////////////
class FFDihedral:
    ATOMTYPES = {}

    # --------------------------------------------------------------------------
    def __init__(self,
        isProper: bool, k1: float, k2: float, k3: float, k4: float,
        periodicity1: int, periodicity2: int, periodicity3: int, periodicity4: int,
        phase1: float, phase2: float, phase3: float, phase4: float,
        types: tuple[str, str, str, str]
    ):
        self.isProper = isProper
        self.k1 = k1
        self.k2 = k2
        self.k3 = k3
        self.k4 = k4
        self.periodicity1 = periodicity1
        self.periodicity2 = periodicity2
        self.periodicity3 = periodicity3
        self.periodicity4 = periodicity4
        self.phase1 = phase1
        self.phase2 = phase2
        self.phase3 = phase3
        self.phase4 = phase4
        self.type1 = self.ATOMTYPES[types[0]] if types[0] else None
        self.type2 = self.ATOMTYPES[types[1]] if types[1] else None
        self.type3 = self.ATOMTYPES[types[2]] if types[2] else None
        self.type4 = self.ATOMTYPES[types[3]] if types[3] else None


    # --------------------------------------------------------------------------
    def __repr__(self):
        return f"FFDihedral({self.type1}, {self.type2}, {self.type3}, {self.type4})"


    # --------------------------------------------------------------------------
    @classmethod
    def register_atomtypes(cls, atomtypes: dict[str, nov.FFAtomType]):
        cls.ATOMTYPES = atomtypes


    # --------------------------------------------------------------------------
    def calc_energy(self, angle: float, contributor: int, proper: bool) -> float:
        match contributor:
            case 1:
                k = self.k1
                phase = self.phase1
                periodicity = self.periodicity1
            case 2:
                k = self.k2
                phase = self.phase2
                periodicity = self.periodicity2
            case 3:
                k = self.k3
                phase = self.phase3
                periodicity = self.periodicity3
            case 4:
                k = self.k4
                phase = self.phase4
                periodicity = self.periodicity4

        if k is None: return 0 # skip irrelevant contributors
        if self.isProper ^ proper: return 0

        ### [ORIGINAL SOURCE] platforms/reference/src/SimTKReference/ReferenceProperDihedralBond.cpp@ReferenceProperDihedralBond::calculateBondIxn
        delta_angle = periodicity * angle - phase # [WIP]
        return k * (1.0 + math.cos(delta_angle)) # [WIP]


# //////////////////////////////////////////////////////////////////////////////

import math

import november as nov

# //////////////////////////////////////////////////////////////////////////////
class FFDihedral(nov.FFInteraction):
    def __init__(self,
        isProper: bool,
        ks: tuple[float, float, float, float],
        periodicities: tuple[int, int, int, int],
        phases: tuple[float, float, float, float],
        type1: "nov.FFAtomType", type2: "nov.FFAtomType",
        type3: "nov.FFAtomType", type4: "nov.FFAtomType",
    ):
        self.isProper = isProper
        self.k1 = ks[0]
        self.k2 = ks[1]
        self.k3 = ks[2]
        self.k4 = ks[3]
        self.periodicity1 = periodicities[0]
        self.periodicity2 = periodicities[1]
        self.periodicity3 = periodicities[2]
        self.periodicity4 = periodicities[3]
        self.phase1 = phases[0]
        self.phase2 = phases[1]
        self.phase3 = phases[2]
        self.phase4 = phases[3]
        self.type1 = type1
        self.type2 = type2
        self.type3 = type3
        self.type4 = type4


    # --------------------------------------------------------------------------
    def __repr__(self):
        return f"FFDihedral({self.type1}, {self.type2}, {self.type3}, {self.type4})"


    # --------------------------------------------------------------------------
    @classmethod
    def register_dihed(cls,
        isProper: bool,
        ks: tuple[float, float, float, float],
        periodicities: tuple[int, int, int, int],
        phases: tuple[float, float, float, float],
        strs_types:   tuple[str|None, str|None, str|None],
        strs_classes: tuple[str|None, str|None, str|None],
    ):
        cls._decide_whether_using_names(strs_types, strs_classes)
        fftypes = cls._get_fftypes(strs_types, strs_classes)
        obj_dihed = cls(isProper, ks, periodicities, phases, *fftypes)
        cls._register_interaction(strs_types, strs_classes, obj_dihed)


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


    # --------------------------------------------------------------------------
    @classmethod
    def get_dihed(cls,
        atom1: nov.FFAtom, atom2: nov.FFAtom,
        atom3: nov.FFAtom, atom4: nov.FFAtom,
        mask: tuple[bool, bool, bool, bool]
    ) -> "nov.FFDihedral | None":
        def _id(a: nov.FFAtom) -> str:
            return a.atom_type.name if cls._USING_NAMES else a.atom_type.atom_class

        key = (
            _id(atom1) if mask[0] else '',
            _id(atom2) if mask[1] else '',
            _id(atom3) if mask[2] else '',
            _id(atom4) if mask[3] else '',
        )
        return cls._get_current_map().get(key)


# //////////////////////////////////////////////////////////////////////////////

import november as nov

# //////////////////////////////////////////////////////////////////////////////
class FFAngle(nov.FFInteraction):
    def __init__(self,
        k: float, angle: float,
        type1: "nov.FFAtomType", type2: "nov.FFAtomType",
        type3: "nov.FFAtomType"
    ):
        self.k = k
        self.angle = angle
        self.type1 = type1
        self.type2 = type2
        self.type3 = type3


    # --------------------------------------------------------------------------
    def __repr__(self):
        return f"FFAngle(k={self.k}, angle={self.angle}, type1={self.type1}, type2={self.type2}, type3={self.type3})"


    # --------------------------------------------------------------------------
    @classmethod
    def register_angle(cls,
        k: float, angle: float,
        strs_types:   tuple[str|None, str|None, str|None],
        strs_classes: tuple[str|None, str|None, str|None],
    ):
        cls._decide_whether_using_names(strs_types, strs_classes)
        fftypes = cls._get_fftypes(strs_types, strs_classes)
        obj_angle = cls(k, angle, *fftypes)
        cls._register_interaction(strs_types, strs_classes, obj_angle)


    # --------------------------------------------------------------------------
    @classmethod
    def get_angle(cls, atom1: nov.FFAtom, atom2: nov.FFAtom, atom3: nov.FFAtom) -> "nov.FFAngle":
        _current_map = cls._get_current_map()
        for key in cls._iter_possible_keys(atom1, atom2, atom3):
            if key in _current_map:
                return _current_map[key]

        raise KeyError(f"{cls.__name__} type not found: {key}")


    # --------------------------------------------------------------------------
    @classmethod
    def _iter_possible_keys(cls, atom1: nov.FFAtom, atom2: nov.FFAtom, atom3: nov.FFAtom):
        """Return possible key tuples for looking up angle parameters."""
        atype1 = atom1.atom_type
        atype2 = atom2.atom_type
        atype3 = atom3.atom_type
        if cls._USING_NAMES:
            yield (atype1.name, atype2.name, atype3.name)
            yield (atype3.name, atype2.name, atype1.name)
        else:
            yield (atype1.atom_class, atype2.atom_class, atype3.atom_class)
            yield (atype3.atom_class, atype2.atom_class, atype1.atom_class)


# //////////////////////////////////////////////////////////////////////////////

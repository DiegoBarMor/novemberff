import november as nov

# //////////////////////////////////////////////////////////////////////////////
class FFBond(nov.FFInteraction):
    def __init__(self,
        k: float, length: float,
        type1: "nov.FFAtomType", type2: "nov.FFAtomType"
    ):
        self.k = k
        self.length = length
        self.type1 = type1
        self.type2 = type2


    # --------------------------------------------------------------------------
    def __repr__(self):
        return f"FFBond(k={self.k}, length={self.length}, type1={self.type1}, type2={self.type2})"


    # --------------------------------------------------------------------------
    @classmethod
    def register_bond(cls,
        k: float, length: float,
        strs_types:   tuple[str|None, str|None],
        strs_classes: tuple[str|None, str|None],
    ):
        cls._decide_whether_using_names(strs_types, strs_classes)
        fftypes = cls._get_fftypes(strs_types, strs_classes)
        obj_bond = cls(k, length, *fftypes)
        cls._register_interaction(strs_types, strs_classes, obj_bond)


    # --------------------------------------------------------------------------
    @classmethod
    def get_bond(cls, atom1: nov.FFAtom, atom2: nov.FFAtom) -> "nov.FFBond":
        _current_map = cls._get_current_map()
        for key in cls._iter_possible_keys(atom1, atom2):
            if key in _current_map:
                return _current_map[key]

        raise KeyError(f"{cls.__name__} type not found: {key}")


    # --------------------------------------------------------------------------
    @classmethod
    def _iter_possible_keys(cls, atom1: nov.FFAtom, atom2: nov.FFAtom):
        """Return possible key tuples for looking up bond parameters."""
        atype1 = atom1.atom_type
        atype2 = atom2.atom_type
        if cls._USING_NAMES:
            yield (atype1.name, atype2.name)
            yield (atype2.name, atype1.name)
        else:
            yield (atype1.atom_class, atype2.atom_class)
            yield (atype2.atom_class, atype1.atom_class)


# //////////////////////////////////////////////////////////////////////////////

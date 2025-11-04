import november as nov

# //////////////////////////////////////////////////////////////////////////////
class FFBond:
    ATOMTYPES_BY_NAME = {}
    ATOMTYPES_BY_CLASS = {}
    USING_NAMES = None # will be set to True/False upon first instantiation

    BONDS = []
    _MAP_BONDS_BY_NAME = {}
    _MAP_BONDS_BY_CLASS = {}

    # --------------------------------------------------------------------------
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
    def register_atomtypes(cls, atomtypes: dict[str, nov.FFAtomType]):
        cls.ATOMTYPES_BY_NAME = atomtypes
        cls.ATOMTYPES_BY_CLASS = {a.atom_class: a for a in atomtypes.values()}


    # --------------------------------------------------------------------------
    @classmethod
    def register_bond(cls,
        k: float, length: float,
        types:   tuple[str|None, str|None],
        classes: tuple[str|None, str|None],
    ):
        t1,t2 = types
        c1,c2 = classes

        if FFBond.USING_NAMES is None:
            FFBond.USING_NAMES = t1 is not None

        if FFBond.USING_NAMES:
            type1 = cls.ATOMTYPES_BY_NAME[t1]
            type2 = cls.ATOMTYPES_BY_NAME[t2]
            bond = cls(k, length, type1, type2)
            cls._MAP_BONDS_BY_NAME[(t1, t2)] = bond
        else:
            type1 = cls.ATOMTYPES_BY_CLASS.get(c1)
            type2 = cls.ATOMTYPES_BY_CLASS.get(c2)
            bond = cls(k, length, type1, type2)
            cls._MAP_BONDS_BY_CLASS[(c1, c2)] = bond

        cls.BONDS.append(bond)


    # --------------------------------------------------------------------------
    @classmethod
    def _iter_possible_keys(cls, atom1: nov.FFAtom, atom2: nov.FFAtom):
        """Return possible key tuples for looking up bond parameters."""
        atype1 = atom1.atom_type
        atype2 = atom2.atom_type
        if FFBond.USING_NAMES:
            yield (atype1.name, atype2.name)
            yield (atype2.name, atype1.name)
        else:
            yield (atype1.atom_class, atype2.atom_class)
            yield (atype2.atom_class, atype1.atom_class)


    # --------------------------------------------------------------------------
    @classmethod
    def get_bond(cls, atom1: nov.FFAtom, atom2: nov.FFAtom):
        _current_map = cls._MAP_BONDS_BY_NAME if FFBond.USING_NAMES else cls._MAP_BONDS_BY_CLASS
        for key in cls._iter_possible_keys(atom1, atom2):
            if key in _current_map:
                return _current_map[key]

        raise KeyError(f"Bond type not found: {key}")


# //////////////////////////////////////////////////////////////////////////////

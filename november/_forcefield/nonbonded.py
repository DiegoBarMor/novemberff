import november as nov

### [AS DEFINED IN] platforms/reference/include/SimTKOpenMMRealType.h
M_PI = 3.14159265358979323846
E_CHARGE = (1.602176634e-19)
AVOGADRO = (6.0221367e23)
EPSILON0 = (1e-6*8.8541878128e-12/(E_CHARGE*E_CHARGE*AVOGADRO)) # (e^2 Na/(kJ nm)) == (e^2/(kJ mol nm))


# //////////////////////////////////////////////////////////////////////////////
class FFNonBonded(nov.FFInteraction):
    ONE_4PI_EPS0 = 1/(4*M_PI*EPSILON0) # 138.935456

    # --------------------------------------------------------------------------
    def __init__(self,
        epsilon: float, sigma: float, charge: float | None,
        atom_type: "nov.FFAtomType"
    ):
        self.epsilon = epsilon
        self.sigma = sigma
        self.charge = charge
        self.type1 = atom_type


    # --------------------------------------------------------------------------
    def __repr__(self):
        return f"FFNonBonded(epsilon={self.epsilon}, sigma={self.sigma}, type={self.type1})"


    # --------------------------------------------------------------------------
    @classmethod
    def register_nonbonded(cls,
        epsilon: float, sigma: float, charge: float | None,
        str_type: str|None,
        str_class: str|None,
    ):
        strs_types   = (str_type,) # nov.FFInteraction expects tuples
        strs_classes = (str_class,)
        cls._decide_whether_using_names(strs_types, strs_classes)
        fftypes = cls._get_fftypes(strs_types, strs_classes)
        obj_nb = cls(epsilon, sigma, charge, *fftypes)
        cls._register_interaction(strs_types, strs_classes, obj_nb)


    # --------------------------------------------------------------------------
    @classmethod
    def get_nonbonded(cls, atom: nov.FFAtom) -> "nov.FFNonBonded":
        _current_map = cls._get_current_map()
        key = (
            atom.atom_type.name if cls._USING_NAMES else atom.atom_type.atom_class,
        ) # nov.FFInteraction expects tuples
        if key in _current_map:
            return _current_map[key]

        raise KeyError(f"{cls.__name__} type not found: {key}")


# //////////////////////////////////////////////////////////////////////////////

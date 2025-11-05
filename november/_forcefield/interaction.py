from abc import ABC

import november as nov

# //////////////////////////////////////////////////////////////////////////////
class FFInteraction(ABC):
    _ATOMTYPES_BY_NAME  = {}
    _ATOMTYPES_BY_CLASS = {}
    _USING_NAMES = None # will be set to True/False upon first interaction registered
    _MAP_INTERACTION_BY_NAME  = {}
    _MAP_INTERACTION_BY_CLASS = {}

    # --------------------------------------------------------------------------
    @classmethod
    def register_atomtypes(cls, atomtypes: dict[str, nov.FFAtomType]):
        cls._ATOMTYPES_BY_NAME = atomtypes
        cls._ATOMTYPES_BY_CLASS = {a.atom_class: a for a in atomtypes.values()}


    # --------------------------------------------------------------------------
    @classmethod
    def _decide_whether_using_names(cls,
        strs_types: tuple[str|None,...], strs_classes: tuple[str|None,...]
    ):
        if cls._USING_NAMES is not None: return

        if not((strs_types[0] is None)^(strs_classes[0] is None)):
            raise ValueError("Forcefield must consistently provide either atom types or atom classes for its interactions.")

        cls._USING_NAMES = strs_classes[0] is None


    # --------------------------------------------------------------------------
    @classmethod
    def _register_interaction(cls,
        strs_types: tuple[str|None,...], strs_classes: tuple[str|None,...],
        obj_interaction: "FFInteraction"
    ):
        if cls._USING_NAMES:
            cls._MAP_INTERACTION_BY_NAME[strs_types] = obj_interaction
            return
        cls._MAP_INTERACTION_BY_CLASS[strs_classes] = obj_interaction


    # --------------------------------------------------------------------------
    @classmethod
    def _get_current_map(cls):
        return cls._MAP_INTERACTION_BY_NAME if cls._USING_NAMES else cls._MAP_INTERACTION_BY_CLASS


    # --------------------------------------------------------------------------
    @classmethod
    def _get_fftypes(cls,
        str_types: tuple[str,...], str_classes: tuple[str,...]
    ) -> list["nov.FFAtomType"]:
        return [
            cls._ATOMTYPES_BY_NAME.get(t) for t in str_types
        ] if cls._USING_NAMES else [
            cls._ATOMTYPES_BY_CLASS.get(c) for c in str_classes
        ]


# //////////////////////////////////////////////////////////////////////////////

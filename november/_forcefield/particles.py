# //////////////////////////////////////////////////////////////////////////////
class FFAtomType:
    def __init__(self, atom_class: str, element: str, mass: float, name: str):
        self.atom_class = atom_class
        self.element = element
        self.mass = mass
        self.name = name

    def __repr__(self): return f"<FFAtomType {self.atom_class}>"


# //////////////////////////////////////////////////////////////////////////////
class FFAtom:
    def __init__(self, name: str, charge: float, atom_type: FFAtomType):
        self.name = name
        self.charge = charge
        self.atom_type = atom_type

    def __repr__(self):
        return f"FFAtom(name={self.name}, charge={self.charge}, type={self.atom_type})"


# //////////////////////////////////////////////////////////////////////////////
class FFResidue:
    def __init__(self, name: str, atoms: dict[str, FFAtom]):
        self.name = name
        self._atoms = atoms


    # --------------------------------------------------------------------------
    def __repr__(self):
        return f"FFResidue(name={self.name}, atoms={' '.join(self._atoms.keys())})"


    # --------------------------------------------------------------------------
    def get_atom_by_name(self, name: str) -> FFAtom | None:
        return self._atoms.get(name)


# //////////////////////////////////////////////////////////////////////////////

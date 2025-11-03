import november as nov

# //////////////////////////////////////////////////////////////////////////////
class FFBond:
    ATOMTYPES = {}

    # --------------------------------------------------------------------------
    def __init__(self,
        k: float, length: float, types: tuple[str, str]
    ):
        self.k = k
        self.length = length
        self.type1 = self.ATOMTYPES[types[0]]
        self.type2 = self.ATOMTYPES[types[1]]


    # --------------------------------------------------------------------------
    def __repr__(self):
        return f"FFBond(k={self.k}, length={self.length}, type1={self.type1}, type2={self.type2})"


    # --------------------------------------------------------------------------
    @classmethod
    def register_atomtypes(cls, atomtypes: dict[str, nov.FFAtomType]):
        cls.ATOMTYPES = atomtypes


# //////////////////////////////////////////////////////////////////////////////

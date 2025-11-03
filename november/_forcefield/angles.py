import november as nov

# //////////////////////////////////////////////////////////////////////////////
class FFAngle:
    ATOMTYPES = {}

    # --------------------------------------------------------------------------
    def __init__(self,
        k: float, angle: float, types: tuple[str, str, str]
    ):
        self.k = k
        self.angle = angle
        self.type1 = self.ATOMTYPES[types[0]]
        self.type2 = self.ATOMTYPES[types[1]]
        self.type3 = self.ATOMTYPES[types[2]]


    # --------------------------------------------------------------------------
    def __repr__(self):
        return f"FFAngle(k={self.k}, angle={self.angle}, type1={self.type1}, type2={self.type2}, type3={self.type3})"


    # --------------------------------------------------------------------------
    @classmethod
    def register_atomtypes(cls, atomtypes: dict[str, nov.FFAtomType]):
        cls.ATOMTYPES = atomtypes


# //////////////////////////////////////////////////////////////////////////////

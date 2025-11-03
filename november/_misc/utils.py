# //////////////////////////////////////////////////////////////////////////////
class Utils:
    # --------------------------------------------------------------------------
    @staticmethod
    def combinations_proper_diheds(a0, a1, a2, a3):
        combinations = [
            (a0, a1, a2, a3),
            (a3, a2, a1, a0)
        ]
        masks = [
            (True,  True,  True,  True ),
            (True,  False, False, True ),
            (False, True,  True,  False),
        ]
        for mask in masks:
            for combo in combinations:
                yield combo, mask


    # --------------------------------------------------------------------------
    @staticmethod
    def combinations_improper_diheds(a0, a1, a2, a3):
        combinations = [
            (a0, a1, a2, a3),
            (a0, a1, a3, a2),
            (a2, a1, a0, a3),
            (a3, a1, a0, a2),
            (a2, a1, a3, a0),
            (a3, a1, a2, a0),
        ]
        masks = [
            (True,  True,  True,  True ),
            (True,  False, False, True ),
            (True,  False, True, True ),
        ]
        for mask in masks:
            for combo in combinations:
                yield combo, mask


# //////////////////////////////////////////////////////////////////////////////

import november as nov

# //////////////////////////////////////////////////////////////////////////////
class ForceField:
    def __init__(self, ff_data: dict):
        self._ffresidues  = ff_data["residues"]
        self._ffangles    = ff_data["angles"]
        self._FFDihedrals    = ff_data["diheds"]
        self._ffnonbonded = ff_data["nonbonded"]

        self.one_4pi_eps0 = nov.FFNonBonded.ONE_4PI_EPS0
        self.lj_14scale = ff_data["LJ14SCALE"]
        self.coulomb_14scale = ff_data["COULOMB14SCALE"]

    # --------------------------------------------------------------------------
    @classmethod
    def from_xml(cls, path_xml):
        xml = nov.FFParserXML(path_xml)
        xml.parse()

        _safe_int   = lambda s: None if s is None else int(s)
        _safe_float = lambda s: None if s is None else float(s)

        node_ff        = xml.root.get_child_by_name("ForceField")
        node_types     = node_ff.get_child_by_name("AtomTypes")
        node_resids    = node_ff.get_child_by_name("Residues")
        node_bonds     = node_ff.get_child_by_name("HarmonicBondForce")
        node_angles    = node_ff.get_child_by_name("HarmonicAngleForce")
        node_diheds    = node_ff.get_child_by_name("PeriodicTorsionForce")
        node_nonbonded = node_ff.get_child_by_name("NonbondedForce")

        ff_data = {
            "types": {}, "residues": {},
            "angles": {}, "diheds": {}, "nonbonded": {},
            "LJ14SCALE":      float(node_nonbonded.get_attr("lj14scale")),
            "COULOMB14SCALE": float(node_nonbonded.get_attr("coulomb14scale")),
        }

        ffatomtypes = {}
        for child in node_types.children:
            name = child.get_attr("name")
            ffatomtypes[name] = nov.FFAtomType(
                atom_class = child.get_attr("class"),
                element    = child.get_attr("element"),
                mass       = float(child.get_attr("mass")),
                name       = name,
            )
        ff_data["types"] = ffatomtypes
        nov.FFBond.register_atomtypes(ffatomtypes)
        nov.FFAngle.register_atomtypes(ffatomtypes)
        nov.FFDihedral.register_atomtypes(ffatomtypes)


        for residue_node in node_resids.children:
            resname = residue_node.get_attr("name")
            atom_nodes = filter(lambda c: c.tag_name == "Atom", residue_node.children)
            atoms = {
                node.get_attr("name") : nov.FFAtom(
                    name      = node.get_attr("name"),
                    charge    = _safe_float(node.get_attr("charge")),
                    atom_type = ffatomtypes[node.get_attr("type")],
                ) for node in atom_nodes
            }
            ff_data["residues"][resname] = nov.FFResidue(resname, atoms)


        for child in node_bonds.children:
            types = (
                child.get_attr("type1"),
                child.get_attr("type2"),
            )
            classes = (
                child.get_attr("class1"),
                child.get_attr("class2"),
            )
            nov.FFBond.register_bond(
                k      = float(child.get_attr("k")),
                length = float(child.get_attr("length")),
                types  = types, classes = classes,
            )

        return cls(ff_data) # [WIP]


        for child in node_angles.children:
            types = (
                child.get_attr("type1"),
                child.get_attr("type2"),
                child.get_attr("type3"),
            )
            classes = (
                child.get_attr("class1"),
                child.get_attr("class2"),
                child.get_attr("class3"),
            )
            ff_data["angles"][types] = nov.FFAngle(
                k     = float(child.get_attr("k")),
                angle = float(child.get_attr("angle")),
                types = types, classes = classes,
            )


        for child in node_diheds.children:
            types = (
                child.get_attr("type1"),
                child.get_attr("type2"),
                child.get_attr("type3"),
                child.get_attr("type4"),
            )
            classes = (
                child.get_attr("class1"),
                child.get_attr("class2"),
                child.get_attr("class3"),
                child.get_attr("class4"),
            )
            ff_data["diheds"][types] = nov.FFDihedral(
                isProper = child.tag_name == "Proper",
                k1 = _safe_float(child.get_attr("k1")),
                k2 = _safe_float(child.get_attr("k2")),
                k3 = _safe_float(child.get_attr("k3")),
                k4 = _safe_float(child.get_attr("k4")),
                periodicity1 = _safe_int(child.get_attr("periodicity1")),
                periodicity2 = _safe_int(child.get_attr("periodicity2")),
                periodicity3 = _safe_int(child.get_attr("periodicity3")),
                periodicity4 = _safe_int(child.get_attr("periodicity4")),
                phase1 = _safe_float(child.get_attr("phase1")),
                phase2 = _safe_float(child.get_attr("phase2")),
                phase3 = _safe_float(child.get_attr("phase3")),
                phase4 = _safe_float(child.get_attr("phase4")),
                types = types, classes = classes,
            )


        for child in node_nonbonded.children:
            if child.tag_name != "Atom": continue
            t1 = child.get_attr("type")
            ff_data["nonbonded"][t1] = nov.FFNonBonded(
                epsilon   = float(child.get_attr("epsilon")),
                sigma     = float(child.get_attr("sigma")),
                atom_type = ffatomtypes[t1],
            )


        return cls(ff_data)


    # --------------------------------------------------------------------------
    def omm2ff(self, atom) -> nov.FFAtom:
        """
        Convert an OpenMM atom to a FF atom.
        """
        ff_residue: nov.FFResidue = self._ffresidues[atom.residue.resname]
        return ff_residue.get_atom_by_name(atom.name)


    # --------------------------------------------------------------------------
    def get_ffbond(self, atom0, atom1) -> nov.FFBond:
        ff_a0 = self.omm2ff(atom0)
        ff_a1 = self.omm2ff(atom1)
        return nov.FFBond.get_bond(ff_a0, ff_a1)


    # --------------------------------------------------------------------------
    def get_ffangle(self, atom0, atom1, atom2) -> nov.FFAngle:
        ff_a0 = self.omm2ff(atom0)
        ff_a1 = self.omm2ff(atom1)
        ff_a2 = self.omm2ff(atom2)

        key = (ff_a0.atom_type.name, ff_a1.atom_type.name, ff_a2.atom_type.name)
        if key in self._ffangles.keys(): return self._ffangles[key]

        key = (ff_a2.atom_type.name, ff_a1.atom_type.name, ff_a0.atom_type.name)
        if key in self._ffangles.keys(): return self._ffangles[key]

        raise KeyError(f"Angle type not found: {ff_a0.atom_type.name}-{ff_a1.atom_type.name}-{ff_a2.atom_type.name}")


    # --------------------------------------------------------------------------
    def iter_ffpropers(self, atoms):
        for ordered_atoms,mask in nov.Utils.combinations_proper_diheds(*atoms):
            a0, a1, a2, a3 = ordered_atoms
            ff_a0 = self.omm2ff(a0)
            ff_a1 = self.omm2ff(a1)
            ff_a2 = self.omm2ff(a2)
            ff_a3 = self.omm2ff(a3)

            key = (
                ff_a0.atom_type.name if mask[0] else '',
                ff_a1.atom_type.name if mask[1] else '',
                ff_a2.atom_type.name if mask[2] else '',
                ff_a3.atom_type.name if mask[3] else '',
            )
            if key in self._FFDihedrals.keys():
                ordered_idxs = (a0.index, a1.index, a2.index, a3.index)
                yield self._FFDihedrals[key], ordered_idxs


    # --------------------------------------------------------------------------
    def iter_ffimpropers(self, atoms):
        ### improper dihedrals in amber should have "a2" at the center (in theory)
        ### https://ambermd.org/doc12/Amber24.pdf#page=279
        ### but in practice, for openmm they have "a1" at the center
        ### graph_bonds.py reports "a0" in the center, so switch a0 and a1

        for ordered_atoms,mask in nov.Utils.combinations_improper_diheds(*atoms):
            a0, a1, a2, a3 = ordered_atoms
            ff_a0 = self.omm2ff(a1) # swap between 0 and 1 is intentional (see explanation above)
            ff_a1 = self.omm2ff(a0)
            ff_a2 = self.omm2ff(a2)
            ff_a3 = self.omm2ff(a3)

            key = (
                ff_a0.atom_type.name if mask[0] else '',
                ff_a1.atom_type.name if mask[1] else '',
                ff_a2.atom_type.name if mask[2] else '',
                ff_a3.atom_type.name if mask[3] else '',
            )
            if key in self._FFDihedrals.keys():
                ordered_idxs = (a0.index, a1.index, a2.index, a3.index)
                yield self._FFDihedrals[key], ordered_idxs


    # --------------------------------------------------------------------------
    def get_ffnonbonded(self, a0, a1):
            ff_a0 = self.omm2ff(a0)
            ff_a1 = self.omm2ff(a1)
            ff_nb0 = self._ffnonbonded[ff_a0.atom_type.name]
            ff_nb1 = self._ffnonbonded[ff_a1.atom_type.name]
            charge0 = ff_a0.charge
            charge1 = ff_a1.charge
            return ff_nb0, ff_nb1, charge0, charge1


# //////////////////////////////////////////////////////////////////////////////

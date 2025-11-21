import novemberff as nov

# //////////////////////////////////////////////////////////////////////////////
class ForceField:
    def __init__(self, ff_data: dict):
        self._ffresidues     = ff_data["residues"]
        self.lj_14scale      = ff_data["LJ14SCALE"]
        self.coulomb_14scale = ff_data["COULOMB14SCALE"]
        self.one_4pi_eps0    = nov.FFNonBonded.ONE_4PI_EPS0

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
        nov.FFNonBonded.register_atomtypes(ffatomtypes)

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
            strs_types = (
                child.get_attr("type1"),
                child.get_attr("type2"),
            )
            strs_classes = (
                child.get_attr("class1"),
                child.get_attr("class2"),
            )
            nov.FFBond.register_bond(
                k      = float(child.get_attr("k")),
                length = float(child.get_attr("length")),
                strs_types = strs_types, strs_classes = strs_classes,
            )

        for child in node_angles.children:
            strs_types = (
                child.get_attr("type1"),
                child.get_attr("type2"),
                child.get_attr("type3"),
            )
            strs_classes = (
                child.get_attr("class1"),
                child.get_attr("class2"),
                child.get_attr("class3"),
            )
            nov.FFAngle.register_angle(
                k     = float(child.get_attr("k")),
                angle = float(child.get_attr("angle")),
                strs_types = strs_types, strs_classes = strs_classes,
            )

        for child in node_diheds.children:
            strs_types = (
                child.get_attr("type1"),
                child.get_attr("type2"),
                child.get_attr("type3"),
                child.get_attr("type4"),
            )
            strs_classes = (
                child.get_attr("class1"),
                child.get_attr("class2"),
                child.get_attr("class3"),
                child.get_attr("class4"),
            )
            nov.FFDihedral.register_dihed(
                isProper = child.tag_name == "Proper",
                ks = [
                    _safe_float(child.get_attr(f"k{i}")) for i in range(1,5)
                ],
                periodicities = [
                    _safe_int(child.get_attr(f"periodicity{i}")) for i in range(1,5)
                ],
                phases = [
                    _safe_float(child.get_attr(f"phase{i}")) for i in range(1,5)
                ],
                strs_types = strs_types, strs_classes = strs_classes,
            )

        for child in node_nonbonded.children:
            if child.tag_name != "Atom": continue
            str_type  = child.get_attr("type")
            str_class = child.get_attr("class")
            nov.FFNonBonded.register_nonbonded(
                epsilon  = float(child.get_attr("epsilon")),
                sigma    = float(child.get_attr("sigma")),
                charge   = _safe_float(child.get_attr("charge")),
                str_type = str_type, str_class = str_class,
            )

        return cls(ff_data)


    # --------------------------------------------------------------------------
    def map_mda2ff(self, atom) -> nov.FFAtom:
        """
        Maps an MDA atom to its correspondant FF atom.
        """
        ff_residue: nov.FFResidue = self._ffresidues[atom.residue.resname]
        ff_atom: nov.FFAtom = ff_residue.get_atom_by_name(atom.name)
        if ff_atom is None:
            raise ValueError(f"Atom '{atom.name}' not found in residue '{atom.residue.resname} (resid: {atom.residue.resid})'")
        return ff_atom


    # --------------------------------------------------------------------------
    def get_ffbond(self, atom0, atom1) -> nov.FFBond:
        ff_a0 = self.map_mda2ff(atom0)
        ff_a1 = self.map_mda2ff(atom1)
        try:
            return nov.FFBond.get_bond(ff_a0, ff_a1)
        except KeyError as e:
            self._err_unkwnown_interaction(e, atom0, atom1)


    # --------------------------------------------------------------------------
    def get_ffangle(self, atom0, atom1, atom2) -> nov.FFAngle:
        ff_a0 = self.map_mda2ff(atom0)
        ff_a1 = self.map_mda2ff(atom1)
        ff_a2 = self.map_mda2ff(atom2)
        try:
            return nov.FFAngle.get_angle(ff_a0, ff_a1, ff_a2)
        except KeyError as e:
            self._err_unkwnown_interaction(e, atom0, atom1, atom2)


    # --------------------------------------------------------------------------
    def iter_ffpropers(self, atoms):
        for ordered_atoms,mask in nov.Utils.combinations_proper_diheds(*atoms):
            a0, a1, a2, a3 = ordered_atoms
            ff_a0 = self.map_mda2ff(a0)
            ff_a1 = self.map_mda2ff(a1)
            ff_a2 = self.map_mda2ff(a2)
            ff_a3 = self.map_mda2ff(a3)
            dihed = nov.FFDihedral.get_dihed(ff_a0, ff_a1, ff_a2, ff_a3, mask)
            if dihed is not None:
                ordered_idxs = (a0.index, a1.index, a2.index, a3.index)
                yield dihed, ordered_idxs


    # --------------------------------------------------------------------------
    def iter_ffimpropers(self, atoms):
        ### improper dihedrals in amber should have "a2" at the center (in theory)
        ### https://ambermd.org/doc12/Amber24.pdf#page=279
        ### but in practice, for openmm they have "a1" at the center
        ### graph_bonds.py reports "a0" in the center, so switch a0 and a1

        for ordered_atoms,mask in nov.Utils.combinations_improper_diheds(*atoms):
            a0, a1, a2, a3 = ordered_atoms
            ff_a0 = self.map_mda2ff(a1) # swap between 0 and 1 is intentional (see explanation above)
            ff_a1 = self.map_mda2ff(a0)
            ff_a2 = self.map_mda2ff(a2)
            ff_a3 = self.map_mda2ff(a3)
            dihed = nov.FFDihedral.get_dihed(ff_a0, ff_a1, ff_a2, ff_a3, mask)
            if dihed is not None:
                ordered_idxs = (a0.index, a1.index, a2.index, a3.index)
                yield dihed, ordered_idxs


    # --------------------------------------------------------------------------
    def get_ffnonbonded(self, a0, a1):
            ff_a0 = self.map_mda2ff(a0)
            ff_a1 = self.map_mda2ff(a1)
            try:
                ff_nb0 = nov.FFNonBonded.get_nonbonded(ff_a0)
                ff_nb1 = nov.FFNonBonded.get_nonbonded(ff_a1)
            except KeyError as e:
                self._err_unkwnown_interaction(e, a0, a1)

            ### charges are sometimes stored in the Residue.Atom XML-nodes (e.g. in RNA.OL3.xml)
            ### other times, they are stored in the Nonbonded XML-nodes (e.g. in amber99sb.xml)
            charge0 = ff_a0.charge if ff_nb0.charge is None else ff_nb0.charge
            charge1 = ff_a1.charge if ff_nb1.charge is None else ff_nb1.charge
            return ff_nb0, ff_nb1, charge0, charge1


    # --------------------------------------------------------------------------
    @staticmethod
    def _err_unkwnown_interaction(exception, *atoms):
        names    = (atom.name for atom in atoms)
        indices  = (str(atom.index) for atom in atoms)
        resnames = (atom.residue.resname for atom in atoms)
        resids   = (str(atom.residue.resid) for atom in atoms)
        raise KeyError(
            f"Unknown interaction between atoms '{'-'.join(names)}' (index: {'-'.join(indices)}) " +
            f"in residues '{'-'.join(resnames)}' (resid: {'-'.join(resids)})'"
        ) from exception


# //////////////////////////////////////////////////////////////////////////////

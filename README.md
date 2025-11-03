# NOVEMBER
**November** (implemeNtatiOn of VErbose aMBER) is a custom implementation of the AMBER forcefield, inspired by [OpenMM](https://github.com/openmm/openmm)'s C++ source. November's aim is not to be performant, but rather to be compact and easy to modify/customize. It also provides by default the option to obtain verbose outputs, i.e. the energies of every singular molecular interaction in arrays (instead of their sum).

* Current outputs supported:
    * Bonded energies
    * Angle energies
    * Dihedral (Torsion) energies
    * LennardJones energies
    * Coulomb energies

* Current forcefield formats supported:
    * XML

* Current forcefields packed along November (at `november/_data`):
    * `RNA.OL3.xml`

<!-- ----------------------------------------------------------------------- -->
# Examples
* display energy sums:
```python3 main.py rna  testdata/1ato.pdb```
```python3 main.py prot testdata/prot.pdb```

* save energy arrays:
```python3 main.py rna  testdata/1ato.pdb testdata/output/1ato```
```python3 main.py prot testdata/prot.pdb testdata/output/prot```


<!-- ----------------------------------------------------------------------- -->
# TODO
* improve CLI.
* add documentation.
* add more forcefields.
* remove OpenMM dependency for PDB parsing.
* improve `bond_graphs`'s graph data structure.
* finish refactoring of parser base class and `forcefield.from_xml` (i.e. should be Parser's responsability to yield the needed data, not Forcefield).
* allow to specify the forcefield to use in `main.py`, either by path or by name.
* add option for saving the molecule's topology features (bonds, angles, etc), as well as their geometry values.

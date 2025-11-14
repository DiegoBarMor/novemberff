# Changelog

## [0.1.0] - 2021-11-10
- Initial upload to PyPI of a preliminary **Novemberff**.

## [0.2.0] - 2021-11-10
- Added possibility to initialize the calculator with a default forcefield name, instead of providing a path to an XML file.

## [0.3.0] - 2021-11-10
- Changed the API of the energy getters.

## [0.4.0] - 2021-11-13
- Fixed missing `package_data` in `setup.py`.
- Added method to `EnergyCalculator` for saving metadata into a CSV file.

## [0.4.1] - 2021-11-14
- Renamed `Forcefield.omm2ff` into `Forcefield.map_mda2ff` and made it assert that the mapped atom is known

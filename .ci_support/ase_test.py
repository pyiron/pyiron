import sys
import pyiron.atomistics.structure.atoms as atoms
import pyiron.atomistics.structure.atom as atom
import pytest

sys.modules["ase.atom"] = atom
sys.modules["ase.atoms"] = atoms
pytest.main(["ase/ase/test/cell"])

import sys
import pyiron.atomistics.structure.atoms as atoms
import pyiron.atomistics.structure.atom as atom
import pytest

import ase

print(ase.__file__)

sys.modules["ase.atom"] = atom
sys.modules["ase.atoms"] = atoms
pytest.main(["ase/ase/test/cell"])

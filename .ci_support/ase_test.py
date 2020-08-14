import sys
import pyiron.atomistics.structure.atoms as atoms
import pyiron.atomistics.structure.atom as atom
import pytest

import ase

test_directory = "/".join(ase.__file__.split("/")[0:-1]) + "/test/"

sys.modules["ase.atom"] = atom
sys.modules["ase.atoms"] = atoms
out = pytest.main([test_directory])

sys.exit(out.value)

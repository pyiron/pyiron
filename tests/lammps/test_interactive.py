import unittest
import numpy as np
import os
from pyiron.base.project.generic import Project
from pyiron.atomistics.structure.atoms import Atoms
from pyiron.base.generic.hdfio import ProjectHDFio
from pyiron.lammps.lammps import Lammps

class InteractiveLibrary(object):
    def __init__(self):
        self._command = []

    def command(self, command_in):
        self._command.append(command_in)

    def scatter_atoms(self, *args):
        self._command.append(' '.join([str(arg) for arg in args]))



class TestLammps(unittest.TestCase):
    def setUp(self):
        self.job._interactive_library = InteractiveLibrary()

    @classmethod
    def setUpClass(cls):
        cls.execution_path = os.path.dirname(os.path.abspath(__file__))
        cls.project = Project(os.path.join(cls.execution_path, 'lammps'))
        cls.job = Lammps(project=ProjectHDFio(project=cls.project, file_name='lammps'), job_name='lammps')
        cls.job.server.run_mode.interactive = True
        cls.job.structure = Atoms(symbols='Fe2', positions=np.outer(np.arange(2), np.ones(3)), cell=2*np.eye(3))

    @classmethod
    def tearDownClass(cls):
        cls.execution_path = os.path.dirname(os.path.abspath(__file__))
        project = Project(os.path.join(cls.execution_path, 'lammps'))
        project.remove_jobs(recursive=True)
        project.remove(enable=True)

    def test_interactive_cells_setter(self):
        self.job.interactive_cells_setter(np.eye(3))
        self.assertEqual(self.job._interactive_library._command[-1],
                         'change_box all x final 0 1.000000 y final 0 1.000000 z final 0 1.000000 units box')

    def test_interactive_positions_setter(self):
        self.job.interactive_positions_setter(np.arange(6).reshape(2,3))
        self.assertTrue(self.job._interactive_library._command[0].startswith('x 1 3'))
        self.assertEqual(self.job._interactive_library._command[1], 'change_box all remap')

    def test_interactive_execute(self):
        self.job._interactive_lammps_input()
        self.assertEqual(self.job._interactive_library._command,
                         ['fix ensemble all nve',
                          'variable dumptime equal 100',
                          'thermo_style custom step temp pe etotal pxx pxy pxz pyy pyz pzz vol',
                          'thermo_modify format float %20.15g',
                          'thermo 100'])


if __name__ == '__main__':
    unittest.main()

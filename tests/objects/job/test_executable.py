import unittest
from pyiron_base.objects.job.executable import Executable


# class TestExecutable(unittest.TestCase):
#     def setUp(self):
#         self.exe_no_existing = Executable(codename='GenericJob', path_binary_codes=['../../../../static/bin_cmmc'],
#                                           overwrite_nt_flag=True)
#         self.exe_default = Executable(codename='Vasp', path_binary_codes=['../../../../static/bin_cmmc'],
#                                       overwrite_nt_flag=True)
#         self.exe_default_mpi = Executable(codename='Vasp', path_binary_codes=['../../../../static/bin_cmmc'],
#                                           overwrite_nt_flag=True)
#         self.exe_default_mpi.mpi = True
#         self.exe_string = Executable(codename='Vasp', path_binary_codes=['../../../../static/bin_cmmc'],
#                                      overwrite_nt_flag=True)
#         self.exe_string.executable_path = '../../../../static/bin_cmmc/vasp/run_vasp_5.3_mpi.sh'
#         self.exe_version_switch = Executable(codename='Vasp', path_binary_codes=['../../../../static/bin_cmmc'],
#                                              overwrite_nt_flag=True)
#         self.exe_version_switch.version = '5.4'
#         self.exe_version_switch_mpi = Executable(codename='Vasp', path_binary_codes=['../../../../static/bin_cmmc'],
#                                                  overwrite_nt_flag=True)
#         self.exe_version_switch_mpi.version = '5.4'
#         self.exe_version_switch_mpi.mpi = True
#
#     def test_version(self):
#         self.assertFalse(self.exe_no_existing.version)
#         self.assertEqual(self.exe_default.version, '5.3')
#         self.assertEqual(self.exe_default_mpi.version, '5.3_mpi')
#         self.assertEqual(self.exe_string.version, '../../../../static/bin_cmmc/vasp/run_vasp_5.3_mpi.sh')
#         self.assertEqual(self.exe_version_switch.version, '5.4')
#         self.assertEqual(self.exe_version_switch_mpi.version, '5.4_mpi')
#
#     def test_mpi(self):
#         self.assertFalse(self.exe_no_existing.mpi)
#         self.assertFalse(self.exe_default.mpi)
#         self.assertTrue(self.exe_default_mpi.mpi)
#         self.assertTrue(self.exe_string.mpi)
#         self.assertFalse(self.exe_version_switch.mpi)
#         self.assertTrue(self.exe_version_switch_mpi.mpi)
#
#     def test_available_versions(self):
#         self.assertFalse(self.exe_no_existing.available_versions)
#         self.assertEqual(self.exe_default.available_versions,
#                          ['5.3', '5.3_col_mpi', '5.3_mpi', '5.4', '5.4.4', '5.4.4_mpi', '5.4_gamma', '5.4_gamma_mpi', '5.4_mpi'])
#         self.assertEqual(self.exe_default_mpi.available_versions,
#                          ['5.3', '5.3_col_mpi', '5.3_mpi', '5.4', '5.4.4', '5.4.4_mpi', '5.4_gamma', '5.4_gamma_mpi', '5.4_mpi'])
#         self.assertEqual(self.exe_string.available_versions,
#                          ['5.3', '5.3_col_mpi', '5.3_mpi', '5.4', '5.4.4', '5.4.4_mpi', '5.4_gamma', '5.4_gamma_mpi', '5.4_mpi'])
#         self.assertEqual(self.exe_version_switch.available_versions,
#                          ['5.3', '5.3_col_mpi', '5.3_mpi', '5.4', '5.4.4', '5.4.4_mpi', '5.4_gamma', '5.4_gamma_mpi', '5.4_mpi'])
#         self.assertEqual(self.exe_version_switch_mpi.available_versions,
#                          ['5.3', '5.3_col_mpi', '5.3_mpi', '5.4', '5.4.4', '5.4.4_mpi', '5.4_gamma', '5.4_gamma_mpi', '5.4_mpi'])
#
#     def test_executable_path(self):
#         self.assertFalse(self.exe_no_existing.executable_path)
#         self.assertEqual(self.exe_default.executable_path,
#                          '../../../../static/bin_cmmc/vasp/run_vasp_5.3.sh')
#         self.assertEqual(self.exe_default_mpi.executable_path,
#                          '../../../../static/bin_cmmc/vasp/run_vasp_5.3_mpi.sh')
#         self.assertEqual(self.exe_string.executable_path,
#                          '../../../../static/bin_cmmc/vasp/run_vasp_5.3_mpi.sh')
#         self.assertEqual(self.exe_version_switch.executable_path,
#                          '../../../../static/bin_cmmc/vasp/run_vasp_5.4.sh')
#         self.assertEqual(self.exe_version_switch_mpi.executable_path,
#                          '../../../../static/bin_cmmc/vasp/run_vasp_5.4_mpi.sh')


if __name__ == '__main__':
    unittest.main()
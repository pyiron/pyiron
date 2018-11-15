import unittest


class TestCommand(unittest.TestCase):
    def test_run(self):
        pass


class TestSshHost(unittest.TestCase):
    def test_path(self):
        pass

    def test_plink(self):
        pass

    def test_run_remote(self):
        pass

    def test_run_python_command(self):
        pass

    def test_run_python_script(self):
        pass

    def test_is_active(self):
        pass

    def test_ping(self):
        pass

    def test_remove_host_directory(self):
        pass

    def test_cp_local_dir_to_host(self):
        pass

    def test_cp_local_file_to_host(self):
        pass


if __name__ == '__main__':
    unittest.main()
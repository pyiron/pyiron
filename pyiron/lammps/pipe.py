from ctypes import c_double, c_int
from multiprocessing import Process, Pipe
import numpy as np

try:
    from lammps import lammps
except ImportError:
    pass


class LammpsLibrary(object):
    def __init__(self):
        lmp_interface = lammps()
        parent_conn, child_conn = Pipe()
        lammps_process = Process(target=self.interactive_run, args=(child_conn, lmp_interface))
        lammps_process.start()
        self._interactive_library = parent_conn

    def command(self, command):
        self._interactive_library.send([self.interactive_lib_command, command])

    def gather_atoms(self, *args):
        self._interactive_library.send([self.interative_gather_atoms] + list(args))
        # self._interactive_library.send([self.interative_gather_atoms, *args])  # Python 3.X only
        return self._interactive_library.recv()

    def scatter_atoms(self, *args):
        self._interactive_library.send([self.interactive_scatter_atoms] + list(args))
        # self._interactive_library.send([self.interactive_scatter_atoms, *args])  # Python 3.X only

    def get_thermo(self, *args):
        self._interactive_library.send([self.interactive_get_thermo] + list(args))
        # self._interactive_library.send([self.interactive_get_thermo, *args])  # Python 3.X only
        return self._interactive_library.recv()

    def extract_compute(self, *args):
        self._interactive_library.send([self.interactive_extract_compute] + list(args))
        # self._interactive_library.send([self.interactive_extract_compute, *args])  # Python 3.X only
        return self._interactive_library.recv()

    def close(self):
        self._interactive_library.send([self.interactive_close])

    @staticmethod
    def interactive_lib_command(conn, job, command):
        job.command(command)

    @staticmethod
    def interative_gather_atoms(conn, job, funct_args):
        return np.array(job.gather_atoms(*funct_args))

    @staticmethod
    def interactive_scatter_atoms(conn, job, funct_args):
        py_vector = funct_args[3]
        if issubclass(type(py_vector[0]), np.integer):
            c_vector = (len(py_vector) * c_int)(*py_vector)
        else:
            c_vector = (len(py_vector) * c_double)(*py_vector)
        job.scatter_atoms(funct_args[0], funct_args[1], funct_args[2], c_vector)

    @staticmethod
    def interactive_get_thermo(conn, job, funct_args):
        return np.array(job.get_thermo(*funct_args))

    @staticmethod
    def interactive_extract_compute(conn, job, funct_args):
        return np.array(job.extract_compute(*funct_args))

    @staticmethod
    def interactive_close(conn, job):
        job.close()
        conn.close()
        return 'exit'

    @staticmethod
    def interactive_run(conn, job):
        while True:
            input_info = conn.recv()
            if isinstance(input_info, list):
                input_function = input_info[0]
                input_args = input_info[1:]
                answer = input_function(conn=conn, job=job, funct_args=input_args)
            else:
                answer = input_info(conn=conn, job=job)
            if isinstance(answer, str) and answer == 'exit':
                break
            elif answer is not None:
                conn.send(answer)

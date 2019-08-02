from ctypes import c_double, c_int
from mpi4py import MPI
import numpy as np
import pickle
import sys

try:
    from lammps import lammps
except ImportError:
    pass


# Lammps executable
job = lammps(cmdargs=['-screen', 'none'])


def extract_compute(funct_args):
    if MPI.COMM_WORLD.rank == 0:
        return np.array(job.extract_compute(*funct_args))


def get_thermo(funct_args):
    if MPI.COMM_WORLD.rank == 0:
        return np.array(job.get_thermo(*funct_args))


def scatter_atoms(funct_args):
    py_vector = funct_args[3]
    if issubclass(type(py_vector[0]), np.integer):
        c_vector = (len(py_vector) * c_int)(*py_vector)
    else:
        c_vector = (len(py_vector) * c_double)(*py_vector)
    job.scatter_atoms(funct_args[0], funct_args[1], funct_args[2], c_vector)


def command(funct_args):
    job.command(funct_args)


def gather_atoms(funct_args):
    return np.array(job.gather_atoms(*funct_args))


def select_cmd(argument):
    """
    Select a lammps command

    Args:
        argument (str): [close, extract_compute, get_thermo, scatter_atoms, command, gather_atoms]

    Returns:
        function: the selected function
    """
    switcher = {f.__name__: f for f in [extract_compute, get_thermo, scatter_atoms, command, gather_atoms]}
    return switcher.get(argument)


if __name__ == '__main__':
    while True:
        if MPI.COMM_WORLD.rank == 0:
            input_dict = pickle.load(sys.stdin.buffer)
            # with open('process.txt', 'a') as file:
            #     print('Input:', input_dict, file=file)
        else:
            input_dict = None
        input_dict = MPI.COMM_WORLD.bcast(input_dict, root=0)
        if input_dict['c'] == 'close':
            job.close()
            break
        output = select_cmd(input_dict['c'])(input_dict['d'])
        if MPI.COMM_WORLD.rank == 0 and output is not None:
            # with open('process.txt', 'a') as file:
            #     print('Output:', output, file=file)
            pickle.dump(output, sys.stdout.buffer)
            sys.stdout.flush()

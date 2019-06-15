from ctypes import c_double, c_int
import numpy as np
import pickle
import sys

try:
    from lammps import lammps
except ImportError:
    pass


# Lammps executable
job = lammps(cmdargs=['-screen', 'lmpscreen.log'])


def extract_compute(funct_args):
    return np.array(job.extract_compute(*funct_args))


def get_thermo(funct_args):
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
        input_dict = pickle.load(sys.stdin.buffer)
        if input_dict['c'] == 'close':
            job.close()
            break
        output = select_cmd(input_dict['c'])(input_dict['d'])
        if output is not None:
            pickle.dump(output, sys.stdout.buffer)
            sys.stdout.flush()

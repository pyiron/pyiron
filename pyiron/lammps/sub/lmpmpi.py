from ctypes import c_double, c_int
from mpi4py import MPI
import numpy as np
import pickle
import sys

try:
    from lammps import lammps
except ImportError:
    pass

# Protocol signals
control_data = bytes([1])
control_stop = bytes([0])

# Lammps executable
job = lammps(comm=MPI.COMM_WORLD, cmdargs=["-screen", "none"])


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
    atoms = job.gather_atoms(*funct_args)
    if MPI.COMM_WORLD.rank == 0:
        return np.array(atoms)


def select_cmd(argument):
    """
    Select a lammps command

    Args:
        argument (str): [close, extract_compute, get_thermo, scatter_atoms, command, gather_atoms]

    Returns:
        function: the selected function
    """
    switcher = {
        f.__name__: f
        for f in [extract_compute, get_thermo, scatter_atoms, command, gather_atoms]
    }
    return switcher.get(argument)


def send_data(arr, buff):
    data_str = pickle.dumps(arr)
    dlen = len(data_str).to_bytes(8, byteorder='big')
    buff.write(control_data)
    buff.write(dlen)
    buff.write(data_str)
    buff.flush()


def recv_data(buff):
    data = buff.read(1)
    if data == control_data:
        data = buff.read(8)
        dlen = int.from_bytes(data, byteorder='big')
        data = buff.read(dlen)
        return pickle.loads(data)
    elif data == control_stop:
        return {'c': 'close', 'b': False}
    else:
        raise ValueError('Unexpected Signal!')


if __name__ == "__main__":
    while True:
        if MPI.COMM_WORLD.rank == 0:
            input_dict = recv_data(buff=sys.stdin.buffer)
            if input_dict['b']:
                with open('process.txt', 'a') as file:
                    if input_dict['c'] == 'command':
                        print('Input:', input_dict, file=file)
                    else:
                        print('Input:', input_dict['c'], file=file)
        else:
            input_dict = None
        input_dict = MPI.COMM_WORLD.bcast(input_dict, root=0)
        if input_dict['c'] == 'close':
            if MPI.COMM_WORLD.rank == 0:
                with open('process.txt', 'a') as file:
                    print('END', file=file)
            # job.close()
            MPI.COMM_WORLD.Barrier()
            MPI.Finalize()
            break
        output = select_cmd(input_dict["c"])(input_dict["d"])
        if MPI.COMM_WORLD.rank == 0 and output is not None:
            if input_dict['b']:
                with open('process.txt', 'a') as file:
                    print('Output:', file=file)
            send_data(arr=output, buff=sys.stdout.buffer)

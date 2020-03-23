#! /usr/bin/python

from molmod.units import *
from yaff import *
import h5py, numpy as np

#Setting up system and force field
system = System.from_file('system.chk')
ff = ForceField.generate(system, 'pars.txt', rcut=15.0*angstrom, alpha_scale=3.2, gcut_scale=1.5, smooth_ei=True)

#Setting up output
f = h5py.File('output.h5', mode='w')
hdf5 = HDF5Writer(f, step=1)
r = h5py.File('restart.h5', mode='w')
restart = RestartWriter(r, step=10000)
hooks = [hdf5, restart]

#Setting up simulation
dof = CartesianDOF(ff, gpos_rms=1e-08, dpos_rms=1e-06)
opt = CGOptimizer(dof, hooks=[hdf5])
opt.run(1000)
system.to_file('opt.chk')


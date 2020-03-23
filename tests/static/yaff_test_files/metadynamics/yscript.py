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

plumed = ForcePartPlumed(ff.system, fn='plumed.dat')
ff.add_part(plumed)
hooks.append(plumed)

temp = 300.0*kelvin
thermo = NHCThermostat(temp, timecon=100.00000000000001*femtosecond)
hooks.append(thermo)

hooks.append(VerletScreenLog(step=1000))
md = VerletIntegrator(ff, 0.5*femtosecond, hooks=hooks)
md.run(1000)


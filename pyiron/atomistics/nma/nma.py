# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import numpy as np
import tamkin

from pyiron.atomistics.job.atomistic import Trajectory

from molmod.units import *
from molmod.constants import *
from molmod.periodic import periodic

import matplotlib.pyplot as pt

class NMA(tamkin.NMA):
    """
    This is a generic module to do a Normal Mode Analysis on a job type,
    which calculates the gradient and the hessian using TAMkin based on a job object.
    With the NMA object you can animate a certain mode and plot the IR spectrum.
    """
    def __init__(self,job):
        if not job['output/generic/hessian'] is None and not job['output/generic/forces'] is None:
            mol = tamkin.Molecule(job['output/structure/numbers'],job['output/structure/positions'],job['output/structure/masses'],
                                  job['output/generic/energy_tot'],job['output/generic/forces']*-1 ,job['output/generic/hessian'])
        else:
            raise ValueError('An NMA calculation requires a gradient and hessian.')
        super(NMA, self).__init__(mol)

    def animate_nma_mode(self,index,amplitude=1.0,frames=24,spacefill=False,particle_size=0.5):
        print("This mode corresponds to a frequency of {} 1/cm".format(self.nma.freqs[index]/lightspeed/(1./centimeter)))
        coordinates = self.nma.coordinates
        symbols = [periodic[n].symbol for n in self.nma.numbers]

        mode = self.nma.modes[:,index]
        if self.nma.masses3 is not None:
            mode /= np.sqrt(self.nma.masses3)
        mode /= np.linalg.norm(mode)

        positions = np.zeros((frames,len(symbols),3))

        for frame in range(frames):
            factor = amplitude*np.sin(2*np.pi*float(frame)/frames)
            positions[frame] = (coordinates + factor*mode.reshape((-1,3)))/angstrom

        try:
            import nglview
        except ImportError:
            raise ImportError("The animate_nma_mode() function requires the package nglview to be installed")

        animation = nglview.show_asetraj(Trajectory(positions,self.structure))
        if spacefill:
            animation.add_spacefill(radius_type='vdw', scale=0.5, radius=particle_size)
            animation.remove_ball_and_stick()
        else:
            animation.add_ball_and_stick()
        return animation

    def _Lorentz(x,p,w):
        """
        Lorentzian line shape function, p is position of max, w is FWHM and x is current frequency
        """
        return 1./(1.+((p-x)/(w/2.))**2)

    def plot_IR_spectrum(self,width=10*lightspeed/centimeter,scale=1.0,intensities=None,charges=None):
        """
        Plot IR spectrum based on Lorentzian width, freqs can be scaled through scale
        Intensities can be provided (e.g. from a Gaussian job) or calculated from the charges
        """
        if not intensities is None:
            assert len(intensities) == (len(self.nma.freqs)-self.nma.zeros)
        if intensities is None and charges is None:
            raise ValueError('This function requires the charges or the intensities to calculate the line shape')
        elif not intensities is None and not charges is None:
            raise ValueError('Please only provide either the intensities or the charges')
        else:
            xr = np.arange(0,5001,1)*lightspeed/centimeter
            alphas = np.zeros(len(xr))

            # Calculate intensities
            amps = self.nma.modes
            for n, (wn, amps) in enumerate(zip(freqs[self.nma.zeros:],amps[self.nma.zeros:])):
                if not charges is None:
                    intensity = 0.0
                    for k in range(3):
                        for i, qi in enumerate(charges):
                            I = 3*(i-1)+k
                            intensity += (qi*An[I])**2
                else:
                    intensity = intensities[n]
                alphas += intensity*_Lorentz(xr,wn,width)
                print('Mode %i:    freq = %.3f 1/cm    IR ampl. = %.3e a.u.' %(n, wn/(lightspeed/centimeter), intensity))


            pt.clf()
            pt.plot(xr/(lightspeed/centimeter),alphas)
            pp.xlabel('Frequency [1/cm]')
            pp.ylabel('Absorption [a.u.]')
            pt.show()

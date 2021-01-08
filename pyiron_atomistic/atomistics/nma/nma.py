# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import numpy as np
import tamkin

from pyiron_atomistic.atomistics.job.atomistic import Trajectory

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
    def __init__(self,job, atomic_units=False):
        self.job = job
        if not job['output/generic/hessian'] is None and not job['output/generic/forces'] is None:
            structure = job.get_structure(-1)
            if atomic_units:
                mol = tamkin.Molecule(structure.get_atomic_numbers(),structure.get_positions(),np.array(structure.get_masses())*amu,
                                  job['output/generic/energy_tot'],job['output/generic/forces']*-1 ,job['output/generic/hessian'])
            else:
                mol = tamkin.Molecule(structure.get_atomic_numbers(),structure.get_positions()*angstrom,np.array(structure.get_masses())*amu,
                                  job['output/generic/energy_tot']*electronvolt,job['output/generic/forces']*-1*electronvolt/angstrom ,job['output/generic/hessian']*electronvolt/angstrom**2)
        else:
            raise ValueError('An NMA calculation requires a gradient and hessian.')
        super(NMA, self).__init__(mol)

    def animate_nma_mode(self,index,amplitude=1.0,frames=24,spacefill=False,particle_size=0.5):
        '''
            Visualize the normal mode corresponding to an index

            **Arguments**

            index       index corresponding to a normal mode

            amplitude   size of the deviation of the normal mode

            frames      number of frames that constitute the full mode (lower means faster movement)

            spacefill   remove atom bonds

            particle size
                        size of the atoms in the structure
        '''
        print("This mode corresponds to a frequency of {} 1/cm".format(self.freqs[index]/lightspeed/(1./centimeter)))
        coordinates = self.coordinates
        symbols = [periodic[n].symbol for n in self.numbers]

        mode = self.modes[:,index]
        if self.masses3 is not None:
            mode /= np.sqrt(self.masses3)
        mode /= np.linalg.norm(mode)

        positions = np.zeros((frames,len(symbols),3))

        for frame in range(frames):
            factor = amplitude*np.sin(2*np.pi*float(frame)/frames)
            positions[frame] = (coordinates + factor*mode.reshape((-1,3)))/angstrom

        try:
            import nglview
        except ImportError:
            raise ImportError("The animate_nma_mode() function requires the package nglview to be installed")

        animation = nglview.show_asetraj(Trajectory(positions,self.job.structure))
        if spacefill:
            animation.add_spacefill(radius_type='vdw', scale=0.5, radius=particle_size)
            animation.remove_ball_and_stick()
        else:
            animation.add_ball_and_stick()
        return animation

    def plot_IR_spectrum(self,width=10*lightspeed/centimeter,scale=1.0,intensities=None,charges=None):
        """
            Plot IR spectrum based on Lorentzian width, freqs can be scaled through scale
            Intensities can be provided (e.g. from a Gaussian job) or calculated from the charges

            **Arguments**

            width       width of the Lorentzian function

            scale       scales the frequencies with this factor

            intensities
                        IR intensities for spectrum, can be read from Gaussian job

            charges     charges to calculate IR intensities, from e.g. Yaff simulation

        """
        if not intensities is None:
            assert len(intensities) == (len(self.freqs)-len(self.zeros))
        if intensities is None and charges is None:
            raise ValueError('This function requires the charges or the intensities to calculate the line shape')
        elif not intensities is None and not charges is None:
            raise ValueError('Please only provide either the intensities or the charges')
        else:
            xr = np.arange(0,5001,1)*lightspeed/centimeter
            alphas = np.zeros(len(xr))

            # Calculate intensities
            amps = self.modes
            freqs = self.freqs * scale
            for n, (wn, ampn) in enumerate(zip(np.delete(freqs,self.zeros),np.delete(amps,self.zeros,axis=0))): #self.zeros contain the indices of the zero frequencies
                if not charges is None:
                    intensity = 0.0
                    for k in range(3):
                        for i, qi in enumerate(charges):
                            idx = 3*i+k
                            intensity += (qi*ampn[idx])**2
                else:
                    intensity = intensities[n]
                alphas += intensity*self._lorentz(xr,wn,width)
                print('Mode %i:    freq = %.3f 1/cm    IR ampl. = %.3e a.u.' %(n, wn/(lightspeed/centimeter), intensity))


            pt.clf()
            pt.plot(xr/(lightspeed/centimeter),alphas)
            pt.xlabel('Frequency [1/cm]')
            pt.ylabel('Absorption [a.u.]')
            pt.show()

    @staticmethod
    def _lorentz(x,p,w):
        """
        Lorentzian line shape function, p is position of max, w is FWHM and x is current frequency
        """
        return 1./(1.+((p-x)/(w/2.))**2)

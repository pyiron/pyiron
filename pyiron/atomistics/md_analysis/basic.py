# coding: utf-8
# Copyright (c) Max-Planck-Institut f?r Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.


import numpy as np
from scipy.stats import chi2
import matplotlib.pyplot as pt
from matplotlib.ticker import MaxNLocator

from molmod.units import *
from molmod.constants import *


def plot_temp_dist(job,temp=None, ndof=None):
    """Plots the distribution of the weighted atomic velocities

       **Arguments:**

       job
            A pyiron job object containing MD data

       **Optional arguments:**

       temp
            The (expected) average temperature used to plot the theoretical
            distributions.

       ndof
            The number of degrees of freedom. If not specified, this is chosen
            to be 3*(number of atoms)

       This type of plot is essential for checking the sanity of a simulation.
       The empirical cumulative distribution is plotted and overlayed with the
       analytical cumulative distribution one would expect if the data were
       taken from an NVT ensemble.

       This type of plot reveals issues with parts that are relatively cold or
       warm compared to the total average temperature. This helps to determine
       (the lack of) thermal equilibrium.
    """
    # Make an array with the weights used to compute the temperature
    structure = job.get_structure(-1)
    weights = np.array(np.array(structure.get_masses())*amu)/boltzmann

    # Load the temperatures from the output file
    temps = job['output/generic/temperature']
    if temp is None:
        temp = temps.mean()

    # System definition
    natom = structure.get_atomic_numbers().shape[0]
    if ndof is None:
        ndof = 3*natom
    sigma = temp*np.sqrt(2.0/ndof)
    temp_step = sigma/5

    # Setup the temperature grid and make the histogram
    temp_grid = np.arange(max(0, temp-3*sigma), temp+5*sigma, temp_step)
    counts = np.histogram(temps.ravel(), bins=temp_grid)[0]
    total = float(len(temps))

    # transform into empirical pdf and cdf
    emp_sys_pdf = counts/total
    emp_sys_cdf = counts.cumsum()/total

    # the analytical form
    rv = chi2(ndof, 0, temp/ndof)
    x_sys = temp_grid[:-1]
    ana_sys_pdf = rv.cdf(temp_grid[1:]) - rv.cdf(temp_grid[:-1])
    ana_sys_cdf = rv.cdf(temp_grid[1:])


    pt.subplot(2, 1, 1)
    pt.title('System (ndof=%i)' % ndof)
    scale = 1/emp_sys_pdf.max()
    pt.plot(x_sys, emp_sys_pdf*scale, 'k-', drawstyle='steps-pre', label='Sim (%.0f)' % (temps.mean()))
    pt.plot(x_sys, ana_sys_pdf*scale, 'r-', drawstyle='steps-pre', label='Exact (%.0f)' % temp)
    pt.axvline(temp, color='k', ls='--')
    pt.ylim(ymin=0)
    pt.xlim(x_sys[0], x_sys[-1])
    pt.ylabel('Rescaled PDF')
    pt.legend(loc=0)
    pt.gca().get_xaxis().set_major_locator(MaxNLocator(nbins=5))

    pt.subplot(2, 1, 2)
    pt.plot(x_sys, emp_sys_cdf, 'k-', drawstyle='steps-pre')
    pt.plot(x_sys, ana_sys_cdf, 'r-', drawstyle='steps-pre')
    pt.axvline(temp, color='k', ls='--')
    pt.xlim(x_sys[0], x_sys[-1])
    pt.ylim(0, 1)
    pt.gca().get_xaxis().set_major_locator(MaxNLocator(nbins=5))
    pt.ylabel('CDF')
    pt.xlabel('Temperature [%s]' % 'K')
    pt.show()

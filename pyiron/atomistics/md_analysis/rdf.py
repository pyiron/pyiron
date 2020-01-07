# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.


import numpy as np
import matplotlib.pyplot as pt

from molmod.units import *
from yaff.pes.ext import Cell

class RDF(object):
    """
    This is a generic module to construct a radial distribution function based on a job object.
    With the RDF object you can plot the rdf and the cdf, and calculate the coordination numbers.
    """
    def __init__(self,job,atom_1,atom_2,rcut=20*angstrom,rspacing=0.01*angstrom,nimage=1,start=0,stop=-1,nf=0,save=False,atomic_units=False):
        ''' Computes RDF for two atoms: atom_1 and atom_2.

            **Arguments:**
                job
                    job object that contains an MD trajectory of the structure

                atom_1,atom_2
                    atom numbers between which the RDF is computed

            **Optional arguments:**

                rcut
                    The cutoff for the RDF analysis. This should be lower than the
                    spacing between the primitive cell planes, multiplied by (1+2*nimage).

                rspacing
                    The width of the bins to build up the RDF.

                nimage
                    number of periodic images taken into account

                start,stop
                    First and last index of the slice that is taken into account

                nf
                    The first nf atoms are not taken into account when constructing the RDF
                    (convenient when studying for instance H2O adsorbed in a framework)

                save
                    Save the RDF data on disk

                atomic_units
                    True if pyiron output is in atomic units
        '''
        self.job = job
        structure = self.job.get_structure(-1)
        numbers = structure.get_atomic_numbers()

        select0 = []
        select1 = []
        for i, num in enumerate(numbers):
            if num == atom_1 and i >= nf:
                select0.append(i)
            elif num == atom_2 and i >= nf:
                select1.append(i)
        if atom_1 == atom_2:
            select1 = None

        # Compute RDF
        self._calculate(rspacing, rcut, start, stop, select0, select1, nimage)
        self.coord_number = self._calc_coord_number()

        if save:
            # Write out a data file containing the radial distance (column 1), the value of the RDF g(r) (column 2)  and the value of the rho g(r) (column 3)
            g = open(self.job.working_directory + 'RDF_' + suffix + '.dat', 'w')
            g.write('# Distance \tRDF \tCRDF \n')
            for i in xrange(len(self.d)):
                g.write(str(self.d[i]/angstrom) + "\t" + str(self.rdf[i]) + "\t" + str(self.cdf[i]*angstrom**3) + "\n")
            g.close()


    def _calculate(self, rspacing, rcut, start, stop, select0, select1, nimage):
        RDFC = RDF_calculator(self.job, rspacing, rcut, start, stop, select0, select1, nimage)
        self.bins = RDFC.bins
        self.d = RDFC.d
        self.rdf = RDFC.rdf
        self.cdf = RDFC.cdf


    @staticmethod
    def _trapezoidal(x, y):
        # Trapezoid rule for integration
        h = (x[-1]-x[0])/(len(x)-1)
        return h*(np.sum(y) - 0.5*y[0] - 0.5*y[-1])


    def _calc_coord_number(self):
        coord_number = np.zeros(len(self.d))
        for i in range(1, len(coord_number)):
            coord_number[i] = self._trapezoidal(self.d[0:i+1], 4*np.pi*self.cdf[0:i+1]*self.d[0:i+1]**2)
        return coord_number


    def plot(self):
        """
            Plot the RDF and coordination number (integral of CDF)
        """
        f, ax = pt.subplots(2, 1, sharex=True)
        ax[0].plot(self.d/angstrom, self.rdf, 'k-', drawstyle='steps-mid')
        ax[0].set_ylabel('RDF')
        ax[0].set_ylim([0,5])
        ax[1].plot(self.d/angstrom,self.coord_number)
        ax[1].set_xlabel(u'Distance [Å]')
        ax[1].set_ylabel('CN')
        ax[1].set_ylim([0,10])
        pt.xlim(self.bins[0]/angstrom, self.bins[-1]/angstrom)
        pt.show()



class RDF_calculator(object):
    """
        This class allows for the calculation of the RDF. This is based on the RDF class of Yaff (+ extension Aran Lamaire for CDF and CN)
    """
    def __init__(self, job, rspacing, rcut, start, stop, select0, select1, nimage):
        # Check arguments
        if select0 is not None:
            if len(select0) != len(set(select0)):
                raise ValueError('No duplicates are allowed in select0')
            if len(select0) == 0:
                raise ValueError('select0 can not be an empty list')
        if select1 is not None:
            if len(select1) != len(set(select1)):
                raise ValueError('No duplicates are allowed in select1')
            if len(select1) == 0:
                raise ValueError('select1 can not be an empty list')
        if select0 is not None and select1 is not None and len(select0) + len(select1) != len(set(select0) | set(select1)):
            raise ValueError('No overlap is allowed between select0 and select1. If you want to compute and RDF within a set of atoms, omit the select1 argument.')
        if select0 is None and select1 is not None:
            raise ValueError('select1 can not be given without select0.')

        # If select0 is None, all atoms are considered

        # Init params
        self.job = job
        self.rcut = rcut
        self.rspacing = rspacing
        self.select0 = select0
        self.select1 = select1

        self.nimage = nimage
        self.nbin = int(self.rcut/self.rspacing)
        self.bins = np.arange(self.nbin+1)*self.rspacing
        self.d = self.bins[:-1] + 0.5*self.rspacing
        self.rdf_sum = np.zeros(self.nbin, float)
        self.cdf_sum = np.zeros(self.nbin, float)

        self.nsample = 0
        self.frames = len(self.job['output/generic/energy_tot'])

        self.pos = self.job['output/generic/positions'] * angstrom # internal units are angstrom
        self.cells = self.job['output/generic/cells'] * angstrom # internal units are angstrom

        # Determine the number of atoms
        structure = job.get_structure(-1)
        numbers = structure.get_atomic_numbers()
        self.natom = numbers.shape[0]

        if self.select0 is None:
            self.natom0 = self.natom
        else:
            self.natom0 = len(self.select0)

        self.pos0 = np.zeros((self.natom0, 3), float)
        # the number of pairs
        if self.select1 is None:
            self.npair = (self.natom0*(self.natom0-1))//2
            self.pos1 = None
        else:
            self.natom1 = len(self.select1)
            self.pos1 = np.zeros((self.natom1, 3), float)
            self.npair = self.natom0*self.natom1
        # multiply the number of pairs by all images
        self.npair *= (1 + 2*self.nimage)**3
        # Prepare the output
        self.work = np.zeros(self.npair, float)

        # Loop over the frames and calculate final results
        if stop < 0:
            stop = self.frames - (stop+1)
        for i in range(start, stop):
            self.read_frame(i)
            self.compute_iteration()
        self.compute_derived()


    # Read frame
    def read_frame(self,i):
        self.get_cell(i)
        if self.select0 is None:
            self.pos0 = self.pos[i]
        else:
            self.pos0 = self.pos[i][self.select0]
        if self.select1 is not None:
            self.pos1 = self.pos[i][self.select1]

    def get_cell(self,i):
        self.cell = Cell(self.cells[i])
        if self.cell.nvec != 3:
            raise ValueError('RDF can only be computed for 3D periodic systems.')
        if (2*self.rcut > self.cell.rspacings*(1+2*self.nimage)).any():
            raise ValueError('The 2*rcut argument ({}) should not exceed any of the cell spacings (max:{}).'.format(2*self.rcut, min(self.cell.rspacings*(1+2*self.nimage))))

    # Compute iteration
    def compute_iteration(self):
        self.cell.compute_distances(self.work, self.pos0, self.pos1, nimage=self.nimage)
        counts = np.histogram(self.work, bins=self.bins)[0]
        normalization = self.npair/(self.cell.volume*(1+2*self.nimage)**3)*(4*np.pi*self.rspacing)*self.d**2
        self.rdf_sum += counts/normalization
        self.cdf_sum += counts/(self.natom0*4*np.pi*self.rspacing*self.d**2)
        self.nsample += 1

    # Compute derived
    def compute_derived(self):
        self.rdf = self.rdf_sum/self.nsample
        if self.select1 is None:
            self.cdf_sum *= (1 + 2*self.nimage)**3*self.natom0**2/self.npair
        self.cdf = self.cdf_sum/self.nsample

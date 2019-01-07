# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from __future__ import print_function
import pandas as pd
import shutil
import os
from pyiron.base.settings.generic import Settings
from pyiron.base.generic.parameters import GenericParameters
from pyiron.atomistics.job.potentials import PotentialAbstract

__author__ = "Joerg Neugebauer, Sudarsan Surendralal, Jan Janssen"
__copyright__ = "Copyright 2017, Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Sudarsan Surendralal"
__email__ = "surendralal@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2017"

s = Settings()


class LammpsPotential(GenericParameters):

    """
    This module helps write commands which help in the control of parameters related to the potential used in LAMMPS
    simulations
    """

    def __init__(self, input_file_name=None):
        super(LammpsPotential, self).__init__(input_file_name=input_file_name,
                                              table_name="potential_inp",
                                              comment_char="#")
        self._potential = None
        self._attributes = {}
        self._df = None

    @property
    def df(self):
        return self._df

    @df.setter
    def df(self, new_dataframe):
        self._df = new_dataframe
        # ToDo: In future lammps should also support more than one potential file - that is currently not implemented.
        try:
            self.load_string(''.join(list(new_dataframe['Config'])[0]))
        except IndexError:
            raise ValueError('Potential not found! '
                             'Validate the potential name by self.potential in self.list_potentials().')

    def remove_structure_block(self):
        self.remove_keys(["units"])
        self.remove_keys(["atom_style"])
        self.remove_keys(["dimension"])

    @property
    def files(self):
        if list(self._df['Filename'])[0]:
            absolute_file_paths = [files for files in list(self._df['Filename'])[0] if os.path.isabs(files)]
            relative_file_paths = [files for files in list(self._df['Filename'])[0] if not os.path.isabs(files)]
            for path in relative_file_paths:
                for resource_path in s.resource_paths:
                    if os.path.exists(os.path.join(resource_path, 'lammps', 'potentials')):
                        resource_path = os.path.join(resource_path, 'lammps', 'potentials')
                    if os.path.exists(os.path.join(resource_path, path)):
                        absolute_file_paths.append(os.path.join(resource_path, path))
                        break
            if len(absolute_file_paths) != len(list(self._df['Filename'])[0]):
                raise ValueError('Was not able to locate the potentials.')
            else:
                return absolute_file_paths

    def copy_pot_files(self, working_directory):
        if self.files is not None:
            _ = [shutil.copy(path_pot, working_directory) for path_pot in self.files]

    def get_element_lst(self):
        return list(self._df['Species'])[0]

    def to_hdf(self, hdf, group_name=None):
        if self._df is not None:
            with hdf.open('potential') as hdf_pot:
                hdf_pot['Config'] = self._df['Config'].values[0]
                hdf_pot['Filename'] = self._df['Filename'].values[0]
                hdf_pot['Name'] = self._df['Name'].values[0]
                hdf_pot['Model'] = self._df['Model'].values[0]
                hdf_pot['Species'] = self._df['Species'].values[0]
        super(LammpsPotential, self).to_hdf(hdf, group_name=group_name)

    def from_hdf(self, hdf, group_name=None):
        with hdf.open('potential') as hdf_pot:
            try:
                self._df = pd.DataFrame({'Config': [hdf_pot['Config']],
                                         'Filename': [hdf_pot['Filename']],
                                         'Name': [hdf_pot['Name']],
                                         'Model': [hdf_pot['Model']],
                                         'Species': [hdf_pot['Species']]})
            except ValueError:
                pass
        super(LammpsPotential, self).from_hdf(hdf, group_name=group_name)


class LammpsPotentialFile(PotentialAbstract):
    """
    The Potential class is derived from the PotentialAbstract class, but instead of loading the potentials from a list,
    the potentials are loaded from a file.

    Args:
        potential_df:
        default_df:
        selected_atoms:
    """

    def __init__(self, potential_df=None, default_df=None, selected_atoms=None):
        if potential_df is None:
            potential_df = self._get_potential_df(plugin_name='lammps',
                                                  file_name_lst={'potentials_lammps.csv'},
                                                  backward_compatibility_name='lammpspotentials')
        super(LammpsPotentialFile, self).__init__(potential_df=potential_df, default_df=default_df,
                                                  selected_atoms=selected_atoms)

    def default(self):
        if self._default_df is not None:
            atoms_str = '_'.join(sorted(self._selected_atoms))
            return self._default_df[(self._default_df['Name'] == self._default_df.loc[atoms_str].values[0])]
        return None

    def find_default(self, element):
        """
        Find the potentials

        Args:
            element (set, str): element or set of elements for which you want the possible LAMMPS potentials
            path (bool): choose whether to return the full path to the potential or just the potential name

        Returns:
            list: of possible potentials for the element or the combination of elements

        """
        if isinstance(element, set):
            element = element
        elif isinstance(element, list):
            element = set(element)
        elif isinstance(element, str):
            element = set([element])
        else:
            raise TypeError('Only, str, list and set supported!')
        element_lst = list(element)
        if self._default_df is not None:
            merged_lst = list(set(self._selected_atoms + element_lst))
            atoms_str = '_'.join(sorted(merged_lst))
            return self._default_df[(self._default_df['Name'] == self._default_df.loc[atoms_str].values[0])]
        return None

    def __getitem__(self, item):
        potential_df = self.find(element=item)
        selected_atoms = self._selected_atoms + [item]
        return LammpsPotentialFile(potential_df=potential_df, default_df=self._default_df, selected_atoms=selected_atoms)


class PairPotential(GenericParameters):
    """
        Create potential file for user defined "pair_style" and potential parameters.
        Input: pair_style string, coefficients datastrame
        
        The input table pair coefficients need to be mapped onto the LAMMPS atom type labels.
    """
    def __init__(self, pair_style, cutoff, coefficient_df, structure):
        self._pair_style = pair_style
        self._cutoff = cutoff
        self._coefficient_df = coefficient_df
        self._structure = structure

        self._config_variables = None


    def _config_file(self):
        self._config_variables = []
        style_str = 'pair_style '+self._pair_style+' '+str(self._cutoff)+ ' \n'
        self._config_variables.append(style_str)
        
        for pair,coeff in zip(self._coefficient_df['pairs'], self._coefficient_df['coefficients']):
            coeff_str = 'pair_coeff '
            coeff_str += ''.join(str(n)+' ' for n in pair)
            coeff_str += ''.join(str(c)+' ' for c in coeff)
            coeff_str += ' \n'
            self._config_variables.append(coeff_str)

        return self._config_variables

    @property
    def elements(self):
        return list(self._structure.get_species_symbols())

    @property
    def potential(self):
        return pd.DataFrame({'Config':[self._config_file()],
                                 'Filename':[[]],
                                 'Model':[self._pair_style],
                                 'Name':['user_defined'],
                                 'Species':[self.elements]})
 

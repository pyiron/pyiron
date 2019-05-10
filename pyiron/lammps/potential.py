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
__copyright__ = "Copyright 2019, Max-Planck-Institut für Eisenforschung GmbH - " \
                "Computational Materials Design (CM) Department"
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

    def set_parameter(self, parameter, elements=[], value='not defined'):
        if self.custom_potential is not None:
            self.custom_potential.set_parameter(parameter=parameter, elements=elements, value=value)
            self._df = self.custom_potential.df

    def to_hdf(self, hdf, group_name=None):
#<<<<<<< HEAD
#        with hdf.open('potential') as hdf_pot:
#            if self.custom_potential is not None:
#                for key in ['Config', 'Filename', 'Name', 'Model', 'Species']:
#                    hdf_pot[key] = self.custom_potential.df[key].values[0]
#                if 'Content' in self.custom_potential.df.columns:
#                    hdf_pot['Content'] = self.custom_potential.df['Content'].values[0]
#            elif self._df is not None:
#                for key in ['Config', 'Filename', 'Name', 'Model', 'Species']:
#                    hdf_pot[key] = self._df[key].values[0]
#                hdf_pot['Content'] = self.potential_content
#=======
        if self._df is not None:
            with hdf.open('potential') as hdf_pot:
                hdf_pot['Config'] = self._df['Config'].values[0]
                hdf_pot['Filename'] = self._df['Filename'].values[0]
                hdf_pot['Name'] = self._df['Name'].values[0]
                hdf_pot['Model'] = self._df['Model'].values[0]
                hdf_pot['Species'] = self._df['Species'].values[0]
#>>>>>>> master
        super(LammpsPotential, self).to_hdf(hdf, group_name=group_name)

    def from_hdf(self, hdf, group_name=None):
        with hdf.open('potential') as hdf_pot:
            try:
#<<<<<<< HEAD
#                keys = ['Config', 'Filename', 'Model', 'Species', 'Content']
#                if "Content" in hdf.list_nodes():
#                    keys.append('Content')
#                    self._potential_content = hdf_pot['Content']
#                self._df = pd.DataFrame({k: [hdf_pot[k]] for k in keys})
#                #if "Content" in hdf.list_nodes():
#                #    self._df = pd.DataFrame({'Config': [hdf_pot['Config']],
#                #                             'Filename': [hdf_pot['Filename']],
#                #                             'Name': [hdf_pot['Name']],
#                #                             'Model': [hdf_pot['Model']],
#                #                             'Species': [hdf_pot['Species']],
#                #                             'Content': [hdf_pot['Content']]})
#                #    self._potential_content = self._df['Content']
#                #else:
#                #    self._df = pd.DataFrame({'Config': [hdf_pot['Config']],
#                #                             'Filename': [hdf_pot['Filename']],
#                #                             'Name': [hdf_pot['Name']],
#                #                             'Model': [hdf_pot['Model']],
#                #                             'Species': [hdf_pot['Species']]})
#=======
                self._df = pd.DataFrame({'Config': [hdf_pot['Config']],
                                         'Filename': [hdf_pot['Filename']],
                                         'Name': [hdf_pot['Name']],
                                         'Model': [hdf_pot['Model']],
                                         'Species': [hdf_pot['Species']]})
#>>>>>>> master
            except ValueError:
                pass
        super(LammpsPotential, self).from_hdf(hdf, group_name=group_name)


#<<<<<<< HEAD
#class CustomPotential(GenericParameters):
#    def __init__(self, structure, pot_type=None, pot_sub_style=None, file_name=None, meam_library=None, eam_combinations=False):
#        super(CustomPotential, self).__init__(comment_char="//",
#                    separator_char="=", end_value_char=';')
#        self._initialized = False
#        self._model = pot_type
#        self._file_name = file_name
#        self._structure = structure
#        self._combinations = []
#        self._value_modified = {}
#        self._element_indices = None
#        self._eam_comb = eam_combinations
#        self.overlay = False
#        self['sub_potential'] = pot_sub_style
#        self._file_eam = []
#        self._output_file_eam = []
#        pair_pots=['lj/cut','morse','buck','mie/cut','yukawa','born','born/coul/long','gauss']
#        if file_name is None:
#            self._initialize(self._model)
#            if self._model not in pair_pots:
#                print('ERROR: Choose a potential type from', pair_pots)
#            self._initialized = True
#        else:
#            if isinstance(pot_type, str) and pot_type=='eam':
#                self._output_file_eam = file_name
#                with open(self._output_file_eam, 'r') as file:
#                    self._file_eam = file.readlines()
#                self._initialize_eam()
#            if isinstance(pot_type, list):
#                self._initialize_hybrid()
#
#    @property
#    def elements(self):
#        if self._element_indices is None:
#            self._element_indices=OrderedDict(sorted(zip(self._structure.get_chemical_symbols(),self._structure.get_chemical_indices()), key=lambda x: x[1]))
#        return list(self._element_indices.keys())
#
#
#    @property
#    def combinations(self):
#        if len(self._combinations)==0:
#            self._combinations=np.array([[self.elements[i], self.elements[j]] for i in range(len(self.elements)) for j in range (i+1)])
#        return self._combinations
#
#    def __setitem__(self, key, value):
#        if key not in self._dataset['Parameter'] and self._initialized:
#            if key.split('_')[0] not in self._dataset['Parameter']:
#                pair_coeff_parameters, pair_style_parameters = self.available_keys(self._model)
#                print('ERROR:\nParameter '+key+' is not defined in '+self._model+' potential.\n'+
#                      'Available parameters are '+', '.join([str(k) for k in pair_coeff_parameters.keys()])+
#                      ', which can be defined globally or for each pair of elements.\n'+
#                      'Furthermore, '+', '.join([str(k) for k in pair_style_parameters.keys()])+' can be defined globally.\n'+
#                      'For more information, please take a look at pair_coeff '+self._model+' on the LAMMPS website.')
#            else:
#                print('ERROR:\nIt seems you chose a pair of elements not given in your structure')
#        else:
#            super(CustomPotential, self).__setitem__(key, value)
#
#    def _initialize(self,pot):
#
#        pair_coeff_parameters, pair_style_parameters = self.available_keys(self._model)
#        for k,v in pair_style_parameters.items():
#            self[k] = v
#        for k,v in pair_coeff_parameters.items():
#            self[k] = v
#            for value in self.combinations:
#                self[k+'_'+str(value[0])+'_'+str(value[1])] = None
#
#        self._value_modified={k:False for k in pair_coeff_parameters.keys()}
#        if 'sub_potential' in list(self._value_modified.keys()):
#            del self._value_modified['sub_potential']
#
#
#    def set_parameter(self, parameter, elements=[], value='not defined'):
#        if isinstance(elements, list) and len(elements)==0:
#            key = parameter
#        else:
#            if not isinstance(elements, list) and not isinstance(elements, str) and value=='not defined':
#                value = elements
#                elements = []
#                key = parameter
#            else:
#                key = parameter+'_'+elements[0]+'_'+elements[1]
#                if parameter+'_'+elements[0]+'_'+elements[1] not in list(self.get_pandas()['Parameter']):
#                    elements[0], elements[1] = elements[1], elements[0]
#                key = parameter+'_'+elements[0]+'_'+elements[1]
#        if value=='not defined':
#            raise ValueError('Value not given. Set a number or None, if the parameter should be removed')
#        self[key]=value
#        self._value_modified[key] = True
#        if len(elements) != 0:
#            pair_coeff_parameters, _ = self.available_keys(self._model)
#            for keys in pair_coeff_parameters.keys():
#                if keys==parameter:
#                    continue
#                if keys+'_'+elements[0]+'_'+elements[1] in list(self.get_pandas()['Parameter']):
#                    if self[keys+'_'+elements[0]+'_'+elements[1]] is None:
#                        self[keys+'_'+elements[0]+'_'+elements[1]] = self[keys]
#
#    def get_parameter(self, parameter, element_one=None, element_two=None):
#        if element_one is None and element_two is None:
#            key = parameter
#        else:
#            key = parameter+'_'+element_one+'_'+element_two
#            if key not in list(self.get_pandas()['Parameter']):
#                key = parameter+'_'+element_two+'_'+element_one
#        if self[key] is None:
#            return ''
#        return self[key]
#
#    def _config_file(self, pot):
#        pair_coeff_parameters, pair_style_parameters = self.available_keys(self._model)
#        config_vars = []
#        pair_style=['pair_style '+self._model+' '+' '.join([str(self[k]) for k in list(pair_style_parameters.keys())])+ ' \n']
#        config_vars.append('pair_coeff * * '+' '.join([str(self.get_parameter(k)) for k in list(pair_coeff_parameters.keys())])+'\n')
#        for i in range(len(self.combinations)):
#            if len(str(''.join([str(self.get_parameter(k, self.combinations[i,0], self.combinations[i,1])) for k in list(pair_coeff_parameters.keys())])))==0 or len(self.elements)==1:
#                continue
#            config_vars.append('pair_coeff '+str(self._element_indices[self.combinations[i,1]]+1)+
#                               ' '+str(self._element_indices[self.combinations[i,0]]+1)+' '+
#                               ' '.join([str(self.get_parameter(k, self.combinations[i,0], self.combinations[i,1])) for k in list(pair_coeff_parameters.keys())])+'\n')
#        config_vars=pair_style+config_vars
#        return config_vars
#
#    def _initialize_eam(self):
#        self._keys_eam = [['No_Elements'],
#                          ['Nrho', 'drho', 'Nr', 'dr', 'cutoff_eam'],
#                          ['atomic_number', 'mass', 'lattice_constant', 'lattice_type']]
#        try:
#            for i in range(len(self._file_eam[3].split())-1):
#                self._keys_eam[0].append('Element'+str(i+1))
#            for i in range(3):
#                for k,v in zip(self._keys_eam[i], self._file_eam[i+3].split()):
#                    self[k] = v
#            if len(self._keys_eam[0]) != self['No_Elements']+1:
#                print('WARNING: Number of elements is not consistent in EAM File')
#        except:
#            print('ERROR: Potential file content not valid')
#            raise
#        self._value_modified={k:False for k in self.get_pandas()['Parameter']}
#        self['model']='eam/alloy'
#
#    def _update_eam_file(self):
#        for i, k in enumerate(self._keys_eam):
#            self._file_eam[i+3] = ' '.join([str(self[kk]) for kk in k])
#            self._file_eam[i+3] += '\n'
#
#    def _config_file_eam(self):
#        config_eam = []
#        for key, value in self._value_modified.items():
#            if not value:
#                print('WARNING: '+key+' is not set. Default value: '+str(self[key]))
#        NULL=[]
#        if len(self.elements) <= self['No_Elements']:
#            self['model'] = 'eam/alloy'
#            config_eam=['pair_style '+self['model']+'\n', 'pair_coeff * * '
#                              +self._output_file_eam.split()[-1]+' '+' '
#                              .join([self[k] for i,k in enumerate(self._keys_eam[0]) if i!=0])+'\n']
#        elif len(self.elements) > self['No_Elements']:
#            self['model'] = 'hybrid'
#            for i in range(len(self.elements)-self['No_Elements']):
#                NULL.append('NULL')
#            config_eam=['pair_style '+self['model']+'\n', 'pair_coeff * * '
#                        +self._output_file_eam.split()[-1]+' '+' '
#                        .join([self[k] for i,k in enumerate(self._keys_eam[0]) if i!=0])+' '+' '
#                        .join(NULL)+'\n']
#        return config_eam
#
#    @property
#    def _config(self):
#        if self._file_name is None:
#            return self._config_file(self._model)
#        if self._model == 'eam':
#            return self._config_file_eam()
#        if self._model == 'hybrid':
#            return self._config_file_hybrid()
#        return None
#
#    def _initialize_hybrid(self):
#        for model in self._model:
#            dict_param, _ = self.available_keys(model)
#            if dict_param is not None:
#                continue
#            if '.lmp' in model:
#                potential_filename = potential_filename.split('.lmp')[0]
#            #potential_db = LammpsPotentialFile()
#            #self.input.potential.df = potential_db.find_by_name(potential_filename)
#        if self['sub_potential'] is None:
#            print ('Error: Substyle for hybrid with EAM is not defined')
#        self._output_file_eam = self._file_name
#        with open(self._output_file_eam, 'r') as file:
#            self._file_eam = file.readlines()
#        self._initialize(self['sub_potential'])
#        self._initialize_eam()
#
#    def _lammps_model(self):
#        if isinstance(self._model, list):
#            if self.overlay:
#                return 'hybrid/overlay'
#            else:
#                return 'hybrid'
#        else:
#            return self._model
#
#    def _config_file_hybrid(self):
#        if self['sub_potential'] is None:
#            print ('Error : Sub Pair Style for hybrid is not defined')
#        eam_pot_config=self._config_file_eam()[1:]
#        for i in range(len(eam_pot_config)):
#            eam_pot_config[i]=eam_pot_config[i].split()
#            eam_pot_config[i].insert(3, self['model'])
#            eam_pot_config[i].append('\n')
#            eam_pot_config[i]=' '.join(eam_pot_config[i])
#        if not self._eam_comb:
#            pair_pot_config=self._config_file(self['sub_potential'])[5:]
#        elif self._eam_comb:
#            pair_pot_config=self._config_file(self['sub_potential'])[2:]
#        for i in range(len(pair_pot_config)):
#            pair_pot_config[i]=pair_pot_config[i].split()
#            pair_pot_config[i].insert(3, self._model)
#            pair_pot_config[i].append('\n')
#            pair_pot_config[i]=' '.join(pair_pot_config[i])
#
#        pair_style=['pair_style hybrid '+self['model']+' '+self._model+' '+str(self['cutoff'])+' \n']
#        pair_pot_config=pair_style+eam_pot_config+pair_pot_config
#
#        return pair_pot_config
#
#    @property
#    def df(self):
#        if self._model == 'eam':
#            self._update_eam_file()
#        return pd.DataFrame({'Config':[self._config],
#                             'Filename':[self._output_file_eam],
#                             'Model':[self._model],
#                             'Name':['custom_potential'],
#                             'Species':[self.elements],
#                             'Content':[self._file_eam]})
#
#    def available_keys(self, pot):
#        '''
#            Return: pairwise parameters, global parameters (with their initial values)
#        '''
#        if pot.startswith('lj'):
#            return OrderedDict([('sigma', 1), ('epsilon', 0)]), OrderedDict([('cutoff', 8.0)])
#        elif pot.startswith('morse'):
#            return OrderedDict([('D0', 0), ('alpha', 1), ('r0', 1)]), OrderedDict([('cutoff', 8.0)])
#        elif pot.startswith('buck'):
#            return OrderedDict([('A', 0), ('rho', 1), ('C', 0)]), OrderedDict([('cutoff', 8.0)])
#        elif pot.startswith('mie'):
#            return OrderedDict([('epsilon', 0), ('sigma', 1), ('gammaR', 12), ('gammaA', 6)]), OrderedDict([('cutoff', 8.0)])
#        elif pot.startswith('yukawa'):
#            return OrderedDict([('A', 0), ('cutoff', 3)]), OrderedDict([('kappa', 1.0), ('cutoff', 8.0)])
#        elif pot.startswith('born'):
#            return OrderedDict([('A', 0), ('rho', 1), ('sigma', 1), ('C', 0), ('D', 0), ('cutoff', 4.0)]), OrderedDict([('cutoff', 4.0)])
#        elif pot.startswith('gauss'):
#            return OrderedDict([('A', 0), ('B', 0), ('cutoff', 3.0)]), OrderedDict([('cutoff', 3.0)])
#        else:
#            None, None
#
#
#=======
#>>>>>>> master
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


class PotentialAvailable(object):
    def __init__(self, list_of_potentials):
        self._list_of_potentials = list_of_potentials

    def __getattr__(self, name):
        if name in self._list_of_potentials:
            return name
        else:
            raise AttributeError

    def __dir__(self):
        return self._list_of_potentials

    def __repr__(self):
        return str(dir(self))

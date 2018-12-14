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
from collections import OrderedDict
import numpy as np
import warnings

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

class CustomPotential(GenericParameters):
    def __init__(self, structure, pot_type=None, pot_sub_style=None, file_name=None, meam_library=None, eam_combinations=False):
        self._initialized = False
        super(CustomPotential, self).__init__(comment_char="//",
                    separator_char="=", end_value_char=';')
        self._model = pot_type
        self._file_name = file_name
        self._structure = structure
        self._combinations = []
        self._Config_vars = []
        self._Config_eam = []
        self._value_modified = {}
        self._element_indices = None
        self._eam_comb = eam_combinations
        self._pair_pots=['lj','morse','buckingham','mie','yukawa','born','gauss']
        self['sub_potential'] = pot_sub_style
        if file_name is None:
            self._initialize(self._model)
            if self._model not in self._pair_pots:
                print('ERROR: Choose a potential type from',self._pair_pots) 
            self._initialized = True
        if file_name is not None:
            self._initialize_hybrid()
            #if pot_type=='eam':
            #    self._output_file_eam = file_name
            #    with open(self._output_file_eam, 'r') as file:
            #        self._file_eam = file.readlines()
            #    self._initialize_eam()
            #if pot_type=='hybrid': 
            #    self._initialize_hybrid()
                
                
    @property
    def elements(self):
        if self._element_indices is None:
            self._element_indices=OrderedDict(sorted(zip(self._structure.get_chemical_symbols(),self._structure.get_chemical_indices()), key=lambda x: x[1]))
        return list(self._element_indices.keys())

    
    @property
    def combinations(self):
        if len(self._combinations)==0:
            self._combinations=np.array([[self.elements[i], self.elements[j]] for i in range(len(self.elements)) for j in range (i+1)])
        return self._combinations
        
    def __setitem__(self, key, value):
        if key not in self._dataset['Parameter'] and self._initialized:
            warnings.warn('Parameter ('+key+') not found in '+self._model+' Potential. Only Parameters '
                          +str(list(self._value_modified.keys()))+' can be set.')
        super(CustomPotential, self).__setitem__(key, value) 
        
    def _initialize(self,pot):
        
        if pot == 'lj':
            for k in ['sigma']+['sigma_'+str(value[0])+'_'+str(value[1]) for value in self.combinations]:
                self[k] = 1
            for k in ['epsilon']+['epsilon_'+str(value[0])+'_'+str(value[1]) for value in self.combinations]:
                self[k] = 0
            self['cutoff']=8.0
            self._model_lammps= 'lj/cut'
            
        if pot == 'morse':
            for k in ['D0']+['D0_'+str(value[0])+'_'+str(value[1]) for value in self.combinations]:
                self[k]=0
            for k in ['alpha']+['alpha_'+str(value[0])+'_'+str(value[1]) for value in self.combinations]:
                self[k]=1
            for k in ['r0']+['r0_'+str(value[0])+'_'+str(value[1]) for value in self.combinations]:
                self[k]=1
            self['cutoff']=8.0
            self._model_lammps= 'morse'
            
        if pot == 'buckingham':
            for k in ['A']+['A_'+str(value[0])+'_'+str(value[1]) for value in self.combinations]:
                self[k]=0
            for k in ['rho']+['rho_'+str(value[0])+'_'+str(value[1]) for value in self.combinations]:
                self[k]=1
            for k in ['C']+['C_'+str(value[0])+'_'+str(value[1]) for value in self.combinations]:
                self[k]=0
            self['cutoff']=8
            self._model_lammps= 'buck'
            
        if pot == 'mie':
            for k in ['epsilon']+['epsilon_'+str(value[0])+'_'+str(value[1]) for value in self.combinations]:
                self[k]=0
            for k in ['sigma']+['sigma_'+str(value[0])+'_'+str(value[1]) for value in self.combinations]:
                self[k]=1
            for k in ['gammaR']+['gammaR_'+str(value[0])+'_'+str(value[1]) for value in self.combinations]:
                self[k]=1
            for k in ['gammaA']+['gammaA_'+str(value[0])+'_'+str(value[1]) for value in self.combinations]:
                self[k]=1
            self['cutoff']=8.0
            self._model_lammps= 'mie/cut'
    
        if pot == 'yukawa':
            for k in ['A']+['A_'+str(value[0])+'_'+str(value[1]) for value in self.combinations]:
                self[k]=1
            for k in ['cutoff']+['cutoff_'+str(value[0])+'_'+str(value[1]) for value in self.combinations]:
                self[k]=3
            self['kappa']=0.0
            self._model_lammps= 'yukawa'

        if pot == 'born':
            for k in ['A']+['A_'+str(value[0])+'_'+str(value[1]) for value in self.combinations]:
                self[k]=0
            for k in ['rho']+['rho_'+str(value[0])+'_'+str(value[1]) for value in self.combinations]:
                self[k]=1
            for k in ['sigma']+['sigma_'+str(value[0])+'_'+str(value[1]) for value in self.combinations]:
                self[k]=1
            for k in ['C']+['C_'+str(value[0])+'_'+str(value[1]) for value in self.combinations]:
                self[k]=0
            for k in ['D']+['D_'+str(value[0])+'_'+str(value[1]) for value in self.combinations]:
                self[k]=0
            for k in ['cutoff_'+str(value[0])+'_'+str(value[1]) for value in self.combinations]:
                self[k] = None
            self['cutoff'] = 4
            self._model_lammps= 'born'
            
        if pot == 'gauss':
            for k in ['A']+['A_'+str(value[0])+'_'+str(value[1]) for value in self.combinations]:
                self[k]=0
            for k in ['B']+['B_'+str(value[0])+'_'+str(value[1]) for value in self.combinations]:
                self[k]=0
            for k in ['cutoff_'+str(value[0])+'_'+str(value[1]) for value in self.combinations]:
                self[k] = None
            self['cutoff'] = 3
            self._model_lammps= 'gauss'
            
        self._value_modified={k:False for k in self.get_pandas()['Parameter']}
        if 'sub_potential' in list(self._value_modified.keys()):
            del self._value_modified['sub_potential']

            
    def set_parameter(self, parameter, element_one=None, element_two=None, value=0):
        if element_one is None and element_two is None:
            key = parameter
        else:
            key = parameter+'_'+element_one+'_'+element_two
            if key not in list(self.get_pandas()['Parameter']):
                key = parameter+'_'+element_two+'_'+element_one
        self[key]=value
        self._value_modified[key] = True

    def get_parameter(self, parameter, element_one=None, element_two=None):
        if element_one is None and element_two is None:
            key = parameter
        else:
            key = parameter+'_'+element_one+'_'+element_two
            if key not in list(self.get_pandas()['Parameter']):
                key = parameter+'_'+element_two+'_'+element_one
        if self[key] is None:
            return ''
        return self[key]
    
    def _config_file(self,pot):
        
        if pot == 'lj':
            for key,value in self._value_modified.items():
                if not value:
                    print('WARNING: '+key+' is not set. Default value: '+str(self[key]))
            pair_style=['pair_style '+self._model_lammps+' '+str(self['cutoff'])+ ' \n']
            self._Config_vars.append('pair_coeff * * '+str(self.get_parameter('epsilon'))+
                                             ' '+str(self.get_parameter('sigma'))+'\n')
            for i in range(len(self.combinations)):  
                if len(self.elements)==1:
                    break
                self._Config_vars.append('pair_coeff '+str(self._element_indices[self.combinations[i,1]]+1)+
                                         ' '+str(self._element_indices[self.combinations[i,0]]+1)+
                                         ' '+str(self.get_parameter('epsilon', self.combinations[i,0], self.combinations[i,1]))+
                                         ' '+str(self.get_parameter('sigma', self.combinations[i,0], self.combinations[i,1]))+'\n')
        
        if pot == 'morse':
            for key,value in self._value_modified.items():
                if not value:
                    print('WARNING: '+key+' is not set. Default value: '+str(self[key]))
            pair_style=['pair_style '+self._model_lammps+' '+str(self['cutoff'])+ ' \n']
            self._Config_vars.append('pair_coeff * * '+str(self.get_parameter('D0'))+
                                             ' '+str(self.get_parameter('alpha'))+' '+str(self.get_parameter('r0'))+'\n')
            for i in range(len(self.combinations)):
                if len(self.elements)==1:
                    break
                self._Config_vars.append('pair_coeff '+str(self._element_indices[self.combinations[i,1]]+1)
                                         +' '+str(self._element_indices[self.combinations[i,0]]+1)
                                         +' '+str(self.get_parameter('D0', self.combinations[i,0], self.combinations[i,1]))
                                         +' '+str(self.get_parameter('alpha', self.combinations[i,0], self.combinations[i,1]))
                                         +' '+str(self.get_parameter('r0', self.combinations[i,0], self.combinations[i,1]))+'\n')

        
        if pot == 'buckingham':
            for key,value in self._value_modified.items():
                if not value:
                    print('WARNING: '+key+' is not set. Default value: '+str(self[key]))
            pair_style=['pair_style '+self._model_lammps+' '+str(self['cutoff'])+ ' \n']
            self._Config_vars.append('pair_coeff * * '+str(self.get_parameter('A'))+
                                             ' '+str(self.get_parameter('rho'))+' '+str(self.get_parameter('C'))+'\n')
            for i in range(len(self.combinations)):
                if len(self.elements)==1:
                    break
                self._Config_vars.append('pair_coeff '+str(self._element_indices[self.combinations[i,1]]+1)
                                         +' '+str(self._element_indices[self.combinations[i,0]]+1)
                                         +' '+str(self.get_parameter('A', self.combinations[i,0], self.combinations[i,1]))
                                         +' '+str(self.get_parameter('rho', self.combinations[i,0], self.combinations[i,1]))
                                         +' '+str(self.get_parameter('C', self.combinations[i,0], self.combinations[i,1]))+'\n')
        
        if pot =='mie':
            for key,value in self._value_modified.items():
                if not value:
                    print('WARNING: '+key+' is not set. Default value: '+str(self[key]))
            pair_style=['pair_style '+self._model_lammps+' '+str(self['cutoff'])+ ' \n']
            self._Config_vars.append('pair_coeff * * '+str(self.get_parameter('epsilon'))
                                             +' '+str(self.get_parameter('sigma'))+' '+str(self.get_parameter('gammaR'))
                                             +' '+str(self.get_parameter('gammaA'))+'\n')
            for i in range(len(self.combinations)):
                if len(self.elements)==1:
                    break
                self._Config_vars.append('pair_coeff '+str(self._element_indices[self.combinations[i,1]]+1)
                                         +' '+str(self._element_indices[self.combinations[i,0]]+1)
                                         +' '+str(self.get_parameter('epsilon', self.combinations[i,0], self.combinations[i,1]))
                                         +' '+str(self.get_parameter('sigma', self.combinations[i,0], self.combinations[i,1]))
                                         +' '+str(self.get_parameter('gammaR', self.combinations[i,0], self.combinations[i,1]))
                                         +' '+str(self.get_parameter('gammaA', self.combinations[i,0], self.combinations[i,1]))+'\n')

            
        
        if pot =='yukawa':
            pair_style=['pair_style '+self._model_lammps+' '+str(self['kappa'])+' '+str(self['cutoff'])+'\n']
            for key,value in self._value_modified.items():
                if not value:
                    print('WARNING: '+key+' is not set. Default value: '+str(self[key]))
            self._Config_vars.append('pair_coeff * * '+str(self.get_parameter('A'))+
                                     ' '+str(self.get_parameter('cutoff'))+'\n')
            for i in range(len(self.combinations)):
                if len(self.elements)==1:
                    break
                self._Config_vars.append('pair_coeff '+str(self._element_indices[self.combinations[i,1]]+1)+
                                         ' '+str(self._element_indices[self.combinations[i,0]]+1)+
                                         ' '+str(self.get_parameter('A', self.combinations[i,0], self.combinations[i,1]))+
                                         ' '+str(self.get_parameter('cutoff', self.combinations[i,0], self.combinations[i,1]))+'\n')
       
        
        if pot == 'born':
            pair_style=['pair_style '+self._model_lammps+' '+str(self['cutoff'])+ ' \n']
            for key,value in self._value_modified.items():
                if not value:
                    print('WARNING: '+key+' is not set. Default value: '+str(self[key]))
            self._Config_vars.append('pair_coeff * * '+str(self.get_parameter('A'))+
                                     ' '+str(self.get_parameter('rho'))+' '+str(self.get_parameter('sigma'))+
                                     ' '+str(self.get_parameter('C'))+' '+str(self.get_parameter('D'))+
                                     ' '+str(self.get_parameter('cutoff'))+'\n')
            for i in range(len(self.combinations)):
                if len(self.elements)==1:
                    break
                self._Config_vars.append('pair_coeff '+str(self._element_indices[self.combinations[i,1]]+1)
                                         +' '+str(self._element_indices[self.combinations[i,0]]+1)
                                         +' '+str(self.get_parameter('A', self.combinations[i,0], self.combinations[i,1]))
                                         +' '+str(self.get_parameter('rho', self.combinations[i,0], self.combinations[i,1]))
                                         +' '+str(self.get_parameter('sigma', self.combinations[i,0], self.combinations[i,1]))
                                         +' '+str(self.get_parameter('C', self.combinations[i,0], self.combinations[i,1]))
                                         +' '+str(self.get_parameter('D', self.combinations[i,0], self.combinations[i,1]))
                                         +' '+str(self.get_parameter('cutoff', self.combinations[i,0], self.combinations[i,1]))+'\n')

        if pot == 'gauss':
            for key,value in self._value_modified.items():
                if not value:
                    print('WARNING: '+key+' is not set. Default value: '+str(self[key]))
            pair_style=['pair_style '+self._model_lammps+' '+str(self['cutoff'])+ ' \n']
            self._Config_vars.append('pair_coeff * * '+str(self.get_parameter('A'))+
                                             ' '+str(self.get_parameter('B'))+' '+str(self.get_parameter('cutoff'))+'\n')
            for i in range(len(self.combinations)):
                if len(self.elements)==1:
                    break
                self._Config_vars.append('pair_coeff '+str(self._element_indices[self.combinations[i,1]]+1)+
                                         ' '+str(self._element_indices[self.combinations[i,0]]+1)+
                                         ' '+str(self.get_parameter('A', self.combinations[i,0], self.combinations[i,1]))+
                                         ' '+str(self.get_parameter('B', self.combinations[i,0], self.combinations[i,1]))+
                                         ' '+str(self.get_parameter('cutoff', self.combinations[i,0], self.combinations[i,1]))+'\n')

        self._Config_vars=pair_style+self._Config_vars
        return self._Config_vars

    def _initialize_eam(self):
        self._keys_eam = [['No_Elements'],
                      ['Nrho', 'drho', 'Nr', 'dr', 'cutoff_eam'],
                      ['atomic_number', 'mass', 'lattice_constant', 'lattice_type']]
        try:
            for i in range(len(self._file_eam[3].split())-1):
                self._keys_eam[0].append('Element'+str(i+1))
            for i in range(3):
                for k,v in zip(self._keys_eam[i], self._file_eam[i+3].split()):
                    self[k] = v
            if len(self._keys_eam[0]) != self['No_Elements']+1:
                print('WARNING: Number of elements is not consistent in EAM File')
        except:
            print('ERROR: Potential file content not valid')
            raise
        self._value_modified={k:False for k in self.get_pandas()['Parameter']}
        self['model']='eam/alloy'
        
    def _update_eam_file(self):
        for i, k in enumerate(self._keys_eam):
            self._file_eam[i+3] = ' '.join([str(self[kk]) for kk in k])
            self._file_eam[i+3] += '\n'
    
    def _config_file_eam(self):
        #if len(self._elements) < self['No_Elements']:
         #   print ('Warning : All elements not defines in structure')
        for key, value in self._value_modified.items():
            if not value:
                print('WARNING: '+key+' is not set. Default value: '+str(self[key]))
        NULL=[]
        if len(self.elements) <= self['No_Elements']:
            self['model'] = 'eam/alloy'
            self._Config_eam=['pair_style '+self['model']+'\n', 'pair_coeff * * '
                              +self._output_file_eam.split()[-1]+' '+' '
                              .join([self[k] for i,k in enumerate(self._keys_eam[0]) if i!=0])+'\n']
        elif len(self.elements) > self['No_Elements']:
            self['model'] = 'hybrid'
            for i in range(len(self.elements)-self['No_Elements']):
                NULL.append('NULL')
            self._Config_eam=['pair_style '+self['model']+'\n', 'pair_coeff * * '
                              +self._output_file_eam.split()[-1]+' '+' '
                              .join([self[k] for i,k in enumerate(self._keys_eam[0]) if i!=0])+' '+' '
                              .join(NULL)+'\n']
        return self._Config_eam

    @property
    def _eam_dataframe(self):
        return pd.DataFrame({
                'Config':[self._config_file_eam()],
                'Filename':[[self._output_file_eam]],
                'Model':[self['model']],
                'Name':['my_potential'],
                'Species':[self.elements]
            })

    def _initialize_hybrid(self):
        if self['sub_potential'] is None: 
            print ('Error: Substyle for hybrid with EAM is not defined')
        self._output_file_eam = self._file_name
        with open(self._output_file_eam, 'r') as file:
            self._file_eam = file.readlines()
        self._initialize(self['sub_potential'])
        self._initialize_eam()

    def _config_file_hybrid(self):
        if self['sub_potential'] is None: 
            print ('Error : Sub Pair Style for hybrid is not defined')
        eam_pot_config=self._config_file_eam()[1:]
        for i in range(len(eam_pot_config)):
            eam_pot_config[i]=eam_pot_config[i].split()
            eam_pot_config[i].insert(3,self['model'])
            eam_pot_config[i].append('\n')
            eam_pot_config[i]=' '.join(eam_pot_config[i])
        if not self._eam_comb:
            pair_pot_config=self._config_file(self['sub_potential'])[5:]
        elif self._eam_comb:
            pair_pot_config=self._config_file(self['sub_potential'])[2:]
        for i in range(len(pair_pot_config)):
            pair_pot_config[i]=pair_pot_config[i].split()
            pair_pot_config[i].insert(3,self._model_lammps)
            pair_pot_config[i].append('\n')
            pair_pot_config[i]=' '.join(pair_pot_config[i])
    
        pair_style=['pair_style hybrid '+self['model']+' '+self._model_lammps+' '+str(self['cutoff'])+' \n']
        pair_pot_config=pair_style+eam_pot_config+pair_pot_config
        
        return pair_pot_config
    
    @property
    def _hybrid_dataframe(self):
        return pd.DataFrame({'Config':[self._config_file_hybrid()],
                                 'Filename':[[]],
                                 'Model':[self._model],
                                 'Name':['my_potential'],
                                 'Species':[self.elements]})

    @property
    def potential(self):
        if self._file_name is None:
            return pd.DataFrame({'Config':[self._config_file(self._model)],
                                 'Filename':[[]],
                                 'Model':[self._model],
                                 'Name':['my_potential'],
                                 'Species':[self.elements]})
        elif self._model == 'eam':
            self._update_eam_file()
            return self._eam_dataframe
        elif self._model == 'hybrid':
            return self._hybrid_dataframe
        else: 
            raise ValueError('Potential not found')
            

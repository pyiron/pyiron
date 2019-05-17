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

class CustomPotential(GenericParameters):
    def __init__(self, structure, pot_type):
        """
            Args:
                structure: pyiron structure. This argument will be erased when this class is integrated in pyiron
                pot_type (list or str): list of potentials to be used.
                file_name (list or str): list of eam potential files. Currently only one file is permitted. Absolute path
                                         or the name of the potential in database
        """
        self._initialized = False # required for setitem
        warnings.simplefilter(action='ignore', category=FutureWarning)
        super(Potential, self).__init__(comment_char="//",
                    separator_char="=", end_value_char=';')
        self._structure = structure
        self._value_modified = {} # checks whether parameters have been modified (rather for debugging)
        if not isinstance(pot_type, list):
            pot_type = [pot_type]
        self._pot_type = pot_type
        assert self._file_name is None or isinstance(self._file_name, str)
        assert structure is not None, 'structure not defined'
        assert isinstance(pot_type[0], str), 'pot_type must be a string or a list of strings'
        assert len(pot_type)<=2, 'Not possible to superpose more than 2 potentials'
        if self._file_name is not None:
            if not os.path.isfile(self._file_name):
                self._file_name = self.find_absolute_path(self._file_name)
            with open(self._file_name, 'r') as ffile:
                self._file_eam = ffile.readlines()
        if self.style=='hybrid':
            self._initialize_hybrid()
        elif self.style=='pairwise':
            self._initialize(self._model_pair_pot[0])
        else:
            try:
                self._initialize_eam()
            except:
                warnings.warn('EAM parsing was not successful; original EAM file is used')
                self._eam_parsing_successful = False
        self._initialized = True

    def __setitem__(self, key, value):
        if hasattr(value, '__len__'):
            value = np.array(value).tolist()
        param=key.split('/')
        if len(param)==1:
            super(Potential, self).__setitem__(key, value)
        else:
            elems = sorted(param[1].split('-'))
            if len(elems)==1:
                key=str(param[0])+'/'+str(elems[0])+'-'+str(elems[0])
            else:
                key = str(param[0])+'/'+str(elems[0])+'-'+str(elems[1])
            super(Potential, self).__setitem__(key, value)
        for kk in self.pair_coeff_parameters:
            if kk == param[0] or kk.endswith('_eam'):
                continue
            if self[key.replace(param[0], kk)] is None and self[key] is not None:
                self[key.replace(param[0], kk)] = self._pair_pot_keys(self._model_pair_pot[0])[0][kk]
            elif self[key.replace(param[0], kk)] is not None and self[key] is None:
                self[key.replace(param[0], kk)] = None

    def __getitem__(self, key):
        param=key.split('/')
        if len(param) != 1:
            elems = sorted(param[1].split('-'))
            if len(elems)==1:
                key = str(param[0])+'/'+'-'+str(elems[0])+str(elems[0])
            else:
                key = str(param[0])+'/'+str(elems[0])+'-'+str(elems[1])
        value = super(Potential, self).__getitem__(key)
        if isinstance(value, list):
            return np.array(value, dtype=float)
        else:
            return value


    @staticmethod
    def find_absolute_path(file_name):
        assert isinstance(file_name, str)
        from pyiron.lammps.potential import LammpsPotential, LammpsPotentialFile
        potential_df = LammpsPotentialFile().find_by_name(file_name)
        assert len(potential_df)>0
        potential = LammpsPotential()
        potential.df = potential_df
        return potential.files[0]

    @property
    def _file_name(self):
        if self._eam_file_name is None:
            pot_diff = list(set(self._pot_type)-set(self._pair_pot_keys().keys()))
            if len(pot_diff)>0:
                self._file_name = pot_diff[0]
            else:
                self._file_name = None
        return self._eam_file_name

    @_file_name.setter
    def _file_name(self, file_name):
        if file_name is not None and not isinstance(file_name, str):
            raise ValueError('file name has to be a string')
        self._eam_file_name = file_name

    @property
    def style(self):
        if self._file_name is not None:
            if len(self._model_pair_pot)>0:
                assert len(self._pot_type)>1
                return 'hybrid'
            else:
                return 'many_body'
        assert len(self._pot_type)==1, 'currently only combinations of eam and pairwise potential are implemented'
        return 'pairwise'

    @property
    def _model_pair_pot(self):
        return list(set(self._pot_type).intersection(self._pair_pot_keys().keys()))

    @property
    def _model_mb(self): # mb stands for 'many-body'
        if self.style=='pairwise':
            return []
        intersect = list(set(self._pot_type).intersection(self._mb_pot_keys().keys()))
        if len(intersect)!=0:
            return intersect
        return ['eam/alloy']

    @property
    def model(self):
        if len(self._pot_type)>1:
            if self.overlay:
                return "hybrid/overlay"
            return "hybrid"
        if len(self._model_mb)!=0:
            return self._model_mb[0]
        return self._model_pair_pot[0]

    @property
    def element_indices(self):
        try:
            return self._element_indices
        except AttributeError:
            self._element_indices=OrderedDict(sorted(zip(self._structure.get_chemical_symbols(), self._structure.get_chemical_indices()), key=lambda x: x[1]))
            return self._element_indices

    @property
    def elements(self):
        return list(self.element_indices.keys())

    @property
    def combinations(self):
        try:
            return self._combinations
        except AttributeError:
            self._combinations=np.array([[elem, self.elements[j]] for i, elem in enumerate(self.elements) for j in range (i+1)])
            return self._combinations

    @property
    def overlay(self):
        """
            whether to overlay two potentials or not. When overlay=True, the two potentials are superposed (default: True)
            For more info: https://lammps.sandia.gov/doc/pair_hybrid.html
        """
        if self.style != 'hybrid':
            return True
        else:
            try:
                return self._overlay
            except AttributeError:
                return True

    @overlay.setter
    def overlay(self, value):
        assert isinstance(value, bool)
        self._overlay = value



    def _initialize(self, pot):
        for k,v in self.pair_style_parameters.items():
            self[k] = v
        for k,v in self.pair_coeff_parameters.items():
            self[k] = v
            for value in self.combinations:
                self[k+'/'+str(value[0])+'-'+str(value[1])] = None
        self._value_modified={k:False for k in self.pair_coeff_parameters.keys()}

    def _initialize_hybrid(self):
        self._initialize_eam()
        self._initialize(self._model_pair_pot[0])


    @property
    def _pair_style_str(self):
        if len(self._pot_type)>1:
            return 'pair_style '+self.model+' '+self._model_mb[0]+' '+self._model_pair_pot[0]+' '+str(self['cutoff'])+'\n'
        else:
            return 'pair_style '+self.model+' '+' '.join([str(self[k]) for k in list(self.pair_style_parameters.keys())])+'\n'

    def _config_file(self, pot):
        config_vars = []
        if self.overlay:
            config_vars = ['pair_coeff * * '+' '.join([str(self[str(k)]) for k in list(self.pair_coeff_parameters.keys())])+'\n']
        for comb in self.combinations:
            if len(str(''.join([str(self[str(k)+'/'+str(comb[0])+'-'+str(comb[1])]) for k in list(self.pair_coeff_parameters.keys())])))==0 or len(self.elements)==1:
                continue
            kk = list(self.pair_coeff_parameters.keys())[0]
            if self.style=='hybrid':
                if self[str(kk)+'/'+str(comb[0])+'-'+str(comb[1])] is None:
                    if len(set(comb)) > len(set(comb).intersection(self._elements_in_eam)):
                        self[str(kk)+'/'+str(comb[0])+'-'+str(comb[1])]= self._pair_pot_keys(self._model_pair_pot[0])[0][kk]
            if self[str(kk)+'/'+str(comb[0])+'-'+str(comb[1])] is None:
                continue
            config_vars.append('pair_coeff '+str(self.element_indices[comb[1]]+1)+
                               ' '+str(self.element_indices[comb[0]]+1)+' '+
                               ' '.join([str(self[str(k)+'/'+str(comb[0])+'-'+str(comb[1])]) for k in list(self.pair_coeff_parameters.keys())])+'\n')
        config_vars = [self._pair_style_str]+config_vars
        return config_vars

    def _config_file_eam(self):
        NULL=[]
        if len(self.elements) > len(set(self.elements).intersection(self._elements_in_eam)):
            assert self.style == 'hybrid', 'Some elements are not defined in the potential file. Choose a pairwise potential to extend it to hybrid'
            for i in range(len(self.elements)-len(set(self.elements).intersection(self._elements_in_eam))):
                NULL.append('NULL')
        config_eam=['pair_coeff * * '+os.getcwd()+'/potential.dat '+
                    ' '.join([self[k] for i,k in enumerate(self._keys_eam[0]) if i!=0])+' '+' '
                    .join(NULL)+'\n']
        return [self._pair_style_str]+config_eam

    def _config_file_hybrid(self):
        eam_pot_config = self._config_file_eam()[1:]
        eam_pot_config = [' '.join(ll.split()[:3]+[self._model_mb[0]]+ll.split()[3:])+'\n' for ll in eam_pot_config]
        pair_pot_config=self._config_file(self._model_pair_pot[0])[1:]
        pair_pot_config = [' '.join(ll.split()[:3]+[self._model_pair_pot[0]]+ll.split()[3:])+'\n' for ll in pair_pot_config]

        pair_pot_config = [self._pair_style_str]+eam_pot_config+pair_pot_config

        return pair_pot_config

    @property
    def _df(self):
        """
            This function returns a dataframe with parameters, Config:contains lammps compatible pair style and pair coeff,
            File_name, Potential Mode.
        """
        if self.style=='hybrid':
            self._update_eam_file()
            config_file = self._config_file_hybrid()
            species = self._elements_in_eam+list(set(self.elements)-set(self._elements_in_eam))
        elif self.style=='pairwise':
            config_file = self._config_file(self._model_pair_pot[0])
            species = self.elements
        else:
            self._update_eam_file()
            config_file = self._config_file_eam()
            species = self._elements_in_eam
        return pandas.DataFrame({'Config':[config_file],
                                 'Filename':[[]], ### Do it for current directory
                                 'Model':[self.model],
                                 'Name':['my_potential'],
                                 'Species':[species]})

    @property
    def available_potentials(self):
        return list(self._pair_pot_keys().keys())+list(self._mb_pot_keys().keys())

    @property
    def pair_coeff_parameters(self):
        if len(self._model_pair_pot)>0:
            return self._pair_pot_keys(self._model_pair_pot[0])[0]
        else:
            return OrderedDict([])

    @property
    def pair_style_parameters(self):
        if len(self._model_pair_pot)>0:
            return self._pair_pot_keys(self._model_pair_pot[0])[1]
        else:
            return OrderedDict([])

    def _pair_pot_keys(self, pot=None):
        """
            This function returns a list of two sets of parameters. The first item contains the values
            for pair_coeff; the second contains the values for pair_style.

            args:
                pot (str): name of the potential.
                           If the name of the potential is not specified,
                           the entire list of parameters is returned
        """
        if pot is not None and not isinstance(pot, str):
            raise ValueError('Potential name must be a string')
        param = {'lj/cut': [OrderedDict([('epsilon', 0), ('sigma', 1)]),
                            OrderedDict([('cutoff', 8.0)])],
                 'morse': [OrderedDict([('D0', 0), ('alpha', 1), ('r0', 1)]),
                           OrderedDict([('cutoff', 8.0)])],
                 'buck': [OrderedDict([('A', 0), ('rho', 1), ('C', 0)]),
                          OrderedDict([('cutoff', 8.0)])],
                 'mie': [OrderedDict([('epsilon', 0), ('sigma', 1), ('gammaR', 1), ('gammaA', 1)]),
                         OrderedDict([('cutoff', 8.0)])],
                 'yukawa': [OrderedDict([('A', 0), ('cutoff', 3)]),
                            OrderedDict([('kappa', 1.0), ('cutoff', 8.0)])],
                 'born': [OrderedDict([('A', 0), ('rho', 1), ('sigma', 1), ('C', 0), ('D', 0), ('cutoff', None)]),
                          OrderedDict([('cutoff', 4.0)])],
                 'gauss': [OrderedDict([('A', 0), ('B', 0), ('cutoff', None)]),
                           OrderedDict([('cutoff', 3.0)])]}
        if pot is None:
            return param
        elif pot in param.keys():
            return param[pot]
        for k,v in param.items():
            if pot.startswith(k):
                return v
        raise NotImplementedError(pot+' is not implemented')


    def analytical_pot(self, pot, points):
        analytic_pot=[]
        if pot.startswith('lj/cut'):
            for i in np.arange(0.01, self['cutoff'], self['cutoff']/points):
                analytic_pot.append(4*self['epsilon']*((self['sigma']/i)**12-(self['sigma']/i)**6))
        return analytic_pot
        # Implement other analytical Potentials

    @property
    def _elements_in_eam(self):
        if self._keys_eam is None or len(self._keys_eam[0])<=1:
            raise ValueError("EAM parsing not performed")
        return [self[k] for k in self._keys_eam[0][1:]]

    def _mb_pot_keys(self, pot=None):
        return {'eam/alloy': None, 'meam': None, 'tersoff': None}


    @property
    def eam_element_comb(self):
        try:
            return self._eam_comb
        except AttributeError:
            eam_list_comb=[]
            self._eam_comb=np.array(self.combinations).tolist()
            for i,k in enumerate(self.combinations):
                if list(set(self.elements)-set(self._elements_in_eam)) in k:
                    warnings.simplefilter(action='ignore', category=FutureWarning)
                    eam_list_comb.append(i)
            for i in sorted(eam_list_comb, reverse=True):
                del self._eam_comb[i]
            return self._eam_comb

    @property
    def _eam_ending(self):
        if self.style == 'hybrid':
            return '_eam'
        else:
            return ''

    @staticmethod
    def _has_non_number(list_of_something):
        if not hasattr(list_of_something, '__len__'):
            list_of_something = [list_of_something]
        for elem in list_of_something:
            try:
                float(elem)
            except ValueError:
                return True
        return False

    def _initialize_eam(self):
        self._keys_eam = [['No_Elements'+self._eam_ending],
                          ['Nrho'+self._eam_ending, 'drho'+self._eam_ending, 'Nr'+self._eam_ending, 'dr'+self._eam_ending, 'cutoff'+self._eam_ending]]
        self._keys_eam[0].extend(['Element_eam_'+str(i+1) for i in range(len(self._file_eam[3].split())-1)])
        for i in range(2):
            for k,v in zip(self._keys_eam[i], self._file_eam[i+3].split()):
                self[k] = v
        for k in ['Nrho'+self._eam_ending, 'Nr'+self._eam_ending]:
            self[k] = int(self[k])
        for k in ['drho'+self._eam_ending, 'dr'+self._eam_ending, 'cutoff'+self._eam_ending]:
            self[k] = float(self[k])
        # self._elem_props=['atomic_number'+self._eam_ending, 'mass'+self._eam_ending, 'lattice_constant'+self._eam_ending, 'lattice_type'+self._eam_ending]
        self._elem_props=['atomic_number'+self._eam_ending, 'mass'+self._eam_ending]
        self['meaningless_item'] = None

        property_lines = list(filter(lambda a: self._has_non_number(a.split()), self._file_eam[5:]))
        if len(self._elements_in_eam) != len(property_lines):
            warnings.warn('EAM parsing might have failed; number of elements defined ('+str(len(self._elements_in_eam))+') != number of element property lines ('+str(len(property_lines))+')')
        for i_elem, elem in enumerate(self._elements_in_eam):
            for i_prop, prop in enumerate(self._elem_props):
                self[str(prop)+'/'+str(elem)+'-'+(elem)] = property_lines[i_elem].split()[i_prop]
        tab_lines = list(filter(lambda a: self._has_non_number(a.split())==False, self._file_eam[5:]))
        tab_lines = ''.join(tab_lines).replace('\n',' ').split()
        start_line = 0
        for elem in self._elements_in_eam:
            self['F'+self._eam_ending+'/'+elem+'-'+elem] = [float(value) for value in tab_lines[start_line:start_line+self['Nrho'+self._eam_ending]]]
            start_line += self['Nrho'+self._eam_ending]
            self['rho'+self._eam_ending+'/'+elem+'-'+elem] = [float(value) for value in tab_lines[start_line:start_line+self['Nr'+self._eam_ending]]]
            start_line += self['Nr'+self._eam_ending]
        for i in range(self['No_Elements'+self._eam_ending]):
            for j in range(i+1):
                self['phi'+self._eam_ending+'/'+self._elements_in_eam[j]+'-'+self._elements_in_eam[i]] = [float(value) for value in tab_lines[start_line:start_line+self['Nr'+self._eam_ending]]]
                start_line += self['Nr'+self._eam_ending]
        if len(self._keys_eam[0]) != self['No_Elements'+self._eam_ending]+1:
            print('WARNING: Number of elements is not consistent in EAM File')
        self._value_modified = {k:False for k in self.get_pandas()['Parameter']}
        self._eam_parsing_successful = True

    @property
    def comment_str(self):
        return ''.join([line for line in self._file_eam[0:3]])

    @property
    def eam_info_str(self):
        first_line = ' '.join([str(self[k]) for k in self._keys_eam[0]])+'\n'
        second_line = ' '.join([str(self[k]) for k in self._keys_eam[1]])+'\n'
        return first_line+second_line


    def element_prop(self, element_1, element_2):
        elem_prop = ([self[str(prop)+'/'+str(element_1)+'-'+str(element_2)] for prop in self._elem_props])
        prop = ' '.join([str(prop) for prop in elem_prop])+'\n'
        return prop


    def F_rho_str(self, element_1, element_2):
        F = np.array(self['F'+self._eam_ending+'/'+str(element_1)+'-'+str(element_2)]).tolist()
        rho = np.array(self['rho'+self._eam_ending+'/'+str(element_1)+'-'+str(element_2)]).tolist()
        F_rho = '\n'.join([str(value) for value in F+rho])+'\n'
        return F_rho


    def phi_str(self, element_1, element_2):
        phi = np.array(self['phi'+self._eam_ending+'/'+str(element_1)+'-'+str(element_2)]).tolist()
        return '\n'.join([str(pp) for pp in phi])+'\n'


    def _update_eam_file(self):
        return_file = self.comment_str+self.eam_info_str
        for elem in self._elements_in_eam:
            return_file += self.element_prop(elem, elem)+self.F_rho_str(elem, elem)
        for elem in self.eam_element_comb:
            return_file += self.phi_str(elem[0], elem[1])
        with open('potential.dat', 'w') as ff:
            for line in return_file:
                ff.write(line)


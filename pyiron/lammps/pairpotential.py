import pandas as pd
from pyiron.base.generic.parameters import GenericParameters


class PairPotential(GenericParameters):
    """
        Create potential file for user defined "pair_style" and potential parameters.
        Input: pair_style string, coefficients dataframe
    """

    def __init__(self, pair_style, cutoff, coefficient_df, structure):
        self._pair_style = pair_style
        self._cutoff = cutoff
        self._coefficient_df = coefficient_df
        self._structure = structure

        self._config_variables = None

    def _config_file(self):
        self._config_variables = []
        style_str = 'pair_style ' + self._pair_style + ' ' + str(self._cutoff) + ' \n'
        self._config_variables.append(style_str)

        for pair, coeff in zip(self._coefficient_df['pairs'], self._coefficient_df['coefficients']):
            coeff_str = 'pair_coeff '
            coeff_str += ''.join(str(n) + ' ' for n in pair)
            coeff_str += ''.join(str(c) + ' ' for c in coeff)
            coeff_str += ' \n'
            self._config_variables.append(coeff_str)

        return self._config_variables

    @property
    def elements(self):
        return list(self._structure.get_species_symbols())

    @property
    def potential(self):
        return pd.DataFrame({'Config': [self._config_file()],
                             'Filename': [[]],
                             'Model': [self._pair_style],
                             'Name': ['user_defined'],
                             'Species': [self.elements]})

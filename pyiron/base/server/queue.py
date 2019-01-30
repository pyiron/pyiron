# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import warnings

"""
Class to define a job queue 
"""

__author__ = "Jan Janssen"
__copyright__ = "Copyright 2019, Max-Planck-Institut für Eisenforschung GmbH - " \
                "Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2017"


class Queue(object):
    """
    Configuration of an individual que

    Args:
        name (str): name of the que
        mini_cores (int): minimum number of cores for an individual job
        maxi_cores (int): maximum number of cores for an individual job
        divisor_list (int, list): single int or list of int of cores the number of cores selected for an individual
                                  job should be divided by. With this we can force the user to always use full nodes
                                  by setting the divisor_list to 20 cores.
        run_time_limit (int): the time limit in seconds for an individual job
    """
    def __init__(self, name, mini_cores, maxi_cores, divisor_list, run_time_limit):
        self.__name__ = name
        self._minimum_number_of_cores = None
        self._maximum_number_of_cores = None
        self._divisors_for_que = None
        self._run_time_limit = None
        self._cores = None
        self._run_time = None
        self.minimum_number_of_cores = mini_cores
        self.maximum_number_of_cores = maxi_cores
        self.divisors_for_que = divisor_list
        self.run_time_limit = run_time_limit

    @property
    def minimum_number_of_cores(self):
        """
        Get minimum number of cores for an individual job

        Returns:
            int: minimum number of cores
        """
        return self._minimum_number_of_cores

    @minimum_number_of_cores.setter
    def minimum_number_of_cores(self, cores):
        """
        Set minimum number of cores for an individual job

        Args:
            cores (int): minimum number of cores
        """
        if self._check_int_cores(cores):
            self._minimum_number_of_cores = cores

    @property
    def maximum_number_of_cores(self):
        """
        Get maximum number of cores for an individual job

        Returns:
            int: minimum number of cores
        """
        return self._maximum_number_of_cores

    @maximum_number_of_cores.setter
    def maximum_number_of_cores(self, cores):
        """
        Set maximum number of cores for an individual job

        Args:
            cores (int): minimum number of cores
        """
        if self._check_int_cores(cores):
            self._maximum_number_of_cores = cores

    @property
    def divisors_for_que(self):
        """
        Force the user to always use full nodes by setting the divisors to the number of cores per node, so the user can
        run a job either on one node or two nodes but not on one node and a single settings from a second node, which is
        insufficient in most cases.

        Returns:
            list: list of divisors
        """
        return self._divisors_for_que

    @divisors_for_que.setter
    def divisors_for_que(self, divisors):
        """
        Force the user to always use full nodes by setting the divisors to the number of cores per node, so the user can
        run a job either on one node or two nodes but not on one node and a single settings from a second node, which is
        insufficient in most cases.

        Args:
            divisors (int, list): single divisor or a list of divisors
        """
        divisor_lst = []
        if isinstance(divisors, list):
            for divisor in divisors:
                if isinstance(divisor, int):
                    divisor_lst.append(divisor)
                else:
                    raise TypeError('The divisors should be either a list of int, or an single integer.')
            self._divisors_for_que = divisor_lst
        elif isinstance(divisors, int):
            divisor_lst.append(divisors)
            self._divisors_for_que = divisor_lst
        else:
            raise TypeError('The divisors should be either a list of int, or an single integer.')

    @property
    def cores(self):
        """
        The number of cores used for the current simulation, this validated that the number of cores selected agrees
        with the requirements of the individual que.

        Returns:
            int: number of cores
        """
        return self._cores

    @cores.setter
    def cores(self, new_cores):
        """
        Set the number of cores used for the current simulation, this validated that the number of cores selected agrees
        with the requirements of the individual que.

        Args:
            new_cores (int): number of cores
        """
        if self._validate_cores(new_cores):
            self._cores = new_cores

    @property
    def run_time_limit(self):
        """
        The run time limit defined by the que.

        Returns:
            int: run time limit
        """
        return self._run_time_limit

    @run_time_limit.setter
    def run_time_limit(self, new_run_time):
        """
        Set the run time limit for a specific que.

        Args:
            new_run_time (int): run time limit
        """
        if self._check_int_runtime(new_run_time):
            self._run_time_limit = new_run_time

    @property
    def run_time(self):
        """
        The run time used for the current simulation, this validated that the run time selected agrees with the
        requirements of the individual que.

        Returns:
            int: run time
        """
        return self._run_time

    @run_time.setter
    def run_time(self, new_run_time):
        """
        Set the run time used for the current simulation, this validated that the run time selected agrees with the
        requirements of the individual que.

        Args:
            new_run_time (int): run time
        """
        if self._check_int_runtime(new_run_time):
            if new_run_time <= self.run_time_limit:
                self._run_time = new_run_time

    @staticmethod
    def _check_int_cores(number):
        """
        internal function to verify that the number of cores always is an int.

        Args:
            number (int): number of cores

        Returns:
            (bool): True - if false a TypeError is raised.
        """
        if not isinstance(number, int):
            raise TypeError('The number of cores should be int.')
        return True

    @staticmethod
    def _check_int_runtime(number):
        """
        internal function to verify that the run time always is an int.

        Args:
            number (int): run time in seconds

        Returns:
            (bool): True - if false a TypeError is raised.
        """
        if not isinstance(number, int):
            raise TypeError('The runtime should be int in Sec.')
        return True

    def _validate_cores(self, cores):
        """
        internal function to verify that the number of cores selected for the current simulation match the requirements
        of the que, in terms of being above the minimum number of cores and below the maximum number while maintaining
        the required steps of increase defined by the divisor list.

        Args:
            cores (int): number of cores

        Returns:
            bool: [True / False]
        """
        if self._check_int_cores(cores):
            if cores > self.maximum_number_of_cores:
                warnings.warn('The maximum number of cores for the Que: ' + self.__name__ + ' is set to ' +
                              str(self.maximum_number_of_cores) + ' so you can not submit a job with ' +
                              str(cores) + ' Cores.', RuntimeWarning)
            elif cores < self.minimum_number_of_cores:
                warnings.warn('The minimum number of cores for the Que: ' + self.__name__ + ' is set to ' +
                              str(self.minimum_number_of_cores) + ' so you can not submit a job with ' +
                              str(cores) + ' Cores.', RuntimeWarning)
            divsor_included = False
            for divsor in self.divisors_for_que:
                if cores % divsor == 0:
                    divsor_included = True
            if not divsor_included:
                warnings.warn('The number of ' + str(cores) + ' Cores, can not be divided by any devisor ' +
                              str(self.divisors_for_que) + ' set for the Que: ' + self.__name__, RuntimeWarning)

            return True

    def __repr__(self):
        """
        Human readable representaion of the que settings.

        Returns:
            str: que settings
        """
        return str({'Name': self.__name__,
                    'Minimum number of cores': self._minimum_number_of_cores,
                    'Maximum number of cores': self._maximum_number_of_cores,
                    'Runtime limit': self._run_time_limit})

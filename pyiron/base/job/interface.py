class JobInterface(object):
    @staticmethod
    def write_input(job):
        """
        Write the input files for the external executable. This method has to be implemented in the individual
        hamiltonians.
        """
        raise NotImplementedError

    @staticmethod
    def collect_output(job):
        """
        Collect the output files of the external executable and store the information in the HDF5 file. This method has
        to be implemented in the individual hamiltonians.
        """
        raise NotImplementedError

    @staticmethod
    def collect_logfiles(job):
        """
        Collect the log files of the external executable and store the information in the HDF5 file. This method has
        to be implemented in the individual hamiltonians.
        """
        pass


class InteractiveInterface(object):
    @staticmethod
    def run_if_interactive(job):
        raise NotImplementedError

    @staticmethod
    def interactive_close(job):
        raise NotImplementedError

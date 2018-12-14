import numpy as np


class FileInterface(object):
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
    def __init__(self):
        self._interactive_library = None
        self._interactive_write_input_files = False
        self._interactive_flush_frequency = 1
        self._interactive_write_frequency = 1
        self.interactive_cache = {}

    @property
    def interactive_flush_frequency(self):
        return self._interactive_flush_frequency

    @interactive_flush_frequency.setter
    def interactive_flush_frequency(self, frequency):
        if not isinstance(frequency, int):
            raise AssertionError()
        self._interactive_flush_frequency = frequency

    @property
    def interactive_write_frequency(self):
        return self._interactive_write_frequency

    @interactive_write_frequency.setter
    def interactive_write_frequency(self, frequency):
        if not isinstance(frequency, int):
            raise AssertionError()
        self._interactive_write_frequency = frequency

    def run_if_interactive(self, job):
        raise NotImplementedError('run_if_interactive() is not implemented!')

    def run_if_interactive_non_modal(self, job):
        raise NotImplementedError('run_if_interactive_non_modal() is not implemented!')

    @staticmethod
    def _extend_hdf(h5, path, key, data):
        if path in h5.list_groups() and key in h5[path].list_nodes():
            current_hdf = h5[path + "/" + key]
            if isinstance(data, list):
                entry = current_hdf.tolist() + data
            else:
                entry = current_hdf.tolist() + data.tolist()
            data = np.array(entry)
        h5[path + "/" + key] = data

    @staticmethod
    def _include_last_step(array, step=1, include_last=False):
        if step == 1:
            return array
        if len(array) > 0:
            if len(array) > step:
                new_array = array[::step]
                index_lst = list(range(len(array)))
                if include_last and index_lst[-1] != index_lst[::step][-1]:
                    new_array.append(array[-1])
                return new_array
            else:
                if include_last:
                    return [array[-1]]
                else:
                    return []
        return []

    def interactive_flush(self, job, path="interactive", include_last_step=False):
        with job.project_hdf5.open("output") as h5:
            for key in self.interactive_cache.keys():
                if len(self.interactive_cache[key]) == 0:
                    continue
                data = self._include_last_step(array=self.interactive_cache[key],
                                               step=self.interactive_write_frequency,
                                               include_last=include_last_step)
                if len(data) > 0 and \
                        isinstance(data[0], list) and \
                        len(np.shape(data)) == 1:
                    self._extend_hdf(h5=h5, path=path, key=key, data=data)
                elif np.array(data).dtype == np.dtype('O'):
                    self._extend_hdf(h5=h5, path=path, key=key, data=data)
                else:
                    self._extend_hdf(h5=h5, path=path, key=key, data=np.array(data))
                self.interactive_cache[key] = []

    def interactive_is_activated(self):
        if self._interactive_library is None:
            return False
        else:
            return True

    def interactive_close(self, job):
        if len(list(self.interactive_cache.keys())) > 0 and \
                len(self.interactive_cache[list(self.interactive_cache.keys())[0]]) != 0:
            self.interactive_flush(path="interactive", include_last_step=True)
        job.project.db.item_update(job.runtime(), job._job_id)
        job.status.finished = True
        self.interactive_cache = {}

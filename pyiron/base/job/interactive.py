import numpy as np
from pyiron.base.job.generic import GenericJob


class InteractiveBase(GenericJob):
    def __init__(self, project, job_name):
        super(InteractiveBase, self).__init__(project, job_name)
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

    def _run_if_running(self):
        if self.server.run_mode.interactive:
            self.run_if_interactive()
        elif self.server.run_mode.interactive_non_modal:
            self.run_if_interactive_non_modal()
        else:
            super(InteractiveBase, self)._run_if_running()

    def _run_if_created(self, que_wait_for=None):
        if self._interactive_write_input_files:
            self.project_hdf5.create_working_directory()
            self.write_input()
            self._copy_restart_files()
            self._write_run_wrapper()
        super(InteractiveBase, self)._run_if_created(que_wait_for=que_wait_for)

    def interactive_is_activated(self):
        if self._interactive_library is None:
            return False
        else:
            return True

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

    def interactive_flush(self, path="interactive", include_last_step=False):
        with self.project_hdf5.open("output") as h5:
            for key in self.interactive_cache.keys():
                if len(self.interactive_cache[key])==0:
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

    def interactive_close(self):
        if len(list(self.interactive_cache.keys())) > 0 and \
                len(self.interactive_cache[list(self.interactive_cache.keys())[0]]) != 0:
            self.interactive_flush(path="interactive", include_last_step=True)
        self.project.db.item_update(self._runtime(), self._job_id)
        self.status.finished = True
        self.interactive_cache = {}

    def interactive_store_in_cache(self, key, value):
        self.interactive_cache[key] = value

    # def __del__(self):
    #     self.interactive_close()

    def run_if_interactive(self):
        raise NotImplementedError('run_if_interactive() is not implemented!')

    def run_if_interactive_non_modal(self):
        raise NotImplementedError('run_if_interactive_non_modal() is not implemented!')

    def to_hdf(self, hdf=None, group_name=None):
        """
        Store the InteractiveBase object in the HDF5 File

        Args:
            hdf (ProjectHDFio): HDF5 group object - optional
            group_name (str): HDF5 subgroup name - optional
        """
        super(InteractiveBase, self).to_hdf(hdf=hdf, group_name=group_name)
        with self.project_hdf5.open("input") as hdf5_input:
            hdf5_input["interactive"] = {"interactive_flush_frequency": self._interactive_flush_frequency,
                                         "interactive_write_frequency": self._interactive_write_frequency}

    def from_hdf(self, hdf=None, group_name=None):
        """
        Restore the InteractiveBase object in the HDF5 File

        Args:
            hdf (ProjectHDFio): HDF5 group object - optional
            group_name (str): HDF5 subgroup name - optional
        """
        super(InteractiveBase, self).from_hdf(hdf=hdf, group_name=group_name)
        with self.project_hdf5.open("input") as hdf5_input:
            if "interactive" in hdf5_input.list_nodes():
                interactive_dict = hdf5_input["interactive"]
                self._interactive_flush_frequency = interactive_dict["interactive_flush_frequency"]
                if "interactive_write_frequency" in interactive_dict.keys():
                    self._interactive_write_frequency = interactive_dict["interactive_write_frequency"]
                else:
                    self._interactive_write_frequency = 1

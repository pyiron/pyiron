from pyiron.base.job.generic import GenericJob
from pyiron.base.job.interface import InteractiveInterface


class InteractiveBase(GenericJob):
    def __init__(self, project, job_name):
        super(InteractiveBase, self).__init__(project, job_name)
        self._interface = InteractiveInterface()

    @property
    def interactive_flush_frequency(self):
        return self._interface.interactive_flush_frequency

    @interactive_flush_frequency.setter
    def interactive_flush_frequency(self, frequency):
        self._interface.interactive_flush_frequency = frequency

    @property
    def interactive_write_frequency(self):
        return self._interface.interactive_write_frequency

    @interactive_write_frequency.setter
    def interactive_write_frequency(self, frequency):
        self._interface.interactive_write_frequency = frequency

    def _run_if_running(self):
        if self.server.run_mode.interactive:
            self.run_if_interactive()
        elif self.server.run_mode.interactive_non_modal:
            self.run_if_interactive_non_modal()
        else:
            super(InteractiveBase, self)._run_if_running()

    def _run_if_created(self, que_wait_for=None):
        if self._interface._interactive_write_input_files:
            self.project_hdf5.create_working_directory()
            self.write_input()
            self._copy_restart_files()
            self._write_run_wrapper()
        super(InteractiveBase, self)._run_if_created(que_wait_for=que_wait_for)

    def interactive_is_activated(self):
        return self._interface.interactive_is_activated()

    def interactive_flush(self, path="interactive", include_last_step=False):
        self._interface.interactive_flush(job=self, path=path, include_last_step=include_last_step)

    def interactive_close(self):
        self._interface.interactive_close(job=self)

    def interactive_store_in_cache(self, key, value):
        self._interface.interactive_cache[key] = value

    # def __del__(self):
    #     self.interactive_close()

    def run_if_interactive(self):
        self._interface.run_if_interactive(job=self)

    def run_if_interactive_non_modal(self):
        self._interface.run_if_interactive_non_modal(job=self)

    def to_hdf(self, hdf=None, group_name=None):
        """
        Store the InteractiveBase object in the HDF5 File

        Args:
            hdf (ProjectHDFio): HDF5 group object - optional
            group_name (str): HDF5 subgroup name - optional
        """
        super(InteractiveBase, self).to_hdf(hdf=hdf, group_name=group_name)
        with self.project_hdf5.open("input") as hdf5_input:
            hdf5_input["interactive"] = {"interactive_flush_frequency": self._interface.interactive_flush_frequency,
                                         "interactive_write_frequency": self._interface.interactive_write_frequency}

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
                self._interface.interactive_flush_frequency = interactive_dict["interactive_flush_frequency"]
                if "interactive_write_frequency" in interactive_dict.keys():
                    self._interface.interactive_write_frequency = interactive_dict["interactive_write_frequency"]
                else:
                    self._interface.interactive_write_frequency = 1

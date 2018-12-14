import numpy as np
import posixpath
from pyiron.base.job.interface import JobInterface, InteractiveInterface
from pyiron.base.pyio.parser import Logstatus


class RandomInterface(JobInterface):
    @staticmethod
    def write_input(job):
        """
        Call routines that generate the codespecifc input files
        """
        job.input.write_file(file_name="input.inp", cwd=job.working_directory)

    def collect_output(self, job):
        """
        Parse the output files of the example job and store the results in the HDF5 File.
        """
        self._collect_output_log(project_hdf5=job.project_hdf5,
                                 working_directory=job.working_directory,
                                 file_name="output.log")
        self._collect_warnings(project_hdf5=job.project_hdf5,
                               working_directory=job.working_directory,
                               file_name="info.log")

    @staticmethod
    def _collect_output_log(project_hdf5, working_directory, file_name="output.log"):
        """
        general purpose routine to extract output from logfile

        Args:
            project_hdf5:
            working_directory:
            file_name (str): output.log - optional
        """
        tag_dict = {"alat": {"arg": "0", "rows": 0},
                    "count": {"arg": "0", "rows": 0},
                    "energy": {"arg": "0", "rows": 0}}
        lf = Logstatus()
        file_name = posixpath.join(working_directory, file_name)
        lf.extract_file(file_name=file_name, tag_dict=tag_dict)
        with project_hdf5.open("output/generic") as h5:
            lf.to_hdf(h5)
            h5["energy_tot"] = np.array(h5["energy"])
            h5["volume"] = np.array(h5["alat"])

    @staticmethod
    def _collect_warnings(project_hdf5, working_directory, file_name="info.log"):
        """
        Collect the warnings if any were written to the info.log file and store them in the HDF5 file
        """
        warnings_lst = []
        with open(posixpath.join(working_directory, file_name), "r") as f:
            lines = f.readlines()
        for line in lines:
            if "WARNING" in line:
                warnings_lst.append(line.split("WARNING"))
                warnings_lst[-1][-1] = warnings_lst[-1][-1].rstrip()
        if len(warnings_lst) > 0:
            warnings_dict = {'Module': [warnings_lst[i][0] for i in range(len(warnings_lst))],
                             'Message': [warnings_lst[i][1] for i in range(len(warnings_lst))]}
            print("module: ", warnings_lst[:][:])
            with project_hdf5.open("output") as hdf_output:
                hdf_output["WARNINGS"] = warnings_dict

    @staticmethod
    def collect_logfiles(job):
        """
        Collect the errors from the info.log file and store them in the HDF5 file
        """
        errors_lst = []
        with open(posixpath.join(job.working_directory, "info.log"), "r") as f:
            lines = f.readlines()
        for line in lines:
            if "ERROR" in line:
                errors_lst.append(line)
        if len(errors_lst) > 0:
            with job.project_hdf5.open("output") as hdf_output:
                hdf_output["ERRORS"] = errors_lst


class InteractiveRandomInterface(RandomInterface, InteractiveInterface):
    def __init__(self):
        super(InteractiveRandomInterface, self).__init__()
        self.interactive_cache = {'alat': [], 'count': [], 'energy': []}

    def run_if_interactive(self, job):
        """
        Run the job as Python library and store the result in the HDF5 File.

        Returns:
            int: job ID
        """
        from pyiron.testing.executable import ExampleExecutable
        job.status.running = True
        alat, count, energy = ExampleExecutable().run_lib(job.input)
        self.interactive_cache['alat'].append(alat)
        self.interactive_cache['count'].append(count)
        self.interactive_cache['energy'].append(energy)

    def interactive_close(self, job):
        job.to_hdf()
        with job.project_hdf5.open("output") as h5:
            h5["generic/energy"] = np.array(self.interactive_cache['energy'])
            h5["generic/volume"] = np.array(self.interactive_cache['alat'])
            h5["generic/alat"] = np.array(self.interactive_cache['alat'])
            h5["generic/count"] = np.array(self.interactive_cache['count'])
            h5["generic/energy_tot"] = np.array(self.interactive_cache['energy'])
        job.project.db.item_update(job.runtime(), job.job_id)
        job.status.finished = True

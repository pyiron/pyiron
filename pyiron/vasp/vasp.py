# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from pyiron.vasp.interactive import VaspInteractive

__author__ = "Sudarsan Surendralal"
__copyright__ = (
    "Copyright 2020, Max-Planck-Institut für Eisenforschung GmbH - "
    "Computational Materials Design (CM) Department"
)
__version__ = "1.0"
__maintainer__ = "Sudarsan Surendralal"
__email__ = "surendralal@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2017"


class Vasp(VaspInteractive):
    """
    Class to setup and run and analyze VASP simulations which is a derivative of pyiron.objects.job.generic.GenericJob.
    The functions in these modules are written in such the function names and attributes are very generic
    (get_structure(), molecular_dynamics(), version) but the functions are written to handle VASP specific input/output.

    Args:
        project (pyiron.project.Project instance):  Specifies the project path among other attributes
        job_name (str): Name of the job

    Examples:
        Let's say you need to run a vasp simulation where you would like to control the input parameters manually. To
        set up a static dft run with Gaussian smearing and a k-point MP mesh of [6, 6, 6]. You would have to set it up
        as shown below:

        >>> ham = Vasp(job_name="trial_job")
        >>> ham.input.incar[IBRION] = -1
        >>> ham.input.incar[ISMEAR] = 0
        >>> ham.input.kpoints.set_kpoints_file(size_of_mesh=[6, 6, 6])

        However, the according to pyiron's philosophy, it is recommended to avoid using code specific tags like IBRION,
        ISMEAR etc. Therefore the recommended way to set this calculation is as follows:

        >>> ham = Vasp(job_name="trial_job")
        >>> ham.calc_static()
        >>> ham.set_occupancy_smearing(smearing="gaussian")
        >>> ham.set_kpoints(mesh=[6, 6, 6])
        The exact same tags as in the first examples are set automatically.

    """

    def __init__(self, project, job_name):
        super(Vasp, self).__init__(project, job_name)
        self.__name__ = "Vasp"
        self.__version__ = (
            None
        )  # Reset the version number to the executable is set automatically
        self._executable_activate(enforce=True)

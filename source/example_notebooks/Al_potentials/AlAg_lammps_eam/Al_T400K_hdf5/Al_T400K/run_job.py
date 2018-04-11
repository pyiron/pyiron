import sys
from pyiron_base.objects.job.wrapper import JobWrapper

debug = False
job = JobWrapper(working_directory='/cmmc/u/mpt13/PyIron_data/website/pyiron_docs/source/example_notebooks/Al_potentials/AlAg_lammps_eam/Al_T400K_hdf5/Al_T400K',
                 job_id= 779484 ,
                 debug=False )
job.run()

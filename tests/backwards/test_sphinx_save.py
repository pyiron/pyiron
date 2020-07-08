from pyiron import Project, __version__
pr = Project(
    "tests/static/backwards/V{}".format(__version__).replace(".", "_")
)
structure = pr.create_ase_bulk("Al")
job = pr.create_job(pr.job_type.Sphinx, "sphinx")
job.structure = structure
job.save()

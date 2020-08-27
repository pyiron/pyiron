from pyiron import Project, __version__
pr = Project(
    "tests/static/backwards/V{}".format(__version__).replace(".", "_")
)
structure = pr.create_ase_bulk("Al")
job = pr.create_job(pr.job_type.StructureContainer, "structurecontainer")
job.structure = structure
job.save()

from pyiron import Project, __version__
pr = Project("tests/static/backwards")
for job in pr.iter_jobs(recursive = True):
    if job.name == "structurecontainer":
        job.run()

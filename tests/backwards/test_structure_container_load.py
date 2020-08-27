from pyiron import Project, __version__
pr = Project("tests/static/backwards")
for job in pr.iter_jobs(recursive = True, convert_to_object = False):
    if job.name == "structure_container":
        job = job.to_object()
        job.run()

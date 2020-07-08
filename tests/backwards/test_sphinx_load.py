from pyiron import Project, __version__
pr = Project("tests/static/backwards/")
for j in pr.iter_jobs(recursive = True):
    if j.name == "sphinx": j.run()

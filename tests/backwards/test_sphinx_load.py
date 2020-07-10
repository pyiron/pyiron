import sys
from pyiron import Project, __version__
pr = Project("tests/static/backwards/")
for job in pr.iter_jobs(recursive = True):
    if job.name == "sphinx": job.run()
    print("job {} loaded from {}".format(job.id, sys.argv[0]))

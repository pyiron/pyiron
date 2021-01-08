import sys
from pyiron_atomistic import Project, __version__
pr = Project("tests/static/backwards/")
for job in pr.iter_jobs(recursive = True, convert_to_object = False):
    if job.name == "sphinx":
        job = job.to_object()
        job.run()
    print("job {} loaded from {}".format(job.id, sys.argv[0]))

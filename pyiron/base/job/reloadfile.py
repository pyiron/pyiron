import os
import sys
import getopt
import h5py
import shutil
from pyiron import Project
from pyiron.base.settings.generic import Settings


s = Settings()


def command_line(argv):
    """
    Parse the command line arguments.

    Args:
    argv: Command line arguments

    """
    input_path = None
    output_path = None
    try:
        opts, args = getopt.getopt(argv, "i:o:h", ["input_path=", "output_path=", "help"])
    except getopt.GetoptError:
        print('python -m pyiron.base.job.reloadfile -i <inputpath> -o <outputpath>')
        sys.exit()
    else:
        for opt, arg in opts:
            if opt in ("-i", "--input_path"):
                input_path = os.path.abspath(arg)
            if opt in ("-o", "--output_path"):
                output_path = os.path.abspath(arg)
        with h5py.File(input_path, 'r') as f:
            job_name = list(f.keys())[0]
        project_path = os.path.join(os.path.dirname(input_path), job_name + '.h5')
        shutil.copy(input_path, project_path)

        file = os.path.basename(project_path)
        job_name = os.path.splitext(file)[0]

        db_protject_path = s.top_path(project_path)
        project = os.path.dirname(project_path)
        db_project = (project + '/').replace(db_protject_path, '')
        job_reload = Project(project).load_from_jobpath(job_id=None,
                                                        db_entry={'id': 1000,
                                                                  'status': '',
                                                                  'chemicalformula': '',
                                                                  'job': job_name,
                                                                  'subjob': '/' + job_name,
                                                                  'projectpath': db_protject_path,
                                                                  'project': db_project,
                                                                  'hamilton': '',
                                                                  'hamversion': '',
                                                                  'parentid': None,
                                                                  'masterid': None},
                                                        convert_to_object=True)
        job_reload.status.initialized = True
        job_reload.server.run_mode.modal = True
        job_reload.run()
        shutil.copy(project_path, output_path)
        sys.exit()


if __name__ == "__main__":
    command_line(sys.argv[1:])

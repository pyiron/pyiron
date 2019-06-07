import sys
import getopt
from pyiron.base.job.wrapper import job_wrapper_function


def command_line(argv):
    """
    Parse the command line arguments.

    Args:
        argv: Command line arguments

    """
    debug = False
    project_path = None
    job_id = None
    try:
        opts, args = getopt.getopt(argv, "dj:p:h", ["debug", "project_path=", "job_id=", "help"])
    except getopt.GetoptError:
        print('cms.py --p <project_path> -j <job_id> <debug>')
        sys.exit()
    else:
        for opt, arg in opts:
            if opt in ("-h", "--help"):
                print('cms.py --p <project_path> -j <job_id> <debug>')
                sys.exit()
            elif opt in ("-d", "--debug"):
                debug = True
            elif opt in ("-j", "--job_id"):
                job_id = arg
            elif opt in ("-p", "--project_path"):
                project_path = arg
        job_wrapper_function(working_directory=project_path, job_id=job_id, debug=debug)
        sys.exit()


if __name__ == "__main__":
    command_line(sys.argv[1:])

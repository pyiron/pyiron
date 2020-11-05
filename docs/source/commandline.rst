.. _commandline:

======================
Command Line Interface
======================

Usage Summary
=============

There's a few command line tools shipped with pyiron to help
administrating and keeping up with your pyiron project as well as some that are
used internally.  All of them are installed by default in the *pyiron* script
that has a few sub commands.

*pyiron install*
    Installs the pyiron resources for the first time, if you don't get them via
    conda.

*pyiron ls*
    list the jobs inside a project and filter them with a few primitives

    Print the run time of all finished jobs
    
        pyiron ls -c job totalcputime -s finished

    Print all jobs with iron
    
        pyiron ls -e Fe

    Print all jobs that successfully finished yesterday and a bit
    
        pyiron ls -s finished -i 1d5h

    Print all jobs that were aborted less than 5 hours ago and match
    "spx.*restart"
    
        pyiron ls -n "spx.*restart" -i 5h -s aborted

*pyiron rm*
    Delete jobs and whole projects from the database and the file system.  If
    you simply *rm* jobs and projects they are still in the database and can
    lead to confusion on pyiron's part.

*pyiron wrapper*
    Runs jobs from the database.  pyiron uses this internally to start jobs on
    the remote cluster nodes, but you can also use it when you set the run mode
    to "manual" or to manually re-run jobs.


Developer Guide
===============

Adding a new sub command is done by adding a new module to ``pyiron.cli``.
This module needs to define a ``register`` and a ``main`` function.  The
former is called with an ``argparse.ArgumentParser`` instance as sole argument
and should define the command line interface in the `usual way
<https://docs.python.org/3/library/argparse.html>`_.  The latter will be called
with the parsed arguments and should just execute whatever it is that utility
should be doing.  Additionally if you need to control the ``formatter_class``
and ``epilog`` keyword arguments when creating the ``argparse.ArgumentParser``
instance you can set the ``formatter`` and ``epilog`` toplevel variables (see
the *ls* sub command for an example).  Finally you must add the module to the
``pyiron.cli.cli_modules`` dict.

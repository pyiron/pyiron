======================
Contributing to pyiron
======================

The following is a set of guidelines for contributing to pyiron, which is
hosted and maintained by the `Max Planck Institut für Eisenforschung`_
on GitHub. These are mostly guidelines to facilitate an efficient
development workflow, and not necessarily rules. Use your best judgment,
and feel free to propose changes even to this document in a pull request.

You can find all the pyiron packages at our `github page`_ .
To create pull requests, you will need to become part of the
pyiron organization. Please email us if you would like to join.

Wait I don't want to read this; I just have a quick question/bugfix!
====================================================================

1. Check out our `FAQ page`_; your question might already be answered there.
2. If your question relates to a bug in pyiron, please briefly search the `issues page`_ and open a new labeled issue if you don't see anything related to your question there.
3. Please feel free just to send one of us a brief, descriptive email with your question, and we'll do our best to get back to you as ASAP as possible.

Table of Contents
=================

`License`_

`What should I know before I get started?`_
  * `pyiron developer meetings`_

..
 * `The structure of pyiron`_
..
 * `The principles of pyiron`_


`How can I contribute?`_
  * `Reporting bugs`_
  * `Suggesting enhancements`_
  * `Your first code contribution`_
  * `Pull requests`_

`Styleguides`_
  * `Git commit messages`_
  * `Python styleguide`_
  * `Documentation styleguide`_

`Additional Notes`_
  * `Issue and pull request labels`_
  * `Build status`_
  * `pyiron releases`_
  
`Debugging`_
  * `My job does not run on the queue`_

License
=======
pyiron is released as an open-source project under the BSD 3-Clause License.
Code contributions should also be considered open-source.

What should I know before I get started?
========================================

.. The structure of pyiron
.. -----------------------

.. The principles of pyiron
.. ------------------------

pyiron developer meetings
-------------------------
If you are interested in discussing pyiron's development, we encourage you to virtually
participate in the weekly pyiron developer meeting at 14:00 german time (GMT+2).
Check the discussion page for details.

How can I contribute?
=====================

Reporting bugs
--------------

    Note: If you find a closed issue that seems like it is the same
    thing that you're experiencing, open a new issue and include a
    link to the original issue in the body of your new one.

**Before Submitting A Bug Report**

Check if you can reproduce the problem in the latest version of pyiron.
Check the `FAQ page`_ for a list of common questions and problems.
Briefly search the issues page for `bugs`_  to see if the problem has already
been reported. If it has and the issue is still open, add a comment
to the existing issue instead of opening a new one.

**How Do I Submit A (Good) Bug Report?**

Bugs are tracked as GitHub issues. You can create an issue on
the pyiron repository by including the following information:

* Use a clear and descriptive title for the issue to identify the problem.
* Describe the exact steps you took so we can reproduce the problem as closely as possible.
* Provide sample code that causes the problem. Include code snippets as markdown code blocks.
* Include information about the environment (OS, python version, how packages were installed) in which you were running pyiron.
* Explain what you expected to happen, and what happened instead.

Suggesting Enhancements
-----------------------

**How Do I Submit A (Good) Enhancement Suggestion?**

Enhancement suggestions are tracked as GitHub issues. You can create an issue on
the pyiron repository by including the following information:

* Use a clear and descriptive title for the issue to identify the suggestion.
* Describe the exact behavior you would expect the suggested feature to produce.
* Provide sample code that you would use to access the feature. If possible, include code for how you think the feature could be built into pyiron's codebase. Include code snippets as markdown code blocks.

Your first code contribution
----------------------------

Unsure where to begin contributing to pyiron? You can start by looking
through these good-first-issue and help-wanted issues:

* `Good first issues`_ - issues which should only require a few lines of code, and a test or two.
* `Help wanted issues`_ - issues which should be a bit more involved than beginner issues.

**Local development**

pyiron can be developed and tested locally. If you are using pyiron to run an
external software package, e.g. VASP or LAMMPS, you might also need to install
those packages locally to run certain integration tests in pyiron.

To get the developmental (git) version of pyiron,

.. code-block::

  git clone https://github.com/pyiron/pyiron.git
  conda env update --name pyiron_dev --file pyiron/.ci_support/environment.yml
  conda activate pyiron_dev
  conda install conda-build
  conda develop pyiron

**Local Testing**

The full test suite is always run automatically when you open a new pull request.  Still it 
sometimes nice to run all or only specific tests on your machine.  To do that run from the repo root, e.g.

.. code-block::

  python -m unittest discover tests
  python -m unittest discover tests/sphinx
  python -m unittest tests/sphinx/test_base.py

Where the first line runs all tests, the second all the sphinx tests and the final line only the tests in that file.
Keep in mind that to run the tests your repository needs to be inside your pyiron project folder and you need to have 
at least the basic resources installed from ``tests/static``.  A neat trick when testing/debugging is to combine the 
pdb and unittest modules like this

.. code-block::

  python -m pdb -m unittest ...
  
This allows you to re-use the sometimes complicated setups for your interactive debugging that might be otherwise
difficult to replicate in a REPL.

Pull requests
-------------

The process described here has several goals:

* Maintain pyiron's quality
* Fix problems that are important to users
* Engage the community in working toward the best possible tools
* Enable a sustainable system for pyiron's maintainers to review contributions

Please follow these steps to have your contribution considered by the maintainers:

* Keep the changes in your pull request as focused as possible- only address one issue per pull request wherever possible.
* Follow the `Styleguides`_
* Assign the appropriate label (see `Issue and pull request labels`_) to your pull request. If you are fixing a specific Github issue, reference the issue directly in the pull request comments.
* If you are aware which maintainer is most closely related to the code you've edited, feel free to request their review.
* After you submit your pull request, verify that all status checks are passing.
* If a status check fails and it seems to be unrelated to your changes, explain why the failure is unrelated as a comment in your pull request.

While the prerequisites above must be satisfied prior to having your
pull request reviewed, the reviewer(s) may ask you to complete
additional design work, tests, or other changes before your pull
request can be ultimately accepted.

Styleguides
===========

Git commit messages
-------------------

* Use the present tense ("Add feature" not "Added feature")
* Use the imperative mood ("Move cursor to..." not "Moves cursor to...")
* Limit the first line to 72 characters or less
* Reference issues and pull requests liberally after the first line
* When only changing documentation, include [ci skip] in the commit title
* Consider starting the commit message with an applicable emoji:

\:art: (``:art:``) improves the format/structure of the code

\:zap: (``:zap:``) improves performance

\:memo: (``:memo:``) adds documentation

\:bug: (``:bug:``) fixes a bug

\:fire: (``:fire:``) removes code or files

\:green_heart: (``:green_heart:``) fixes the CI build

\:white_check_mark: (``:white_check_mark:``) adds tests

Managing git commits is much easier using an IDE (we recommend PyCharm).

Python styleguide
-----------------

Please follow `PEP8 conventions`_ for all python code added to pyiron. Pull
requests will be checked for PEP8 plus a few other security issues with
`Codacy`_, and will be rejected if they do not meet the specified
formatting criteria.

Any new features should include coverage with a unit test, such that
your pull request does not decrease pyiron's overall coverage. This
will be automatically tested within the ci test suite and `Coveralls`_.

Deprecation warning template
----------------------------
*XXX is deprecated as of vers. A.B.C. It is not guaranteed to be in service in vers. D.E.F*

Documentation styleguide
------------------------

All new/modified functions should include a docstring that follows
the `Google Python Docstring format`_.

Documentation is built automatically with `Sphinx`_; any manually created
documentation should be added as a restructured text (.rst) file
under pyiron/docs/source.

Notebooks created to exemplify features in pyiron are very useful, and
can even be used as integration tests. If you have added a major feature,
consider creating a notebook to show its usage under pyiron/notebooks/.
See the other examples that are already there.

Additional notes
================

Issue and pull request labels
-----------------------------

We use the following tags to organize pyiron Github issues
and pull requests:

* bug: something isn't working
* duplicate: this issue/pull request already existed
* enhancement: new feature or request
* good first issue: easy fix for beginners
* help wanted: extra attention is needed
* invalid: this doesn't seem right
* question: further information is requested
* wontfix: this will not be worked on
* stale: inactive after 2 weeks

Build status
------------

The build status for pyiron and all sub packages are given below

.. image:: https://coveralls.io/repos/github/pyiron/pyiron/badge.svg?branch=master
    :target: https://coveralls.io/github/pyiron/pyiron?branch=master
    :alt: Coverage Status

.. image:: https://api.codacy.com/project/badge/Grade/c513254f10004df5a1f5c76425c6584b
    :target: https://app.codacy.com/app/pyiron-runner/pyiron?utm_source=github.com&utm_medium=referral&utm_content=pyiron/pyiron&utm_campaign=Badge_Grade_Settings
    :alt: Codacy Badge

.. image:: https://anaconda.org/conda-forge/pyiron/badges/latest_release_date.svg
    :target: https://anaconda.org/conda-forge/pyiron/
    :alt: Release_Date

.. image:: https://travis-ci.org/pyiron/pyiron.svg?branch=master
    :target: https://travis-ci.org/pyiron/pyiron
    :alt: Build Status

.. image:: https://ci.appveyor.com/api/projects/status/wfdgqkxca1i19xcq/branch/master?svg=true
    :target: https://ci.appveyor.com/project/pyiron-runner/pyiron/branch/master
    :alt: Build status

pyiron releases
---------------

.. image:: https://anaconda.org/conda-forge/pyiron/badges/downloads.svg
    :target: https://anaconda.org/conda-forge/pyiron/
    :alt: Downloads

For the pyiron release management we use git tags::

   https://git-scm.com/book/en/v2/Git-Basics-Tagging

The tag format consists of a tag_prefix (<package name>-) and the release version, for example::

   pyiron-0.2.0

For the automated versioning we use::

   https://github.com/warner/python-versioneer/

So the configuration of the release is included in setup.cfg::

   https://github.com/pyiron/pyiron_base/blob/master/setup.cfg

As the pyiron packages are pure python packages – we use only the Linux Python 3.7 job to build the packages, as defined in the .travis.yml file::

   https://github.com/pyiron/pyiron_base/blob/master/.travis.yml

The python 3.7 linux tests therefore takes more time, compared to the other tests on travis.

Just like each other commit to the master branch the tagged releases are pushed to pypi.org and anaconda.org::

   https://pypi.org/project/pyiron-base/#history
   https://anaconda.org/pyiron/pyiron_base

The major difference for pypi (pip) is that tagged releases are the default for pip while installing prerelease versions using pip requires the `--pre` flag.
`pip install --pre pyiron`

Those pre-release versions are named `<version_number>.post0.dev<release number>` ::

   0.2.0.post0.dev1

For anaconda the prereleases are pushed to the pyiron channel and can be installed using:
conda install -c pyiron pyiron

On the other hand the tagged releases are available through conda-forge, as soon as the corresponding packages are merged::

   https://github.com/conda-forge/pyiron-feedstock
   conda install -c conda-forge pyiron

So for both conda and pip both the prereleases as well as the official releases are available.

.. _Max Planck Institut für Eisenforschung: https://mpie.de
.. _github page: https://github.com/pyiron
.. _issues page: https://github.com/pyiron/pyiron/issues
.. _FAQ page: https://github.com/pyiron/pyiron/docs/source/faq.html
.. _bugs: https://github.com/pyiron/pyiron/issues?q=is%3Aopen+is%3Aissue+label%3A%22bug%22
.. _Good first issues: https://github.com/pyiron/pyiron/issues?q=is%3Aopen+is%3Aissue+label%3A%22good+first+issue%22
.. _Help wanted issues: https://github.com/pyiron/pyiron/issues?q=is%3Aissue+is%3Aopen+label%3A%22help+wanted%22
.. _PEP8 conventions: https://www.python.org/dev/peps/pep-0008/
.. _Codacy: https://www.codacy.com/
.. _Coveralls: https://coveralls.io/
.. _Google Python Docstring format: http://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html
.. _Sphinx: https://www.sphinx-doc.org/en/master/

Debugging
================
My job does not run on the queue
-----------------------------

In case a job runs properly while executing it locally (or on the head node), but not when you submit it to a queue,

**1. Check if the job class is available in the project:**

In this example, we want a custom job class ``ProtoMD`` from the module ``pyiron_contrib``:

.. code-block::

  from pyiron import Project
  import pyiron_contrib  # only if importing a custom job class

  pr = Project("debug")
  dir(pr.job_type)

This should output:

.. code-block::

  >>> ['AtomisticExampleJob',
   'Atoms',
   'ConvEncutParallel',
   ...
   ...
   'ProtoMD']

If you see your job class in the list, proceed to step 3. If not, 


**2. Check if the job class in initialized in ``__init__.py`` of the module**

Make sure that the ``__init__.py`` of your module (here, ``pyiron_contrib``) initializes the job class in the following format:

.. code-block::

  from pyiron import Project
  from pyiron.base.job.jobtype import JOB_CLASS_DICT

  # Make classes available for new pyiron version
  JOB_CLASS_DICT['ProtoMD'] = 'pyiron_contrib.protocol.compound.md'  # the path of your job class
  
  
**3. Confirm that the job class can be instantiatied**

Create a new job, but instead of running it, save it:

.. code-block::

  job = pr.create_job(job_type = pr.job_type.ProtoMD, job_name = 'job')
  ...  # input parameters that the job requires
  ...
  job.save()

  >>> 98  # this is the job id of the saved job

Note down the job id, then run the following line:

.. code-block::

  job["TYPE"]

This should output an instance of the job class:

.. code-block::

  >>> "<class 'pyiron_contrib.protocol.compound.md.ProtoMD'>"

Now we know that the job class is indeed available in the project and can be instantiated.

**4. Debug using a second notebook**

Submitting and running a job on the queue, is essentially the same as saving a job in one notebook, but loading and executing it in another notebook.

In **a new notebook** , load the job that you just saved, using its job id. You may or may not import the module (here, ``pyiron_conntirb``):

.. code-block::

  from pyiron import Project
  # we do not import pyiron_contrib here, becasue it should not be necessary

  pr = Project("second_notebook")
  reloaded_job = pr.load(98)  # 98 is the job id of the previously saved job
  reloaded_job.run(run_again=True)

If the job loads and runs properly, the job should also run properly on the queue. This also means that there may be a bug in your custom job class. Debug the job class, and repeat steps 3 and 4 till you no longer get an error in step 4.

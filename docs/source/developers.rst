.. _developers:


==============
For developers
==============

If you are a pyiron developer or would like to contribute to pyiron, this page is for you.


Find us on GitHub
=================

You can find all the pyiron packages at out `github page`_ . Please contact us if you interested in joining
the organization.


Documentation
=============

Individual functions and modules are documented in the `Google Style Python Docstrings format`_. The webpages are
designed and built using the Restructured text (rst) format with the `sphinx document generator`_ which generates rst
files from the docstrings and jupyter notebooks. The rst files are available in the `docs`_ repository. The
example notebooks from which the tutorial pages are generated are in the `notebooks`_ directory.

Testing and CI
==============

Ideally, every module has unittests which have to pass for all platforms. In github we sue the TracisCI for Linux and
MacOS and Appveyor for Windows builds


Build Status
============

The build status for pyiron and all sub packages are given below

pyiron
------
.. image:: https://anaconda.org/pyiron/pyiron/badges/latest_release_date.svg

.. image:: https://travis-ci.org/pyiron/pyiron.svg?branch=master

.. image:: https://anaconda.org/pyiron/pyiron/badges/downloads.svg

.. image:: https://ci.appveyor.com/api/projects/status/wfdgqkxca1i19xcq/branch/master?svg=true



pyiron release
==============

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

.. toctree::
   :maxdepth:3

.. _github page: https://github.com/pyiron
.. _Google Style Python Docstrings format: http://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html
.. _sphinx document generator: http://www.sphinx-doc.org/en/stable/index.html
.. _pyiron_docs: https://github.com/pyiron/pyiron_docs
.. _examples: https://github.com/pyiron/examples

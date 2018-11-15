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
files from the docstrings and jupyter notebooks. The rst files are available in the `pyiron_docs`_ repository. The
example notebooks from which the tutorial pages are generated are in the `examples`_ directory.

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

pyiron_base
-----------

.. image:: https://travis-ci.org/pyiron/pyiron_base.svg?branch=master

.. image:: https://ci.appveyor.com/api/projects/status/c9w3tjyffnw1d47x/branch/master?svg=true


pyiron_atomistics
-----------------

.. image:: https://travis-ci.org/pyiron/pyiron_atomistics.svg?branch=master

.. image:: https://ci.appveyor.com/api/projects/status/57f61ea4t01l1rqg/branch/master?svg=true


pyiron_dft
----------

.. image:: https://travis-ci.org/pyiron/pyiron_dft.svg?branch=master

.. image:: https://ci.appveyor.com/api/projects/status/tu2owtwrmjsh98yr/branch/master?svg=true


pyiron_vasp
----------

.. image:: https://travis-ci.org/pyiron/pyiron_vasp.svg?branch=master

.. image:: https://ci.appveyor.com/api/projects/status/h7w6b1m3pd7hc4n9/branch/master?svg=true



pyiron_lammps
-------------

.. image:: https://travis-ci.org/pyiron/pyiron_lammps.svg?branch=master

.. image:: https://ci.appveyor.com/api/projects/status/hhcy3dmjy6ffdy53/branch/master?svg=true


pyiron_example_job
------------------

.. image:: https://travis-ci.org/pyiron/pyiron_example_job.svg?branch=master

.. image:: https://ci.appveyor.com/api/projects/status/4gs490vgif1bl0v5/branch/master?svg=true


pyiron release
==============

For the pyiron release management we use git tags::

   https://git-scm.com/book/en/v2/Git-Basics-Tagging

The tag format consists of a tag_prefix (<package name>-) and the release version, for example::

   pyiron_base-0.1.6

For the automated versioning we use::

   https://github.com/warner/python-versioneer/

So the configuration of the release is included in setup.cfg::

   https://github.com/pyiron/pyiron_base/blob/master/setup.cfg

As the pyiron packages are pure python packages – we use only the Linux Python 3.6 job to build the packages, as defined in the .travis.yml file::

   https://github.com/pyiron/pyiron_base/blob/master/.travis.yml

The python 3.6 linux tests therefore takes more time, compared to the other tests on travis.

Just like each other commit to the master branch the tagged releases are pushed to pypi.org and anaconda.org::

   https://pypi.org/project/pyiron-base/#history
   https://anaconda.org/pyiron/pyiron_base

The major difference for pypi (pip) is that tagged releases are the default for pip while installing prerelease versions using pip requires the `--pre` flag.
`pip install --pre pyiron_base`

Those pre-release versions are named `<version_number>.post0.dev<release number>` ::

   0.1.6.post0.dev11

For anaconda the prereleases are pushed to the pyiron channel and can be installed using:
conda install -c pyiron pyiron_base

On the other hand the tagged releases are available through conda-forge, as soon as the corresponding packages are merged::


   https://github.com/conda-forge/pyiron_base-feedstock
   conda install -c conda-forge pyiron_base

So for both conda and pip both the prereleases as well as the official releases are available.

.. toctree::
   :maxdepth:3

.. _github page: https://github.com/pyiron
.. _Google Style Python Docstrings format: http://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html
.. _sphinx document generator: http://www.sphinx-doc.org/en/stable/index.html
.. _pyiron_docs: https://github.com/pyiron/pyiron_docs
.. _examples: https://github.com/pyiron/examples

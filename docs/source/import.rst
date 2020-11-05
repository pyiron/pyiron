.. _import:

====================================
Importing finished VASP calculations
====================================

Finished VASP calculations that were created outside of pyiron can be imported using the following script::

    from pyiron.project import Project
    pr = Project('imported_jobs')
    # Searches and imports vasp jobs from 'vasp_directory'
    path_to_import = "vasp_directory"
    pr.import_from_path(path=path_to_import, recursive=True)

The calculations are imported into the project 'imported_jobs'. The recursive function imports vasp directories within
each vasp directory if present.

.. note::  This functionality best works when both the vasprun.xml and OUTCAR files are present in the directories. The
   import would work only id the vasprun.xml file exists too. If the vasprun.xml file does not exist, the OUTCAR and CONTCAR
   files must be present

.. toctree::
   :maxdepth:4

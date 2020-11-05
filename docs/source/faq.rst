.. _faq:

===
FAQ
===

How to cite pyiron?
===================
To cite pyiron and the corresponding codes, please follow the instructions on the `publication page <citation>`_.

What units does pyiron use?
===========================
* mass = atomic mass units
* distance = Angstroms
* time = femtoseconds
* energy = eV
* velocity = Angstroms/femtoseconds
* force = eV/Angstrom
* temperature = Kelvin
* pressure = GPa
* charge = multiple of electron charge (1.0 is a proton)

How to import existing calculation?
===================================

How to import structures from files or existing databases?
==========================================================

How to install pyiron?
======================
pyiron is designed to be installed as centralized service on your local computer cluster, rather than a local installation
on each individual workstation. To test pyiron online or with a local installation, please follow the instructions on the
`installation page <installation>`_.

How to use a custom Pseudo potential in VASP?
=============================================

How to use VASP tags which are not supported by pyiron?
=======================================================

How to use a custom potential in LAMMPS?
========================================
A custom empirical potential (here, a hybrid potential) can be defined in the following format:

.. code-block::

  custom_potential = pd.DataFrame({    
    'Name': ['SrTiO3_Pedone'],
    'Filename': [[]],
    'Model': ['Custom'],
    'Species': [['O', 'Sr', 'Ti']],
    'Config': [['atom_style full\n',  # I use 'full' here as atom_style 'charge' gives the same result
                '## create groups ###\n', 
                'group O type 1\n', 
                'group Sr type 2\n', 
                'group Ti type 3\n', 
                '\n', 
                '## set charges - beside manually ###\n', 
                'set group O charge -1.2000\n',  
                'set group Sr charge 1.2000\n',
                'set group Ti charge 2.4000\n',
                '\n', 
                'pair_style hybrid/overlay morse 15.0 mie/cut 15.0 coul/long 15.0 beck 15.0\n', 
                'pair_coeff * * coul/long\n', 
                'pair_coeff 1 2 beck 3.0 0 0 0 0\n', 
                'pair_coeff 1 3 beck 1.0 0 0 0 0\n', 
                'pair_coeff 1 1 beck 22.0 0 0 0 0\n', 
                'pair_coeff 1 2 mie/cut 3.0 1.0 12.0 0\n', 
                'pair_coeff 1 3 mie/cut 1.0 1.0 12.0 0\n', 
                'pair_coeff 1 1 mie/cut 22.0 1.0 12.0 0\n', 
                'pair_coeff 1 2 morse 0.019623 1.8860 3.32833\n', 
                'pair_coeff 1 3 morse 0.024235 2.2547 2.708943\n', 
                'pair_coeff 1 1 morse 0.042395 1.3793 3.618701\n', 
                'kspace_style ewald 1.0e-8\n']]
  })
  
The lines in ``Config`` will be written to the LAMMPS ``potential.inp`` file. Make sure that the arrangement of the species in ``Species`` is the same as the group types ``create groups`` within ``Config``. Otherwise, a mixup or the species may occur in the LAMMPS ``structure.inp`` file.

The potential can then be used by assigning ``job.potential = custom_potential``.

How to extend the potential database inside pyiron?
===================================================

How to link your own executable?
================================

How to send a calculation to the background ?
=============================================

How to submit a calculation to the queuing system?
==================================================

How to setup spin constraint calculation?
=========================================

What is the meaning of the name - pyiron?
=========================================
pyiron is the combination of **py** + **iron** connecting Python, the programming language with iron as pyiron was
initially developed at the Max Planck Institut für Eisenforschung (iron research).

Which output quantities are stored in pyiron?
=============================================

.. include:: outputformats.rst


.. toctree::
   :maxdepth:2

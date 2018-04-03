"""
Setuptools based setup module
"""
from setuptools import setup, find_packages


setup(
    name='pyiron',
    version='0.0.8',
    description='pyiron IDE base',
    long_description='http://pyiron.org',

    url='https://github.com/jan-janssen/pyirondist',
    author='Jan Janssen (MPIE)',
    author_email='janssen@mpie.de',
    license='BSD',

    classifiers=[

        'Development Status :: 4 - Beta',
        # Indicate who your project is intended for
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Physics',

        # Pick your license as you wish (should match "license" above)
        'License :: OSI Approved :: BSD License',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 3.5',
    ],

    keywords='pyiron',
    packages=find_packages(exclude=["tests"]),
    install_requires=['pyiron_base',
                      'pyiron_atomistics',
                      'pyiron_dft',
                      'pyiron_example_job',
                      'pyiron_lammps',
                      'pyiron_vasp']
    )

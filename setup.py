"""
Setuptools based setup module
"""
from setuptools import setup, find_packages
import versioneer

setup(
    name='pyiron',
    version=versioneer.get_version(),
    description='pyiron - an integrated development environment (IDE) for computational materials science.',
    long_description='http://pyiron.org',

    url='https://github.com/pyiron/pyiron',
    author='Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department',
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
    packages=find_packages(exclude=["*tests*", "*docs*"]),
    install_requires=['ase',
                      'h5io>=0.1.1',
                      'h5py',
                      'matplotlib',
                      'numpy',
                      'pandas',
                      'pathlib2',
                      'phonopy',
                      'psutil',
                      'six',
                      'sqlalchemy',
                      'tables',
                      'spglib'],
    cmdclass=versioneer.get_cmdclass(),
    )

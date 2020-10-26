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

    classifiers=['Development Status :: 5 - Production/Stable',
                 'Topic :: Scientific/Engineering :: Physics',
                 'License :: OSI Approved :: BSD License',
                 'Intended Audience :: Science/Research',
                 'Operating System :: OS Independent',
                 'Programming Language :: Python :: 3',
                 'Programming Language :: Python :: 3.4',
                 'Programming Language :: Python :: 3.5',
                 'Programming Language :: Python :: 3.6',
                 'Programming Language :: Python :: 3.7',
                 'Programming Language :: Python :: 3.8'],

    keywords='pyiron',
    packages=find_packages(exclude=["*tests*", "*docs*", "*binder*", "*conda*", "*notebooks*", "*.ci_support*"]),
    install_requires=[
        'ase>=3.19',
        'defusedxml>=0.6.0',
        'dill>=0.3.1.1',
        'future>=0.18.2',
        'gitpython>=3.1.0',
        'h5io>=0.1.1',
        'h5py>=2.10.0',
        'matplotlib>=3.2.0',
        'mendeleev>=0.5.2',
        'molmod>=1.4.5',
        'numpy>=1.18.1',
        'pandas>=1.0.1',
        'pathlib2>=2.3.5',
        'phonopy>=2.4.2',
        'psutil>=5.7.0',
        'pyfileindex>=0.0.4',
        'pyiron_base>=0.1.23',
        'pysqa>=0.0.11',
        'quickff>=2.2.4',
        'scipy>=1.4.1',
        'seekpath>=1.9.4',
        'six>=1.14.0',
        'scikit-learn>=0.22',
        'spglib>=1.14.1',
        'sqlalchemy>=1.3.14',
        'tables>=3.6.1',
        'tamkin>=1.2.6',
        'tqdm>=4.43.0',
        'yaff>=1.4.2'
    ],
    cmdclass=versioneer.get_cmdclass(),

    )

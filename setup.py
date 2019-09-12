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
                 'Programming Language :: Python :: 2.7',
                 'Programming Language :: Python :: 3',
                 'Programming Language :: Python :: 3.4',
                 'Programming Language :: Python :: 3.5',
                 'Programming Language :: Python :: 3.6',
                 'Programming Language :: Python :: 3.7'],

    keywords='pyiron',
    packages=find_packages(exclude=["*tests*", "*docs*", "*binder*", "*conda*", "*notebooks*", "*.ci_support*"]),
    install_requires=[
        'ase>=3.17',
        'defusedxml>=0.5.0',
        'dill>=0.3.0',
        'future>=0.17.1',
        'h5io>=0.1.1',
        'h5py>=2.10.0',
        'matplotlib>=2.2.4',
        'mendeleev>=0.5.1',
        'numpy>=1.16.4',
        'pandas>=0.24.2',
        'pathlib2>=2.3.4',
        'phonopy>=2.3.2',
        'psutil>=5.6.3',
        'pysqa>=0.0.4',
        'scipy>=1.2.1',
        'six>=1.12.0',
        'spglib>=1.14.1',
        'sqlalchemy>=1.3.8',
        'tables>=3.5.2',
        'tqdm>=4.35.0'
    ],
    cmdclass=versioneer.get_cmdclass(),
    )

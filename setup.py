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
                 'Programming Language :: Python :: 3.8',
                 'Programming Language :: Python :: 3.9'],

    keywords='pyiron',
    packages=find_packages(exclude=["*tests*", "*docs*", "*binder*", "*conda*", "*notebooks*", "*.ci_support*"]),
    install_requires=[
        'ase>=3.20.1',
        'defusedxml>=0.6.0',
        'future>=0.18.2',
        'h5py>=3.1.0',
        'matplotlib>=3.3.3',
        'mendeleev>=0.6.1',
        'molmod>=1.4.8',
        'numpy>=1.19.4',
        'pandas>=1.2.0',
        'phonopy>=2.8.1',
        'pyiron_base>=0.1.42',
        'quickff>=2.2.4',
        'scipy>=1.6.0',
        'seekpath>=1.9.5',
        'six>=1.15.0',
        'scikit-learn>=0.24.0',
        'spglib>=1.16.0',
        'tables>=3.6.1',
        'tamkin>=1.2.6',
        'yaff>=1.4.2'
    ],
    cmdclass=versioneer.get_cmdclass(),

    )

"""
Setuptools based setup module
"""
from setuptools import setup, find_packages


setup(
    name='pyiron_atomistics',
    version='0.0.1',
    description='pyiron IDE plugin for atomistic calculations',
    long_description='http://pyiron.org',

    url='https://github.com/jan-janssen/pyiron_atomistics',
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
    packages=find_packages(),
    install_requires=['ase',
                      'pyiron_base',
                      'spglib']
    )

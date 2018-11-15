from testing.randomatomistic import AtomisticExampleJob
from testing.randomatomistic import ExampleJob
__all__ = ['AtomisticExampleJob', 'ExampleJob']

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions

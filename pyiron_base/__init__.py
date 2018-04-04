__version__ = '0.1'
__all__ = ["hamilton",
           "objects",
           "hamutilities",
           "structure",
           "utilities",
           "gui",
           "masterjobs",
           "workbench.py",
           "workbench_gui.py",
           "setup.py",
           "update.py"]

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions

# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.


def fix_ipython_autocomplete(enable=True):
    """Change autocomplete behavior for IPython > 6.x
        
    Parameter
    ---------
    enable : bool (default True)
        Is use the trick.

    Notes
    -----
    Since IPython > 6.x the ``jedi`` package is using for autocomplete by default.
    But in some cases, the autocomplete doesn't work correctly wrong (see e.g.
    `here <https://github.com/ipython/ipython/issues/11653>`_).
    
    To set the correct behaviour we should use in IPython environment::

        %config Completer.use_jedi = False

    or add to IPython config (``<HOME>\.ipython\profile_default\ipython_config.py``)::

        c.Completer.use_jedi = False
    """

    try:
        __IPYTHON__
    except NameError:
        pass
    else:
        from IPython import __version__
        major = int(__version__.split('.')[0])
        if major >= 6:
            from IPython import get_ipython
            get_ipython().Completer.use_jedi = not enable

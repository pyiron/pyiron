@ECHO OFF

REM Command file for automatic scanning of the source code and generation of the rst files
REM add the following changes to conf.py
REM
REM sys.path.insert(0, os.path.abspath('../src'))
REM extensions:     'sphinx.ext.autodoc', 'sphinx.ext.viewcode',
REM

sphinx-apidoc -f -e -o ./apidoc ../PyIron
make html
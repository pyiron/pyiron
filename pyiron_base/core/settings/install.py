import os
from zipfile import ZipFile
from shutil import copytree, rmtree
import tempfile
import sys

if sys.version_info >= (3,):
    import urllib.request as urllib2
else:
    import urllib as urllib2


def _download_resources(zip_file="resources.zip",
                        resource_directory="~/pyiron/resources",
                        giturl_for_zip_file="https://github.com/pyiron/pyiron-resources/archive/master.zip",
                        git_folder_name="pyiron-resources-master"):
    user_directory = os.path.normpath(os.path.abspath(os.path.expanduser(resource_directory)))
    temp_directory = tempfile.gettempdir()
    temp_zip_file = os.path.join(temp_directory, zip_file)
    temp_extract_folder = os.path.join(temp_directory, git_folder_name)
    urllib2.urlretrieve(giturl_for_zip_file, temp_zip_file)
    if os.path.exists(user_directory):
        raise ValueError('The resource directory exists already, therefore it can not be created: ', user_directory)
    with ZipFile(temp_zip_file) as zip_file_object:
        zip_file_object.extractall(temp_directory)
    copytree(temp_extract_folder, user_directory)
    os.remove(temp_zip_file)
    rmtree(temp_extract_folder)


def _write_config_file(file_name='~/.pyiron'):
    config_file = os.path.normpath(os.path.abspath(os.path.expanduser(file_name)))
    if not os.path.isfile(config_file):
        with open(config_file, 'w') as cf:
            cf.writelines(['[DEFAULT]\n',
                           'PROJECT_PATHS = ~/pyiron/projects\n',
                           'RESOURCE_PATHS = ~/pyiron/resources\n'])
        project_path = os.path.normpath(os.path.abspath(os.path.expanduser('~/pyiron/projects')))
        if not os.path.exists(project_path):
            os.makedirs(project_path)


def install_pyiron(config_file_name='~/.pyiron',
                   zip_file="resources.zip",
                   resource_directory="~/pyiron/resources",
                   giturl_for_zip_file="https://github.com/pyiron/pyiron-resources/archive/master.zip",
                   git_folder_name="pyiron-resources-master"):
    _write_config_file(file_name=config_file_name)
    _download_resources(zip_file=zip_file,
                        resource_directory=resource_directory,
                        giturl_for_zip_file=giturl_for_zip_file,
                        git_folder_name=git_folder_name)

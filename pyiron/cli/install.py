# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.
"""
Install pyiron config and resources for the first time.
"""

from pyiron.base.settings.install import install_pyiron

__author__ = "Marvin Poul"
__copyright__ = (
    "Copyright 2020, Max-Planck-Institut für Eisenforschung GmbH - "
    "Computational Materials Design (CM) Department"
)
__version__ = "0.1"
__maintainer__ = "Marvin Poul"
__email__ = "poul@mpie.de"
__status__ = "development"
__date__ = "Jun 26, 2020"

def register(parser):
    parser.add_argument(
            "-c", "--config", default = "~/.pyiron",
            help = "path where config file should be written"
    )
    parser.add_argument(
            "-r", "--resources", default = "~/pyiron/resources",
            help = "path where resources should be installed"
    )
    parser.add_argument(
            "-u", "--url",
            default = "https://github.com/pyiron/pyiron-resources/archive/master.zip",
            help = "url to download zipped resources"
    )
    parser.add_argument(
            "-p", "--project", default = "~/pyiron/projects",
            help = "path where pyiron should expect projects to run"
    )

def main(args):
    install_pyiron(
        config_file_name=args.config,
        project_path=args.project,
        resource_directory=args.resources,
        giturl_for_zip_file=args.url
    )

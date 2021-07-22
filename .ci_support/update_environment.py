# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import json
import re
import sys

import yaml

environment_file = '.ci_support/environment.yml'
name_mapping_file = '.ci_support/pypi_vs_conda_names.json'


class EnvironmentUpdater:
    def __init__(self, package_name, from_version, to_version):
        """
        Updates the version of a package in the conda environment file.

        Parameters:
            package_name: Name of the package to update as available on PyPI
            from_version: Version the package is before the update
            to_version: Version to which the package should be updated
        """
        self.from_version = from_version
        self.to_version = to_version
        with open(name_mapping_file, 'r') as f:
            self._name_conversion_dict = json.load(f)

        with open(environment_file, 'r') as f:
            self.environment = yaml.safe_load(f)

        self.package_name = self._convert_package_name(package_name)

    def _convert_package_name(self, name):
        if name in self._name_conversion_dict.keys():
            result = self._name_conversion_dict[name]
        else:
            result = name
        return result

    def _update_dependencies(self):
        updated_dependencies = []

        for dep in self.environment['dependencies']:
            updated_dependencies.append(re.sub(
                r'(' + self.package_name + '.*)' + self.from_version,
                r'\g<1>' + self.to_version,
                dep
            ))

        self.environment['dependencies'] = updated_dependencies

    def _write(self):
        with open(environment_file, 'w') as f:
            yaml.safe_dump(self.environment, f)

    def update_dependencies(self):
        """Update the version of the requested dependency in the environment file"""
        self._update_dependencies()
        self._write()


if len(sys.argv) != 7 or not (sys.argv[1] == 'Bump' and sys.argv[3] == 'from' and sys.argv[5] == 'to'):
    raise ValueError(f"Title of a dependabot PR 'Bump <package> from <version> to <version>' expected, "
                     f"but got {' '.join(sys.argv[1:])}")
package_to_update = sys.argv[2]
from_version = sys.argv[4]
to_version = sys.argv[6]

updater = EnvironmentUpdater(package_to_update, from_version, to_version)
updater.update_dependencies()

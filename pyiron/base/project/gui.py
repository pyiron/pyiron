# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import numpy as np
import os
import ipywidgets as widgets
from IPython.display import display
import html
import pandas

__author__ = "Joerg Neugebauer, Jan Janssen"
__copyright__ = (
    "Copyright 2020, Max-Planck-Institut für Eisenforschung GmbH - "
    "Computational Materials Design (CM) Department"
)
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "development"
__date__ = "Sep 1, 2017"


def is_h5_file(name):
    return "h5" == name.split(".")[-1]


def is_h5_dir(name):
    return "hdf5" == name.split("_")[-1]


class ProjectGUI:
    """
    Gui to quickly browse through projects and objects
    Note: requires "%matplotlib notebook" at the beginning of your notebook (to refresh plots)
    """

    def __init__(self, project):
        self._ref_path = project.path
        self.project = project.copy()
        self.project._inspect_mode = True
        self.parent = None
        self.name = None
        self.fig, self.ax = None, None
        self.w_group = None
        self.w_node = None
        self.w_file = None
        self.w_text = None
        self.w_tab = None
        self.w_path = None
        self.w_type = None
        #         self.fig, self.ax = plt.subplots()

        self.create_widgets()
        self.connect_widgets()
        self.display()

    def create_widgets(self):
        """

        Returns:

        """
        select = widgets.Select  # Multiple
        self.w_group = select(description="Groups:")
        self.w_node = select(description="Nodes:")
        self.w_file = select(description="Files:")
        self.w_text = widgets.HTML()
        # self.w_text = widgets.Textarea()
        self.w_text.layout.height = "330px"
        self.w_text.layout.width = "580px"
        self.w_text.disabled = True
        #         self.w_plot = self.fig
        w_list = widgets.VBox([self.w_group, self.w_node, self.w_file])
        self.w_tab = widgets.HBox([w_list, self.w_text])
        #         tab = widgets.Tab(children=[self.w_group, self.w_node, self.w_file, self.w_text])
        #         [tab.set_title(num, name) for num, name in enumerate(['groups', 'nodes', 'files', 'text'])]
        #         self.w_tab = tab

        self.w_path = widgets.Text(name="Path: ")
        self.w_path.layout.width = "680px"
        self.w_type = widgets.Text(name="Type: ")

        self.refresh_view()

    def connect_widgets(self):
        """

        Returns:

        """
        self.w_group.observe(self.on_value_change, names="value")
        self.w_node.observe(self.on_value_change, names="value")
        self.w_file.observe(self.on_value_change, names="value")

    def display(self):
        """

        Returns:

        """
        w_txt = widgets.HBox([self.w_path, self.w_type])
        display(widgets.VBox([self.w_tab, w_txt]))

    def plot_array(self, val):
        """

        Args:
            val:

        Returns:

        """
        try:
            import pylab as plt
        except ImportError:
            import matplotlib.pyplot as plt

        plt.ioff()
        if self.fig is None:
            self.fig, self.ax = plt.subplots()
        else:
            self.ax.clear()

        # self.ax.set_title(self.name)
        if val.ndim == 1:
            self.ax.plot(val)
        elif val.ndim == 2:
            if len(val) == 1:
                self.ax.plot(val[0])
            else:
                self.ax.plot(val)
        elif val.ndim == 3:
            self.ax.plot(val[:, :, 0])
        # self.fig.canvas.draw()
        self.w_text.value = self.plot_to_html()
        plt.close()

    def plot_to_html(self):
        import base64
        import io

        # write image data to a string buffer and get the PNG image bytes
        buf = io.BytesIO()
        self.fig.set_size_inches(8, 4)
        self.fig.savefig(buf, format="png")
        buf.seek(0)
        return """<img src='data:image/png;base64,{}'/>""".format(
            base64.b64encode(buf.getvalue()).decode("ascii")
        )

    def on_value_change(self, change):
        """

        Args:
            change:

        Returns:

        """
        name = change["new"]
        self.w_text.value = ""
        if name == "..":
            self.move_up()
            self.w_group.value = "."
        else:
            if name is not None:
                if isinstance(name, str):
                    self.update_project(name)

    def get_rel_path(self, path):
        """

        Args:
            path:

        Returns:

        """
        return os.path.relpath(path, self._ref_path).replace("\\", "/")

    @staticmethod
    def dict_to_str(my_dict):
        """

        Args:
            my_dict:

        Returns:

        """
        # eol = os.linesep
        eol = "<br>"
        if "Parameter" in my_dict.keys():
            key = html.escape(my_dict["Parameter"])
            val = html.escape(my_dict["Value"])
            com = html.escape(my_dict["Comment"])
            table = [
                "{}: {} {} {}".format(key, val, com, eol)
                for key, val, com in zip(key, val, com)
            ]
        else:
            table = ["{}: {} {}".format(key, val, eol) for key, val in my_dict.items()]
        return "".join(table)

    def refresh_view(self):
        """

        Returns:

        """
        eol = os.linesep
        self.w_type.value = str(type(self.project))
        if isinstance(self.project, str):
            self.w_text.value = html.escape(self.project)
            self._move_up()
            # self.w_path.value = self.get_rel_path(self.parent.path + '/' + self.name)
        elif isinstance(self.project, dict):
            self.w_text.value = self.dict_to_str(self.project)
            self._move_up()
            # self.w_path.value = self.get_rel_path(self.parent.path + '/' + self.name)
        elif isinstance(self.project, (int, float)):
            self.w_text.value = html.escape(str(self.project))
            self._move_up()
            # self.w_path.value = self.get_rel_path(self.parent.path + '/' + self.name)
        elif isinstance(self.project, list):
            max_length = 2000  # performance of widget above is extremely poor
            if len(self.project) < max_length:
                self.w_text.value = "<br>".join(self.project)
            else:
                self.w_text.value = (
                    "<br>".join(self.project[:max_length])
                    + eol
                    + " .... file too long: skipped ...."
                )

            self.w_type.value = "list: {} lines".format(len(self.project))
            self._move_up()
            # self.w_path.value = self.get_rel_path(self.parent.path + '/' + self.name)
        elif isinstance(self.project, np.ndarray):
            self.plot_array(self.project)
            self._move_up()
            # self.w_path.value = self.get_rel_path(self.parent.path + '/' + self.name)
        elif "data_dict" in self.project.list_nodes():
            self.w_text.value = pandas.DataFrame(self.project["data_dict"]).to_html()
            self._move_up()
            # self.w_path.value = self.get_rel_path(self.parent.path + '/' + self.name)
        elif self.project is None:
            raise ValueError(
                "project is None: {}".format(type(self.project)), self.parent
            )
        else:
            self.w_group.options = [".", ".."]
            self.w_node.options = ["."]
            self.w_file.options = ["."]

            if hasattr(self.project, "path"):
                self.w_path.value = self.get_rel_path(self.project.path)
            else:
                print("path missing: ", type(self.project))

            groups = sorted(
                [el for el in self.project.list_groups() if not is_h5_dir(el)]
            )
            self.w_group.options = list(self.w_group.options) + groups

            nodes = sorted([el for el in self.project.list_nodes()])
            self.w_node.options = list(self.w_node.options) + nodes

            if hasattr(self.project, "list_files"):
                files = sorted(
                    [el for el in self.project.list_files() if not is_h5_file(el)]
                )
                self.w_file.options = list(self.w_file.options) + files
            else:
                self.w_file.options = []

    def update_project(self, name):
        """

        Args:
            name:

        Returns:

        """
        if name == ".":
            return
        self.name = name
        self.parent = self.project.copy()
        self.project = self.project[name]
        self.refresh_view()

    def _move_up(self):
        """

        Returns:

        """
        if hasattr(self.project, "path"):
            self.project = self.project[".."]
        else:
            self.project = self.parent

        if self.parent is None:
            self.w_path.value = "/".join(self.w_path.value.split("/")[:-1])
        else:
            self.w_path.value = self.get_rel_path(self.parent.path + "/" + self.name)

    def move_up(self):
        """

        Returns:

        """
        self._move_up()
        self.refresh_view()

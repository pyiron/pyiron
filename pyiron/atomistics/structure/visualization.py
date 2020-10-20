# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from ase.atoms import Atoms as ASEAtoms, get_distances as ase_get_distances, Atom as ASEAtom
import numpy as np
import warnings
from matplotlib.colors import rgb2hex
from pyiron_base import Settings
from scipy.interpolate import interp1d

__author__ = "Joerg Neugebauer, Sudarsan Surendralal"
__copyright__ = (
    "Copyright 2020, Max-Planck-Institut für Eisenforschung GmbH - "
    "Computational Materials Design (CM) Department"
)
__version__ = "1.0"
__maintainer__ = "Sudarsan Surendralal"
__email__ = "surendralal@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2017"

s = Settings()


class Visualization(object):
    """
    Visualization class. Available functions are:

        plot3d: NGLView based structure plot tool
        plot3d_ase: ASE version of plot3d
        plot_plotly: plotly based structure plot tool

    Read DocStrings of individual functions for more info

    """

    def __init__(self, structure):
        self.ref_structure = structure

    @staticmethod
    def _ngl_write_cell(a1, a2, a3, f1=90, f2=90, f3=90):
        """
        Writes a PDB-formatted line to represent the simulation cell.

        Args:
            a1, a2, a3 (float): Lengths of the cell vectors.
            f1, f2, f3 (float): Angles between the cell vectors (which angles exactly?) (in degrees).

        Returns:
            (str): The line defining the cell in PDB format.
        """
        return "CRYST1 {:8.3f} {:8.3f} {:8.3f} {:6.2f} {:6.2f} {:6.2f} P 1\n".format(
            a1, a2, a3, f1, f2, f3
        )

    @staticmethod
    def _ngl_write_atom(
        num,
        species,
        x,
        y,
        z,
        group=None,
        num2=None,
        occupancy=1.0,
        temperature_factor=0.0,
    ):
        """
        Writes a PDB-formatted line to represent an atom.

        Args:
            num (int): Atomic index.
            species (str): Elemental species.
            x, y, z (float): Cartesian coordinates of the atom.
            group (str): A...group name? (Default is None, repeat elemental species.)
            num2 (int): An "alternate" index. (Don't ask me...) (Default is None, repeat first number.)
            occupancy (float): PDB occupancy parameter. (Default is 1.)
            temperature_factor (float): PDB temperature factor parameter. (Default is 0.

        Returns:
            (str): The line defining an atom in PDB format

        Warnings:
            * The [PDB docs](https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html) indicate that
                the xyz coordinates might need to be in some sort of orthogonal basis. If you have weird behaviour,
                this might be a good place to investigate.
        """
        if group is None:
            group = species
        if num2 is None:
            num2 = num
        return "ATOM {:>6} {:>4} {:>4} {:>5} {:10.3f} {:7.3f} {:7.3f} {:5.2f} {:5.2f} {:>11} \n".format(
            num, species, group, num2, x, y, z, occupancy, temperature_factor, species
        )

    def _ngl_write_structure(self, elements, positions, cell):
        """
        Turns structure information into a NGLView-readable protein-database-formatted string.

        Args:
            elements (numpy.ndarray/list): Element symbol for each atom.
            positions (numpy.ndarray/list): Vector of Cartesian atom positions.
            cell (numpy.ndarray/list): Simulation cell Bravais matrix.

        Returns:
            (str): The PDB-formatted representation of the structure.
        """
        from ase.geometry import cell_to_cellpar, cellpar_to_cell
        if cell is None or any(np.max(cell, axis=0) < 1e-2):
            # Define a dummy cell if it doesn't exist (eg. for clusters)
            max_pos = np.max(positions, axis=0)
            max_pos[np.abs(max_pos) < 1e-2] = 10
            cell = np.eye(3) * max_pos
        cellpar = cell_to_cellpar(cell)
        exportedcell = cellpar_to_cell(cellpar)
        rotation = np.linalg.solve(cell, exportedcell)

        pdb_str = self._ngl_write_cell(*cellpar)
        pdb_str += "MODEL     1\n"

        if rotation is not None:
            positions = np.array(positions).dot(rotation)

        for i, p in enumerate(positions):
            pdb_str += self._ngl_write_atom(i, elements[i], *p)

        pdb_str += "ENDMDL \n"
        return pdb_str

    def _atomic_number_to_radius(self, atomic_number, shift=0.2, slope=0.1, scale=1.0):
        """
        Give the atomic radius for plotting, which scales like the root of the atomic number.

        Args:
            atomic_number (int/float): The atomic number.
            shift (float): A constant addition to the radius. (Default is 0.2.)
            slope (float): A multiplier for the root of the atomic number. (Default is 0.1)
            scale (float): How much to rescale the whole thing by.

        Returns:
            (float): The radius. (Not physical, just for visualization!)
        """
        return (shift + slope * np.sqrt(atomic_number)) * scale

    def _add_colorscheme_spacefill(
        self, view, elements, atomic_numbers, particle_size, scheme="element"
    ):
        """
        Set NGLView spacefill parameters according to a color-scheme.

        Args:
            view (NGLWidget): The widget to work on.
            elements (numpy.ndarray/list): Elemental symbols.
            atomic_numbers (numpy.ndarray/list): Integer atomic numbers for determining atomic size.
            particle_size (float): A scale factor for the atomic size.
            scheme (str): The scheme to use. (Default is "element".)

            Possible NGLView color schemes:
              " ", "picking", "random", "uniform", "atomindex", "residueindex",
              "chainindex", "modelindex", "sstruc", "element", "resname", "bfactor",
              "hydrophobicity", "value", "volume", "occupancy"

        Returns:
            (nglview.NGLWidget): The modified widget.
        """
        for elem, num in set(list(zip(elements, atomic_numbers))):
            view.add_spacefill(
                selection="#" + elem,
                radius_type="vdw",
                radius=self._atomic_number_to_radius(num, scale=particle_size),
                color_scheme=scheme,
            )
        return view

    def _add_custom_color_spacefill(self, view, atomic_numbers, particle_size, colors):
        """
        Set NGLView spacefill parameters according to per-atom colors.

        Args:
            view (NGLWidget): The widget to work on.
            atomic_numbers (numpy.ndarray/list): Integer atomic numbers for determining atomic size.
            particle_size (float): A scale factor for the atomic size.
            colors (numpy.ndarray/list): A per-atom list of HTML or hex color codes.

        Returns:
            (nglview.NGLWidget): The modified widget.
        """
        for n, num in enumerate(atomic_numbers):
            view.add_spacefill(
                selection=[n],
                radius_type="vdw",
                radius=self._atomic_number_to_radius(num, scale=particle_size),
                color=colors[n],
            )
        return view

    @staticmethod
    def _scalars_to_hex_colors(scalar_field, start=None, end=None, cmap=None):
        """
        Convert scalar values to hex codes using a colormap.

        Args:
            scalar_field (numpy.ndarray/list): Scalars to convert.
            start (float): Scalar value to map to the bottom of the colormap (values below are clipped). (Default is
                None, use the minimal scalar value.)
            end (float): Scalar value to map to the top of the colormap (values above are clipped).  (Default is
                None, use the maximal scalar value.)
            cmap (matplotlib.cm): The colormap to use. (Default is None, which gives a blue-red divergent map.)

        Returns:
            (list): The corresponding hex codes for each scalar value passed in.
        """
        if start is None:
            start = np.amin(scalar_field)
        if end is None:
            end = np.amax(scalar_field)
        interp = interp1d([start, end], [0, 1])
        remapped_field = interp(
            np.clip(scalar_field, start, end)
        )  # Map field onto [0,1]

        if cmap is None:
            try:
                from seaborn import diverging_palette
            except ImportError:
                print(
                    "The package seaborn needs to be installed for the plot3d() function!"
                )
            cmap = diverging_palette(245, 15, as_cmap=True)  # A nice blue-red palette

        return [
            rgb2hex(cmap(scalar)[:3]) for scalar in remapped_field
        ]  # The slice gets RGB but leaves alpha

    @staticmethod
    def _get_orientation(view_plane):
        """
        A helper method to plot3d, which generates a rotation matrix from the input `view_plane`, and returns a
        flattened list of len = 16. This flattened list becomes the input argument to `view.contol.orient`.

        Args:
            view_plane (numpy.ndarray/list): A Nx3-array/list (N = 1,2,3); the first 3d-component of the array
                specifies which plane of the system to view (for example, [1, 0, 0], [1, 1, 0] or the [1, 1, 1] planes),
                the second 3d-component (if specified, otherwise [1, 0, 0]) gives the horizontal direction, and the
                third component (if specified) is the vertical component, which is ignored and calculated internally.
                The orthonormality of the orientation is internally ensured, and therefore is not required in the
                function call.

        Returns:
            (list): orientation tensor
        """
        if len(np.array(view_plane).flatten()) % 3 != 0:
            raise ValueError("The shape of view plane should be (N, 3), where N = 1, 2 or 3. Refer docs for more info.")
        view_plane = np.array(view_plane).reshape(-1, 3)
        rotation_matrix = np.roll(np.eye(3), -1, axis=0)
        rotation_matrix[:len(view_plane)] = view_plane
        rotation_matrix /= np.linalg.norm(rotation_matrix, axis=-1)[:, np.newaxis]
        rotation_matrix[1] -= np.dot(rotation_matrix[0], rotation_matrix[1]) * rotation_matrix[0]  # Gran-Schmidt
        rotation_matrix[2] = np.cross(rotation_matrix[0], rotation_matrix[1])  # Specify third axis
        if np.isclose(np.linalg.det(rotation_matrix), 0):
            return np.eye(3)  # view_plane = [0,0,1] is the default view of NGLview, so we do not modify it
        return np.roll(rotation_matrix / np.linalg.norm(rotation_matrix, axis=-1)[:, np.newaxis], 2, axis=0).T

    def _get_flattened_orientation(self, view_plane, distance_from_camera):
        """
        A helper method to plot3d, which generates a rotation matrix from the input `view_plane`, and returns a
        flattened list of len = 16. This flattened list becomes the input argument to `view.contol.orient`.

        Args:
            view_plane (numpy.ndarray/list): A Nx3-array/list (N = 1,2,3); the first 3d-component of the array
                specifies which plane of the system to view (for example, [1, 0, 0], [1, 1, 0] or the [1, 1, 1] planes),
                the second 3d-component (if specified, otherwise [1, 0, 0]) gives the horizontal direction, and the
                third component (if specified) is the vertical component, which is ignored and calculated internally.
                The orthonormality of the orientation is internally ensured, and therefore is not required in the
                function call.
            distance_from_camera (float): Distance of the camera from the structure. Higher = farther away.

        Returns:
            (list): Flattened list of len = 16, which is the input argument to `view.contol.orient`
        """
        if distance_from_camera <= 0:
            raise ValueError("´distance_from_camera´ must be a positive float!")
        flattened_orientation = np.eye(4)
        flattened_orientation[:3, :3] = self._get_orientation(view_plane)

        return (distance_from_camera * flattened_orientation).ravel().tolist()

    def plot3d_plotly(
        self,
        scalar_field=None,
        particle_size=1.0,
        camera="orthographic",
        view_plane=np.array([1, 1, 1]),
        distance_from_camera=1.25,
        opacity=1,
    ):
        """
        Make a 3D plot of the atomic structure.

        Args:
            camera (str): 'perspective' or 'orthographic'. (Default is 'perspective'.)
            particle_size (float): Size of the particles. (Default is 1.)
            scalar_field (numpy.ndarray): Color each atom according to the array value (Default is None, use coloring
                scheme.)
            view_plane (numpy.ndarray): A Nx3-array (N = 1,2,3); the first 3d-component of the array specifies
                which plane of the system to view (for example, [1, 0, 0], [1, 1, 0] or the [1, 1, 1] planes), the
                second 3d-component (if specified, otherwise [1, 0, 0]) gives the horizontal direction, and the third
                component (if specified) is the vertical component, which is ignored and calculated internally. The
                orthonormality of the orientation is internally ensured, and therefore is not required in the function
                call. (Default is np.array([0, 0, 1]), which is view normal to the x-y plane.)
            distance_from_camera (float): Distance of the camera from the structure. Higher = farther away.
                (Default is 14, which also seems to be the NGLView default value.)
            opacity (float): opacity

        Returns:
            (plotly.express): The NGLView widget itself, which can be operated on further or viewed as-is.

        """
        try:
            import plotly.express as px
        except ImportError:
            print("plotly not installed")
            return
        parent_basis = self.ref_structure.get_parent_basis()
        elements = parent_basis.get_chemical_symbols()
        atomic_numbers = parent_basis.get_atomic_numbers()
        if scalar_field is None:
            scalar_field = elements
        fig = px.scatter_3d(x=self.ref_structure.positions[:,0],
                            y=self.ref_structure.positions[:,1],
                            z=self.ref_structure.positions[:,2],
                            color=scalar_field,
                            opacity=opacity,
                            size=self._atomic_number_to_radius(atomic_numbers, scale=particle_size/(0.1*self.ref_structure.get_volume()**(1/3))))
        fig.layout.scene.camera.projection.type = camera
        rot = self._get_orientation(view_plane).T
        rot[0,:] *= distance_from_camera
        angle = dict(
            up=dict(x=rot[2,0], y=rot[2,1], z=rot[2,2]),
            eye=dict(x=rot[0,0], y=rot[0,1], z=rot[0,2])
        )
        fig.update_layout(scene_camera=angle)
        fig.update_traces(marker=dict(line=dict(width=0.1, color='DarkSlateGrey')))
        return fig

    def plot3d(
        self,
        show_cell=True,
        show_axes=True,
        camera="orthographic",
        spacefill=True,
        particle_size=1.0,
        select_atoms=None,
        background="white",
        color_scheme=None,
        colors=None,
        scalar_field=None,
        scalar_start=None,
        scalar_end=None,
        scalar_cmap=None,
        vector_field=None,
        vector_color=None,
        magnetic_moments=False,
        custom_array=None,
        custom_3darray=None,
        view_plane=np.array([0, 0, 1]),
        distance_from_camera=14.0
    ):
        """
        Plot3d relies on NGLView to visualize atomic structures. Here, we construct a string in the "protein database"
        ("pdb") format, then turn it into an NGLView "structure". PDB is a white-space sensitive format, so the
        string snippets are carefully formatted.

        The final widget is returned. If it is assigned to a variable, the visualization is suppressed until that
        variable is evaluated, and in the meantime more NGL operations can be applied to it to modify the visualization.

        Args:
            show_cell (bool): Whether or not to show the frame. (Default is True.)
            show_axes (bool): Whether or not to show xyz axes. (Default is True.)
            camera (str): 'perspective' or 'orthographic'. (Default is 'perspective'.)
            spacefill (bool): Whether to use a space-filling or ball-and-stick representation. (Default is True, use
                space-filling atoms.)
            particle_size (float): Size of the particles. (Default is 1.)
            select_atoms (numpy.ndarray): Indices of atoms to show, either as integers or a boolean array mask.
                (Default is None, show all atoms.)
            background (str): Background color. (Default is 'white'.)
            color_scheme (str): NGLView color scheme to use. (Default is None, color by element.)
            colors (numpy.ndarray): A per-atom array of HTML color names or hex color codes to use for atomic colors.
                (Default is None, use coloring scheme.)
            scalar_field (numpy.ndarray): Color each atom according to the array value (Default is None, use coloring
                scheme.)
            scalar_start (float): The scalar value to be mapped onto the low end of the color map (lower values are
                clipped). (Default is None, use the minimum value in `scalar_field`.)
            scalar_end (float): The scalar value to be mapped onto the high end of the color map (higher values are
                clipped). (Default is None, use the maximum value in `scalar_field`.)
            scalar_cmap (matplotlib.cm): The colormap to use. (Default is None, giving a blue-red divergent map.)
            vector_field (numpy.ndarray): Add vectors (3 values) originating at each atom. (Default is None, no
                vectors.)
            vector_color (numpy.ndarray): Colors for the vectors (only available with vector_field). (Default is None,
                vectors are colored by their direction.)
            magnetic_moments (bool): Plot magnetic moments as 'scalar_field' or 'vector_field'.
            view_plane (numpy.ndarray): A Nx3-array (N = 1,2,3); the first 3d-component of the array specifies
                which plane of the system to view (for example, [1, 0, 0], [1, 1, 0] or the [1, 1, 1] planes), the
                second 3d-component (if specified, otherwise [1, 0, 0]) gives the horizontal direction, and the third
                component (if specified) is the vertical component, which is ignored and calculated internally. The
                orthonormality of the orientation is internally ensured, and therefore is not required in the function
                call. (Default is np.array([0, 0, 1]), which is view normal to the x-y plane.)
            distance_from_camera (float): Distance of the camera from the structure. Higher = farther away.
                (Default is 14, which also seems to be the NGLView default value.)

            Possible NGLView color schemes:
              " ", "picking", "random", "uniform", "atomindex", "residueindex",
              "chainindex", "modelindex", "sstruc", "element", "resname", "bfactor",
              "hydrophobicity", "value", "volume", "occupancy"

        Returns:
            (nglview.NGLWidget): The NGLView widget itself, which can be operated on further or viewed as-is.

        Warnings:
            * Many features only work with space-filling atoms (e.g. coloring by a scalar field).
            * The colour interpretation of some hex codes is weird, e.g. 'green'.
        """
        try:  # If the graphical packages are not available, the GUI will not work.
            import nglview
        except ImportError:
            raise ImportError(
                "The package nglview needs to be installed for the plot3d() function!"
            )

        if custom_array is not None:
            warnings.warn(
                "custom_array is deprecated. Use scalar_field instead",
                DeprecationWarning,
            )
            scalar_field = custom_array

        if custom_3darray is not None:
            warnings.warn(
                "custom_3darray is deprecated. Use vector_field instead",
                DeprecationWarning,
            )
            vector_field = custom_3darray

        if magnetic_moments is True and hasattr(self.ref_structure, 'spin'):
            if len(self.ref_structure.get_initial_magnetic_moments().shape) == 1:
                scalar_field = self.ref_structure.get_initial_magnetic_moments()
            else:
                vector_field = self.ref_structure.get_initial_magnetic_moments()

        parent_basis = self.ref_structure.get_parent_basis()
        elements = parent_basis.get_chemical_symbols()
        atomic_numbers = parent_basis.get_atomic_numbers()
        positions = self.ref_structure.positions

        # If `select_atoms` was given, visualize only a subset of the `parent_basis`
        if select_atoms is not None:
            select_atoms = np.array(select_atoms, dtype=int)
            elements = elements[select_atoms]
            atomic_numbers = atomic_numbers[select_atoms]
            positions = positions[select_atoms]
            if colors is not None:
                colors = np.array(colors)
                colors = colors[select_atoms]
            if scalar_field is not None:
                scalar_field = np.array(scalar_field)
                scalar_field = scalar_field[select_atoms]
            if vector_field is not None:
                vector_field = np.array(vector_field)
                vector_field = vector_field[select_atoms]
            if vector_color is not None:
                vector_color = np.array(vector_color)
                vector_color = vector_color[select_atoms]

        # Write the nglview protein-database-formatted string
        struct = nglview.TextStructure(
            self._ngl_write_structure(elements, positions, self.ref_structure.cell)
        )

        # Parse the string into the displayable widget
        view = nglview.NGLWidget(struct)

        if spacefill:
            # Color by scheme
            if color_scheme is not None:
                if colors is not None:
                    warnings.warn("`color_scheme` is overriding `colors`")
                if scalar_field is not None:
                    warnings.warn("`color_scheme` is overriding `scalar_field`")
                view = self._add_colorscheme_spacefill(
                    view, elements, atomic_numbers, particle_size, color_scheme
                )
            # Color by per-atom colors
            elif colors is not None:
                if scalar_field is not None:
                    warnings.warn("`colors` is overriding `scalar_field`")
                view = self._add_custom_color_spacefill(
                    view, atomic_numbers, particle_size, colors
                )
            # Color by per-atom scalars
            elif scalar_field is not None:  # Color by per-atom scalars
                colors = self._scalars_to_hex_colors(
                    scalar_field, scalar_start, scalar_end, scalar_cmap
                )
                view = self._add_custom_color_spacefill(
                    view, atomic_numbers, particle_size, colors
                )
            # Color by element
            else:
                view = self._add_colorscheme_spacefill(
                    view, elements, atomic_numbers, particle_size
                )
            view.remove_ball_and_stick()
        else:
            view.add_ball_and_stick()

        if show_cell:
            if parent_basis.cell is not None:
                if all(np.max(parent_basis.cell, axis=0) > 1e-2):
                    view.add_unitcell()

        if vector_color is None and vector_field is not None:
            vector_color = (
                0.5
                * np.array(vector_field)
                / np.linalg.norm(vector_field, axis=-1)[:, np.newaxis]
                + 0.5
            )
        elif (
            vector_field is not None and vector_field is not None
        ):  # WARNING: There must be a bug here...
            try:
                if vector_color.shape != np.ones((len(self.ref_structure), 3)).shape:
                    vector_color = np.outer(
                        np.ones(len(self.ref_structure)), vector_color / np.linalg.norm(vector_color)
                    )
            except AttributeError:
                vector_color = np.ones((len(self.ref_structure), 3)) * vector_color

        if vector_field is not None:
            for arr, pos, col in zip(vector_field, positions, vector_color):
                view.shape.add_arrow(list(pos), list(pos + arr), list(col), 0.2)

        if show_axes:  # Add axes
            axes_origin = -np.ones(3)
            arrow_radius = 0.1
            text_size = 1
            text_color = [0, 0, 0]
            arrow_names = ["x", "y", "z"]

            for n in [0, 1, 2]:
                start = list(axes_origin)
                shift = np.zeros(3)
                shift[n] = 1
                end = list(start + shift)
                color = list(shift)
                # We cast as list to avoid JSON warnings
                view.shape.add_arrow(start, end, color, arrow_radius)
                view.shape.add_text(end, text_color, text_size, arrow_names[n])

        if camera != "perspective" and camera != "orthographic":
            warnings.warn(
                "Only perspective or orthographic is (likely to be) permitted for camera"
            )

        view.camera = camera
        view.background = background

        orientation = self._get_flattened_orientation(view_plane=view_plane,
                                                      distance_from_camera=distance_from_camera)
        view.control.orient(orientation)

        return view

    def plot3d_ase(
        self,
        spacefill=True,
        show_cell=True,
        camera="perspective",
        particle_size=0.5,
        background="white",
        color_scheme="element",
        show_axes=True,
    ):
        """
        Possible color schemes:
          " ", "picking", "random", "uniform", "atomindex", "residueindex",
          "chainindex", "modelindex", "sstruc", "element", "resname", "bfactor",
          "hydrophobicity", "value", "volume", "occupancy"
        Returns:
        """
        try:  # If the graphical packages are not available, the GUI will not work.
            import nglview
        except ImportError:
            raise ImportError(
                "The package nglview needs to be installed for the plot3d() function!"
            )
        # Always visualize the parent basis
        parent_basis = self.ref_structure.get_parent_basis()
        view = nglview.show_ase(parent_basis)
        if spacefill:
            view.add_spacefill(
                radius_type="vdw", color_scheme=color_scheme, radius=particle_size
            )
            # view.add_spacefill(radius=1.0)
            view.remove_ball_and_stick()
        else:
            view.add_ball_and_stick()
        if show_cell:
            if parent_basis.cell is not None:
                if all(np.max(parent_basis.cell, axis=0) > 1e-2):
                    view.add_unitcell()
        if show_axes:
            view.shape.add_arrow([-2, -2, -2], [2, -2, -2], [1, 0, 0], 0.5)
            view.shape.add_arrow([-2, -2, -2], [-2, 2, -2], [0, 1, 0], 0.5)
            view.shape.add_arrow([-2, -2, -2], [-2, -2, 2], [0, 0, 1], 0.5)
        if camera != "perspective" and camera != "orthographic":
            print("Only perspective or orthographic is permitted")
            return None
        view.camera = camera
        view.background = background
        return view



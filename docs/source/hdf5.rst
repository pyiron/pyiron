===============================
HDF5 Serialization Architecture
===============================

---------
Structure
---------

Each hierachical object lives under its own group in the hdf, i.e. objects that
are attributes of another must have their own sub-group in that larger objects
group.  In its group each object must store
    - 'TYPE' equal to `str(type(self))`
      this provides the module path and class name from which pyiron will load
      a class
    - 'NAME' equal to `type(self).__name__`
      the unqualified class name, informational only
They may also store
    - 'HDF_VERSION' equal to a version string with format MAJOR.MINOR.PATCH
      the version of the structure of the type *in HDF5*; all classes
      must be able to read from HDF5 with at least the same MAJOR release, but
      explicit breaking behaviour should be very rare
    - 'VERSION' equal to a version string with format MAJOR.MINOR.PATCH
      the version of the functionality of the class; higher version must not
      change the HDF5 structure unless they also change HDF_VERSION

For example a class defined like this
.. code-block::

    class Foo:
        def __init__(self, parameter):
            self.bar = Bar()
            self.baz = Baz()
            self.parameter = parameter

should be serialized as
.. code-block::
    foo/
    foo/TYPE
    foo/NAME
    foo/VERSION
    foo/HDF_VERSION
    foo/parameter
    foo/bar/
    foo/bar/TYPE
    foo/bar/NAME
    foo/bar/VERSION
    foo/bar/HDF_VERSION
    foo/baz/
    foo/baz/TYPE
    foo/baz/NAME
    foo/baz/VERSION
    foo/baz/HDF_VERSION

---------------
Writing to HDF5
---------------

Each type must define a `to_hdf(self, hdf, group_name = None)` method that
takes the given `hdf` object, creates a subgroup called `group_name` in it (if
given) and then serializes itself to this group.  Some objects may keep a
default `ProjectHDFio` object during their lifetime (e.g. jobs), in this case
`hdf` maybe an optional parameter.

-----------------
Reading from HDF5
-----------------

Each type must define a `from_hdf(self, hdf, group_name = None)` method and may
define a `from_hdf_args(cls, hdf)`.  `from_hdf()` restores the state of the
already initialized object from the information stored in the HDF5 file.
`from_hdf_args()` reads the required parameters to instantiate the object from
HDF5 and returns them in a `dict`.

To read an object from a given `ProjectHDFio` path, call the `to_object()`
method.  This will first call `import_class` to read the class object, then
`make_from_hdf()` to instantiate it, if the class defines `from_hd_args()` it
will be called to supply the correct init parameters.  `to_object()` can also
be supplied with additional paraters to overrride the ones written to HDF5, in
particular it will always provide `job_name` and `project`.  However only those
parameters that are needed (i.e. declared by that classes' `__init__()`) will
be passed.

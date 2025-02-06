.. _install:

Installing GRGL
===============

The simplest way to install GRGL is via pip:

::

	pip install pygrgl

The resulting installation has two parts:

  1. The ``grg`` command, used for creating or summarizing GRGs. See ``grg --help``.
  2. The Python API, accessed via ``import pygrgl``. The API can be used to create
     or access GRGs - see the `Python API reference <python_api.html>`_ for more
     details.

For usage as a C++ library you'll need to clone (``git clone --recursive``)
the `GRGL repository <https://github.com/aprilweilab/grgl>`_. We use `CMake <https://cmake.org/>`_
as our build system, so the easiest way to incorporate GRGL into your project
is via CMake's `add_subdirectory <https://cmake.org/cmake/help/latest/command/add_subdirectory.html>`_ method.

Advanced: Installing from source
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You will want to install from source if you need support for BGEN file input, but otherwise
the simpler pip-based installation is likely sufficient. It is generally best
practice to use a `virtual environment <https://docs.python.org/3/library/venv.html>`_
when installing Python packages.

Installing from source:

::

  python setup.py bdist_wheel               # Compiles C++, builds a wheel in the dist/ directory
  pip install --force-reinstall dist/*.whl  # Install from wheel

Installing from source *and enabling BGEN support*:

::

  GRGL_BGEN=1 python setup.py bdist_wheel   # Compiles C++, builds a wheel in the dist/ directory
  pip install --force-reinstall dist/*.whl  # Install from wheel

There are other environment variables that control the behavior of the GRGL build as well:

- ``GRGL_GSL=1``: Enable GNU scientific library for computing p-values with GWAS.
- ``GRGL_DEBUG=1``: Build the C++ code in debug mode.
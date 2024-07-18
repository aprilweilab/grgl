from pprint import pprint
from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
import os
import shutil
import subprocess
import sys

C_MODULE_NAME = "_grgl"
ARG_DEBUG = "--debug-build"
ARG_BGEN = "--bgen"
ARG_COPYBINS = "--copy-bins"
ARG_GSL = "--gsl"

THISDIR = os.path.realpath(os.path.dirname(__file__))

copy_bins = False # Copy executables to the top-level directory?
extra_cmake_args = []
build_type = "Release"
for arg in sys.argv[1:]:
    if arg == ARG_DEBUG:
        build_type = "Debug"
        sys.argv.remove(ARG_DEBUG)
    elif arg == ARG_BGEN:
        extra_cmake_args.append("-DENABLE_BGEN=ON")
        sys.argv.remove(ARG_BGEN)
    elif arg == ARG_COPYBINS:
        copy_bins = True
        sys.argv.remove(ARG_COPYBINS)
    elif arg == ARG_GSL:
        extra_cmake_args.append("-DENABLE_GSL=ON")
        sys.argv.remove(ARG_GSL)

class CMakeExtension(Extension):
    def __init__(self, name, cmake_lists_dir=".", sources=[], extra_executables=[], **kwa):
        Extension.__init__(self, name, sources=sources, **kwa)
        self.cmake_lists_dir = os.path.abspath(cmake_lists_dir)
        self.extra_executables = extra_executables

class CMakeBuild(build_ext):
    def build_extensions(self):
        try:
            out = subprocess.check_output(["cmake", "--version"])
        except OSError:
            raise RuntimeError("Cannot find CMake executable")

        for ext in self.extensions:
            extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
            cmake_args = [
                "-DPYTHON_SUPPORT=ON",
                "-DCMAKE_BUILD_TYPE=%s" % build_type,
                "-DENABLE_CHECKERS=OFF",
                "-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}".format(build_type.upper(), extdir),
                "-DCMAKE_ARCHIVE_OUTPUT_DIRECTORY_{}={}".format(build_type.upper(), self.build_temp),
                "-DPYTHON_EXECUTABLE={}".format(sys.executable),
            ] + extra_cmake_args
            pprint(cmake_args)

            if not os.path.exists(self.build_temp):
                os.makedirs(self.build_temp)

            # Config and build the extension
            subprocess.check_call(["cmake", ext.cmake_lists_dir] + cmake_args,
                                  cwd=self.build_temp, stdout=sys.stdout)
            subprocess.check_call(["cmake", "--build", ".", "--config", build_type, "--", "-j"],
                                  cwd=self.build_temp, stdout=sys.stdout)

            for executable in ext.extra_executables:
                shutil.copy2(os.path.join(self.build_temp, executable),
                             os.path.join(extdir, executable))
                if copy_bins:
                    shutil.copy2(os.path.join(self.build_temp, executable),
                                 os.path.join(THISDIR, executable))

PACKAGE_NAME = "pygrgl"
EXECUTABLES = ["grgl", "grg-merge", "grgp", "gconverter", "gindexer"]

with open(os.path.join(THISDIR, "include", "grgl", "version.h")) as vf:
    for line in vf:
        line = line.strip()
        if line.startswith("#define GRGL_MAJOR_VERSION"):
            major_version = int(line.split(" ")[-1])
        if line.startswith("#define GRGL_MINOR_VERSION"):
            minor_version = int(line.split(" ")[-1])
version = f"{major_version}.{minor_version}"

setup(name=PACKAGE_NAME,
      packages=find_packages(),
      version=version,
      description="Genotype Representation Graph Library",
      author="Drew DeHaas, April Wei",
      author_email="",
      url="https://aprilweilab.github.io/",
      ext_modules=[CMakeExtension(C_MODULE_NAME, extra_executables=EXECUTABLES)],
      cmdclass={"build_ext": CMakeBuild},
      zip_safe=False,
      classifiers=[
          "Programming Language :: Python :: 3",
          "Operating System :: MacOS",
          "Operating System :: POSIX :: Linux",
      ],
      entry_points = {
        "console_scripts": ["grg=pygrgl.cli:main"],
      }
)

import subprocess, os
from cgitb import html
import os
import sys
sys.path.insert(0, os.path.abspath('../'))

# Useful for RTD Doxygen integration: https://breathe.readthedocs.io/en/latest/readthedocs.html
DOC_OUTPUT = "doc"
orig_dir = os.getcwd()
try:
    if not os.path.exists(DOC_OUTPUT):
        os.mkdir(DOC_OUTPUT)
    use_dir = os.environ.get("DOC_BUILD_DIR")
    if use_dir is not None:
        os.chdir(use_dir)
    subprocess.check_call("doxygen ../Doxyfile.in", shell=True)
finally:
    os.chdir(orig_dir)


extensions = ["breathe", "sphinx.ext.autodoc", "sphinx.ext.autosummary"]

html_theme = "sphinx_rtd_theme"

# Breathe configuration
breathe_projects = {
    "grgl": DOC_OUTPUT + "/xml/",
}
breathe_default_project = "grgl"


autosummary_generate = True  # Turn on sphinx.ext.autosummary

project = "GRGL"

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

extensions = ["breathe", "sphinx.ext.autodoc", "sphinx.ext.autosummary", "IPython.sphinxext.ipython_console_highlighting"]

html_theme = "pydata_sphinx_theme"
html_static_path = ["_static"]
html_js_files = [
    ("custom-icons.js", {"defer": "defer"}),
]
html_theme_options = {
    "external_links": [
        {
            "url": "https://aprilweilab.github.io",
            "name": "Wei Lab",
        }
    ],
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/aprilweilab/grgl",
            "icon": "fa-brands fa-github",
        },
        {
            "name": "PyPI",
            "url": "https://pypi.org/project/pygrgl",
            "icon": "fa-custom fa-pypi",
        }
    ],
}

# Breathe configuration
breathe_projects = {
    "grgl": DOC_OUTPUT + "/xml/",
}
breathe_default_project = "grgl"

autosummary_generate = True  # Turn on sphinx.ext.autosummary

project = "GRGL"

exclude_patterns = ['_build', '**.ipynb_checkpoints']

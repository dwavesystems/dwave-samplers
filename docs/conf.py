# -*- coding: utf-8 -*-

# This file contains function linkcode_resolve, based on
# https://github.com/numpy/numpy/blob/main/doc/source/conf.py,
# which is licensed under the BSD 3-Clause "New" or "Revised"
# license: ./licenses/numpy.rst

import os
import sys


# root_directory = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# sys.path.insert(len(sys.path), os.path.abspath("."))
# sys.path.insert(len(sys.path), root_directory)

from dwave.samplers import __version__ as version
from dwave.samplers import __version__ as release

# -- Project information - these are special values used by sphinx. -------

copyright = "D-Wave Systems Inc."

project = "D-Wave Samplers"

# -- General configuration ------------------------------------------------

extensions = [
    "sphinx.ext.autosummary",
    "sphinx.ext.autodoc",
    "sphinx.ext.coverage",
    "sphinx.ext.doctest",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx.ext.todo",
    "sphinx.ext.githubpages",
    "sphinx.ext.ifconfig",
    'reno.sphinxext',
]

autosummary_generate = True

source_suffix = [".rst", ".md"]

root_doc = "index"  # before Sphinx 4.0, named master_doc

add_module_names = False

exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

pygments_style = "sphinx"

todo_include_todos = True

doctest_global_setup = ""


# -- Options for HTML output ----------------------------------------------

html_theme = "pydata_sphinx_theme"

html_theme_options = {
    "collapse_navigation": True,
    "show_prev_next": False,
}
html_sidebars = {"**": ["search-field", "sidebar-nav-bs"]}  # remove ads
#html_static_path = ["_static"]


# -- Panels ---------------------------------------------------------------
panels_add_bootstrap_css = False

# -- Intersphinx ----------------------------------------------------------
intersphinx_mapping = {
    "python": ("https://docs.python.org/3/", None),
    "numpy": ("http://docs.scipy.org/doc/numpy/", None),
    "oceandocs": ("https://docs.ocean.dwavesys.com/en/stable/", None),
}

# breathe, for pulling in doxygen generated C++
breathe_projects = {"cqmsolver": "xml"}
breathe_default_project = "cqmsolver"


# -- linkcheck -------------------------------------------------------------

# this can be removed once we've deployed the package
linkcheck_ignore = [r'https://pypi.python.org/pypi/dwave-samplers']

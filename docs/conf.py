# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys

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
    "sphinx_design",
    "reno.sphinxext",
]

autosummary_generate = True

source_suffix = [".rst", ".md"]

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

# TODO: replace oceandocs & sysdocs_gettingstarted
intersphinx_mapping = {
    "python": ("https://docs.python.org/3/", None),
    "numpy": ("http://docs.scipy.org/doc/numpy/", None),
    "oceandocs": ("https://docs.ocean.dwavesys.com/en/stable/", None),
}

# breathe, for pulling in doxygen generated C++
breathe_projects = {"cqmsolver": "xml"}
breathe_default_project = "cqmsolver"
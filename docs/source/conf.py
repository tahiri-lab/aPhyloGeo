# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
import os
import sys
sys.path.insert(0, os.path.abspath('../'))

project = 'aPhyloGeo'
copyright = '2024, Tahiri Lab'
author = 'Tahiri Lab'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
        'sphinx.ext.duration',
        'sphinx.ext.doctest',
        'sphinx.ext.autodoc',
        'sphinx.ext.autosummary',
        'sphinx.ext.coverage',
        'sphinx.ext.napoleon',
        'sphinx.ext.viewcode',
        'sphinx.ext.githubpages',
        'sphinx_rtd_theme',
        'sphinx.ext.intersphinx',
        ]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
autodoc_member_order = 'bysource'


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_context = {
    'display_github': True,
    'github_user': 'tahiri-lab',
    'github_repo': 'aPhiloGeo',
    'github_version': 'main/',
    'conf_py_path': 'docs/',
    'source_suffix': '.rst',
    'commit': False,
}

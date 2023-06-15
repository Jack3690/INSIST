""" Configuration file for the Sphinx documentation builder."""

# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import sys

sys.path.insert(0, os.path.abspath('.'))
sys.path.insert(0, os.path.abspath("../../src"))
package_path = os.path.abspath('../../src')
os.environ['PYTHONPATH'] = ':'.join((package_path, os.environ.get('PYTHONPATH', '')))

project = 'insist-pista'
copyright = '2023, Avinash CK and INSIST Science Team'
author = 'Avinash CK'
release = '1.0.27'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.napoleon',
	    'sphinx.ext.autodoc',
	    'sphinx.ext.viewcode',
	    'sphinx.ext.coverage',
	    'sphinx.ext.autosectionlabel',
	    'jupyter_sphinx',
]
templates_path = ['_templates']
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']


# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".

#Stuff below is scavenged from the web
html_logo = "../logo.png"
html_theme_options = {
    "navigation_with_keys": True,
}

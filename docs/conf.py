# -*- coding: utf-8 -*-

import sys
import os
import re

sys.path.insert(0, os.path.abspath("../../src"))
package_path = os.path.abspath('../../src')
os.environ['PYTHONPATH'] = ':'.join((package_path, os.environ.get('PYTHONPATH', '')))

from sphinx_rtd_theme import __version__ as theme_version
from sphinx_rtd_theme import __version_full__ as theme_version_full
from sphinx.locale import _

sys.path.insert(0, os.path.abspath('.'))
sys.path.insert(0, os.path.abspath("../../src"))
package_path = os.path.abspath('../../src')
os.environ['PYTHONPATH'] = ':'.join((package_path, os.environ.get('PYTHONPATH', '')))

author = u'Avinash CK and INSIST Science team at the Indian Institute of Astrophysics '
release = 1.0.0
copyright = author
language = 'en'

extensions = [
    'sphinx.ext.intersphinx',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
    'sphinx_rtd_theme',
]

templates_path = ['_templates']
source_suffix = '.rst'
exclude_patterns = []

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static/']

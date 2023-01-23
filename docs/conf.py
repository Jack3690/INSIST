# -*- coding: utf-8 -*-

import sys
import os
import re

sys.path.insert(0, os.path.abspath("../../src"))
package_path = os.path.abspath('../../src')
os.environ['PYTHONPATH'] = ':'.join((package_path, os.environ.get('PYTHONPATH', '')))


sys.path.insert(0, os.path.abspath('.'))
sys.path.insert(0, os.path.abspath("../../src"))
package_path = os.path.abspath('../../src')
os.environ['PYTHONPATH'] = ':'.join((package_path, os.environ.get('PYTHONPATH', '')))

author = u'Avinash CK and INSIST Science team at the Indian Institute of Astrophysics '
release = '1.0.0'
copyright = author
language = 'en'

extensions = ['sphinx.ext.napoleon',
    'sphinx.ext.autodoc',
    'sphinx.ext.viewcode',
    'sphinx.ext.coverage',
    'sphinx.ext.autosectionlabel',
    'jupyter_sphinx',
]

templates_path = ['_templates']
source_suffix = '.rst'
exclude_patterns = []

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static/']

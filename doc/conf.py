# -*- coding: utf-8 -*-
#
# Fabber documentation build configuration file
#
# This file is execfile()d with the current directory set to its
# containing dir.

# -- General configuration ------------------------------------------------
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.viewcode',
    'sphinx.ext.mathjax'
]

project = u'CEST models for Fabber'
copyright = u'2018, Martin Craig'
author = u'Martin Craig'
build_dir = u"_build"

version = u''
release = u''
language = None
templates_path = ['_templates']
source_suffix = '.rst'
master_doc = 'index'
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# -- Options for HTML output ----------------------------------------------

import sphinx_rtd_theme
html_theme = 'sphinx_rtd_theme'
html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]

# -- Options for LaTeX output ---------------------------------------------

latex_elements = {
}

latex_documents = [
    (master_doc, 'fabber_cest.tex', u'Fabber CEST documentation',
     u'Martin Craig', 'manual'),
]

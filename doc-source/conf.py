#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# This file is managed by `git_helper`. Don't edit it directly

import os
import re
import sys

sys.path.append(os.path.abspath('.'))
sys.path.append(os.path.abspath('..'))

from sphinx.locale import _

from __pkginfo__ import __version__

import ipynb2rst
nitpicky = True
autodoc_default_options = {'autosummary': True}

github_url = f"https://github.com/domdfcoding/PyMassSpec"

rst_prolog = f""".. |pkgname| replace:: PyMassSpec
.. |pkgname2| replace:: ``PyMassSpec``
.. |browse_github| replace:: `Browse the GitHub Repository <{github_url}>`__
.. |ghurl| replace:: {github_url}
"""

author = "PyMassSpec Authors"
project = "PyMassSpec"
slug = re.sub(r'\W+', '-', project.lower())
release = version = __version__
copyright = "Copyright 2019-2020 Dominic Davis-Foster"
language = 'en'
package_root = "pyms"

extensions = [
		'sphinx.ext.intersphinx',
		'sphinx.ext.autodoc',
		'sphinx.ext.mathjax',
		'sphinx.ext.viewcode',
		'sphinxcontrib.httpdomain',
		"sphinxcontrib.extras_require",
		"sphinx.ext.todo",
		"sphinxemoji.sphinxemoji",
		'autodocsumm',
'nbsphinx',

		]

sphinxemoji_style = 'twemoji'
todo_include_todos = bool(os.environ.get("SHOW_TODOS", False))

templates_path = ['_templates']
html_static_path = ['_static']
source_suffix = '.rst'
exclude_patterns = []

master_doc = 'index'
suppress_warnings = ['image.nonlocal_uri']
pygments_style = 'default'

intersphinx_mapping = {
		'rtd': ('https://docs.readthedocs.io/en/latest/', None),
		'sphinx': ('http://www.sphinx-doc.org/en/stable/', None),
		'python': ('https://docs.python.org/3/', None),
		"NumPy": ('https://numpy.org/doc/stable/', None),
		"SciPy": ('https://docs.scipy.org/doc/scipy/reference', None),
		"matplotlib": ('https://matplotlib.org', None),
		"h5py": ('https://docs.h5py.org/en/latest/', None),
		"Sphinx": ('https://www.sphinx-doc.org/en/stable/', None),
		"Django": ('https://docs.djangoproject.com/en/dev/', 'https://docs.djangoproject.com/en/dev/_objects/'),
		"sarge": ('https://sarge.readthedocs.io/en/latest/', None),
		"attrs": ('https://www.attrs.org/en/stable/', None),

		}

html_theme = 'sphinx_rtd_theme'
html_theme_options = {
		'logo_only': False,  # True will show just the logo
		'includehidden': False,

		}
html_theme_path = ["../.."]
# html_logo = "logo/pyms.png"
html_show_sourcelink = False  # True will show link to source

html_context = {
		# Github Settings
		"display_github": True,  # Integrate GitHub
		"github_user": "domdfcoding",  # Username
		"github_repo": "PyMassSpec",  # Repo name
		"github_version": "master",  # Version
		"conf_py_path": "/",  # Path in the checkout to the docs root
		}

htmlhelp_basename = slug

latex_documents = [
		('index', '{0}.tex'.format(slug), project, author, 'manual'),
		]

man_pages = [
		('index', slug, project, [author], 1)
		]

texinfo_documents = [
		('index', slug, project, author, slug, project, 'Miscellaneous'),
		]


# Extensions to theme docs
def setup(app):
	from sphinx.domains.python import PyField
	from sphinx.util.docfields import Field

	app.add_object_type(
			'confval',
			'confval',
			objname='configuration value',
			indextemplate='pair: %s; configuration value',
			doc_field_types=[
					PyField(
							'type',
							label=_('Type'),
							has_arg=False,
							names=('type',),
							bodyrolename='class'
							),
					Field(
							'default',
							label=_('Default'),
							has_arg=False,
							names=('default',),
							),
					]
			)

# Configuration file for the Sphinx documentation builder.

# -- Project information

project = 'pykegg'
copyright = '2024, Noriaki Sato'
author = 'Noriaki Sato'

release = '0.1'
version = '0.0.5'

# -- General configuration
import os
import sys
# sys.path.insert(0, os.path.abspath('../../'))

extensions = [
    'nbsphinx',
    'sphinx.ext.duration',
    'sphinx.ext.napoleon',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
]

autosummary_generate = True
intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']

# -- Options for HTML output

html_theme = 'sphinx_book_theme'
html_use_opensearch = "http://pykegg.readthedocs.io"
html_search_language = "en"


# -- Options for EPUB output
epub_show_urls = 'footnote'

source_suffix = '.rst'

# html_static_path = ['_static']
# html_style = 'custom.css'
# html_css_files = ["custom.css"]

def run_apidoc(_):
    from sphinx.ext.apidoc import main
    import os
    import sys
    sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
    cur_dir = os.path.abspath(os.path.dirname(__file__))
    module = os.path.join(cur_dir,"..","src", "pykegg")
    main(['--force', '-T', '-E',
        '-o', cur_dir, module])

def setup(app):
    app.connect('builder-inited', run_apidoc)
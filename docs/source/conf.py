author = 'KavrakiLab'

release = '0.1'
version = '0.0.5'

# -- General configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    'sphinx_rtd_theme',
    'nbsphinx'
]

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']

# -- Options for HTML output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_logo = "_static/logo.svg"
html_theme_options = {
    'logo_only': True,
    'display_version': False,
}

# -- Options for EPUB output
epub_show_urls = 'footnote'

# -- Copy the notebooks to the right place, from https://github.com/spatialaudio/nbsphinx/issues/170

import shutil
import os

print("Copy example notebooks into docs/_examples")

def all_but_ipynb(dir, contents):
    result = []
    for c in contents:
        if os.path.isfile(os.path.join(dir,c)) and (not c.endswith(".ipynb")):
            result += [c]
    return result

project_root = os.path.dirname(os.path.abspath(__file__))
shutil.rmtree(os.path.join(project_root, "./_notebooks"), ignore_errors=True)
shutil.copytree(os.path.join(project_root, "../../notebooks"),
                os.path.join(project_root, "./_notebooks"),
                ignore=all_but_ipynb)

# -- Setup engens environment

os.chdir('..')
os.system('./linux_setup.sh')

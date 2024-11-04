# -*- coding: utf-8 -*-
#
import sys, os

source_parsers = {
   '.md': 'recommonmark.parser.CommonMarkParser',
}
source_suffix = '.rst'

master_doc = 'index'

project = u'Appositionizer'

# The short X.Y version.
version = ''  # Is set by calling `setup.py docs`
# The full version, including alpha/beta/rc tags.
release = ''  # Is set by calling `setup.py docs`

pygments_style = 'sphinx'

html_theme = 'sphinx-bluebrain-theme'
html_static_path = ['_static']

# Any file that's being included from another file, e.g.,
# :include: ...
# Must be excluded here to prevent duplication of labels.
exclude_patterns = []

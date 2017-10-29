#!/usr/bin/python -S
"""
highlight.py
"""

__author__ = 'Andy Chu'


import cgi
from ._jsontemplate import Template, SplitMeta, MakeTokenRegex


_TEMPLATE = None

COMMENT, DIRECTIVE, SUBSTITUTION = range(3)


def AsHtml(template_str, meta='{}', format_char='|'):

  global _TEMPLATE
  if not _TEMPLATE:
    _TEMPLATE = Template(
        '<span style="color: {color|htmltag};">{token|html}</span>')

  meta_left, meta_right = SplitMeta(meta)
  token_re = MakeTokenRegex(meta_left, meta_right)
  tokens = token_re.split(template_str)

  html = []

  for i, token in enumerate(tokens):

    # By the definition of re.split, even tokens are literal strings, and odd
    # tokens are directives.
    if i % 2 == 0:
      html.append(cgi.escape(token))
    else:
      # Because of the regex, the token should be at least 2 characters long
      c = token[1]

      if c == '#':
        token_type = COMMENT
      elif c == '.':
        token_type = DIRECTIVE
      else:
        token_type = SUBSTITUTION

      # TODO: Use classes, and make comments italic
      color = {
          COMMENT: 'red',
          DIRECTIVE: 'blue',
          SUBSTITUTION: 'green'
          }[token_type]
      html.append(_TEMPLATE.expand({'color': color, 'token': token}))

  # Without <pre>, we would have to turn newlines into line breaks, etc.
  return '<pre>%s</pre>' % ''.join(html)


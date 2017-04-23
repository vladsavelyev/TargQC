#!/usr/bin/python -S
"""
formatters.py

This module should implement the standard list of formatters.

It also provides a method LookupChain for *composing lookup chains* for
formatters.

Formatter lookup chaining is not to be confused with plain formatter chaining,
e.g.:

  {variable|html|json}

If anyone has any better names for the two types of chaining, let the mailing
list know.
"""

__author__ = 'Andy Chu'


import os
import sys

import _jsontemplate as jsontemplate  # For TemplateFileInclude


class Error(Exception):
  """Base class for all exceptions raised by this module."""


def LookupChain(lookup_func_list):
  """Returns a *function* suitable for passing as the more_formatters argument
  to Template.

  NOTE: In Java, this would be implemented using the 'Composite' pattern.  A
  *list* of formatter lookup function behaves the same as a *single* formatter
  lookup funcion.

  Note the distinction between formatter *lookup* functions and formatter
  functions here.
  """
  def MoreFormatters(formatter_name):
    for lookup_func in lookup_func_list:
      formatter_func = lookup_func(formatter_name)
      if formatter_func is not None:
        return formatter_func

  return MoreFormatters


def PythonPercentFormat(format_str):
  """Use Python % format strings as template format specifiers."""

  if format_str.startswith('printf '):
    fmt = format_str[len('printf '):]
    return lambda value: fmt % value
  else:
    return None


# Seam for testing
_open = open

# Cache of compiled templates.  In Java, this might need to be a
# ConcurrentHashMap like the tokenization regex cache.
_compiled_template_cache = {}


class TemplateFileInclude(object):
  """Template include mechanism.

  The relative path is specified as an argument to the template.
  """

  def __init__(self, root_dir):
    self.root_dir = root_dir

  def __call__(self, format_str):
    """Returns a formatter function."""

    if format_str.startswith('template-file '):
      relative_path = format_str[len('template-file '):]
      full_path = os.path.join(self.root_dir, relative_path)

      if full_path not in _compiled_template_cache:
        f = _open(full_path)
        _compiled_template_cache[full_path] = jsontemplate.FromFile(f)
        f.close()

      return _compiled_template_cache[full_path].expand  # a 'bound method'

    else:
      return None  # this lookup is not applicable


class Json(object):
  """Format arbitrary nodes as JSON.
  
  It takes a function which converts JSON structures to strings as a parameter.

  All this does is relieve the user of having to remember the standard names
  'json' and 'js-string'.  Just pass your program's JSON serializer in here.
  """

  def __init__(self, json_func):
    self.json_func = json_func

  def __call__(self, format_str):
    """Returns a formatter function."""
    if format_str in ('json', 'js-string'):
      return self.json_func

    else:
      return None  # this lookup is not applicable


def Plural(format_str):
  """Returns whether the value should be considered a plural value.

  Integers greater than 1 are plural, and lists with length greater than one are
  too.
  """
  if format_str.startswith('plural?'):
    i = len('plural?')

    try:
      splitchar = format_str[i]  # Usually a space, but could be something else
      _, plural_val, singular_val = format_str.split(splitchar)
    except IndexError:
      raise Error('plural? must have exactly 2 arguments')

    def Formatter(value):
      plural = False
      if isinstance(value, int) and value > 1:
        plural = True
      if isinstance(value, list) and len(value) > 1:
        plural = True

      if plural:
        return plural_val
      else:
        return singular_val

    return Formatter

  else:
    return None  # this lookup is not applicable


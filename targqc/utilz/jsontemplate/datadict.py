#!/usr/bin/python -S
"""
datadict.py

This module demonstrates (recursive) mutations of the data dictionary.  You may
want to use a similar technique in various circumstances.

This may be combined with modifying the Template.render() or expand() functions
so that they do the mutation first.
"""

__author__ = 'Andy Chu'


def AddIndex(node, index=None):
  """
  Recursively add the current index (with respect to a repeated section) in all
  data dictionaries.
  """
  if isinstance(node, list):
    for i, item in enumerate(node):
      AddIndex(item, index=i)
  elif isinstance(node, dict):
    if index is not None:
      node['index'] = index
    for key in node:
      AddIndex(node[key])

#--------------------------------------------------------------------------
#     This file is part of OASA - a free chemical python library
#     Copyright (C) 2003-2008 Beda Kosata <beda@zirael.org>

#     This program is free software; you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation; either version 2 of the License, or
#     (at your option) any later version.

#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.

#     Complete text of GNU GPL can be found in the file gpl.txt in the
#     main directory of the program

#--------------------------------------------------------------------------

"""Vertex class.

"""

import copy


class vertex(object):
  """Simple vertex class.

  Normaly would not be needed but it can speed up many analytical tasks
  to store data directly in vertex and not get them from the graph
  connectivity matrix.
  Vertex has a value attribute used to store arbitrary objects.
  """
  attrs_to_copy = ("value",)

  def __init__(self):
    self.properties_ = {} # used to store intermediate properties such as distances etc.
    self.value = None  # used to store any object associated with the vertex
    self._neighbors = {} # set of all neighbors
    self._clean_cache()


  def __str__(self):
    return ("vertex, value=%s, degree=%d, " % (str(self.value), self.degree)) + str(self.properties_)


  def _clean_cache(self):
    self._cache = {}


  def copy(self):
    other = self.__class__()
    for attr in self.attrs_to_copy:
      setattr(other, attr, copy.copy(getattr(self, attr)))
    return other


  def add_neighbor(self, v, e):
    """adds a neighbor connected via e"""
    self._clean_cache()
    self._neighbors[e] = v


  def remove_neighbor(self, v):
    self._clean_cache()
    to_del = None
    for k, vv in self._neighbors.items():
      if v == vv:
        to_del = k
        break
    if to_del:
      del self._neighbors[to_del]
    else:
      raise Exception("Cannot remove non-existing neighbor")


  def remove_edge_and_neighbor(self, e):
    self._clean_cache()
    if e in self._neighbors.keys():
      del self._neighbors[e]
    else:
      raise Exception("Cannot remove non-existing edge", e)


  @property
  def neighbors(self):
    """Neighboring vertices.

    """
    return [v for (e,v) in self._neighbors.items()
                  if not e.disconnected]


  def get_neighbor_connected_via( self, e):
    return self._neighbors[e]


  def get_edge_leading_to(self, a):
    for b, at in self._neighbors.items():
      if a == at:
        return b
    return None


  @property
  def degree(self):
    """Degree of the vertex.

    """
    return len(self.neighbors)


  def get_neighbors_with_distance( self, d):
    ret = []
    for v in self.neighbors:
      if 'd' in v.properties_ and v.properties_['d'] == d:
        ret.append(v)
    return ret


  def get_neighbor_edge_pairs(self):
    for e, v in self._neighbors.items():
      if not e.disconnected:
        yield e,v


  @property
  def neighbor_edges(self):
    """Neighboring edges.

    """
    return [e for e in self._neighbors.keys()
                if not e.disconnected]


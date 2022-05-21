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

import copy

from warnings import warn



class edge(object):

  attrs_to_copy = ("disconnected",)

  def __init__(self, vs=[]):
    self._vertices = []
    self.set_vertices(vs)
    self.properties_ = {}
    self.disconnected = False


  def __str__(self):
    return "edge between %s %s" % tuple(map(str, self.vertices))


  def copy(self):
    other = self.__class__()
    for attr in self.attrs_to_copy:
      setattr(other, attr, copy.copy(getattr(self, attr)))
    return other


  def set_vertices(self, vs=[]):
    # Ring perception algorithm relies on allowing both vertices to be the same
    if vs and len(vs) == 2:
      self._vertices = list(vs)


  def get_vertices(self):
    return self._vertices


  @property
  def neighbor_edges(self):
    v1, v2 = self.vertices
    out1 = [e for e in v1.neighbor_edges if e != self]
    out2 = [e for e in v2.neighbor_edges if e != self]
    return out1 + out2


  def get_neighbor_edges2(self):
    """Return 2 lists of neighbor edges (one for one side, one for the other).

    """
    v1, v2 = self.vertices
    out1 = [e for e in v1.neighbor_edges if e != self]
    out2 = [e for e in v2.neighbor_edges if e != self]
    return out1, out2


  @property
  def disconnected(self):
    return self._disconnected


  @disconnected.setter
  def disconnected(self, d):
    self._disconnected = d


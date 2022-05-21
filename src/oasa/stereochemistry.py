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

from . import oasa_exceptions



class stereochemistry(object):

  def __init__( self, center=None, references=None, value=None):
    self.center = center
    if not references:
      self.references = []
    else:
      self.references = references
    self.value = value


  @property
  def value(self):
    return self._value


  @value.setter
  def value(self, value):
    self._value = value


  @property
  def references(self):
    return self._references


  @references.setter
  def references(self, references):
    self._references = references


  @property
  def center(self):
    return self._center


  @center.setter
  def center(self, center):
    self._center = center



class cis_trans_stereochemistry(stereochemistry):

  UNDEFINED = 0
  OPPOSITE_SIDE = 1
  SAME_SIDE = 2

  # Override value
  @property
  def value(self):
    return stereochemistry.value.__get__(self)


  @value.setter
  def value(self, value):
    if value not in (self.SAME_SIDE, self.OPPOSITE_SIDE, self.UNDEFINED):
      raise oasa_exceptions.oasa_stereochemistry_error( "invalid stereochemistry identifier '%s'" % value)
    stereochemistry.value.__set__(self, value)


  # Override references
  @property
  def references(self):
    return stereochemistry.references.__get__(self)


  @references.setter
  def references(self, references):
    if len( references) != 4:
      raise oasa_exceptions.oasa_stereochemistry_error( "wrong number of references in stereochemistry specification '%s'" % len( references))
    stereochemistry.references.__set__(self, references)


  def get_other_end( self, ref):
    if not ref in self.references:
      raise ValueError("Submitted object is not referenced in this stereochemistry object.")
    ref1, _r1, _r2, ref2 = self.references
    return ref is ref1 and ref2 or ref1



class tetrahedral_stereochemistry(stereochemistry):

  UNDEFINED = 0
  # as in CIP - put the last atom behind the reference and observer the first 3
  # this should create the same result as SMILES notation - put firts in front and observe the rest
  CLOCKWISE = 1
  ANTICLOCKWISE = 2

  # Override value
  @property
  def value(self):
    return stereochemistry.value.__get__(self)


  @value.setter
  def value(self, value):
    if value not in (self.CLOCKWISE, self.ANTICLOCKWISE, self.UNDEFINED):
      raise oasa_exceptions.oasa_stereochemistry_error( "invalid stereochemistry identifier '%s'" % value)
    stereochemistry.value.__set__(self, value)


  # Override references
  @property
  def references(self):
    return stereochemistry.references.__get__(self)


  @references.setter
  def references(self, references):
    if len( references) != 4:
      raise oasa_exceptions.oasa_stereochemistry_error( "wrong number of references in stereochemistry specification '%s'" % len( references))
    stereochemistry.references.__set__(self, references)



class explicit_hydrogen(object):
  """Placeholder for explicit hydrogen in stereochemistry references.

  """
  def __eq__(self, other):
    if isinstance(other, explicit_hydrogen):
      return True
    else:
      return False


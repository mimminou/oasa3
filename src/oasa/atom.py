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

import sys
import copy
import itertools

from warnings import warn
from operator import itemgetter

from . import graph
from . import periodic_table as PT
from .chem_vertex import chem_vertex
from .common import is_uniquely_sorted
from .oasa_exceptions import oasa_invalid_atom_symbol



class atom(chem_vertex):
  ## ("value","charge","x","y","z","multiplicity","valency","charge","free_sites")
  attrs_to_copy = chem_vertex.attrs_to_copy + ("symbol", "isotope","explicit_hydrogens")

  def __init__( self, symbol='C', charge=0, coords=None):
    chem_vertex.__init__( self, coords=coords)
    self.symbol = symbol
    self.charge = charge
    self.isotope = None
    self.explicit_hydrogens = 0


  def matches( self, other):
    if not isinstance( other, atom):
      return False
    if self.symbol == other.symbol and self.valency == other.valency and self.multiplicity == other.multiplicity:
      # check charge only if other has it set to something non-zero
      if other.charge and self.charge != other.charge:
        return False
      return True
    return False


  @property
  def symbol(self):
    """Atom symbol.

    """
    return self._symbol


  @symbol.setter
  def symbol(self, symbol):
    self._clean_cache()
    try:
      self.valency = PT.periodic_table[ symbol]['valency'][0]
      self.symbol_number = PT.periodic_table[ symbol]['ord']
    except KeyError:
      raise oasa_invalid_atom_symbol( "invalid atom symbol", symbol)
    self._symbol = symbol


  # Overrides chem_vertex occupied_valency
  @property
  def occupied_valency(self):
    """Atoms occupied valency.

    """
    try:
      return self._cache['occupied_valency']
    except KeyError:
      pass

    bonds_alternating_aromatic = 0
    bonds_single_aromatic = 0
    odd_aromatic = False
    for b in self._neighbors.keys():
      order = b.order
      if order == 4:
        if not odd_aromatic:
          order = 1
        else:
          order = 2
        bonds_single_aromatic += 1
        bonds_alternating_aromatic += order
        odd_aromatic = not odd_aromatic
      else:
        bonds_single_aromatic += order
        bonds_alternating_aromatic += order

    if self.charge:
      # now we deal with charge
      if abs( self.charge) > 1:
        # charges higher than one should always decrease valency
        charge = abs( self.charge)
      elif (self.symbol in PT.accept_cation) and (self.charge == 1) and (self.valency <= PT.accept_cation[self.symbol]):
        # elements that can accept cations to increase their valency (NH4+)
        charge = -1
      elif (self.symbol in PT.accept_anion) and (self.charge == -1) and (self.valency <= PT.accept_anion[self.symbol]):
        # elements that can accept anions to increase their valency (BH4-)
        charge = -1
      else:
        # otherwise charge reduces valency
        charge = abs( self.charge)
    else:
      charge = 0

    x = bonds_alternating_aromatic+charge+self.multiplicity-1+self.explicit_hydrogens
    if x > self.valency:
      # we have computed occupied_valency using alternating single and double
      # bonds for aromatic bonds, however it led to occupied_valency being higher
      # than valency - lets replace it using single bonds for aromatic
      # (this fixed thiophene where occupied_valency of S would be computed to be 3
      #  and valency raise would be triggered)
      x = bonds_single_aromatic+charge+self.multiplicity-1+self.explicit_hydrogens
    #self._cache['occupied_valency'] = x
    return x


  # Overrides chem_vertex multiplicity
  @property
  def multiplicity(self):
    """Atom multiplicity.

    """
    return self._multiplicity


  @multiplicity.setter
  def multiplicity(self, multiplicity):
    chem_vertex.multiplicity.__set__(self, multiplicity)
    if self.free_valency < 0:
      self.raise_valency_to_senseful_value()


  @property
  def free_sites(self):
    """Atoms free_sites.

    """
    if self._free_sites > self.free_valency:
      return self.free_valency
    return self._free_sites


  @free_sites.setter
  def free_sites(self, free_sites):
    self._free_sites = free_sites


  @property
  def isotope(self):
    """Isotope.

    """
    return self._isotope


  @isotope.setter
  def isotope(self, isotope):
    if isotope is not None and not isinstance(isotope, int):
      # isotope must be a number or None
      raise oasa_exceptions.oasa_invalid_value( "isotope", isotope)
    self._isotope = isotope


  @property
  def electronegativity(self):
    """Atom's electronegativity.

    """
    try:
      return PT.periodic_table[self.symbol]['en']
    except KeyError:
      return None


  @property
  def oxidation_number(self):
    """Atom's oxidation number.

    """
    en = self.charge
    for e, n in self.get_neighbor_edge_pairs():
      if isinstance( n, atom) and n.symbol != self.symbol:
        en += e.order * (n.electronegativity > self.electronegativity and 1 or -1)
    hen = PT.periodic_table['H']['en']
    en += self.free_valency * (hen > self.electronegativity and 1 or -1)
    return en


  @property
  def electron_pairs(self):
    """Number of electron pairs on the atom.

    """
    return (PT.periodic_table[self.symbol]['els'] - sum([b.order for b in self.neighbor_edges]) - self.charge - self.free_valency - self.multiplicity + 1) / 2.0


  def __str__( self):
    return "atom '%s'" % str( self.symbol)


  def get_formula_dict(self):
    """Return formula as dictionary.

    That can be passed to functions in periodic_table.
    """
    ret = PT.formula_dict( self.symbol)
    if self.free_valency + self.explicit_hydrogens > 0:
      ret['H'] = self.free_valency + self.explicit_hydrogens
    return ret


  def raise_valency_to_senseful_value(self):
    """Set atom valency to the lowest possible.

    So that free_valency if non-negative (when possible) or highest possible,
    does not lower valency when set to higher than necessary value.
    """
    while self.free_valency < 0:
      if not self.raise_valency():
        return


  def raise_valency(self):
    """Used if valency < occupied_valency to try to find a higher one.

    """
    for v in PT.periodic_table[self.symbol]['valency']:
      if v > self.valency:
        self.valency = v
        #self._clean_cache()
        return True
    return False


  def get_hydrogen_count(self):
    return self.explicit_hydrogens + self.free_valency


  def is_chiral(self):
    """This code is CIP (Cahn-Ingold-Prelog) based.

    And therefore not necessarily the fastest for this job,
    however in newer versions it will be able to take care
    of the chirality of other centres and therefor will be universal,
    for now it takes only care of the connectivity.
    """
    if len( self._neighbors) < 4:
      # this is not true for N,P and similar !!!
      return False
    cips = []
    for a in self.neighbors:
      cips.append( [[], a, a.gen_CIP_sequence( came_from = self)])
    while not is_uniquely_sorted( cips, cip_sorting_function):
      for cip in cips:
        try:
          cip[0].append(next(cip[2]).symbol_number)
        except StopIteration:
          cip[0].append( None)
        except AttributeError:
          cip[0].append( None)  # None was yielded and it has no symbol_number :)
      # to test if we have already reached the end
      for i in range( len( cips)-1):
        if cips[i][0] == cips[i+1][0] and cips[i][0][-3:-1] == [None, None]:
          return False
      cips.sort(key=itemgetter(0), reverse=True)
    return True


  def get_neighbors_CIP_sorted(self):
    """Return neighbors sorted according to the CIP rules.

    """
    cips = []
    for a in self.neighbors:
      cips.append( [[], a, a.gen_CIP_sequence( came_from = self)])
    while not is_uniquely_sorted( cips, cip_sorting_function):
      for cip in cips:
        try:
          cip[0].append(next(cip[2]).symbol_number)
        except StopIteration:
          cip[0].append( None)
        except AttributeError:
          cip[0].append( None)  # None was yielded and it has no symbol_number :)
      # to test if we have already reached the end
      all_finished = len( cip[0]) > 1
      for cip in cips:
        if len( cip[0]) > 1:
          if cip[0][-3:-1] != [None, None]:
            all_finished = 0
      if all_finished:
        break # it can't be uniquely sorted, we are out of atoms
      cips.sort(key=itemgetter(0), reverse=True)
    return [cip[1] for cip in cips]


  def gen_CIP_sequence(self, iter_over=None, came_from=None):
    """Generate the CIP (Cahn-Ingold-Prelog) stream of atoms.

    Suitable for comparison in searches for chiral centres,
    their configuration etc.
    The values in different layers (with raising distance from self) are
    separated by Nones.
    """
    yield self
    yield None
    neighs = self.neighbors
    if came_from:
      assert came_from in neighs
      neighs.remove( came_from)
    # end if there are no neighs
    if not neighs:
      raise StopIteration
    # generators of the neighs
    cips = []
    for a in neighs:
      cips.append( [[], [], a.gen_CIP_sequence( came_from = self), 0])
      #         number_list, atom_list, cip_generator, last yielded index
    while cips:
      # this returns the already generated cip of the new first atom in later cycles
      to_remove = []
      for cip in cips:
        try:
          cs = itertools.takewhile(lambda x: x is not None, cip[2])
        except StopIteration:
          to_remove.append( cip)
          continue
        cip[1] += cs
        cip[0] += [c.symbol_number for c in cs]
      [cips.remove( cip) for cip in to_remove]
      cips.sort(key=itemgetter(0), reverse=True)
      for cip in cips:
        for i in range( cip[3], len( cip[1])):
          yield cip[1][i]
        cip[3] = len( cip[1])
      yield None


  def get_highest_possible_free_valency(self):
    """Used in case of aromatic bonds.

    Takes all aromatic bonds as single, thus giving the maximum free valency
    that would be possible if all these localized to single.
    """
    return self.valency - chem_vertex.occupied_valency.__get__(self)



def cip_sorting_function(a, b):
  if a[0] == b[0]:
    return 0
  elif a[0] < b[0]:
    return 1
  elif a[0] > b[0]:
    return -1



##################################################
# TODO

# chirality for those possesing a free electron pair

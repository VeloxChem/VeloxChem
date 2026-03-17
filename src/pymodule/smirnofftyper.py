#
#                                   VELOXCHEM
#              ----------------------------------------------------
#                          An Electronic Structure Code
#
#  SPDX-License-Identifier: BSD-3-Clause
#
#  Copyright 2018-2025 VeloxChem developers
#
#  Redistribution and use in source and binary forms, with or without modification,
#  are permitted provided that the following conditions are met:
#
#  1. Redistributions of source code must retain the above copyright notice, this
#     list of conditions and the following disclaimer.
#  2. Redistributions in binary form must reproduce the above copyright notice,
#     this list of conditions and the following disclaimer in the documentation
#     and/or other materials provided with the distribution.
#  3. Neither the name of the copyright holder nor the names of its contributors
#     may be used to endorse or promote products derived from this software without
#     specific prior written permission.
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
#  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
#  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
#  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
#  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
#  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
#  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
#  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
#  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

"""
smirnofftyper.py — Self-contained SMIRNOFF parameter applier
=============================================================
Applies OpenFF-format (.offxml) force fields to molecules using only RDKit
and Python stdlib (xml.etree, urllib). 

Supports  : Bonds, Angles, ProperTorsions, ImproperTorsions, vdW.
No support: virtual sites, fractional bond-order interpolation,
            library charges (charges handled externally by VeloxChem RESP).

Assignment rule (SMIRNOFF spec §2 — "last match wins"):
    Iterate all patterns in file order; the LAST match for each interaction
    wins.  More-specific patterns appear later in the file and thus override
    generic ones.

Output unit conventions (matching VeloxChem / GROMACS internals):
    lengths     nm
    angles      degrees
    bond k      kJ mol⁻¹ nm⁻²
    angle k     kJ mol⁻¹ rad⁻²
    dihedral k  kJ mol⁻¹   (idivf already applied)
    epsilon     kJ mol⁻¹
    sigma       nm

"""

import re
import math
import xml.etree.ElementTree as ET
import urllib.request


# Unit conversion constants 

_KCAL_TO_KJ        = 4.184          # 1 kcal/mol → kJ/mol
_ANG_TO_NM         = 0.1            # 1 Å → nm
# 1 kcal mol⁻¹ Å⁻²  →  kJ mol⁻¹ nm⁻²   (÷ 0.1² = ×100, then ×4.184)
_KCAL_ANG2_TO_KJ_NM2 = _KCAL_TO_KJ / (_ANG_TO_NM ** 2)
# 1 kcal mol⁻¹ rad⁻²  →  kJ mol⁻¹ rad⁻²
_KCAL_RAD2_TO_KJ_RAD2 = _KCAL_TO_KJ


# Quantity parsing 

def _parse_quantity(text):
    """
    Parse a SMIRNOFF unit string such as
      '1.520 * angstrom'
      '531.0 * angstrom**-2 * mole**-1 * kilocalorie'
      '116.5 * degree'
    Returns (float_value, unit_str_lower).
    """
    text = text.strip()
    parts = text.split('*', 1)
    value = float(parts[0].strip())
    unit  = parts[1].strip().lower() if len(parts) > 1 else ''
    return value, unit


def _to_nm(value, unit):
    if 'angstrom' in unit:
        return value * _ANG_TO_NM
    if 'nanometer' in unit or unit == 'nm':
        return value
    raise ValueError(f'Unknown length unit: {unit!r}')


def _to_degree(value, unit):
    if 'degree' in unit:
        return value
    if 'radian' in unit:
        return math.degrees(value)
    raise ValueError(f'Unknown angle unit: {unit!r}')


def _bond_k_to_kj_nm2(value, unit):
    """kcal mol⁻¹ Å⁻²  →  kJ mol⁻¹ nm⁻²"""
    if 'kilocalorie' in unit or 'kcal' in unit:
        if 'angstrom' in unit:
            return value * _KCAL_ANG2_TO_KJ_NM2
        if 'nanometer' in unit or 'nm' in unit:
            return value * _KCAL_TO_KJ
    if 'kilojoule' in unit or 'kj' in unit:
        if 'angstrom' in unit:
            return value / (_ANG_TO_NM ** 2)
        return value
    raise ValueError(f'Unknown bond k unit: {unit!r}')


def _angle_k_to_kj_rad2(value, unit):
    """kcal mol⁻¹ rad⁻²  →  kJ mol⁻¹ rad⁻²"""
    if 'kilocalorie' in unit or 'kcal' in unit:
        return value * _KCAL_RAD2_TO_KJ_RAD2
    if 'kilojoule' in unit or 'kj' in unit:
        return value
    raise ValueError(f'Unknown angle k unit: {unit!r}')


def _torsion_k_to_kj(value, unit):
    """kcal mol⁻¹  →  kJ mol⁻¹"""
    if 'kilocalorie' in unit or 'kcal' in unit:
        return value * _KCAL_TO_KJ
    if 'kilojoule' in unit or 'kj' in unit:
        return value
    raise ValueError(f'Unknown torsion k unit: {unit!r}')


def _epsilon_to_kj(value, unit):
    """kcal mol⁻¹  →  kJ mol⁻¹"""
    return _torsion_k_to_kj(value, unit)


def _sigma_to_nm(value, unit):
    return _to_nm(value, unit)


# .offxml fetching 

# Commit-pinned raw URLs for MIT-licensed force field files.
_OFFXML_URLS = {
    'openff-2.0.0': (
        'https://raw.githubusercontent.com/openforcefield/openff-forcefields/'
        'main/openforcefields/offxml/openff-2.0.0.offxml'
    ),
    'openff-2.1.0': (
        'https://raw.githubusercontent.com/openforcefield/openff-forcefields/'
        'main/openforcefields/offxml/openff-2.1.0.offxml'
    ),
    'openff-2.2.0': (
        'https://raw.githubusercontent.com/openforcefield/openff-forcefields/'
        'main/openforcefields/offxml/openff-2.2.0.offxml'
    ),
}


def fetch_offxml(ff_name):
    """
    Download and return the parsed XML root for the named force field.
    Raises RuntimeError with a helpful message on network failure.
    """
    if ff_name not in _OFFXML_URLS:
        raise ValueError(
            f'SmirnoffTyper: unknown force field {ff_name!r}. '
            f'Available: {sorted(_OFFXML_URLS)}'
        )
    url = _OFFXML_URLS[ff_name]
    try:
        with urllib.request.urlopen(url, timeout=30) as resp:
            xml_bytes = resp.read()
    except Exception as exc:
        raise RuntimeError(
            f'SmirnoffTyper: could not download {ff_name} from\n  {url}\n'
            f'Error: {exc}\n'
            'Check your internet connection, or supply the .offxml path via '
            'SmirnoffTyper(offxml_path=...).'
        ) from exc
    return ET.fromstring(xml_bytes.decode('utf-8'))


def load_offxml_file(path):
    """Parse a local .offxml file and return the XML root element."""
    return ET.parse(path).getroot()


# SMIRKS matching helpers 

def _smirks_to_query(smirks):
    """
    Convert a SMIRKS pattern string to an RDKit query mol.
    SMIRKS is a superset of SMARTS so MolFromSmarts works directly.
    Returns None if the pattern cannot be parsed.
    """
    from rdkit import Chem
    return Chem.MolFromSmarts(smirks)


def _tagged_indices(smirks, match_tuple, query):
    """
    Given an RDKit substructure match and the query mol, return a dict
    {tag_number: mol_atom_index} for every atom map number > 0.
    """
    tag_to_idx = {}
    for q_idx, mol_idx in enumerate(match_tuple):
        tag = query.GetAtomWithIdx(q_idx).GetAtomMapNum()
        if tag > 0:
            tag_to_idx[tag] = mol_idx
    return tag_to_idx


# Molecule builder 

def build_rdkit_mol(connectivity_matrix, elem_ids):
    """
    Build and sanitize an RDKit Mol from VeloxChem connectivity data.

    :param connectivity_matrix: 2-D numpy array; integer bond orders (0/1/2/3)
    :param elem_ids:            array-like of atomic numbers (int)
    :returns:                   sanitized RDKit Mol
    """
    from rdkit import Chem
    n  = len(elem_ids)
    rw = Chem.RWMol()

    for z in elem_ids:
        rw.AddAtom(Chem.Atom(int(z)))

    for i in range(n):
        for j in range(i + 1, n):
            bo = int(round(connectivity_matrix[i, j]))
            if bo == 0:
                continue
            bt = {
                1: Chem.BondType.SINGLE,
                2: Chem.BondType.DOUBLE,
                3: Chem.BondType.TRIPLE,
            }.get(bo, Chem.BondType.SINGLE)
            rw.AddBond(i, j, bt)

    mol = rw.GetMol()
    Chem.SanitizeMol(mol)   # perceives aromaticity, ring info, valence
    return mol


# Core typer 

class SmirnoffTyper:
    """
    Applies a SMIRNOFF (.offxml) force field to a molecule using only RDKit.

    Usage::

        # Download from GitHub (MIT licensed):
        typer = SmirnoffTyper('openff-2.0.0')

        # Or use a local file:
        typer = SmirnoffTyper(offxml_path='/path/to/openff-2.0.0.offxml')

        rdmol  = build_rdkit_mol(connectivity_matrix, elem_ids)
        params = typer.assign_all(rdmol, bond_indices,
                                   angle_indices, dihedral_indices)

    params keys:
        'atoms'     {i: {'sigma': nm, 'epsilon': kJ/mol}}
        'bonds'     {(i,j): {'length': nm, 'k': kJ/mol/nm²}}
        'angles'    {(i,j,k): {'angle': deg, 'k': kJ/mol/rad²}}
        'propers'   {(i,j,k,l): [{'k', 'phase', 'periodicity'}]}
        'impropers' {(i,j,k,l): [{'k', 'phase', 'periodicity'}]}
    """

    def __init__(self, ff_name=None, offxml_path=None):
        if offxml_path is not None:
            root = load_offxml_file(offxml_path)
        elif ff_name is not None:
            root = fetch_offxml(ff_name)
        else:
            raise ValueError('Provide ff_name or offxml_path.')
        self._root = root
        # Pre-compile all patterns once at construction time
        self._bonds     = self._parse_bonds()
        self._angles    = self._parse_angles()
        self._propers   = self._parse_propers()
        self._impropers = self._parse_impropers()
        self._vdw       = self._parse_vdw()

    # Parsers 

    def _parse_bonds(self):
        """
        Returns list of (smirks, query, params_dict) in file order.
        params_dict: {'length': nm, 'k': kJ/mol/nm², 'id': str}
        """
        out = []
        section = self._root.find('Bonds')
        if section is None:
            return out
        for el in section.findall('Bond'):
            smirks = el.get('smirks', '')
            lv, lu = _parse_quantity(el.get('length', '0 * angstrom'))
            kv, ku = _parse_quantity(el.get('k', '0 * angstrom**-2 * mole**-1 * kilocalorie'))
            q = _smirks_to_query(smirks)
            if q is None:
                continue
            out.append((smirks, q, {
                'length': _to_nm(lv, lu),
                'k':      _bond_k_to_kj_nm2(kv, ku),
                'id':     el.get('id', ''),
            }))
        return out

    def _parse_angles(self):
        """
        Returns list of (smirks, query, params_dict) in file order.
        params_dict: {'angle': deg, 'k': kJ/mol/rad², 'id': str}
        """
        out = []
        section = self._root.find('Angles')
        if section is None:
            return out
        for el in section.findall('Angle'):
            smirks = el.get('smirks', '')
            av, au = _parse_quantity(el.get('angle', '109.5 * degree'))
            kv, ku = _parse_quantity(el.get('k', '0 * mole**-1 * radian**-2 * kilocalorie'))
            q = _smirks_to_query(smirks)
            if q is None:
                continue
            out.append((smirks, q, {
                'angle': _to_degree(av, au),
                'k':     _angle_k_to_kj_rad2(kv, ku),
                'id':    el.get('id', ''),
            }))
        return out

    def _parse_propers(self):
        """
        Returns list of (smirks, query, terms_list) in file order.
        terms_list: [{'k': kJ/mol, 'phase': deg, 'periodicity': int}]
        idivf is folded into k during parsing.
        """
        out = []
        section = self._root.find('ProperTorsions')
        if section is None:
            return out

        try:
            default_idivf = float(section.get('default_idivf', '1'))
        except ValueError:
            default_idivf = 1.0

        for el in section.findall('Proper'):
            smirks = el.get('smirks', '')
            terms = []
            n = 1
            while el.get(f'k{n}') is not None:
                kv, ku      = _parse_quantity(el.get(f'k{n}'))
                pv, pu      = _parse_quantity(el.get(f'phase{n}', '0.0 * degree'))
                periodicity = int(el.get(f'periodicity{n}', '1'))
                try:
                    idivf = float(el.get(f'idivf{n}', str(default_idivf)))
                except ValueError:
                    idivf = 1.0   
                terms.append({
                    'k':           _torsion_k_to_kj(kv, ku) / idivf,
                    'phase':       _to_degree(pv, pu),
                    'periodicity': periodicity,
                })
                n += 1
            if not terms:
                continue
            q = _smirks_to_query(smirks)
            if q is None:
                continue
            out.append((smirks, q, terms))
        return out

    def _parse_impropers(self):
        """
        Returns list of (smirks, query, terms_list) in file order.
        SMIRNOFF convention: tag :2 is the central atom.
        idivf defaults to 3 for impropers (trefoil), folded into k.
        """
        out = []
        section = self._root.find('ImproperTorsions')
        if section is None:
            return out
        try:
            default_idivf = float(section.get('default_idivf', '3'))
        except ValueError:
            default_idivf = 3.0

        for el in section.findall('Improper'):
            smirks = el.get('smirks', '')
            terms = []
            n = 1
            while el.get(f'k{n}') is not None:
                kv, ku      = _parse_quantity(el.get(f'k{n}'))
                pv, pu      = _parse_quantity(el.get(f'phase{n}', '180.0 * degree'))
                periodicity = int(el.get(f'periodicity{n}', '2'))
                try:
                    idivf = float(el.get(f'idivf{n}', str(default_idivf)))
                except ValueError:
                    idivf = default_idivf
                terms.append({
                    'k':           _torsion_k_to_kj(kv, ku) / idivf,
                    'phase':       _to_degree(pv, pu),
                    'periodicity': periodicity,
                })
                n += 1
            if not terms:
                continue
            q = _smirks_to_query(smirks)
            if q is None:
                continue
            out.append((smirks, q, terms))
        return out

    def _parse_vdw(self):
        """
        Returns list of (smirks, query, params_dict) in file order.
        params_dict: {'sigma': nm, 'epsilon': kJ/mol, 'id': str}
        Handles both sigma and rmin_half attributes.
        """
        out = []
        section = self._root.find('vdW')
        if section is None:
            return out
        for el in section.findall('Atom'):
            smirks  = el.get('smirks', '')
            eps_str = el.get('epsilon')
            sig_str = el.get('sigma')
            rmin_str = el.get('rmin_half')
            if not (smirks and eps_str):
                continue
            ev, eu = _parse_quantity(eps_str)
            if sig_str:
                sv, su = _parse_quantity(sig_str)
                sigma = _sigma_to_nm(sv, su)
            elif rmin_str:
                rv, ru = _parse_quantity(rmin_str)
                # sigma = rmin_half * 2 / 2^(1/6)
                sigma = _to_nm(rv, ru) * 2.0 / (2.0 ** (1.0 / 6.0))
            else:
                continue
            q = _smirks_to_query(smirks)
            if q is None:
                continue
            out.append((smirks, q, {
                'sigma':   sigma,
                'epsilon': _epsilon_to_kj(ev, eu),
                'id':      el.get('id', ''),
            }))
        return out

    # Last-match-wins assignment 

    def _last_match(self, entry_list, mol, key_atoms):
        """
        Walk entry_list in order; return the params of the LAST entry whose
        SMIRKS matches mol such that the tagged atoms (:1,:2,...) correspond
        to key_atoms in order (or reversed, for symmetric interactions).

        Returns None if no pattern matches.
        """
        result = None
        n_tags = len(key_atoms)
        for _smirks, query, params in entry_list:
            for match in mol.GetSubstructMatches(query, useChirality=False):
                tagged = _tagged_indices(_smirks, match, query)
                ordered = tuple(tagged.get(t) for t in range(1, n_tags + 1))
                if ordered == key_atoms or ordered == key_atoms[::-1]:
                    result = params
                    break   # found a match for this pattern; continue to next
        return result

    # Public assignment methods 

    def assign_atoms(self, mol):
        """
        Assign vdW parameters to each atom.
        Returns {atom_idx: {'sigma': nm, 'epsilon': kJ/mol}}.
        """
        n_atoms = mol.GetNumAtoms()
        result  = {}
        for i in range(n_atoms):
            params = self._last_match(self._vdw, mol, (i,))
            result[i] = {
                'sigma':   params['sigma']   if params else 0.0,
                'epsilon': params['epsilon'] if params else 0.0,
            }
        return result

    def assign_bonds(self, mol, bond_indices):
        """
        Assign bond parameters.
        Returns {(i,j): {'length': nm, 'k': kJ/mol/nm², 'id': str}}.
        """
        result = {}
        for ij in bond_indices:
            params = self._last_match(self._bonds, mol, tuple(ij))
            if params is not None:
                result[tuple(ij)] = dict(params)
        return result

    def assign_angles(self, mol, angle_indices):
        """
        Assign angle parameters.
        Returns {(i,j,k): {'angle': deg, 'k': kJ/mol/rad², 'id': str}}.
        """
        result = {}
        for ijk in angle_indices:
            params = self._last_match(self._angles, mol, tuple(ijk))
            if params is not None:
                result[tuple(ijk)] = dict(params)
        return result

    def assign_propers(self, mol, dihedral_indices):
        """
        Assign proper torsion parameters.
        A single dihedral may yield multiple terms (multi-periodicity).
        Returns {(i,j,k,l): [{'k', 'phase', 'periodicity'}]}.
        """
        result = {}
        for ijkl in dihedral_indices:
            terms = self._last_match(self._propers, mol, tuple(ijkl))
            if terms is not None:
                result[tuple(ijkl)] = list(terms)
        return result

    def assign_impropers(self, mol):
        """
        Assign improper torsion parameters (SMIRNOFF trefoil convention).
        Tag :2 is the centre; the three peripheral atoms are sorted to give a
        canonical ordering, avoiding triple-counting.

        Returns {(a, centre, b, c): [{'k', 'phase', 'periodicity'}]}
        where a < b < c are the peripheral atoms.
        """
        result = {}
        for _smirks, query, terms in self._impropers:
            for match in mol.GetSubstructMatches(query, useChirality=False):
                tagged  = _tagged_indices(_smirks, match, query)
                centre  = tagged.get(2)
                periph  = [tagged.get(t) for t in (1, 3, 4)]
                if None in (centre, *periph):
                    continue
                periph_sorted = tuple(sorted(periph))
                key = (periph_sorted[0], centre,
                       periph_sorted[1], periph_sorted[2])
                # last match wins
                result[key] = list(terms)
        return result

    def assign_all(self, mol, bond_indices, angle_indices, dihedral_indices):
        """
        Convenience: run all assignment methods and return a combined dict
        with keys 'atoms', 'bonds', 'angles', 'propers', 'impropers'.
        """
        return {
            'atoms':     self.assign_atoms(mol),
            'bonds':     self.assign_bonds(mol, bond_indices),
            'angles':    self.assign_angles(mol, angle_indices),
            'propers':   self.assign_propers(mol, dihedral_indices),
            'impropers': self.assign_impropers(mol),
        }


#  Module-level convenience function

def parametrize_molecule(
    connectivity_matrix,
    elem_ids,
    bond_indices,
    angle_indices,
    dihedral_indices,
    ff_name='openff-2.0.0',
    offxml_path=None,
):
    """
    One-shot parametrization: build RDKit mol, load force field, assign all.

    :param connectivity_matrix: numpy (n×n) integer bond-order matrix
    :param elem_ids:            array-like of atomic numbers
    :param bond_indices:        list of (i,j) tuples
    :param angle_indices:       list of (i,j,k) tuples
    :param dihedral_indices:    list of (i,j,k,l) tuples
    :param ff_name:             'openff-2.0.0', 'openff-2.1.0', 'openff-2.2.0'
    :param offxml_path:         path to a local .offxml file (overrides ff_name)
    :returns:                   dict with keys 'atoms','bonds','angles',
                                'propers','impropers'
    """
    mol   = build_rdkit_mol(connectivity_matrix, elem_ids)
    typer = SmirnoffTyper(ff_name=ff_name, offxml_path=offxml_path)
    return typer.assign_all(mol, bond_indices, angle_indices, dihedral_indices)

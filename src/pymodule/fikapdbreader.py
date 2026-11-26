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
from mpi4py import MPI
from pathlib import Path
import numpy as np

from .veloxchemlib import (mpi_master, bohr_in_angstrom,
                           chemical_element_identifier, chemical_element_mass,
                           FikaClassicalSystem)
from .molecule import Molecule
from .fikaforcefield import FikaForceFieldReader


class FikaPdbResidue:
    """
    A residue of a PDB file.

    Instance variables
        - name: The residue name as in the file, e.g. 'MOL' or 'HOH'.
        - chain: The chain identifier (one character, may be blank).
        - number: The residue sequence number.
        - insertion_code: The insertion code (one character, may be blank).
        - atom_names: The atom names, in file order.
        - elements: The atomic numbers, in file order.
        - coordinates: The coordinates (bohr), shape (atoms, 3).
    """

    def __init__(self, name, chain, number, insertion_code):
        """
        Initializes an empty residue.
        """

        self.name = name
        self.chain = chain
        self.number = number
        self.insertion_code = insertion_code
        self.atom_names = []
        self.elements = []
        self.coordinates = []

    def key(self):
        """
        Gets the identity of the residue: name, chain, number and insertion
        code.
        """

        return (self.name, self.chain, self.number, self.insertion_code)

    def centre_of_mass(self):
        """
        Gets the centre of mass (bohr) with VeloxChem's atomic masses.
        """

        masses = np.array([chemical_element_mass(z) for z in self.elements])
        return masses @ self.coordinates / np.sum(masses)


class FikaPdbReader:
    """
    Reads the ATOM and HETATM records of a single-model PDB file into
    residues, and builds from them the QM solute (a VeloxChem molecule) and
    the classical (MM) system of the solvent.

    Fixed columns of PDB v3.3 are read; coordinates are in angstrom and kept
    in bohr. Consecutive atoms with equal residue name, chain, number and
    insertion code form a residue. Reading stops at END; other records are
    ignored. The element symbol (columns 77-78) is required; alternate
    locations, MODEL records and residues split by another residue are
    rejected. Errors raise ValueError naming the source and line.

    :param comm:
        The MPI communicator: the master rank reads files and broadcasts
        their text.

    Instance variables
        - comm: The MPI communicator.
        - rank: The rank of this process.
    """

    def __init__(self, comm=None):
        """
        Initializes the PDB reader.
        """

        if comm is None:
            comm = MPI.COMM_WORLD

        self.comm = comm
        self.rank = comm.Get_rank()

    def parse(self, text, source='<string>'):
        """
        Parses a PDB document.

        :param text:
            The text of the PDB document.
        :param source:
            The name of the text in error messages.

        :return:
            The list of FikaPdbResidue objects, in file order.
        """

        to_bohr = 1.0 / bohr_in_angstrom()
        residues = []
        finished = set()

        for number, line in enumerate(text.splitlines(), start=1):

            def fail(reason):
                raise ValueError(f'{source}: line {number}: {reason}')

            record = line[0:6].strip()
            if record == 'END':
                break
            if record == 'MODEL':
                fail('MODEL records are not supported (single model only)')
            if record not in ('ATOM', 'HETATM'):
                continue

            if line[16:17] not in ('', ' '):
                fail(f"alternate location '{line[16]}' is not supported")

            symbol = line[76:78].strip()
            if not symbol:
                fail('missing element symbol (columns 77-78)')
            element = chemical_element_identifier(symbol.upper())
            if element <= 0:
                fail(f"unknown element symbol '{symbol}'")

            position = []
            for first in (31, 39, 47):
                field = line[first - 1:first + 7].strip()
                try:
                    value = float(field)
                except ValueError:
                    value = None
                if value is None or not np.isfinite(value):
                    fail(f"invalid coordinate '{field}' in columns " +
                         f'{first}-{first + 7}')
                position.append(value * to_bohr)

            name = line[17:20].strip()
            if not name:
                fail('missing residue name (columns 18-20)')
            field = line[22:26].strip()
            try:
                sequence = int(field)
            except ValueError:
                fail(f"invalid residue number '{field}'")
            chain = line[21:22] or ' '
            insertion = line[26:27] or ' '
            key = (name, chain, sequence, insertion)

            if not residues or residues[-1].key() != key:
                if residues:
                    finished.add(residues[-1].key())
                if key in finished:
                    fail(f"residue {name} {sequence} of chain '{chain}' is " +
                         'split by another residue')
                residues.append(FikaPdbResidue(*key))

            residue = residues[-1]
            residue.atom_names.append(line[12:16].strip())
            residue.elements.append(element)
            residue.coordinates.append(position)

        if not residues:
            raise ValueError(f'{source}: no ATOM or HETATM records')
        for residue in residues:
            residue.coordinates = np.array(residue.coordinates)
        return residues

    def read(self, path):
        """
        Reads a PDB file.

        :param path:
            The path of the file.

        :return:
            The list of FikaPdbResidue objects, in file order.
        """

        path = Path(path)
        text = None
        if self.rank == mpi_master():
            try:
                text = path.read_text()
            except OSError:
                text = None
        text = self.comm.bcast(text, root=mpi_master())
        if text is None:
            raise FileNotFoundError(f'cannot open PDB file {path}')
        return self.parse(text, str(path))

    @staticmethod
    def _solutes(residues, solute):

        return [r for r in residues if r.name.lower() == solute.lower()]

    def get_solute(self, residues, solute='MOL'):
        """
        Gets the solute as a VeloxChem molecule.

        :param residues:
            The residues of a PDB file.
        :param solute:
            The residue name of the solute (case insensitive).

        :return:
            The molecule of the single residue named `solute`, atoms in file
            order.
        """

        solutes = self._solutes(residues, solute)
        if len(solutes) != 1:
            raise ValueError(f'expected one {solute} residue, found ' +
                             f'{len(solutes)}')
        return Molecule(solutes[0].elements, solutes[0].coordinates, 'bohr')

    def get_classical_system(self,
                             residues,
                             solvent,
                             radius=None,
                             shells=None,
                             solute='MOL',
                             unit='angstrom'):
        """
        Builds the classical system of the solvent residues (all residues but
        the solute) with force fields from FikaForceFieldReader, in one of
        three modes:

        - all solvent residues: solvent = {residue name: force-field label},
          each residue in the region of its force field's kind;
        - radius: as above, only residues whose centre of mass is within
          `radius` of the solute's;
        - shells = (inner, outer): solvent = {'polarizable': {...},
          'nonpolarizable': {...}}, residues within `inner` of the solute's
          centre of mass in the polarizable region (polarizable force
          fields), those between `inner` and `outer` in the nonpolarizable
          region (nonpolarizable force fields); an empty mapping skips its
          region.

        Residues are added in file order; their indices count the residues of
        each name within their region from 0. The radius and shell modes need
        exactly one solute residue, the first mode at most one.

        :param residues:
            The residues of a PDB file.
        :param solvent:
            The force fields of the solvent residues (see above).
        :param radius:
            The radius of the radius mode.
        :param shells:
            The (inner, outer) radii of the shell mode.
        :param solute:
            The residue name of the solute (case insensitive).
        :param unit:
            The unit of the radii, 'angstrom' or 'bohr'.

        :return:
            The FikaClassicalSystem and the identities of its residues:
            {'polarizable': [...], 'nonpolarizable': [...]} with (name,
            chain, number, insertion code) of each residue in region order.
        """

        if radius is not None and shells is not None:
            raise ValueError('give either a radius or shells, not both')
        if unit.lower() not in ('angstrom', 'bohr'):
            raise ValueError(f"unknown unit '{unit}'")
        scale = (1.0 / bohr_in_angstrom()) if unit.lower() == 'angstrom' else 1.0

        def radius_in_bohr(value, what):
            if not (value > 0.0 and np.isfinite(value)):
                raise ValueError(f'{what} must be positive and finite')
            return value * scale

        reader = FikaForceFieldReader(self.comm)

        def load(mapping):
            fields = {}
            for name, label in mapping.items():
                if name.lower() in fields:
                    raise ValueError(f'residue {name} is listed twice')
                fields[name.lower()] = reader.load(label, name)
            return fields

        solutes = self._solutes(residues, solute)
        if len(solutes) > 1:
            raise ValueError(f'more than one {solute} residue')

        centre = None
        if shells is None:
            solvent_fields = load(solvent)
            cutoff = (np.inf if radius is None else radius_in_bohr(
                radius, 'the radius'))

            def select(distance):
                return solvent_fields if distance <= cutoff else None
        else:
            inner = radius_in_bohr(shells[0], 'the inner radius')
            outer = radius_in_bohr(shells[1], 'the outer radius')
            if inner > outer:
                raise ValueError('the inner radius exceeds the outer radius')
            unknown = set(solvent) - {'polarizable', 'nonpolarizable'}
            if unknown:
                raise ValueError('shell force fields are given as ' +
                                 "{'polarizable': ..., 'nonpolarizable': ...}")
            regions = {
                kind: load(solvent.get(kind, {}))
                for kind in ('polarizable', 'nonpolarizable')
            }
            if not regions['polarizable'] and not regions['nonpolarizable']:
                raise ValueError('no solvent force fields')
            for kind, fields in regions.items():
                for name, field in fields.items():
                    if field.polarizable != (kind == 'polarizable'):
                        raise ValueError(
                            f'force field {field.label} of residue {name} is ' +
                            ('' if field.polarizable else 'not ') +
                            f'polarizable, but assigned to the {kind} region')

            def select(distance):
                if distance <= inner:
                    return regions['polarizable'] or None
                if distance <= outer:
                    return regions['nonpolarizable'] or None
                return None

        if radius is not None or shells is not None:
            if not solutes:
                raise ValueError(f'no {solute} residue')
            centre = solutes[0].centre_of_mass()

        system = FikaClassicalSystem()
        identities = {'polarizable': [], 'nonpolarizable': []}
        added = set()

        # Consecutive residues of one name and force field go in one call.
        run = None

        def flush():
            if run is not None:
                field, members = run
                system.add_residues(members[0].name, field.label,
                                    members[0].elements,
                                    np.array([m.coordinates for m in members]),
                                    field.polarizable)

        for residue in residues:
            if residue.name.lower() == solute.lower():
                continue
            distance = (np.linalg.norm(residue.centre_of_mass() - centre)
                        if centre is not None else 0.0)
            selected = select(distance)
            if selected is None:
                continue
            describe = (f"{residue.name} {residue.number} of chain " +
                        f"'{residue.chain}'")
            field = selected.get(residue.name.lower())
            if field is None:
                raise ValueError(f'no force field for residue {describe}')
            if len(residue.elements) != field.atom_count:
                raise ValueError(f'residue {describe} has ' +
                                 f'{len(residue.elements)} atoms, force field ' +
                                 f'{field.label} {field.atom_count}')
            if (field.label, field.residue) not in added:
                system.add_force_field(field)
                added.add((field.label, field.residue))
            identities['polarizable' if field.polarizable else
                       'nonpolarizable'].append(residue.key())

            if (run is not None and run[0] is field and
                    run[1][0].elements == residue.elements):
                run[1].append(residue)
            else:
                flush()
                run = (field, [residue])
        flush()

        return system, identities

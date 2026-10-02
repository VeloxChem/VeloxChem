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
from os import environ
from pathlib import Path

from .veloxchemlib import mpi_master, FikaForceField
from .environment import get_data_path


class FikaForceFieldReader:
    """
    Reads force fields in the fika .fff format into FikaForceField objects,
    one per residue. The format is line oriented; blank lines and text after
    '#' are ignored, atom indices are 1-based and parameters in atomic units:

        fika-force-field 1 <label> <polarizable|nonpolarizable>
        residue <name> <atoms>
        charges                              (required; atoms 1..n in order)
        <atom> <q>
        polarizabilities isotropic <n>       (required in a polarizable force
        <atom> <alpha>                        field, absent otherwise; atoms
        end                                   in increasing order)

    Dipoles, quadrupoles, anisotropic polarizabilities and reference
    geometries are part of the format but not supported yet; they are
    rejected. Force fields are looked up by label (case insensitive) as
    <label>.fff in the directories of VLXFORCEFIELDPATH (separated by ':'),
    then in the force-field library of VeloxChem.

    Errors in a force field raise ValueError naming the source and line; a
    force field that is not found raises FileNotFoundError.

    :param comm:
        The MPI communicator: the master rank reads the files and broadcasts
        their text.

    Instance variables
        - comm: The MPI communicator.
        - rank: The rank of this process.
    """

    _unsupported = ('dipoles', 'quadrupoles', 'geometry')

    def __init__(self, comm=None):
        """
        Initializes the force-field reader.
        """

        if comm is None:
            comm = MPI.COMM_WORLD

        self.comm = comm
        self.rank = comm.Get_rank()

    @staticmethod
    def get_library_path():
        """
        Gets the location of the force-field library of VeloxChem.

        :return:
            The directory of the shipped .fff files.
        """

        return get_data_path() / 'force_fields'

    def get_search_path(self):
        """
        Gets the directories searched for force fields, in order.

        :return:
            The list of directories.
        """

        directories = []
        for entry in environ.get('VLXFORCEFIELDPATH', '').split(':'):
            if entry:
                directories.append(Path(entry))
        directories.append(self.get_library_path())
        return directories

    def parse(self, text, source='<string>'):
        """
        Parses a force field.

        :param text:
            The text of a .fff document.
        :param source:
            The name of the text in error messages.

        :return:
            The list of FikaForceField objects, one per residue.
        """

        lines = _Lines(text, source)

        fields = lines.expect('the fika-force-field header')
        if fields[0] != 'fika-force-field':
            lines.fail("expected 'fika-force-field' header")
        lines.check(fields, 4, 'header')
        if lines.count(fields[1], 'format version') != 1:
            lines.fail(f'unsupported format version {fields[1]}')
        label = fields[2]
        flag = fields[3].lower()
        if flag not in ('polarizable', 'nonpolarizable'):
            lines.fail("expected 'polarizable' or 'nonpolarizable', " +
                       f"found '{fields[3]}'")
        polarizable = (flag == 'polarizable')

        force_fields = []
        while True:
            fields = lines.next()
            if fields is None:
                break
            if fields[0] != 'residue':
                lines.fail(f"expected 'residue', found '{fields[0]}'")
            start = lines.number
            residue, charges, alphas = self._read_residue(lines, fields)

            def fail_at_start(reason):
                raise ValueError(
                    f'{source}: line {start}: residue {residue}: {reason}')

            if (alphas is not None) != polarizable:
                fail_at_start('a polarizable force field needs polarizabilities'
                              if polarizable else
                              'a nonpolarizable force field has no ' +
                              'polarizabilities')
            if any(f.residue == residue.lower() for f in force_fields):
                fail_at_start('repeated residue')
            try:
                force_fields.append(
                    FikaForceField(label, residue, charges, alphas))
            except ValueError as error:
                fail_at_start(str(error))

        if not force_fields:
            raise ValueError(f'{source}: force field {label} has no residues')
        return force_fields

    def read(self, path):
        """
        Reads a force field from a .fff file.

        :param path:
            The path of the file.

        :return:
            The list of FikaForceField objects, one per residue.
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
            raise FileNotFoundError(f'cannot open force-field file {path}')
        return self.parse(text, str(path))

    def load(self, label, residue=None):
        """
        Loads a force field by label from the search path.

        :param label:
            The force-field label (case insensitive), e.g. 'TIP3P'.
        :param residue:
            The residue name (case insensitive), or None for all residues.

        :return:
            The list of FikaForceField objects, or the FikaForceField of the
            residue.
        """

        name = f'{label.lower()}.fff'
        path = None
        if self.rank == mpi_master():
            for directory in self.get_search_path():
                if (directory / name).is_file():
                    path = str(directory / name)
                    break
        path = self.comm.bcast(path, root=mpi_master())
        if path is None:
            raise FileNotFoundError(
                f'force field {label} not found (searched ' + ', '.join(
                    str(d) for d in self.get_search_path()) + ')')

        force_fields = self.read(path)
        if residue is None:
            return force_fields
        for force_field in force_fields:
            if force_field.residue == residue.lower():
                return force_field
        raise ValueError(f'force field {label} has no residue {residue}')

    def _read_residue(self, lines, header):
        """
        Reads one residue block after its header line.

        :return:
            The residue name, its charges and its polarizabilities
            ({0-based atom: alpha}, or None).
        """

        lines.check(header, 3, 'residue header')
        residue = header[1]
        atoms = lines.count(header[2], 'atom count')
        if atoms == 0:
            lines.fail('a residue needs at least one atom')

        charges = None
        alphas = None
        while True:
            fields = lines.expect(f"'end' of residue {residue}")
            section = fields[0]
            if section == 'end':
                lines.check(fields, 1, 'end')
                break
            elif section == 'charges':
                if charges is not None:
                    lines.fail("repeated section 'charges'")
                lines.check(fields, 1, 'charges header')
                charges = []
                for a in range(atoms):
                    fields = lines.expect('charge line')
                    lines.check(fields, 2, 'charge line')
                    if lines.atom(fields[0], atoms) != a:
                        lines.fail(f'charges must list atoms 1..{atoms} in order')
                    charges.append(lines.real(fields[1]))
            elif section == 'polarizabilities':
                if alphas is not None:
                    lines.fail("repeated section 'polarizabilities'")
                lines.check(fields, 3, 'polarizabilities header')
                kind = fields[1].lower()
                if kind == 'anisotropic':
                    lines.fail('anisotropic polarizabilities are not supported yet')
                if kind != 'isotropic':
                    lines.fail(f"unknown polarizability kind '{fields[1]}'")
                count = lines.count(fields[2], 'polarizability count')
                alphas = {}
                previous = -1
                for k in range(count):
                    fields = lines.expect('polarizability line')
                    lines.check(fields, 2, 'polarizability line')
                    atom = lines.atom(fields[0], atoms)
                    if atom <= previous:
                        lines.fail('polarizabilities must list atoms in ' +
                                   'increasing order')
                    previous = atom
                    alphas[atom] = lines.real(fields[1])
            elif section in self._unsupported:
                lines.fail(f"section '{section}' is not supported yet")
            else:
                lines.fail(f"unknown section '{section}'")

        if charges is None:
            lines.fail(f'residue {residue} has no charges')
        return residue, charges, alphas


class _Lines:
    """
    Significant lines of a .fff document (comments after '#' removed, blank
    lines skipped) with their numbers.
    """

    def __init__(self, text, source):

        self._lines = text.splitlines()
        self._source = source
        self.number = 0

    def next(self):

        while self.number < len(self._lines):
            line = self._lines[self.number].split('#', 1)[0]
            self.number += 1
            fields = line.split()
            if fields:
                return fields
        return None

    def expect(self, what):

        fields = self.next()
        if fields is None:
            raise ValueError(
                f'{self._source}: unexpected end of file: expected {what}')
        return fields

    def fail(self, reason):

        raise ValueError(f'{self._source}: line {self.number}: {reason}')

    def check(self, fields, expected, what):

        if len(fields) != expected:
            self.fail(f'{what} needs {expected} fields, found {len(fields)}')

    def count(self, field, what):

        if not field.isdigit():
            self.fail(f"invalid {what} '{field}'")
        return int(field)

    def real(self, field):

        try:
            value = float(field)
        except ValueError:
            value = None
        if value is None or value != value or value in (float('inf'),
                                                        float('-inf')):
            self.fail(f"invalid number '{field}'")
        return value

    def atom(self, field, atoms):

        index = self.count(field, 'atom index')
        if index < 1 or index > atoms:
            self.fail(f'atom index {field} outside 1..{atoms}')
        return index - 1

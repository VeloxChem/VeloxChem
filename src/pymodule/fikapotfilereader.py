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
                           chemical_element_identifier, FikaClassicalSystem,
                           FikaForceField, FikaResidue)
from .inputparser import InputParser


class FikaPotfileReader:
    """
    Reads a potential file with an ``@environment`` section into a fika
    classical system.

    The environment lists the atoms as

        element x y z residue resid atom_name

    and may be followed by ``@charges`` and ``@polarizabilities`` sections
    holding one template per residue name. A residue name ending in ``_pe``
    is polarizable, one ending in ``_npe`` is nonpolarizable; without these
    suffixes, the presence of a polarizability template decides. Coordinates
    are in the unit given by the ``units`` entry (angstrom by default) and
    are converted to bohr. Anisotropic polarizabilities are reduced to their
    isotropic part (trace/3), as in the fika embedding of a PDB file.

    :param comm:
        The MPI communicator: the master rank reads the file and broadcasts
        its contents.

    Instance variables
        - comm: The MPI communicator.
        - rank: The rank of this process.
    """

    def __init__(self, comm=None):
        """
        Initializes the potential-file reader.
        """

        if comm is None:
            comm = MPI.COMM_WORLD

        self.comm = comm
        self.rank = comm.Get_rank()

    def read(self, path):
        """
        Reads a potential file.

        :param path:
            The path of the file.

        :return:
            The FikaClassicalSystem of the environment.
        """

        path = Path(path)
        data = None
        if self.rank == mpi_master():
            try:
                data = InputParser(path).input_dict
            except OSError:
                data = None
        data = self.comm.bcast(data, root=mpi_master())
        if data is None:
            raise FileNotFoundError(f'cannot open potential file {path}')
        return self._build(data, str(path))

    def _build(self, data, source):

        environment = data.get('environment')
        if not isinstance(environment, dict) or 'xyz' not in environment:
            raise ValueError(
                f'{source}: missing @environment section with xyz')

        units = str(environment.get('units', 'angstrom')).lower()
        if len(units) >= 3 and units == 'angstrom'[:len(units)]:
            scale = 1.0 / bohr_in_angstrom()
        elif units in ('au', 'bohr'):
            scale = 1.0
        else:
            raise ValueError(f'{source}: unknown units {units}')

        xyz = environment['xyz']
        if isinstance(xyz, str):
            xyz = [xyz]

        residues = {'order': [], 'names': []}
        last_key = None
        for number, line in enumerate(xyz, start=1):
            fields = line.split()
            if len(fields) not in (6, 7):
                raise ValueError(
                    f'{source}: xyz line {number}: needs element, three ' +
                    'coordinates, residue, resid and atom name')
            element = chemical_element_identifier(fields[0].upper())
            if element <= 0:
                raise ValueError(
                    f'{source}: xyz line {number}: unknown element ' +
                    f"'{fields[0]}'")
            try:
                coordinates = [float(x) for x in fields[1:4]]
            except ValueError:
                raise ValueError(
                    f'{source}: xyz line {number}: invalid coordinates')
            if not all(np.isfinite(coordinates)):
                raise ValueError(
                    f'{source}: xyz line {number}: invalid coordinates')
            try:
                resid = int(fields[5])
            except ValueError:
                raise ValueError(
                    f'{source}: xyz line {number}: invalid resid ' +
                    f"'{fields[5]}'")
            name = fields[4]
            key = resid
            if key not in residues:
                residues[key] = {
                    'name': name,
                    'elements': [],
                    'coordinates': [],
                }
                residues['order'].append(key)
                if name not in residues['names']:
                    residues['names'].append(name)
            elif last_key != key:
                raise ValueError(
                    f'{source}: xyz line {number}: residue {name} {resid} ' +
                    'is split')
            residues[key]['elements'].append(element)
            residues[key]['coordinates'].append(
                [x * scale for x in coordinates])
            last_key = key

        for key in residues['order']:
            residues[key]['coordinates'] = np.array(
                residues[key]['coordinates'])

        charges = self._parse_charges(data.get('charges', []), source)
        polarizabilities = self._parse_polarizabilities(
            data.get('polarizabilities', []), source)

        fields = {}
        for name in residues['names']:
            first = next(key for key in residues['order']
                         if residues[key]['name'] == name)
            count = len(residues[first]['elements'])
            q = charges.get(name)
            if q is None:
                q = [0.0] * count
            elif len(q) != count:
                raise ValueError(
                    f'{source}: residue {name} has {count} atoms, but ' +
                    f'{len(q)} charges')
            if name.endswith('_npe'):
                polarizable = False
            elif name.endswith('_pe'):
                polarizable = True
            else:
                polarizable = name in polarizabilities
            alphas = None
            if polarizable:
                pol = polarizabilities.get(name)
                if pol is None:
                    raise ValueError(
                        f'{source}: polarizable residue {name} has no ' +
                        'polarizabilities')
                if len(pol) != count:
                    raise ValueError(
                        f'{source}: residue {name} has {count} atoms, but ' +
                        f'{len(pol)} polarizabilities')
                alphas = {
                    i: (p[0] + p[3] + p[5]) / 3.0
                    for i, p in enumerate(pol)
                }
            fields[name] = (f'potfile_{len(fields)}', polarizable, q, alphas)

        system = FikaClassicalSystem()
        for name, (label, _, q, alphas) in fields.items():
            system.add_force_field(FikaForceField(label, name, q, alphas))

        index_of = {}
        for key in residues['order']:
            residue = residues[key]
            name = residue['name']
            label, polarizable, q, _ = fields[name]
            if len(residue['elements']) != len(q):
                raise ValueError(
                    f'{source}: residue {name} {key} has ' +
                    f'{len(residue["elements"])} atoms, but {len(q)} charges')
            index = index_of.get(name, 0)
            index_of[name] = index + 1
            system.add_residue(
                FikaResidue(name, index, label, residue['elements'],
                            residue['coordinates'], polarizable))

        return system

    def _parse_charges(self, lines, source):

        charges = {}
        for number, line in enumerate(lines, start=1):
            fields = line.split()
            if len(fields) != 3:
                raise ValueError(
                    f'{source}: charge line {number}: needs element, ' +
                    'charge and residue')
            try:
                charge = float(fields[1])
            except ValueError:
                raise ValueError(
                    f'{source}: charge line {number}: invalid charge ' +
                    f"'{fields[1]}'")
            if not np.isfinite(charge):
                raise ValueError(
                    f'{source}: charge line {number}: invalid charge ' +
                    f"'{fields[1]}'")
            charges.setdefault(fields[2], []).append(charge)
        return charges

    def _parse_polarizabilities(self, lines, source):

        polarizabilities = {}
        for number, line in enumerate(lines, start=1):
            fields = line.split()
            if len(fields) != 8:
                raise ValueError(
                    f'{source}: polarizability line {number}: needs ' +
                    'element, six tensor components and residue')
            try:
                values = [float(x) for x in fields[1:7]]
            except ValueError:
                raise ValueError(
                    f'{source}: polarizability line {number}: invalid ' +
                    'polarizability')
            if not all(np.isfinite(values)):
                raise ValueError(
                    f'{source}: polarizability line {number}: invalid ' +
                    'polarizability')
            polarizabilities.setdefault(fields[7], []).append(values)
        return polarizabilities

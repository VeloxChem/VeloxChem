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

from collections import OrderedDict
from pathlib import Path

from .errorhandler import assert_msg_critical
from .molecule import Molecule


class SimdPolarizableEmbeddingDriver:
    """
    Reads a solvated structure and separates it into the solute and the
    molecules of the environment.

    The first molecule of the file is the solute and every other one belongs to
    the environment. Each environment molecule is given a label of its residue
    name and its number within that residue name, `S01-1`, `S01-2` and so on,
    numbered from one and starting again for each kind of molecule: the residue
    name says which species it is and the number says which of them.

    :param ostream:
        The output stream.

    Instance variables
        - solute: The molecule the environment surrounds.
        - solvent: The molecules of the environment, as (label, molecule) pairs.
        - box: The edges of the periodic box in angstrom, or None.
    """

    def __init__(self, ostream=None):
        """
        Initializes the driver.
        """

        self.ostream = ostream

        self.solute = None
        self.solvent = []
        self.box = None

        self.solute_label = None

    def read_pdb(self, pdb_file):
        """
        Reads a solvated structure and separates it.

        :param pdb_file:
            The name of the PDB file.

        :return:
            The driver, so that a caller may read and use it in one expression.
        """

        path = Path(pdb_file)

        assert_msg_critical(
            path.is_file(),
            f'SimdPolarizableEmbeddingDriver: {pdb_file} is not a file')

        residues, box = self._read_residues(path)

        assert_msg_critical(
            len(residues) > 0,
            f'SimdPolarizableEmbeddingDriver: {pdb_file} holds no atoms')

        self.box = box

        # NOTE: **the first residue of the file and not the one which is named
        # like a solute.** What a solvated structure calls its solute is whatever
        # the program which wrote it chose, and this way the reader does not have
        # to know. An ordered dictionary is what keeps that first.

        keys = list(residues.keys())

        self.solute_label = self._label_of(keys[0], 1)
        self.solute = self._molecule_of(residues[keys[0]])

        # NOTE: the number restarts for each residue name, so a mixed environment
        # reads as S01-1 ... S01-n and S02-1 ... S02-m rather than as one run of
        # numbers which says nothing about what each molecule is.

        counters = {}
        self.solvent = []

        for key in keys[1:]:
            name = key[0]
            counters[name] = counters.get(name, 0) + 1
            self.solvent.append(
                (self._label_of(key, counters[name]),
                 self._molecule_of(residues[key])))

        return self

    def number_of_solvent_molecules(self):
        """
        Gets the number of molecules in the environment.

        :return:
            The number of molecules.
        """

        return len(self.solvent)

    def solvent_species(self):
        """
        Gets how many molecules of each residue name the environment holds.

        :return:
            A dictionary of the residue name and the count.
        """

        species = OrderedDict()

        for label, _ in self.solvent:
            name = label.rsplit('-', 1)[0]
            species[name] = species.get(name, 0) + 1

        return species

    @staticmethod
    def _label_of(key, number):
        """
        Builds the label of a molecule from its residue name and its number.

        :param key:
            The residue name and the residue number of the molecule.
        :param number:
            The number of the molecule within its residue name.

        :return:
            The label.
        """

        return f'{key[0]}-{number}'

    @staticmethod
    def _molecule_of(residue):
        """
        Builds a molecule from the atoms of one residue.

        :param residue:
            The labels and the coordinates of the atoms.

        :return:
            The molecule.
        """

        return Molecule(residue['labels'], residue['coordinates'], 'angstrom')

    @staticmethod
    def _read_residues(path):
        """
        Reads the atoms of a PDB file, gathered by residue.

        :param path:
            The path of the PDB file.

        :return:
            The residues, keyed by residue name and number and in the order the
            file gives them, and the edges of the box or None.
        """

        residues = OrderedDict()
        box = None

        with open(path, 'r') as handle:
            for line in handle:
                if line.startswith('CRYST1'):
                    box = tuple(float(line[i:i + 9]) for i in (6, 15, 24))
                    continue

                if not line.startswith(('ATOM', 'HETATM')):
                    continue

                # NOTE: the element as the file gives it, and the first letter of
                # the atom name only where it does not. The residue name is three
                # columns wide and is read as three: read as two, S01 comes back
                # as S0 and every molecule of the environment answers to the same
                # wrong name.

                label = line[76:78].strip()

                if not label:
                    label = line[12:16].strip()[:1]

                key = (line[17:20].strip(), int(line[22:26]))

                if key not in residues:
                    residues[key] = {'labels': [], 'coordinates': []}

                residues[key]['labels'].append(label)
                residues[key]['coordinates'].append(
                    [float(line[i:i + 8]) for i in range(30, 54, 8)])

        return residues, box

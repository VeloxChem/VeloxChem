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

import numpy as np

from .errorhandler import assert_msg_critical
from .veloxchemlib import mpi_master
from .molecule import Molecule
from .oneeints import compute_simd_nuclear_potential_integrals
from .veloxchemlib import bohr_in_angstrom
from .veloxchemlib import PolarizableEmbedding
from .veloxchemlib import PolarizableForceField

# NOTE: how many numbers a multipole of each order is written with, in the order
# the sites keep them: a charge, then x y z, then xx xy xz yy yz zz.
_MULTIPOLE_SIZES = {0: 1, 1: 3, 2: 6}

# NOTE: what each order is called, for the messages. The permanent multipoles of
# an order the environment carries but this driver cannot yet contract are named
# rather than passed over, because passing them over is an energy which is wrong
# by a little and looks like an energy which is right.
_MULTIPOLE_NAMES = {0: 'charges', 1: 'dipoles', 2: 'quadrupoles'}


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
        - force_fields: The force fields of a potential file, keyed by name.
        - embedding: The environment the force fields and the structure make up.
        - shells: How many molecules each shell took and how many were dropped.
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

        self.force_fields = OrderedDict()

        self.embedding = None

        self.shells = None

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

    def read_potential(self, potential_file):
        """
        Reads the force fields of a potential file.

        The force fields are kept as they are read and are not matched to the
        molecules of the environment: a potential file says what a kind of
        molecule is made of and says nothing about which molecules there are.

        :param potential_file:
            The path of the potential file.

        :return:
            The driver.
        """

        self.force_fields = read_potential_file(potential_file)

        return self

    def build_embedding(self):
        """
        Builds the environment of the structure and the force fields.

        Each molecule of the environment is given the force field of its
        species, which is the part of its label before the number: a molecule
        labelled `WAT-17` is described by the force field named `WAT`.

        Which region a molecule goes to follows from its force field, and is
        not a choice: one carrying a polarizability goes to the polarizable
        region and one without goes to the nonpolarizable region. The regions
        hold to that themselves, the nonpolarizable one refusing a force field
        which polarizes, so a rule which said anything else would produce
        refusals nobody asked for.

        The solute is not part of it. The environment is what surrounds the
        molecule the wave function is of.

        :return:
            The environment.
        """

        assert_msg_critical(
            self.solvent or self.solute is not None,
            'build_embedding: no structure has been read')

        assert_msg_critical(
            len(self.force_fields) > 0,
            'build_embedding: no force field has been read')

        embedding = PolarizableEmbedding()

        for label, molecule in self.solvent:
            force_field = self._force_field_for(label, molecule,
                                                self.force_fields,
                                                'build_embedding')

            region = (embedding.polarizable_region()
                      if force_field.is_polarizable() else
                      embedding.nonpolarizable_region())

            region.add_molecule(molecule, force_field)

        self.embedding = embedding

        return embedding

    def build_shells(self,
                     pdb_file,
                     polarizable_radius,
                     nonpolarizable_radius,
                     polarizable_potential,
                     nonpolarizable_potential,
                     units='angstrom'):
        """
        Builds the environment as two shells about the solute.

        The distance of a molecule of the environment is that between its own
        centre of mass and the centre of mass of the solute. A molecule within
        the polarizable radius goes to the polarizable region, one beyond it
        and within the nonpolarizable radius goes to the nonpolarizable region,
        and one beyond that is left out. How many went each way is kept in
        `shells`, so a radius which drops most of the box says so.

        The two shells are given separate potential files because they describe
        the same species differently: the water of the inner shell polarizes
        and the water of the outer shell does not, and the two force fields are
        both named for that species. They do not collide, since each region
        holds its own. The outer file must carry no polarizability, which the
        nonpolarizable region sees to itself.

        NOTE: a sphere about one point suits a solute whose atoms sit near
        their own centre of mass. A long solute is cut by distance from its
        middle instead of by distance from itself, so the molecules in contact
        with its ends fall farther out than those in contact with its middle.

        :param pdb_file:
            The path of the PDB file of the solvated structure.
        :param polarizable_radius:
            The radius the polarizable shell reaches to.
        :param nonpolarizable_radius:
            The radius the nonpolarizable shell reaches to.
        :param polarizable_potential:
            The path of the potential file of the polarizable shell.
        :param nonpolarizable_potential:
            The path of the potential file of the nonpolarizable shell.
        :param units:
            The units the two radii are given in, angstrom or au.

        :return:
            The environment.
        """

        assert_msg_critical(
            str(units).lower() in ('angstrom', 'au', 'bohr'),
            f'build_shells: the radii are given in angstrom or in au, and '
            f'these are in {units}')

        assert_msg_critical(
            polarizable_radius >= 0.0,
            f'build_shells: the polarizable radius {polarizable_radius} is '
            f'negative')

        assert_msg_critical(
            nonpolarizable_radius >= polarizable_radius,
            f'build_shells: the nonpolarizable shell reaches to '
            f'{nonpolarizable_radius} and the polarizable one to '
            f'{polarizable_radius}, which is outside it')

        to_bohr = (1.0 / bohr_in_angstrom()
                   if str(units).lower() == 'angstrom' else 1.0)

        inner = polarizable_radius * to_bohr
        outer = nonpolarizable_radius * to_bohr

        self.read_pdb(pdb_file)

        force_fields = {
            'polarizable': read_potential_file(polarizable_potential),
            'nonpolarizable': read_potential_file(nonpolarizable_potential),
        }

        # NOTE: checked here and not left to the region, which would refuse on
        # the first molecule of the outer shell instead of before any of them.
        # The two files are the caller's choice and this is the one way of
        # getting them the wrong way round.

        polarizing = [
            name for name, force_field in
            force_fields['nonpolarizable'].items()
            if force_field.is_polarizable()
        ]

        assert_msg_critical(
            not polarizing,
            f'build_shells: the outer shell does not polarize and '
            f'{Path(nonpolarizable_potential).name} carries a polarizability '
            f'for {", ".join(polarizing)}')

        centre = np.array(self.solute.center_of_mass_in_bohr())

        embedding = PolarizableEmbedding()

        regions = {
            'polarizable': embedding.polarizable_region(),
            'nonpolarizable': embedding.nonpolarizable_region(),
        }

        shells = OrderedDict([('polarizable', 0), ('nonpolarizable', 0),
                              ('dropped', 0)])

        for label, molecule in self.solvent:
            distance = np.linalg.norm(
                np.array(molecule.center_of_mass_in_bohr()) - centre)

            if distance <= inner:
                shell = 'polarizable'
            elif distance <= outer:
                shell = 'nonpolarizable'
            else:
                shells['dropped'] += 1
                continue

            force_field = self._force_field_for(
                label, molecule, force_fields[shell],
                f'build_shells: in the {shell} shell')

            regions[shell].add_molecule(molecule, force_field)

            shells[shell] += 1

        self.embedding = embedding
        self.shells = shells

        return embedding

    def gather_multipoles(self):
        """
        Gets the permanent multipoles of the environment.

        The walk itself is the embedding's own: it pairs each force field with
        the molecules it describes, site against atom, and holds the flat
        arrays the integral drivers take. It is done once and again only after
        the environment has been added to, so asking for it costs nothing and
        the induced dipole solve may ask once an iteration.

        Both regions are walked. The permanent multipoles of the whole
        environment act on the wave function, and what tells the two regions
        apart is the induction which is added on top of this, not this.

        NOTE: a charge is a number and belongs wherever its atom is. A dipole
        is a vector and a quadrupole is a tensor, and both are written in the
        frame the force field was made in. An instance of the molecule turned
        some other way needs them turned with it, and nothing turns anything.
        Whoever contracts the orders above zero has to rotate them first, and
        nothing downstream can tell that they were not: the charges stay right
        while the higher moments point the wrong way, which reads as a small
        error in the energy rather than as a mistake.

        :return:
            The multipoles by order, each an array of coordinates in bohr and
            an array of values. An order no site carries is absent.
        """

        assert_msg_critical(
            self.embedding is not None,
            'gather_multipoles: no environment has been built')

        multipoles = OrderedDict()

        for order in sorted(_MULTIPOLE_NAMES):
            places, values = self.embedding.permanent_multipoles(order)

            if len(values):
                multipoles[order] = (places, values)

        return multipoles

    def compute_permanent_fock(self, qm_molecule, basis, comm=None):
        """
        Computes what the permanent multipoles of the environment add to the
        Fock matrix of the quantum region.

        Only the charges are contracted, with the SIMD nuclear attraction
        driver. An environment carrying dipoles or quadrupoles is refused and
        named rather than computed without them.

        :param qm_molecule:
            The molecule of the quantum region.
        :param basis:
            The molecular basis set of the quantum region.
        :param comm:
            The MPI communicator to divide the sites over, or None to compute
            every site here.

        :return:
            The matrix to add to the Fock matrix on the master rank, and None
            on every other rank. It is the **whole** of the contribution and
            not a share of it: a caller which reduces it again gets it as many
            times over as there are ranks.
        """

        multipoles = self.gather_multipoles()

        self._refuse_what_is_not_implemented(multipoles)

        places, charges = multipoles.get(
            0, (np.zeros((0, 3)), np.zeros(0)))

        if comm is not None and comm.Get_size() > 1:
            start, end = self._share_of(len(charges), comm.Get_rank(),
                                        comm.Get_size())
            places, charges = places[start:end, :], charges[start:end]

        # NOTE: an empty share is a matrix of zeros of the right shape, so a
        # rank with no sites needs no special case.

        matrix = compute_simd_nuclear_potential_integrals(
            qm_molecule, basis, charges, places)

        if comm is not None and comm.Get_size() > 1:

            # NOTE: gathered on the master rank alone, which is where the Fock
            # matrix is built. The other ranks answer None rather than a zero
            # matrix or a share of one, so that a rank which was not meant to
            # use this cannot quietly add the wrong thing.

            matrix = comm.reduce(matrix, root=mpi_master())

        return matrix

    def compute_permanent_energy(self, qm_molecule, basis):
        """
        Computes what the permanent multipoles of the environment do to the
        nuclei of the quantum region.

        This is the other half of the permanent electrostatics. The Fock matrix
        carries what the environment does to the electrons, and this is what it
        does to the nuclei; a total energy without it is wrong by the whole of
        this term.

        The sum itself is the embedding's own and is threaded over the sites,
        which for a real environment are the millions while the nuclei are the
        tens.

        :param qm_molecule:
            The molecule of the quantum region.
        :param basis:
            The molecular basis set of the quantum region, which says what each
            nuclear charge is once its core is described by a potential.

        :return:
            The energy.
        """

        self._refuse_what_is_not_implemented(self.gather_multipoles())

        return self.embedding.permanent_nuclear_energy(qm_molecule, basis)

    def _refuse_what_is_not_implemented(self, multipoles):
        """
        Refuses an environment whose multipoles cannot yet be contracted.

        :param multipoles:
            The multipoles by order.
        """

        for order in sorted(multipoles):
            if order == 0:
                continue

            carriers = self._force_fields_carrying(order)

            assert_msg_critical(
                False,
                f'The permanent {_MULTIPOLE_NAMES[order]} of the environment '
                f'are not implemented yet, and {", ".join(carriers)} '
                f'carries them')

    def _force_fields_carrying(self, order):
        """
        Finds which force fields carry a multipole of an order.

        :param order:
            The order.

        :return:
            The names of the force fields, in the order the regions hold them.
        """

        names = []

        for region in (self.embedding.get_polarizable_region(),
                       self.embedding.get_nonpolarizable_region()):

            for force_field in region.get_force_fields():
                carries = any(site.get_order() >= order
                              for site in force_field.get_sites())

                if carries and force_field.get_name() not in names:
                    names.append(force_field.get_name())

        return names

    @staticmethod
    def _share_of(count, rank, nodes):
        """
        Divides a count into a contiguous share for each rank.

        :param count:
            How many there are.
        :param rank:
            The rank to answer for.
        :param nodes:
            How many ranks there are.

        :return:
            The first and the last of this rank's share.
        """

        ave, res = divmod(count, nodes)

        counts = [ave + 1 if p < res else ave for p in range(nodes)]

        return sum(counts[:rank]), sum(counts[:rank + 1])

    @staticmethod
    def _force_field_for(label, molecule, force_fields, what):
        """
        Finds the force field of the species a molecule belongs to.

        :param label:
            The label of the molecule, whose species is the part of it before
            the number.
        :param molecule:
            The molecule.
        :param force_fields:
            The force fields to look in, keyed by name.
        :param what:
            What is being built, for the message if there is no force field.

        :return:
            The force field.
        """

        species = label.rsplit('-', 1)[0]

        assert_msg_critical(
            species in force_fields,
            f'{what}: there is no force field named {species}, which {label} '
            f'is made of. The potential file holds '
            f'{", ".join(force_fields) if force_fields else "none"}')

        force_field = force_fields[species]

        # NOTE: the regions check this as well, and name the force field and
        # the two counts. What they cannot name is which molecule of the
        # environment was being added, and with five hundred of them that is
        # the part worth having.

        assert_msg_critical(
            len(force_field) == molecule.number_of_atoms(),
            f'{what}: the force field {species} has {len(force_field)} '
            f'site(s) and {label} has {molecule.number_of_atoms()} atom(s)')

        return force_field

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


def read_potential_file(potential_file):
    """
    Reads the force fields of a potential file.

    A potential file holds force field data and nothing else: the parameters of
    each kind of molecule, with no coordinates, no molecules and no exclusion
    lists. Which molecules there are is what a structure file says, and the two
    are read separately.

    :param potential_file:
        The path of the potential file.

    :return:
        The force fields, keyed by name and in the order the file gives them.
    """

    return _PotentialFileReader(potential_file).read()


def _tensor_of(values):
    """
    Builds a symmetric tensor of its six components.

    :param values:
        The components, xx xy xz yy yz zz.

    :return:
        The tensor, three by three.
    """

    xx, xy, xz, yy, yz, zz = values

    return np.array([[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]])


class _PotentialFileReader:
    """
    Reads a potential file.

    The file is a sequence of `@FORCEFIELD` blocks. Each opens with a `name:`
    and a `sites:` count and is followed by the sections which give the
    parameters of those sites, indexed from one:

        `@MULTIPOLES`, holding `ORDER 0`, `ORDER 1` and `ORDER 2` subsections.
        A site carries the highest moment it is listed with, so a site absent
        from `ORDER 2` has no quadrupole rather than one which is zero.

        `@POLARIZABILITIES`, holding either `ISOTROPIC` or `ORDER 1 1`. Which
        of the two is said and not worked out from the numbers: an isotropic
        polarizability is a multiplication where an anisotropic one is a solve,
        and a parameter which happens to be diagonal is not the former.

        `@WIDTHS`, which gives a site the two widths of a Gaussian form. A site
        it does not name is a point.

    Everything after a `#` is a comment, and blank lines are passed over.

    :param potential_file:
        The path of the potential file.
    """

    def __init__(self, potential_file):
        """
        Initializes the reader.
        """

        self.path = Path(potential_file)

        self.lines = []
        self.position = 0

    def read(self):
        """
        Reads the force fields.

        :return:
            The force fields, keyed by name and in the order the file gives them.
        """

        self._load()

        force_fields = OrderedDict()

        while not self._at_end():
            line_number, text = self._peek()

            self._require(
                text.upper().startswith('@FORCEFIELD'), line_number,
                f'expected a @FORCEFIELD section and found {text.split()[0]}')

            name, force_field = self._read_force_field()

            self._require(name not in force_fields, line_number,
                          f'the file gives a force field named {name} twice')

            force_fields[name] = force_field

        return force_fields

    def _load(self):
        """
        Reads the lines of the file which carry something.
        """

        assert_msg_critical(self.path.is_file(),
                            f'read_potential_file: no such file {self.path}')

        with open(self.path, 'r') as handle:
            for number, raw in enumerate(handle, start=1):
                text = raw.split('#', 1)[0].strip()

                if text:
                    self.lines.append((number, text))

        assert_msg_critical(
            len(self.lines) > 0,
            f'read_potential_file: {self.path.name} holds no force field')

    def _read_force_field(self):
        """
        Reads one @FORCEFIELD block.

        :return:
            The name of the force field and the force field.
        """

        opening, _ = self._next()

        name, nsites = self._read_header(opening)

        charges, dipoles, quadrupoles = {}, {}, {}
        isotropic, anisotropic, widths = {}, {}, {}

        while not self._at_end():
            line_number, text = self._peek()

            section = text.upper()

            if section.startswith('@FORCEFIELD'):
                break

            self._next()

            if section.startswith('@MULTIPOLES'):
                self._read_multipoles(nsites, line_number, charges, dipoles,
                                      quadrupoles)

            elif section.startswith('@POLARIZABILITIES'):
                self._read_polarizabilities(nsites, line_number, isotropic,
                                            anisotropic)

            elif section.startswith('@WIDTHS'):
                self._read_rows(nsites, 2, widths, '@WIDTHS')

            else:
                self._refuse(line_number,
                             f'there is no section named {text.split()[0]}')

        return name, self._force_field_of(name, nsites, charges, dipoles,
                                          quadrupoles, isotropic, anisotropic,
                                          widths)

    def _read_header(self, opening):
        """
        Reads the keys which open a @FORCEFIELD block.

        :param opening:
            The number of the line the block opens on.

        :return:
            The name of the force field and the number of sites.
        """

        name = None
        nsites = None

        while not self._at_end() and not self._peek()[1].startswith('@'):
            line_number, text = self._next()

            self._require(':' in text, line_number,
                          'expected a key and a value, as "name: S01"')

            key, _, value = text.partition(':')

            key = key.strip().lower()
            value = value.strip()

            if key == 'name':
                self._require(name is None, line_number,
                              'the force field is named twice')
                self._require(value != '', line_number, 'the name is empty')

                name = value

            elif key == 'sites':
                self._require(nsites is None, line_number,
                              'the number of sites is given twice')

                nsites = self._integer(value, line_number, 'the number of sites')

                self._require(nsites > 0, line_number,
                              f'a force field of {nsites} sites is not one')

            elif key == 'units':

                # NOTE: refused rather than converted. The parameters of this
                # layer are atomic units throughout, and a file written in
                # Debye or in cubic angstrom read as though it were not is a
                # wrong energy rather than a failure.

                self._require(
                    value.lower() in ('au', 'a.u.', 'atomic'), line_number,
                    f'the parameters are read in atomic units and this file '
                    f'says {value}')

            else:
                self._refuse(line_number,
                             f'a force field has no key named {key}')

        self._require(name is not None, opening, 'the force field has no name')
        self._require(nsites is not None, opening,
                      'the force field does not say how many sites it has')

        return name, nsites

    def _read_multipoles(self, nsites, section, charges, dipoles, quadrupoles):
        """
        Reads a @MULTIPOLES section.

        :param nsites:
            The number of sites of the force field.
        :param section:
            The number of the line the section opens on.
        :param charges:
            The charges read so far, which this adds to.
        :param dipoles:
            The dipoles read so far, which this adds to.
        :param quadrupoles:
            The quadrupoles read so far, which this adds to.
        """

        targets = {0: charges, 1: dipoles, 2: quadrupoles}

        orders = 0

        while not self._at_end() and self._peek()[1].upper().startswith('ORDER'):
            line_number, text = self._next()

            fields = text.split()

            self._require(len(fields) == 2, line_number,
                          'expected "ORDER k", with k one of 0, 1 or 2')

            order = self._integer(fields[1], line_number, 'the order')

            self._require(order in targets, line_number,
                          f'the order {order} is not one of 0, 1 or 2')

            self._read_rows(nsites, _MULTIPOLE_SIZES[order], targets[order],
                            f'ORDER {order}')

            orders += 1

        self._require(orders > 0, section,
                      '@MULTIPOLES holds no ORDER subsection')

    def _read_polarizabilities(self, nsites, section, isotropic, anisotropic):
        """
        Reads a @POLARIZABILITIES section.

        :param nsites:
            The number of sites of the force field.
        :param section:
            The number of the line the section opens on.
        :param isotropic:
            The isotropic polarizabilities read so far, which this adds to.
        :param anisotropic:
            The anisotropic polarizabilities read so far, which this adds to.
        """

        self._require(not self._at_end(), section,
                      '@POLARIZABILITIES says nothing')

        line_number, text = self._next()

        kind = text.upper()

        if kind.startswith('ISOTROPIC'):
            self._read_rows(nsites, 1, isotropic, 'ISOTROPIC')

        elif kind.startswith('ORDER'):
            fields = text.split()

            self._require(fields[1:] == ['1', '1'], line_number,
                          'the only polarizability is "ORDER 1 1"')

            self._read_rows(nsites, 6, anisotropic, 'ORDER 1 1')

        else:
            self._refuse(
                line_number,
                '@POLARIZABILITIES takes ISOTROPIC or ORDER 1 1 and this is '
                f'{text.split()[0]}')

    def _read_rows(self, nsites, nvalues, target, what):
        """
        Reads a count and that many rows of a site index and its values.

        :param nsites:
            The number of sites of the force field.
        :param nvalues:
            How many values a row carries.
        :param target:
            The rows read so far, which this adds to.
        :param what:
            The name of the subsection, for the message if a row is not one.
        """

        line_number, text = self._next()

        count = self._integer(text, line_number, f'the number of rows of {what}')

        self._require(count >= 0, line_number,
                      f'{what} says it has {count} rows')

        for _ in range(count):
            self._require(
                not self._at_end() and not self._peek()[1].startswith('@'),
                line_number,
                f'{what} says it has {count} rows and it ends before that many')

            row_number, row = self._next()

            fields = row.split()

            self._require(
                len(fields) == nvalues + 1, row_number,
                f'{what} takes a site index and {nvalues} value(s), and this '
                f'row has {len(fields)} field(s)')

            index = self._integer(fields[0], row_number, 'the site index')

            self._require(
                1 <= index <= nsites, row_number,
                f'the site index {index} is not between 1 and {nsites}')

            self._require(index not in target, row_number,
                          f'{what} gives the site {index} twice')

            target[index] = (tuple(
                self._number(field, row_number, 'a value')
                for field in fields[1:]), row_number)

    def _force_field_of(self, name, nsites, charges, dipoles, quadrupoles,
                        isotropic, anisotropic, widths):
        """
        Builds the force field of what the sections gave.

        :param name:
            The name of the force field.
        :param nsites:
            The number of sites.
        :param charges:
            The charges.
        :param dipoles:
            The dipoles.
        :param quadrupoles:
            The quadrupoles.
        :param isotropic:
            The isotropic polarizabilities.
        :param anisotropic:
            The anisotropic polarizabilities.
        :param widths:
            The widths of the Gaussian sites.

        :return:
            The force field.
        """

        sites = []

        for index in range(1, nsites + 1):
            site = {}

            if index in charges:
                site['charge'] = charges[index][0][0]

            if index in dipoles:
                site['dipole'] = np.array(dipoles[index][0])

            if index in quadrupoles:
                site['quadrupole'] = _tensor_of(quadrupoles[index][0])

            if index in isotropic:
                self._require(
                    index not in anisotropic, anisotropic[index][1]
                    if index in anisotropic else 0,
                    f'the site {index} is given a polarizability twice')

                site['polarizability'] = isotropic[index][0][0]

            elif index in anisotropic:
                site['polarizability'] = _tensor_of(anisotropic[index][0])

            if index in widths:
                values, row_number = widths[index]

                self._require(
                    min(values) > 0.0, row_number,
                    f'the width of the site {index} is not positive')

                site['multipole_width'] = values[0]
                site['polarizability_width'] = values[1]

            sites.append(site)

        return PolarizableForceField(name=name, sites=sites)

    def _at_end(self):
        """
        Checks whether every line has been read.

        :return:
            True if it has.
        """

        return self.position >= len(self.lines)

    def _peek(self):
        """
        Gets the line which comes next, without taking it.

        :return:
            Its number and its text.
        """

        return self.lines[self.position]

    def _next(self):
        """
        Takes the line which comes next.

        :return:
            Its number and its text.
        """

        line = self.lines[self.position]

        self.position += 1

        return line

    def _integer(self, text, line_number, what):
        """
        Reads a whole number.

        :param text:
            The text of it.
        :param line_number:
            The number of the line it is on.
        :param what:
            What it is, for the message if it is not a number.

        :return:
            The number.
        """

        try:
            return int(text)
        except ValueError:
            self._refuse(line_number, f'{what} is not a whole number: {text}')

    def _number(self, text, line_number, what):
        """
        Reads a number.

        :param text:
            The text of it.
        :param line_number:
            The number of the line it is on.
        :param what:
            What it is, for the message if it is not a number.

        :return:
            The number.
        """

        try:
            return float(text)
        except ValueError:
            self._refuse(line_number, f'{what} is not a number: {text}')

    def _require(self, condition, line_number, message):
        """
        Refuses the file unless the condition holds.

        :param condition:
            The condition.
        :param line_number:
            The number of the line at fault.
        :param message:
            The reason, if it does not hold.
        """

        if not condition:
            self._refuse(line_number, message)

    def _refuse(self, line_number, message):
        """
        Refuses the file.

        :param line_number:
            The number of the line at fault.
        :param message:
            The reason.
        """

        assert_msg_critical(
            False, f'{self.path.name}, line {line_number}: {message}')

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
import json
import time
import numpy as np

from .veloxchemlib import (mpi_master, FikaChargeSummation,
                           FikaClassicalSystem, FikaQmmmEmbeddingDriver,
                           FikaQmmmEmbeddingOptions, FikaQmmmSources,
                           FikaTholeDamping)
from .errorhandler import assert_msg_critical


def is_fika_embedding(embedding):
    """
    Checks whether an embedding dictionary selects the fika embedding.

    :param embedding:
        The embedding dictionary (or None).

    :return:
        True if settings.embedding_method is 'fika' (case insensitive).
    """

    if not isinstance(embedding, dict):
        return False
    method = embedding.get('settings', {}).get('embedding_method', '')
    return isinstance(method, str) and method.lower() == 'fika'


def fika_embedding_text(embedding):
    """
    Gets a canonical text of a fika embedding dictionary (its settings and
    inputs; objects by type), e.g. to compare checkpoints.

    :param embedding:
        The embedding dictionary.

    :return:
        The text.
    """

    return json.dumps(embedding, sort_keys=True,
                      default=lambda obj: f'<{type(obj).__name__}>')


def fika_embedding_sanity_check(embedding):
    """
    Checks a fika embedding dictionary:

        settings:
            - embedding_method: 'fika'
            - damping: 'thole' (default) or None
            - thole_a: the Thole parameter (default 2.1304)
            - induced_dipoles: tolerance, field_accuracy, max_iterations,
              summation ('automatic', 'direct' or 'multipole')
            - fock_threshold, fock_summation
        inputs: one of
            - pdb_file, solvent, and optionally radius or shells, unit
              ('angstrom' or 'bohr') and solute (residue name, 'MOL'): the QM
              molecule is the solute of the PDB file;
            - objects: {'classical_system': FikaClassicalSystem}.

    :param embedding:
        The embedding dictionary.
    """

    assert_msg_critical(is_fika_embedding(embedding),
                        'fika embedding: embedding_method must be fika')
    settings = embedding['settings']
    inputs = embedding.get('inputs')
    assert_msg_critical(isinstance(inputs, dict),
                        "fika embedding: missing 'inputs' dictionary")

    known = {
        'embedding_method', 'damping', 'thole_a', 'induced_dipoles',
        'fock_threshold', 'fock_summation'
    }
    unknown = set(settings) - known
    assert_msg_critical(
        not unknown,
        f'fika embedding: unknown settings {sorted(unknown)}')
    damping = settings.get('damping', 'thole')
    assert_msg_critical(
        damping is None or str(damping).lower() in ('thole', 'none'),
        "fika embedding: damping must be 'thole' or None")
    dipoles = settings.get('induced_dipoles', {})
    assert_msg_critical(
        isinstance(dipoles, dict) and not (set(dipoles) - {
            'tolerance', 'field_accuracy', 'max_iterations', 'summation'
        }), 'fika embedding: induced_dipoles takes tolerance, ' +
        'field_accuracy, max_iterations and summation')
    for key in ('summation',):
        if key in dipoles:
            assert_msg_critical(
                str(dipoles[key]).lower() in ('automatic', 'direct',
                                              'multipole'),
                'fika embedding: summation must be automatic, direct or ' +
                'multipole')

    has_pdb = 'pdb_file' in inputs
    has_objects = 'objects' in inputs
    assert_msg_critical(
        has_pdb != has_objects,
        "fika embedding: inputs need either 'pdb_file' or 'objects'")
    if has_pdb:
        assert_msg_critical('solvent' in inputs,
                            "fika embedding: inputs need 'solvent'")
        unknown = set(inputs) - {
            'pdb_file', 'solvent', 'radius', 'shells', 'unit', 'solute'
        }
    else:
        assert_msg_critical(
            isinstance(inputs['objects'], dict) and isinstance(
                inputs['objects'].get('classical_system'),
                FikaClassicalSystem),
            "fika embedding: objects need a 'classical_system' " +
            '(FikaClassicalSystem)')
        unknown = set(inputs) - {'objects'}
    assert_msg_critical(not unknown,
                        f'fika embedding: unknown inputs {sorted(unknown)}')


class FikaEmbedding:
    """
    Polarizable embedding through fika: the classical system (from a PDB
    file with the fika force-field library, or given as an object) and the
    options of FikaQmmmEmbeddingDriver. The embedding is computed on the
    master rank (with its OpenMP threads); the other ranks get no
    contribution, as the SCF and response drivers use it there only.

    Tolerances follow the convergence threshold of the calling driver: the
    induced-dipole tolerance min(1e-8, conv_thresh / 10), the field accuracy
    a tenth of it and the Fock screening threshold min(1e-12, conv_thresh *
    1e-4), unless set.

    :param molecule:
        The QM molecule; with a PDB file it must be the solute.
    :param ao_basis:
        The AO basis set.
    :param options:
        The embedding dictionary (see fika_embedding_sanity_check).
    :param comm:
        The MPI communicator.
    :param conv_thresh:
        The convergence threshold of the calling driver.
    """

    def __init__(self, molecule, ao_basis, options, comm=None, conv_thresh=1.0e-6):

        if comm is None:
            comm = MPI.COMM_WORLD

        fika_embedding_sanity_check(options)

        self.comm = comm
        self.rank = comm.Get_rank()
        self.molecule = molecule
        self.basis = ao_basis
        self.options = options

        settings = options['settings']
        inputs = options['inputs']

        if 'pdb_file' in inputs:
            from .fikapdbreader import FikaPdbReader

            reader = FikaPdbReader(comm)
            residues = reader.read(inputs['pdb_file'])
            solute = inputs.get('solute', 'MOL')
            self.classical_system, self.residues = reader.get_classical_system(
                residues,
                inputs['solvent'],
                radius=inputs.get('radius'),
                shells=inputs.get('shells'),
                solute=solute,
                unit=inputs.get('unit', 'angstrom'))
            self._check_molecule(reader.get_solute(residues, solute))
        else:
            self.classical_system = inputs['objects']['classical_system']
            self.residues = None

        damping = settings.get('damping', 'thole')
        if damping is None or str(damping).lower() == 'none':
            self.damping = None
        else:
            self.damping = FikaTholeDamping(settings.get('thole_a', 2.1304))

        dipoles = settings.get('induced_dipoles', {})
        self.driver_options = FikaQmmmEmbeddingOptions()
        induced = self.driver_options.induced
        induced.tolerance = dipoles.get('tolerance',
                                        min(1.0e-8, 0.1 * conv_thresh))
        induced.field_accuracy = dipoles.get('field_accuracy',
                                             0.1 * induced.tolerance)
        if 'max_iterations' in dipoles:
            induced.max_iterations = dipoles['max_iterations']
        induced.summation = self._summation(dipoles.get('summation',
                                                        'automatic'))
        self.driver_options.induced = induced
        self.driver_options.fock_threshold = settings.get(
            'fock_threshold', min(1.0e-12, 1.0e-4 * conv_thresh))
        self.driver_options.fock_summation = self._summation(
            settings.get('fock_summation', 'automatic'))

        self.result = None
        self.timing = 0.0

    @staticmethod
    def _summation(value):

        return {
            'automatic': FikaChargeSummation.automatic,
            'direct': FikaChargeSummation.direct,
            'multipole': FikaChargeSummation.multipole,
        }[str(value).lower()]

    def _check_molecule(self, solute):
        """
        Checks that the QM molecule is the solute of the PDB file.
        """

        same = np.array_equal(self.molecule.get_element_ids(),
                              solute.get_element_ids())
        if same:
            same = np.max(
                np.abs(self.molecule.get_coordinates_in_bohr() -
                       solute.get_coordinates_in_bohr()),
                initial=0.0) < 1.0e-8
        assert_msg_critical(
            same, 'fika embedding: the molecule is not the solute of the ' +
            'PDB file (take it from FikaPdbReader.get_solute)')

    def _compute(self, density_matrix, sources):
        """
        Computes the embedding of a density on the master rank.
        """

        t0 = time.time()
        options = self.driver_options
        options.sources = sources
        result = FikaQmmmEmbeddingDriver(options).compute(
            self.molecule, self.basis, np.asarray(density_matrix),
            self.classical_system, self.damping)
        assert_msg_critical(
            result.induced.converged,
            'fika embedding: induced dipoles did not converge in ' +
            f'{result.induced.iterations} iterations (residual ' +
            f'{result.induced.residual:.2e})')
        self.timing += time.time() - t0
        return result

    def get_info(self):
        """
        Gets lines describing the embedding (for the output header).
        """

        system = self.classical_system
        lines = ['Polarizable embedding: fika']
        if 'pdb_file' in self.options['inputs']:
            inputs = self.options['inputs']
            lines.append(f"- pdb_file        : {inputs['pdb_file']}")
            lines.append(f"- solvent         : {inputs['solvent']}")
            for key in ('radius', 'shells'):
                if key in inputs:
                    lines.append(f'- {key:<15s} : {inputs[key]} ' +
                                 inputs.get('unit', 'angstrom'))
        lines.append('- residues        : ' +
                     f'{system.number_of_residues(True)} polarizable, ' +
                     f'{system.number_of_residues(False)} nonpolarizable')
        lines.append('- damping         : ' + (
            'none' if self.damping is None else f'Thole (a = {self.damping.a})'))
        induced = self.driver_options.induced
        lines.append(f'- tolerance       : {induced.tolerance:.1e} ' +
                     f'(field accuracy {induced.field_accuracy:.1e}, ' +
                     f'Fock threshold {self.driver_options.fock_threshold:.1e})')
        return lines


class FikaEmbeddingSCF(FikaEmbedding):
    """
    The fika embedding of the SCF density: energy and Fock contribution, the
    induced dipoles of each iteration starting from those of the previous
    one.
    """

    def __init__(self, molecule, ao_basis, options, comm=None, conv_thresh=1.0e-6):

        super().__init__(molecule, ao_basis, options, comm, conv_thresh)
        self.iterations = []

    def compute_pe_contributions(self, density_matrix):
        """
        Computes the embedding energy and Fock contribution of the total
        density (master rank).

        :param density_matrix:
            The total (alpha + beta) AO density matrix.

        :return:
            The embedding energy and Fock contribution (0 and None on other
            ranks).
        """

        if self.rank != mpi_master():
            return 0.0, None

        if self.result is not None:
            self.driver_options.induced.initial_guess = self.result.induced.dipoles
        self.result = self._compute(density_matrix, FikaQmmmSources.all)
        self.iterations.append(self.result.induced.iterations)
        return self.result.energy(), self.result.fock()

    def get_induced_dipoles(self):
        """
        Gets the induced dipoles of the last evaluation (master rank), in the
        order of fika_polarizable_sites.
        """

        return None if self.result is None else self.result.induced.dipoles

    def get_pe_summary(self):
        """
        Gets lines summarizing the embedding energy (master rank).
        """

        if self.result is None:
            return []
        r = self.result
        return [
            'Polarizable Embedding (fika) Energy Contributions',
            '------------------------------------',
            'Electrostatic Contribution         :' +
            f'{r.electron_permanent_energy + r.nuclear_permanent_energy:20.10f} a.u.',
            'Induced Contribution               :' +
            f'{r.polarization_energy:20.10f} a.u.',
            '------------------------------------',
            f'Induced dipoles: {len(r.induced.dipoles)} sites, ' +
            f'{sum(self.iterations)} CG iterations in {len(self.iterations)} ' +
            f'evaluations, {self.timing:.2f} sec.',
        ]


class FikaEmbeddingLRS(FikaEmbedding):
    """
    The fika embedding in linear response: the Fock contribution of the
    dipoles induced by the electron field of a perturbed density alone.
    """

    def compute_pe_contributions(self, density_matrix):
        """
        Computes the Fock contribution of a perturbed density (master rank).

        :param density_matrix:
            The perturbed total AO density matrix (any symmetry; its
            symmetric part acts).

        :return:
            The Fock contribution (None on other ranks).
        """

        if self.rank != mpi_master():
            return None

        return self._compute(density_matrix,
                             FikaQmmmSources.electrons_only).induced_fock

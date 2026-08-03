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

import sys
import numpy as np

from .veloxchemlib import mpi_master, hartree_in_ev
from .outputstream import OutputStream
from .errorhandler import assert_msg_critical
from .openqpscfdriver import OpenQPScfDriver


class OpenQPExcitedStatesDriver:
    """
    Implements the OpenQP excited-states driver for MRSF-/SF-/TD-DFT.

    OpenQP performs the reference SCF and the excited-state response
    internally.  This driver configures the OpenQP calculation
    (method='tdhf'), runs the reference and excitation through the shared
    OpenQP engine, and returns the excited-state energies.

    :param openqp_scf_drv:
        The OpenQP SCF driver providing the OpenQP engine.

    Instance variables
        - tddft_type: OpenQP TD type (`mrsf`, `sf`, `tda`, or `rpa`).
        - nstates: Number of excited states (OpenQP roots) to compute.
        - tddft_multiplicity: TD spin multiplicity (1 => singlet).
    """

    def __init__(self, openqp_scf_drv):
        """
        Initializes the OpenQP excited-states driver.
        """

        errmsg = 'OpenQPExcitedStatesDriver: invalid OpenQP SCF driver.'
        assert_msg_critical(isinstance(openqp_scf_drv, OpenQPScfDriver), errmsg)

        self.openqp_driver = openqp_scf_drv
        self.comm = openqp_scf_drv.comm
        self.rank = openqp_scf_drv.rank
        self.nodes = openqp_scf_drv.nodes
        self.ostream = openqp_scf_drv.ostream

        self.tddft_type = 'mrsf'
        self.nstates = 10
        self.tddft_multiplicity = 1

        self.flag = 'OpenQP Excited States Driver'

    def set_tddft_type(self, tddft_type):
        """
        Sets the OpenQP TD type (`mrsf`, `sf`, `tda`, or `rpa`).
        """

        ttype = str(tddft_type).strip().lower()
        assert_msg_critical(
            ttype in ('mrsf', 'sf', 'tda', 'rpa', 'umrsf'),
            f'OpenQPExcitedStatesDriver: invalid tddft_type "{tddft_type}"')
        self.tddft_type = ttype

    def set_nstates(self, nstates):
        """
        Sets the number of excited-state roots to compute.
        """

        self.nstates = int(nstates)

    def set_multiplicity(self, multiplicity):
        """
        Sets the TD spin multiplicity (1 => singlet, 3 => triplet).
        """

        self.tddft_multiplicity = int(multiplicity)

    def update_settings(self, rsp_dict=None, method_dict=None):
        """
        Updates settings in the OpenQP excited-states driver.

        :param rsp_dict:
            The input dictionary of response settings.
        :param method_dict:
            The input dictionary of method settings.
        """

        if rsp_dict is None:
            rsp_dict = {}

        if 'tddft_type' in rsp_dict:
            self.set_tddft_type(rsp_dict['tddft_type'])
        if 'method' in rsp_dict:
            self.set_tddft_type(rsp_dict['method'])
        if 'nstates' in rsp_dict:
            self.set_nstates(rsp_dict['nstates'])
        if 'nstate' in rsp_dict:
            self.set_nstates(rsp_dict['nstate'])
        if 'multiplicity' in rsp_dict:
            self.set_multiplicity(rsp_dict['multiplicity'])

    def compute(self, molecule):
        """
        Computes the OpenQP excited-state energies.

        :param molecule:
            The molecule.

        :return:
            A results dictionary on the master rank, otherwise None.
            Keys: 'reference_energy', 'ground_state_energy', 'state_energies'
            (absolute energies of the OpenQP roots S0, S1, ...),
            'excitation_energies' (relative to the ground-state root, in
            Hartree), and 'excitation_energies_ev' (in eV).
        """

        errmsg = 'OpenQPExcitedStatesDriver: oqp is not available. '
        errmsg += 'Please install OpenQP (pip install openqp) in this '
        errmsg += 'environment.'
        assert_msg_critical(OpenQPScfDriver.is_available(), errmsg)

        if self.rank == mpi_master():
            self.print_title()
            results = self._compute_master(molecule)
        else:
            results = None

        if self.rank == mpi_master():
            self.print_results(results)
            self.ostream.print_blank()
            self.ostream.flush()
            return results
        return None

    def _compute_master(self, molecule):
        drv = self.openqp_driver

        scf_type = 'rohf' if self.tddft_type in ('mrsf', 'sf', 'umrsf') else \
            drv._effective_scf_type(molecule)
        multiplicity = 3 if self.tddft_type in ('mrsf', 'sf', 'umrsf') else \
            drv._effective_multiplicity(molecule)
        functional = drv.functional if drv.method == 'dft' else ''

        mol = drv.run_oqp(
            molecule,
            method='tdhf',
            functional=functional,
            scf_type=scf_type,
            multiplicity=multiplicity,
            exc_multiplicity=self.tddft_multiplicity,
            tdhf_type=self.tddft_type,
            nstate=self.nstates,
        )

        energies = np.array(mol.energies, dtype=float)

        # OpenQP energy layout: index 0 is the reference (SCF) energy; indices
        # 1.. are the physical roots (S0, S1, ...).  Excitation energies are
        # taken relative to the lowest root (the ground state, index 1).
        reference_energy = float(energies[0])
        state_energies = energies[1:].copy()
        ground_state_energy = float(state_energies[0])
        excitation_energies = state_energies - ground_state_energy

        drv._energy = reference_energy

        return {
            'reference_energy': reference_energy,
            'ground_state_energy': ground_state_energy,
            'state_energies': state_energies,
            'excitation_energies': excitation_energies,
            'excitation_energies_ev': excitation_energies * hartree_in_ev(),
        }

    def print_title(self):
        """
        Prints the title for the OpenQP excited-states calculation.
        """

        self.ostream.print_blank()
        self.ostream.print_header('OpenQP Excited States Driver')
        self.ostream.print_header(28 * '=')
        self.ostream.print_blank()
        self.ostream.print_info(
            f'TD type            : {self.tddft_type.upper()}')
        self.ostream.print_info(f'Number of states   : {self.nstates}')
        self.ostream.flush()

    def print_results(self, results):
        """
        Prints the OpenQP excited-state energies.
        """

        if results is None:
            return

        self.ostream.print_blank()
        self.ostream.print_info(
            f'Reference energy   : {results["reference_energy"]:18.10f} a.u.')
        self.ostream.print_info(
            f'Ground-state energy: {results["ground_state_energy"]:18.10f} a.u.')
        self.ostream.print_blank()
        self.ostream.print_info('  State        E (a.u.)        dE (eV)')
        for i, (e_abs, e_exc_ev) in enumerate(zip(
                results['state_energies'], results['excitation_energies_ev'])):
            self.ostream.print_info(
                f'  S{i:<3d}   {e_abs:18.10f}   {e_exc_ev:12.4f}')
        self.ostream.flush()

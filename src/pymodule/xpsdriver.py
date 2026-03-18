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
import numpy as np
import copy
import sys

from .veloxchemlib import hartree_in_ev, mpi_master
from .outputstream import OutputStream
from .scfunrestdriver import ScfUnrestrictedDriver
from .errorhandler import assert_msg_critical
from .spectrumplot import plot_xps_spectrum


class XPSDriver:
    """
    Implements XPS (X-ray Photoelectron Spectroscopy) calculations
    using the full core-hole (FCH) approximation.

    # vlxtag: RHF, XPS
    # vlxtag: RKS, XPS

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - energy_ranges: Dictionary of core orbital energy ranges for different elements (in Hartree).
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes the XPS driver.
        """

        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(sys.stdout)
            else:
                ostream = OutputStream(None)

        self.comm = comm
        self.ostream = ostream
        self.rank = self.comm.Get_rank()

        # Core orbital energy ranges (in Hartree) for different elements
        # These are approximate ranges for identifying 1s core orbitals
        self.energy_ranges = {
            'C': (-10.4, -10.1),
            'O': (-19.4, -19.1),
            'S': (-9.0, -8.7),
            'N': (-14.5, -14.2),
        }

    def update_settings(self, xps_dict, method_dict=None):
        """
        Updates settings in XPS driver.

        :param xps_dict:
            The dictionary of XPS input.
        :param method_dict:
            The dictionary of method settings.
        """

        if method_dict is None:
            method_dict = {}

        xps_keywords = {}

        for key in xps_dict:
            if key in xps_keywords:
                setattr(self, key, xps_dict[key])

    def _find_core_orbital_indices(self, scf_results, element_symbol):
        """
        Finds the indices of core orbitals for a given element.

        :param scf_results:
            The SCF results dictionary containing orbital energies.
        :param element_symbol:
            The chemical symbol of the element (e.g., 'C', 'O', 'N', 'S').

        :return:
            List of orbital indices corresponding to core orbitals of the element.
        """

        if 'E_alpha' in scf_results:
            orbital_energies = scf_results['E_alpha']
        elif 'orbital_energies_alpha' in scf_results:
            orbital_energies = scf_results['orbital_energies_alpha']
        else:
            assert_msg_critical(
                False,
                f'XPSDriver._find_core_orbital_indices: Cannot find orbital energies in scf_results')

        assert_msg_critical(
            element_symbol in self.energy_ranges,
            f'XPSDriver._find_core_orbital_indices: Element {element_symbol} not supported. '
            f'Supported elements: {list(self.energy_ranges.keys())}')

        emin, emax = self.energy_ranges[element_symbol]
        indices = [
            index for index, energy in enumerate(orbital_energies)
            if emin <= energy <= emax
        ]

        assert_msg_critical(
            len(indices) > 0,
            f'XPSDriver._find_core_orbital_indices: No core orbital found for element '
            f'{element_symbol} in the energy range [{emin}, {emax}] Hartree.')

        return indices

    def compute(self, molecule, basis, scf_driver, element):
        """
        Computes core ionization energies for all core orbitals of the
        specified element using the full core-hole (FCH) approximation.

        :param molecule:
            The molecule object (neutral ground state).
        :param basis:
            The molecular basis set.
        :param scf_driver:
            The converged SCF driver object (e.g., ScfRestrictedDriver).
        :param element:
            The chemical symbol of the element (e.g., 'C', 'O', 'N', 'S').

        :return:
            Dictionary with element as key and list of (orbital_index, ionization_energy)
            tuples as value. Ionization energies are in eV.
        """

        if self.rank == mpi_master():
            self.ostream.print_blank()
            self.ostream.print_header('XPS Driver')
            self.ostream.print_header('=' * 80)
            self.ostream.print_blank()
            self.ostream.print_info(f'Computing core ionization energies for element: {element}')
            self.ostream.print_blank()
            self.ostream.flush()

        # Get SCF results
        scf_results = scf_driver.scf_tensors

        # Get ground state energy
        gs_energy = scf_driver.get_scf_energy()

        # Get molecular orbitals
        orbs = scf_driver.molecular_orbitals

        # Get XC functional from ground state calculation
        xcfun = scf_driver.xcfun

        # Find core orbital indices for the specified element
        if self.rank == mpi_master():
            core_orbital_indices = self._find_core_orbital_indices(scf_results, element)
            self.ostream.print_info(
                f'Found {len(core_orbital_indices)} core orbital(s) for {element}: '
                f'{core_orbital_indices}')
            self.ostream.print_blank()
            self.ostream.flush()
        else:
            core_orbital_indices = None

        core_orbital_indices = self.comm.bcast(core_orbital_indices, root=mpi_master())

        # Create molecular ion (+1 charge, doublet spin state)
        molecular_ion = copy.deepcopy(molecule)
        molecular_ion.set_charge(molecule.get_charge() + 1)
        molecular_ion.set_multiplicity(2)

        # Number of alpha electrons in neutral molecule
        nalpha = molecule.number_of_alpha_electrons()

        # Calculate ionization energies for each core orbital
        ionization_energies = []

        for mo_index in core_orbital_indices:
            if self.rank == mpi_master():
                self.ostream.print_info(f'Computing FCH for orbital index {mo_index}...')
                self.ostream.flush()

            # Create occupation lists
            # For FCH: remove one beta electron from the core orbital
            occa = np.arange(0, nalpha, 1).tolist()
            occb = np.arange(0, nalpha, 1).tolist()
            occb.remove(mo_index)

            # Setup unrestricted SCF calculation with FCH
            scf_ion = ScfUnrestrictedDriver(self.comm, self.ostream)
            scf_ion.xcfun = xcfun  # Use same functional as ground state

            # Apply maximum overlap method to maintain core hole
            scf_ion.maximum_overlap(molecular_ion, basis, orbs, occa, occb)

            # Compute FCH state
            fch_results = scf_ion.compute(molecular_ion, basis)
            fch_energy = scf_ion.get_scf_energy()

            # Calculate core ionization energy (in eV)
            ie = (fch_energy - gs_energy) * hartree_in_ev()

            ionization_energies.append((mo_index, ie))

            if self.rank == mpi_master():
                self.ostream.print_info(
                    f'  Orbital {mo_index}: Ionization Energy = {ie:.2f} eV')
                self.ostream.print_blank()
                self.ostream.flush()

        # Create results dictionary
        results = {element: ionization_energies}

        if self.rank == mpi_master():
            self.ostream.print_blank()
            self.ostream.print_header('XPS Calculation Completed')
            self.ostream.print_header('=' * 80)
            self.ostream.print_blank()
            self.ostream.flush()

        return results

    def print_results(self, results):
        """
        Prints XPS results in a formatted table.

        :param results:
            The results dictionary from compute method.
        """

        if self.rank != mpi_master():
            return

        self.ostream.print_blank()
        self.ostream.print_header('XPS Core Ionization Energies')
        self.ostream.print_header('=' * 60)
        self.ostream.print_blank()

        for element, ionization_data in results.items():
            self.ostream.print_info(f'Element: {element}')
            self.ostream.print_info('-' * 60)
            self.ostream.print_info(f'{"Orbital Index":<20} {"Ionization Energy (eV)":<30}')
            self.ostream.print_info('-' * 60)

            for orbital_idx, ie in ionization_data:
                self.ostream.print_info(f'{orbital_idx:<20} {ie:>25.2f}')

            self.ostream.print_blank()

        self.ostream.print_header('=' * 60)
        self.ostream.print_blank()
        self.ostream.flush()

    def plot_spectrum(self,
                      results,
                      broadening_type="lorentzian",
                      broadening_value=0.5,
                      ax=None):
        """
        Plot the XPS spectrum from the computed results.

        :param results:
            The results dictionary from compute method.
        :param broadening_type:
            The type of broadening to use. Either 'lorentzian' or 'gaussian'.
        :param broadening_value:
            The broadening value (FWHM) in eV.
        :param ax:
            The matplotlib axis to plot on. If None, a new figure is created.

        :return:
            The matplotlib axis object.
        """

        return plot_xps_spectrum(results,
                                broadening_type=broadening_type,
                                broadening_value=broadening_value,
                                ax=ax)

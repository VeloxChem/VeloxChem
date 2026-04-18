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
from .sanitychecks import scf_results_sanity_check

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

        # TODO: can HF and DFT be combined in the same range?
        # Core orbital energy ranges (in Hartree) for different elements and methods
        # Based on typical Hartree-Fock calculations
        self.energy_ranges_hf = {
            'C': (-11.5, -11.0),
            'N': (-15.7, -15.6),
            'O': (-20.7, -20.5),
            'F': (-26.4, -26.3),
            'S': (-92.0, -91.0),
        }

        # Based on typical DFT (B3LYP/PBE) calculations
        self.energy_ranges_dft = {
            'C': (-10.5, -9.2),
            'N': (-14.6, -13.8),
            'O': (-19.4, -18.0),
            'F': (-24.9, -24.4),
            'S': (-89.0, -87.0),
        }

        # Default to DFT ranges for backward compatibility
        self.energy_ranges = self.energy_ranges_dft

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

    @staticmethod
    def _get_xcfun_label(xcfun):
        """
        Get the label of a xcfun object.

        :param xcfun:
            The xcfun object.
        :return:
            The xcfun label.
        """

        if xcfun is None:
            return 'HF'

        if isinstance(xcfun, str):
            return xcfun

        return xcfun.get_func_label()

    def _detect_method_type(self, xcfun):
        """
        Detects whether the calculation is HF or DFT based on xcfun.

        :param xcfun:
            The exchange-correlation functional name.

        :return:
            'hf' for Hartree-Fock, 'dft' for DFT methods.
        """

        # HF-like functionals (no correlation or exchange-only)
        hf_like_funcs = ['hf', 'slater']

        xcfun_label = self._get_xcfun_label(xcfun)

        if xcfun_label.lower() in hf_like_funcs:
            return 'hf'
        else:
            return 'dft'

    def _get_energy_ranges(self, method_type):
        """
        Gets energy ranges for the specified method type.

        :param method_type:
            Either 'hf' or 'dft'.

        :return:
            Dictionary of energy ranges for each element.
        """

        if method_type == 'hf':
            return self.energy_ranges_hf
        else:
            return self.energy_ranges_dft

    def _get_ao_to_atom_mapping(self, molecule, basis):
        """
        Creates a mapping from AO index to atom index.

        :param molecule:
            The molecule object.
        :param basis:
            The molecular basis set.

        :return:
            List where ao_to_atom[i] gives the atom index for AO i.
        """

        n_atoms = molecule.number_of_atoms()
        max_angl = basis.max_angular_momentum()

        ao_to_atom = []

        # Loop over angular momentum shells
        for angl in range(max_angl + 1):
            n_sph = angl * 2 + 1  # Number of spherical components

            # Loop over spherical components
            for isph in range(n_sph):

                # Loop over atoms
                for atom_idx in range(n_atoms):
                    n_ao = basis.number_of_basis_functions([atom_idx], angl)

                    # All AOs of this type on this atom
                    for iao in range(n_ao):
                        ao_to_atom.append(atom_idx)

        return ao_to_atom

    def _assign_orbital_to_atom(self, mo_index, mo_coefficients, molecule, basis, element):
        """
        Assigns a core molecular orbital to a specific atom based on MO coefficient analysis.

        Uses Mulliken-like population analysis: the atom with the largest sum of squared
        MO coefficients is assigned as the parent atom of the core orbital.

        :param mo_index:
            Index of the molecular orbital.
        :param mo_coefficients:
            MO coefficient matrix (C_alpha), shape (n_aos, n_mos).
        :param molecule:
            The molecule object.
        :param basis:
            The molecular basis set.
        :param element:
            Element symbol (e.g., 'C', 'O', 'N').

        :return:
            Tuple of (atom_index, contribution) where atom_index is the 0-based index
            and contribution is the fractional contribution (0-1).
        """

        # Get MO coefficients for this orbital
        mo_coeff = mo_coefficients[:, mo_index]

        # Get atom labels
        atom_labels = molecule.get_labels()
        n_atoms = molecule.number_of_atoms()

        # Map each AO to its parent atom
        ao_to_atom_map = self._get_ao_to_atom_mapping(molecule, basis)

        # Calculate contribution from each atom (Mulliken-like analysis)
        atom_contributions = np.zeros(n_atoms)

        for ao_idx, atom_idx in enumerate(ao_to_atom_map):
            # Sum of squared coefficients for this atom
            atom_contributions[atom_idx] += mo_coeff[ao_idx]**2

        # Find atoms of the target element
        target_atom_indices = [i for i, label in enumerate(atom_labels)
                               if label == element]

        # Find which target atom has the largest contribution
        max_contrib = 0.0
        assigned_atom = None

        for atom_idx in target_atom_indices:
            if atom_contributions[atom_idx] > max_contrib:
                max_contrib = atom_contributions[atom_idx]
                assigned_atom = atom_idx

        return assigned_atom, max_contrib

    def _find_core_orbital_indices(self, molecule, basis, scf_results, element_symbol, method_type):
        """
        Finds the indices of core orbitals for a given element and assigns them to specific atoms.

        :param molecule:
            The molecule object.
        :param basis:
            The molecular basis set.
        :param scf_results:
            The SCF results dictionary containing orbital energies and coefficients.
        :param element_symbol:
            The chemical symbol of the element (e.g., 'C', 'N', 'O', 'F', 'S').
        :param method_type:
            The method type ('hf' or 'dft') for selecting appropriate energy ranges.

        :return:
            List of tuples (mo_index, atom_index, contribution) where:
            - mo_index: orbital index
            - atom_index: 0-based atom index
            - contribution: fractional contribution (0-1) of the atom to the MO
        """

        if 'E_alpha' in scf_results:
            orbital_energies = scf_results['E_alpha']
        elif 'orbital_energies_alpha' in scf_results:
            orbital_energies = scf_results['orbital_energies_alpha']
        else:
            assert_msg_critical(
                False,
                f'XPSDriver._find_core_orbital_indices: Cannot find orbital energies in scf_results')

        # Get appropriate energy ranges based on method type
        energy_ranges = self._get_energy_ranges(method_type)

        assert_msg_critical(
            element_symbol in energy_ranges,
            f'XPSDriver._find_core_orbital_indices: Element {element_symbol} not supported. '
            f'Supported elements: {list(energy_ranges.keys())}')

        # Count atoms of this element in the molecule
        elem_labels = molecule.get_labels()
        expected_count = elem_labels.count(element_symbol)

        assert_msg_critical(
            expected_count > 0,
            'XpsDriver: No {element_symbol} atoms found in molecule.')

        # Start with method-specific energy range
        emin, emax = energy_ranges[element_symbol]
        expansion_factor = 0.0
        max_expansions = 5
        
        for expansion in range(max_expansions):
            # Expand range if needed
            if expansion > 0:
                expansion_factor = 0.1 * expansion
                e_range = emax - emin
                emin_expanded = emin - e_range * expansion_factor
                emax_expanded = emax + e_range * expansion_factor
            else:
                emin_expanded = emin
                emax_expanded = emax

            indices = [
                index for index, energy in enumerate(orbital_energies)
                if emin_expanded <= energy <= emax_expanded
            ]

            if len(indices) == expected_count:
                if expansion > 0:
                    self.ostream.print_info(
                        f'Found {len(indices)} core orbital(s) for {element_symbol} '
                        f'with expanded range [{emin_expanded:.3f}, {emax_expanded:.3f}] Ha')
                break
            elif len(indices) > expected_count:
                assert_msg_critical(
                    False,
                    f'XPSDriver._find_core_orbital_indices: Found {len(indices)} core orbitals '
                    f'but expected {expected_count} for {element_symbol}. '
                    'Core orbital identification is ambiguous; adjust energy_ranges.')

        assert_msg_critical(
            len(indices) > 0,
            f'XPSDriver._find_core_orbital_indices: No core orbital found for element '
            f'{element_symbol}. Expected {expected_count} based on molecule composition.')

        if len(indices) < expected_count and self.rank == mpi_master():
            self.ostream.print_info(
                f'Warning: Found only {len(indices)} core orbital(s) but expected {expected_count} '
                f'for {element_symbol}. Consider adjusting energy_ranges.')

        # Get MO coefficients for atom assignment
        assert_msg_critical(
            'C_alpha' in scf_results,
            f'XPSDriver._find_core_orbital_indices: Cannot find MO coefficients in scf_results')

        mo_coefficients = scf_results['C_alpha']

        # Assign each core orbital to an atom
        assignments = []
        for mo_idx in indices:
            atom_idx, contribution = self._assign_orbital_to_atom(
                mo_idx, mo_coefficients, molecule, basis, element_symbol)
            assignments.append((mo_idx, atom_idx, contribution))

        return assignments

    def compute(self, molecule, basis, scf_driver, element=None, elements=None):
        """
        Computes core ionization energies for specified element(s) using 
        the full core-hole (FCH) approximation.

        :param molecule:
            The molecule object (neutral ground state).
        :param basis:
            The molecular basis set.
        :param scf_driver:
            The converged SCF driver object (e.g., ScfRestrictedDriver).
        :param element:
            Single element as string (e.g., 'C', 'O', 'N', 'S').
        :param elements:
            List of element symbols (e.g., ['C', 'O']).

        :return:
            Dictionary with element as key and list of (orbital_index, ionization_energy)
            tuples as value. Ionization energies are in eV.
        """

        # Handle both element (single string) and elements (list)
        if element is not None and elements is not None:
            assert_msg_critical(
                False,
                'XPSDriver.compute: Specify either element or elements, not both.')

        if element is None and elements is None:
            assert_msg_critical(
                False,
                'XPSDriver.compute: Must specify either element or elements parameter.')
        
        if element is not None:
            assert_msg_critical(
                isinstance(element, str),
                'XPSDriver.compute: Invalid element input.')
            elements_list = [element]
        elif elements is not None:
            assert_msg_critical(
                isinstance(elements, (list, tuple)),
                'XPSDriver.compute: Invalid elements input.')
            elements_list = list(elements)

        # Detect method type (HF or DFT)
        method_type = self._detect_method_type(scf_driver.xcfun)

        energy_ranges = self._get_energy_ranges(method_type)

        # The current FCH setup assumes a closed-shell restricted reference.
        assert_msg_critical(
            scf_driver.scf_type == 'restricted' and molecule.get_multiplicity() == 1,
            'XPSDriver.compute: Only closed-shell restricted SCF references are supported.')

        # Validate all elements are supported
        for elem in elements_list:
            assert_msg_critical(
                elem in energy_ranges,
                f'XPSDriver.compute: Element {elem} not supported. '
                f'Supported elements: {list(energy_ranges.keys())}')

        if self.rank == mpi_master():
            # TODO: this could be moved to a self.print_header?
            self.ostream.print_blank()
            self.ostream.print_header('XPS Driver')
            self.ostream.print_header('=' * 12)
            self.ostream.print_blank()
            elem_str = ', '.join(elements_list)
            self.ostream.print_info(f'Computing core ionization energies for element(s): {elem_str}')
            self.ostream.print_info(f'Method type: {method_type.upper()}')
            xcfun_label = self._get_xcfun_label(scf_driver.xcfun)
            if xcfun_label.lower() != 'hf':
                self.ostream.print_info(f'XC functional: {xcfun_label}')
            self.ostream.print_info(f'Using {method_type.upper()} energy ranges for core orbital identification')
            self.ostream.print_blank()
            self.ostream.flush()

        # Get SCF results
        scf_results = scf_driver.scf_tensors

        # Get ground state energy
        gs_energy = scf_driver.get_scf_energy()

        # Get molecular orbitals
        orbs = scf_driver.molecular_orbitals

        # Create molecular ion (+1 charge, doublet spin state)
        molecular_ion = copy.deepcopy(molecule)
        molecular_ion.set_charge(molecule.get_charge() + 1)
        molecular_ion.set_multiplicity(2)

        # Number of alpha electrons in neutral molecule
        nalpha = molecule.number_of_alpha_electrons()

        # Results dictionary
        results = {}

        # Process each element
        for elem in elements_list:
            if self.rank == mpi_master():
                self.ostream.print_blank()
                self.ostream.print_info(f'Processing element: {elem}')
                self.ostream.print_info('-' * 60)

            # Find core orbital indices for this element and assign to atoms
            if self.rank == mpi_master():
                core_orbital_assignments = self._find_core_orbital_indices(
                    molecule, basis, scf_results, elem, method_type)
                if len(core_orbital_assignments) > 0:
                    self.ostream.print_info(
                        f'Found {len(core_orbital_assignments)} core orbital(s) for {elem}:')
                    delocalized = []
                    for mo_idx, atom_idx, contrib in core_orbital_assignments:
                        self.ostream.print_info(
                            f'  MO {mo_idx} -> Atom {atom_idx+1} ({elem}) with {contrib:.1%} contribution')
                        if contrib < 0.75:
                            delocalized.append(mo_idx)
                    if delocalized:
                        # TODO: Localize orbitals. If multiple edges, core orbitals must be
                        # localized separately for each edge.
                        n_deloc = len(delocalized)
                        warning = f' The molecule contains {n_deloc} delocalized core orbitals.'
                        self.ostream.print_blank()
                        self.ostream.print_warning(warning)
                        warning = f' For correct XPS spectra, localize these orbitals first: '
                        for mo_idx in delocalized:
                            warning += f' MO {mo_idx} '
                        self.ostream.print_warning(warning)
                self.ostream.print_blank()
                self.ostream.flush()
            else:
                core_orbital_assignments = None

            core_orbital_assignments = self.comm.bcast(core_orbital_assignments, root=mpi_master())

            # Skip if no orbitals found
            if len(core_orbital_assignments) == 0:
                continue

            # Calculate ionization energies for each core orbital
            ionization_energies = []

            for mo_index, atom_index, contribution in core_orbital_assignments:
                if self.rank == mpi_master():
                    self.ostream.print_info(
                        f'Computing FCH for {elem} orbital {mo_index} '
                        f'(Atom {atom_index+1}, contribution: {contribution:.1%})...')
                    self.ostream.flush()

                # Create occupation lists
                # For FCH: remove one beta electron from the core orbital
                occa = np.arange(0, nalpha, 1).tolist()
                occb = np.arange(0, nalpha, 1).tolist()
                occb.remove(mo_index)

                # Setup unrestricted SCF calculation with FCH
                scf_ion = ScfUnrestrictedDriver(self.comm, self.ostream)

                # Scf sanity check
                scf_results_sanity_check(scf_ion, scf_results)

                # Maximum number of iterations (not included in scf_results)
                scf_ion.max_iter = scf_driver.max_iter

                # Apply maximum overlap method to maintain core hole
                scf_ion.maximum_overlap(molecular_ion, basis, orbs, occa, occb)

                # Compute FCH state
                fch_results = scf_ion.compute(molecular_ion, basis)
                fch_energy = scf_ion.get_scf_energy()

                # Calculate core ionization energy (in eV)
                ie = (fch_energy - gs_energy) * hartree_in_ev()

                ionization_energies.append((mo_index, atom_index, ie, contribution))

                if self.rank == mpi_master():
                    self.ostream.print_info(
                        f'  Orbital {mo_index} (Atom {atom_index+1}): Ionization Energy = {ie:.2f} eV')
                    self.ostream.print_blank()
                    self.ostream.flush()

            results[elem] = ionization_energies

        if self.rank == mpi_master():
            self.ostream.print_blank()
            self.ostream.print_header('*** XPS Calculation Completed.'.ljust(92))
            self.ostream.print_blank()
            self.ostream.flush()

        return results

    def print_results(self, results):
        """
        Prints XPS results in a formatted table with atom assignments.

        :param results:
            The results dictionary from compute method.
        """

        if self.rank != mpi_master():
            return

        self.ostream.print_blank()
        self.ostream.print_header('XPS Core Ionization Energies')
        self.ostream.print_header('=' * 80)
        self.ostream.print_blank()

        for element, ionization_data in results.items():
            self.ostream.print_info(f'Element: {element}')
            self.ostream.print_info('-' * 80)
            self.ostream.print_info(
                f'{"Atom":<8} {"MO Index":<12} {"IE (eV)":<15} {"Localization":<15}')
            self.ostream.print_info('-' * 80)

            # Sort by atom index for better readability
            sorted_data = sorted(ionization_data, key=lambda x: x[1])

            for mo_idx, atom_idx, ie, contribution in sorted_data:
                self.ostream.print_info(
                    f'{atom_idx+1:<8} {mo_idx:<12} {ie:>12.2f}   {contribution:>13.1%}')

            self.ostream.print_blank()

        self.ostream.print_header('=' * 80)
        self.ostream.print_blank()
        self.ostream.flush()

    def plot_spectrum(self,
                      results,
                      element=None,
                      broadening_type="lorentzian",
                      broadening_value=0.5,
                      color='vlx',
                      show_atom_labels=True,
                      color_by_atom=False,
                      ax=None):
        """
        Plot the XPS spectrum for a single element from the computed results.

        :param results:
            The results dictionary from compute method.
        :param element:
            Element symbol to plot (e.g., 'C', 'O', 'N', 'F', 'S').
            If None and results contains only one element, that element is plotted.
            If None and results contains multiple elements, an error is raised.
        :param broadening_type:
            The type of broadening to use. Either 'lorentzian' or 'gaussian'.
        :param broadening_value:
            The broadening value (FWHM) in eV.
        :param color:
            Color scheme for plotting. Either 'vlx' for VeloxChem default color (darkcyan)
            or 'cpk' for CPK coloring (element-specific colors). Default is 'vlx'.
        :param show_atom_labels:
            If True, display atom indices as labels above peaks. Default is True.
        :param color_by_atom:
            If True, color each peak according to its atom index instead of using 
            element color. Default is False.
        :param ax:
            The matplotlib axis to plot on. If None, a new figure is created.

        :return:
            The matplotlib axis object.
        """

        plot_xps_spectrum(results,
                          element=element,
                          broadening_type=broadening_type,
                          broadening_value=broadening_value,
                          color=color,
                          show_atom_labels=show_atom_labels,
                          color_by_atom=color_by_atom,
                          ax=ax)

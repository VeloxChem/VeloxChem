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

import numpy as np
import time as tm
import math
import copy

from .veloxchemlib import OverlapGeom200Driver
from .veloxchemlib import OverlapGeom101Driver
from .veloxchemlib import OverlapGeom100Driver
from .veloxchemlib import KineticEnergyGeom200Driver
from .veloxchemlib import KineticEnergyGeom101Driver
from .veloxchemlib import NuclearPotentialGeom200Driver
from .veloxchemlib import NuclearPotentialGeom020Driver
from .veloxchemlib import NuclearPotentialGeom110Driver
from .veloxchemlib import NuclearPotentialGeom101Driver
from .veloxchemlib import ElectricDipoleMomentGeom100Driver
from .veloxchemlib import FockGeom2000Driver
from .veloxchemlib import FockGeom1100Driver
from .veloxchemlib import FockGeom1010Driver
from .veloxchemlib import XCMolecularHessian
from .veloxchemlib import T4CScreener
from .veloxchemlib import make_matrix, mat_t
from .veloxchemlib import mpi_master, bohr_in_angstrom, hartree_in_kjpermol
from .molecule import Molecule
from .molecularbasis import MolecularBasis
from .griddriver import GridDriver
from .scfgradientdriver import ScfGradientDriver
from .hessiandriver import HessianDriver
from .firstorderprop import FirstOrderProperties
from .hessianorbitalresponse import HessianOrbitalResponse
from .profiler import Profiler
from .matrices import Matrices
from .dftutils import get_default_grid_level
from .errorhandler import assert_msg_critical
from .oneeints import compute_electric_dipole_integrals
from .sanitychecks import (molecule_sanity_check, scf_results_sanity_check,
                           dft_sanity_check, pe_sanity_check)


class ScfHessianDriver(HessianDriver):
    """
    Implements SCF Hessian driver.

    :param scf_drv:
        The SCF driver.

    Instance variables
        - hessian: The Hessian in Hartree per Bohr**2.
        - flag: The type of Hessian driver.
        - perturbed_density: The perturbed density
    """

    def __init__(self, scf_drv):
        """
        Initializes SCF Hessian driver.
        """

        super().__init__(scf_drv.comm, scf_drv.ostream)

        self.flag = 'SCF Hessian Driver'
        self.scf_driver = scf_drv

        self.perturbed_density = None

        # TODO TEMPORARY FLAG
        # Only run orbital response for performance testing
        self.orbrsp_only = False

        # flag for printing the Hessian
        self.do_print_hessian = False

        self.atom_pairs = None

        # TODO: determine _block_size_factor for SCF Hessian driver
        # self._block_size_factor = 4

        self._xcfun_ldstaging = scf_drv._xcfun_ldstaging

        self.use_subcomms = False

        # option dictionaries from input
        # TODO: cleanup
        self.method_dict = {}
        self.cphf_dict = {}

        self._input_keywords['hessian'].update({
            'orbrsp_only':
            ('bool', 'whether to only run CPHF orbital response'),
            'use_subcomms':
            ('bool', 'whether to use subcommunicators in orbital response'),
        })

    def update_settings(self, method_dict, hess_dict=None, cphf_dict=None):
        """
        Updates settings in ScfHessianDriver.

        :param method_dict:
            The input dictionary of method settings group.
        :param hess_dict:
            The input dictionary of Hessian settings group.
        :param cphf_dict:
            The input dictionary of CPHF (orbital response) settings.
        :param rsp_dict:
            The input dictionary for linear response settings
            (needed to compute the polarizability gradient).
        """

        super().update_settings(method_dict, hess_dict)

        if cphf_dict is None:
            cphf_dict = {}

        self.cphf_dict = dict(cphf_dict)

    def compute(self, molecule, ao_basis):
        """
        Computes the analytical or numerical nuclear Hessian.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        """

        if self.rank == mpi_master():
            self.print_header()

        start_time = tm.time()

        profiler = Profiler({
            'timing': self.timing,
            'profiling': self.profiling,
            'memory_profiling': self.memory_profiling,
            'memory_tracing': self.memory_tracing,
        })
        self.profiler = profiler
        self.ostream.unmute()
        self.ostream.print_info("Creating profiler for Hessian driver...")
        self.ostream.mute()

        # Save the electronic energy
        self.elec_energy = self.scf_driver.get_scf_energy()

        # TODO TEMPORARY
        if self.orbrsp_only:
            self.ostream.print_header(
                '*** WARNING only computing Hessian orbital response!')
            self.compute_orbital_response(molecule, ao_basis)
            self.ostream.print_header(
                '*** Hessian orbital response only: DONE  ***')
            self.ostream.flush()
            return

        if self.numerical:
            self.compute_numerical(molecule, ao_basis)
        else:
            self.compute_analytical(molecule, ao_basis, profiler,
                                    self.atom_pairs)

        if self.rank == mpi_master():
            # print Hessian
            if self.do_print_hessian:
                self.print_geometry(molecule)
                self.ostream.print_blank()
                self.print_hessian(molecule)

            valstr = '*** Time spent in Hessian calculation: '
            valstr += '{:.2f} sec ***'.format(tm.time() - start_time)
            self.ostream.print_header(valstr)
            self.ostream.print_blank()
            self.ostream.print_blank()
            self.ostream.flush()

    def compute_numerical(self, molecule, ao_basis):
        """
        Performs the calculation of a numerical Hessian based only
        on the energy.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        """

        self.ostream.mute()

        # atom labels and atom basis labels
        labels = molecule.get_labels()
        atom_basis_labels = molecule.get_atom_basis_labels()

        # main basis label
        basis_label = ao_basis.get_label()

        # number of atoms
        natm = molecule.number_of_atoms()

        # atom coordinates (nx3)
        coords = molecule.get_coordinates_in_bohr()

        # charge and spin multiplicity
        charge = molecule.get_charge()
        multiplicity = molecule.get_multiplicity()

        # Hessian
        hessian = np.zeros((natm, 3, natm, 3))

        # First-order properties for gradient of dipole moment
        prop = FirstOrderProperties(self.comm, self.ostream)
        # numerical gradient (3 dipole components, no. atoms x 3 atom coords)
        if self.rank == mpi_master():
            self.dipole_gradient = np.zeros((3, 3 * natm))

        grad_drv = ScfGradientDriver(self.scf_driver)

        for i in range(natm):

            self.ostream.unmute()
            self.ostream.print_info(f'Processing atom {i + 1}/{natm}...')
            self.ostream.flush()
            self.ostream.mute()

            for x in range(3):

                coords[i, x] += self.delta_h
                new_mol = Molecule(labels, coords, 'au', atom_basis_labels)
                new_mol.set_charge(charge)
                new_mol.set_multiplicity(multiplicity)
                new_bas = MolecularBasis.read(new_mol, basis_label)
                scf_results = self.scf_driver.compute(new_mol, new_bas)
                grad_drv.compute(new_mol, new_bas, scf_results)
                grad_plus = grad_drv.get_gradient().copy()

                prop.compute_scf_prop(new_mol, new_bas, scf_results)
                if self.rank == mpi_master():
                    mu_plus = prop.get_property('dipole moment')

                coords[i, x] -= 2.0 * self.delta_h
                new_mol = Molecule(labels, coords, 'au', atom_basis_labels)
                new_mol.set_charge(charge)
                new_mol.set_multiplicity(multiplicity)
                new_bas = MolecularBasis.read(new_mol, basis_label)
                scf_results = self.scf_driver.compute(new_mol, new_bas)
                grad_drv.compute(new_mol, new_bas, scf_results)
                grad_minus = grad_drv.get_gradient().copy()

                prop.compute_scf_prop(new_mol, new_bas, scf_results)
                if self.rank == mpi_master():
                    mu_minus = prop.get_property('dipole moment')

                coords[i, x] += self.delta_h
                hessian[i, x, :, :] = ((grad_plus - grad_minus) /
                                       (2.0 * self.delta_h))

                if self.rank == mpi_master():
                    for c in range(3):
                        self.dipole_gradient[c, 3 * i +
                                             x] = ((mu_plus[c] - mu_minus[c]) /
                                                   (2.0 * self.delta_h))

        # reshaped Hessian as member variable
        self.hessian = hessian.reshape(natm * 3, natm * 3)

        # restore scf_drv to initial state
        scf_results = self.scf_driver.compute(molecule, ao_basis)
        assert_msg_critical(self.scf_driver.is_converged,
                            'ScfHessianDriver: SCF did not converge')
        self.ostream.unmute()

    def compute_analytical(self, molecule, ao_basis, profiler, atom_pairs=None):
        """
        Computes the analytical nuclear Hessian.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param profiler:
            The profiler.
        :param atom_pairs:
            The atom pairs to compute the Hessian for.
        """

        profiler.timing = True
        profiler.set_timing_key('atom pairs')

        assert_msg_critical(
            self.scf_driver.scf_type == 'restricted',
            'ScfHessianDriver: Analytical Hessian only implemented ' +
            'for restricted case')

        assert_msg_critical(
            self.scf_driver.solvation_model is None,
            'ScfHessianDriver: Solvation model not implemented')

        # sanity checks
        molecule_sanity_check(molecule)
        scf_results_sanity_check(self, self.scf_driver.scf_tensors)
        dft_sanity_check(self, 'compute')

        # use determine_xc_hessian_grid_level here to ensure early exit for
        # unsupported cases
        if self._dft:
            self.determine_xc_hessian_grid_level(
                molecule, get_default_grid_level(self.scf_driver.xcfun))

        self.ostream.print_info('Computing analytical Hessian...')
        self.ostream.print_blank()
        hess_ref = 'P. Deglmann, F. Furche, R. Ahlrichs,'
        hess_ref += ' Chem. Phys. Lett. 2002, 362, 511-518.'
        self.ostream.print_reference('Reference: ' + hess_ref)
        self.ostream.print_blank()
        self.ostream.flush()

        # Preparation

        natm = molecule.number_of_atoms()
        scf_tensors = self.scf_driver.scf_tensors

        if self.rank == mpi_master():
            density = scf_tensors['D_alpha']
            mo = scf_tensors['C_alpha']
            nocc = molecule.number_of_alpha_electrons()
            mo_occ = mo[:, :nocc]
            mo_energies = scf_tensors['E_alpha']
            eocc = mo_energies[:nocc]
            omega_ao = -np.linalg.multi_dot([mo_occ, np.diag(eocc), mo_occ.T])

        else:
            density = None
            omega_ao = None

        density = self.comm.bcast(density, root=mpi_master())
        omega_ao = self.comm.bcast(omega_ao, root=mpi_master())

        # CPHF equations

        cphf_solver = HessianOrbitalResponse(self.comm, self.ostream)
        cphf_solver.update_settings(self.cphf_dict, self.method_dict)

        # TODO: double check analytical Hessian with PE
        assert_msg_critical(
            not self.scf_driver._pe,
            'ScfHessianDriver: Analytical Hessian with ' +
            'polarizable embedding (PE) not yet available')

        if self.scf_driver._pe:

            # TODO: double check
            cphf_solver.embedding = self.scf_driver.embedding
            pe_sanity_check(cphf_solver, molecule=molecule)

            from .embedding import PolarizableEmbeddingHess
            cphf_solver._embedding_hess_drv = PolarizableEmbeddingHess(
                molecule=molecule,
                ao_basis=ao_basis,
                options=self.scf_driver.embedding,
                comm=self.comm,
                density=density * 2)

        # TODO: double check propagation of cphf settings
        cphf_keywords = {
            'timing', 'profiling', 'memory_profiling', 'memory_tracing',
            'use_subcomms', 'filename'
        }
        for key in cphf_keywords:
            setattr(cphf_solver, key, getattr(self, key))

        # todo add atom pair option to cphf solver
        profiler.start_timer('total')
        profiler.start_timer('cphf')
        cphf_solver.compute(molecule, ao_basis, scf_tensors, atom_pairs)
        profiler.stop_timer('cphf')
        cphf_solution_dict = cphf_solver.cphf_results
        dist_cphf_ov = cphf_solution_dict['dist_cphf_ov']
        dist_cphf_rhs = cphf_solution_dict['dist_cphf_rhs']

        hessian_first_integral_derivatives = cphf_solution_dict[
            'hessian_first_integral_derivatives']
        hessian_eri_overlap = cphf_solution_dict['hessian_eri_overlap']

        # First-order contributions
        profiler.start_timer('1e')
        t1 = tm.time()

        # RHS contracted with CPHF coefficients (ov)
        hessian_cphf_coeff_rhs = np.zeros((natm, 3, natm, 3))

        atoms = []
        if atom_pairs is not None:
            for i, j in atom_pairs:
                if i not in atoms:
                    atoms.append(i)
                if j not in atoms:
                    atoms.append(j)
        else:
            atoms = range(natm)

        # TODO: use alternative way to partition atoms
        local_atoms = atoms[self.rank::self.nodes]

        for i in atoms:
            for x in range(3):
                dist_cphf_ov_ix_data = dist_cphf_ov[i * 3 + x].data

                for j in atoms:
                    if atom_pairs is not None:
                        if i != j and (i, j) not in atom_pairs and (
                                j, i) not in atom_pairs:
                            continue
                    for y in range(3):
                        hess_ijxy = 4.0 * (np.dot(
                            dist_cphf_ov_ix_data,
                            dist_cphf_rhs[j * 3 + y].data))
                        hessian_cphf_coeff_rhs[i, x, j, y] += hess_ijxy
                        if i != j:
                            hessian_cphf_coeff_rhs[j, y, i, x] += hess_ijxy

        hessian_cphf_coeff_rhs = self.comm.reduce(hessian_cphf_coeff_rhs,
                                                  root=mpi_master())

        if self.rank == mpi_master():
            hessian_first_order_derivatives = (
                hessian_cphf_coeff_rhs + hessian_first_integral_derivatives +
                hessian_eri_overlap)
        else:
            hessian_first_order_derivatives = None

        self.ostream.print_info('First order derivative contributions' +
                                ' to the Hessian computed in' +
                                ' {:.2f} sec.'.format(tm.time() - t1))
        self.ostream.print_blank()
        self.ostream.flush()
        profiler.stop_timer('1e')

        # Second-order contributions
        profiler.start_timer('2e pure')
        t2 = tm.time()

        ovlp_hess_200_drv = OverlapGeom200Driver()
        ovlp_hess_101_drv = OverlapGeom101Driver()

        kin_hess_200_drv = KineticEnergyGeom200Driver()
        kin_hess_101_drv = KineticEnergyGeom101Driver()

        npot_hess_200_drv = NuclearPotentialGeom200Driver()
        npot_hess_020_drv = NuclearPotentialGeom020Driver()
        npot_hess_110_drv = NuclearPotentialGeom110Driver()
        npot_hess_101_drv = NuclearPotentialGeom101Driver()

        fock_hess_2000_drv = FockGeom2000Driver()
        fock_hess_1100_drv = FockGeom1100Driver()
        fock_hess_1010_drv = FockGeom1010Driver()

        # determine fock_type and exchange_scaling_factor
        if self._dft:
            xcfun = self.scf_driver.xcfun
            if xcfun.is_hybrid():
                fock_type = '2jkx'
                exchange_scaling_factor = xcfun.get_frac_exact_exchange()
                fock_factor = 1.0
            else:
                fock_type = 'j'
                exchange_scaling_factor = 0.0
                fock_factor = 2.0
        else:
            fock_type = '2jk'
            exchange_scaling_factor = 1.0
            fock_factor = 1.0

        # TODO: range-separated Fock
        need_omega = (self._dft and self.scf_driver.xcfun.is_range_separated())
        if need_omega:
            assert_msg_critical(
                False, 'ScfHessianDriver: Not implemented for' +
                ' range-separated functional')

        den_mat_for_fock = make_matrix(ao_basis, mat_t.symmetric)
        den_mat_for_fock.set_values(density)

        den_mat_for_fock2 = make_matrix(ao_basis, mat_t.general)
        den_mat_for_fock2.set_values(density)

        screener = T4CScreener()
        screener.partition(ao_basis, molecule, 'eri')

        thresh_int = int(-math.log10(self.scf_driver.eri_thresh))

        if self.scf_driver.point_charges is not None:
            mm_coords = []
            mm_charges = []
            npoints = self.scf_driver.point_charges.shape[1]
            for p in range(npoints):
                xyz_p = self.scf_driver.point_charges[:3, p]
                chg_p = self.scf_driver.point_charges[3, p]
                mm_coords.append(xyz_p.copy())
                mm_charges.append(chg_p)
        else:
            mm_coords = None
            mm_charges = None

        # Parts related to second-order integral derivatives
        hessian_2nd_order_derivatives = np.zeros((natm, natm, 3, 3))

        for i in local_atoms:

            ovlp_hess_200_mats = ovlp_hess_200_drv.compute(
                molecule, ao_basis, i)

            for x, label_x in enumerate('XYZ'):
                for y, label_y in enumerate('XYZ'):
                    ovlp_label = label_x + label_y if x <= y else label_y + label_x
                    ovlp_iixy = ovlp_hess_200_mats.matrix_to_numpy(ovlp_label)
                    hessian_2nd_order_derivatives[i, i, x, y] += 2.0 * (np.sum(
                        omega_ao * (ovlp_iixy + ovlp_iixy.T)))

            ovlp_hess_200_mats = Matrices()

            kin_hess_200_mats = kin_hess_200_drv.compute(molecule, ao_basis, i)

            for x, label_x in enumerate('XYZ'):
                for y, label_y in enumerate('XYZ'):
                    kin_label = label_x + label_y if x <= y else label_y + label_x
                    kin_200_iixy = kin_hess_200_mats.matrix_to_numpy(kin_label)
                    hessian_2nd_order_derivatives[i, i, x, y] += 2.0 * (np.sum(
                        density * (kin_200_iixy + kin_200_iixy.T)))

            kin_hess_200_mats = Matrices()

            npot_hess_200_mats = npot_hess_200_drv.compute(
                molecule, ao_basis, i)
            npot_hess_020_mats = npot_hess_020_drv.compute(
                molecule, ao_basis, i)

            for x, label_x in enumerate('XYZ'):
                for y, label_y in enumerate('XYZ'):
                    npot_label = label_x + label_y if x <= y else label_y + label_x
                    npot_200_iixy = npot_hess_200_mats.matrix_to_numpy(
                        npot_label)
                    npot_020_iixy = npot_hess_020_mats.matrix_to_numpy(
                        npot_label)
                    # TODO: move minus sign into function call (such as in oneints)
                    hessian_2nd_order_derivatives[i, i, x, y] += -2.0 * (np.sum(
                        density *
                        (npot_200_iixy + npot_200_iixy.T + npot_020_iixy)))

            npot_hess_200_mats = Matrices()
            npot_hess_020_mats = Matrices()

            if self.scf_driver.point_charges is not None:
                hmats_200 = npot_hess_200_drv.compute(molecule, ao_basis, i,
                                                      mm_coords, mm_charges)

                for x, label_x in enumerate('XYZ'):
                    for y, label_y in enumerate('XYZ'):
                        npot_label = label_x + label_y if x <= y else label_y + label_x
                        npot_200_iixy = hmats_200.matrix_to_numpy(npot_label)
                        # TODO: move minus sign into function call (such as in oneints)
                        hessian_2nd_order_derivatives[i, i, x, y] += -2.0 * (
                            np.sum(density * (npot_200_iixy + npot_200_iixy.T)))

                hmats_200 = Matrices()

            screener_atom_i = T4CScreener()
            screener_atom_i.partition_atom(ao_basis, molecule, 'eri', i)

            fock_hess_2000 = fock_hess_2000_drv.compute(
                ao_basis, screener_atom_i, screener, den_mat_for_fock,
                den_mat_for_fock2, i, fock_type, exchange_scaling_factor, 0.0,
                thresh_int)

            # 'XX', 'XY', 'XZ', 'YY', 'YZ', 'ZZ'
            xy_pairs_upper_triang = [(x, y) for x in range(3)
                                     for y in range(x, 3)]

            for idx, (x, y) in enumerate(xy_pairs_upper_triang):
                hess_val = fock_factor * fock_hess_2000[idx]
                hessian_2nd_order_derivatives[i, i, x, y] += hess_val
                if x != y:
                    hessian_2nd_order_derivatives[i, i, y, x] += hess_val
        profiler.stop_timer('2e pure')
        # do only upper triangular matrix
        profiler.start_timer('2e mixed')
        # TODO: use alternative way to partition atom pairs
        if atom_pairs is None:
            all_atom_pairs = [(i, j) for i in range(natm)
                              for j in range(i, natm)]
            local_atom_pairs = all_atom_pairs[self.rank::self.nodes]
        else:
            all_atom_pairs = copy.copy(atom_pairs)
            for i in atoms:
                all_atom_pairs.append((i, i))
            local_atom_pairs = all_atom_pairs[self.rank::self.nodes]

        for i, j in local_atom_pairs:

            ovlp_hess_101_mats = ovlp_hess_101_drv.compute(
                molecule, ao_basis, i, j)

            for x, label_x in enumerate('XYZ'):
                for y, label_y in enumerate('XYZ'):
                    ovlp_label = f'{label_x}_{label_y}'
                    ovlp_ijxy = ovlp_hess_101_mats.matrix_to_numpy(ovlp_label)
                    hessian_2nd_order_derivatives[i, j, x, y] += 2.0 * (np.sum(
                        omega_ao * (ovlp_ijxy + ovlp_ijxy.T)))

            ovlp_hess_101_mats = Matrices()

            kin_hess_101_mats = kin_hess_101_drv.compute(
                molecule, ao_basis, i, j)

            for x, label_x in enumerate('XYZ'):
                for y, label_y in enumerate('XYZ'):
                    kin_label = f'{label_x}_{label_y}'
                    kin_101_ijxy = kin_hess_101_mats.matrix_to_numpy(kin_label)
                    hessian_2nd_order_derivatives[i, j, x, y] += 2.0 * (np.sum(
                        density * (kin_101_ijxy + kin_101_ijxy.T)))

            kin_hess_101_mats = Matrices()

            npot_hess_110_mats_ij = npot_hess_110_drv.compute(
                molecule, ao_basis, i, j)
            npot_hess_110_mats_ji = npot_hess_110_drv.compute(
                molecule, ao_basis, j, i)
            npot_hess_101_mats = npot_hess_101_drv.compute(
                molecule, ao_basis, i, j)

            for x, label_x in enumerate('XYZ'):
                for y, label_y in enumerate('XYZ'):
                    npot_xy_label = f'{label_x}_{label_y}'
                    npot_yx_label = f'{label_y}_{label_x}'
                    npot_110_ijxy = (
                        npot_hess_110_mats_ij.matrix_to_numpy(npot_xy_label) +
                        npot_hess_110_mats_ji.matrix_to_numpy(npot_yx_label))
                    npot_101_ijxy = npot_hess_101_mats.matrix_to_numpy(
                        npot_xy_label)
                    # TODO: move minus sign into function call (such as in oneints)
                    hessian_2nd_order_derivatives[i, j, x, y] += -2.0 * (np.sum(
                        density * (npot_110_ijxy + npot_110_ijxy.T +
                                   npot_101_ijxy + npot_101_ijxy.T)))

            npot_hess_110_mats_ij = Matrices()
            npot_hess_110_mats_ji = Matrices()
            npot_hess_101_mats = Matrices()

            if self.scf_driver.point_charges is not None:
                hmats_101 = npot_hess_101_drv.compute(molecule, ao_basis, i, j,
                                                      mm_coords, mm_charges)

                for x, label_x in enumerate('XYZ'):
                    for y, label_y in enumerate('XYZ'):
                        npot_xy_label = f'{label_x}_{label_y}'
                        npot_101_ijxy = hmats_101.matrix_to_numpy(npot_xy_label)
                        # TODO: move minus sign into function call (such as in oneints)
                        hessian_2nd_order_derivatives[i, j, x, y] += -2.0 * (
                            np.sum(density * (npot_101_ijxy + npot_101_ijxy.T)))

                hmats_101 = Matrices()

            screener_atom_pair = T4CScreener()
            screener_atom_pair.partition_atom_pair(ao_basis, molecule, 'eri', i,
                                                   j)

            fock_hess_1100 = fock_hess_1100_drv.compute(
                ao_basis, screener_atom_pair, screener, den_mat_for_fock,
                den_mat_for_fock2, i, j, fock_type, exchange_scaling_factor,
                0.0, thresh_int)

            screener_atom_i = T4CScreener()
            screener_atom_i.partition_atom(ao_basis, molecule, 'eri', i)

            screener_atom_j = T4CScreener()
            screener_atom_j.partition_atom(ao_basis, molecule, 'eri', j)

            # Note: use general matrix on both sides
            fock_hess_1010 = fock_hess_1010_drv.compute(
                ao_basis, screener_atom_i, screener_atom_j, den_mat_for_fock2,
                den_mat_for_fock2, i, j, fock_type, exchange_scaling_factor,
                0.0, thresh_int)

            # 'X_X', 'X_Y', 'X_Z', 'Y_X', 'Y_Y', 'Y_Z', 'Z_X', 'Z_Y', 'Z_Z'
            xy_pairs = [(x, y) for x in range(3) for y in range(3)]

            for idx, (x, y) in enumerate(xy_pairs):
                hessian_2nd_order_derivatives[i, j, x, y] += fock_factor * (
                    fock_hess_1100[idx] + fock_hess_1010[idx])

            # lower triangle is transpose of the upper part
            if i != j:
                hessian_2nd_order_derivatives[j, i] += (
                    hessian_2nd_order_derivatives[i, j].T)

        hessian_2nd_order_derivatives = self.comm.reduce(
            hessian_2nd_order_derivatives, root=mpi_master())
        profiler.stop_timer('2e mixed')
        # DFT:
        if self._dft:
            profiler.start_timer('DFT')
            grid_drv = GridDriver(self.comm)
            grid_level = (get_default_grid_level(self.scf_driver.xcfun)
                          if self.scf_driver.grid_level is None else
                          self.scf_driver.grid_level)

            # determine grid level for XC Hessian
            grid_level = self.determine_xc_hessian_grid_level(
                molecule, grid_level)

            grid_drv.set_level(grid_level)
            mol_grid = grid_drv.generate(molecule)

            xc_mol_hess = XCMolecularHessian()
            hessian_dft_xc = xc_mol_hess.integrate_exc_hessian(
                molecule, ao_basis, [density], mol_grid,
                self.scf_driver.xcfun.get_func_label())
            hessian_dft_xc = self.comm.reduce(hessian_dft_xc, root=mpi_master())
            profiler.stop_timer('DFT')

        # nuclei-point charges contribution
        if self.scf_driver.point_charges is not None:

            hessian_point_charges = np.zeros((natm, natm, 3, 3))

            qm_coords = molecule.get_coordinates_in_bohr()
            nuclear_charges = molecule.get_element_ids()

            for i in range(self.rank, natm, self.nodes):
                if atom_pairs is not None:
                    if i not in atoms:
                        continue

                q_i = nuclear_charges[i]
                xyz_i = qm_coords[i]

                for j in range(len(mm_charges)):
                    q_j = mm_charges[j]
                    xyz_j = mm_coords[j]

                    vec_ij = xyz_j - xyz_i
                    rij = np.linalg.norm(vec_ij)

                    hess_ii = q_i * q_j * (3.0 * np.outer(vec_ij, vec_ij) /
                                           rij**5 - np.eye(3) / rij**3)

                    hessian_point_charges[i, i, :, :] += hess_ii

            if self.scf_driver.qm_vdw_params is not None:

                for i in range(self.rank, natm, self.nodes):
                    if atom_pairs is not None:
                        if i not in atoms:
                            continue

                    xyz_i = qm_coords[i]
                    sigma_i = self.scf_driver.qm_vdw_params[i, 0]
                    epsilon_i = self.scf_driver.qm_vdw_params[i, 1]

                    for j in range(len(mm_charges)):
                        xyz_j = self.scf_driver.point_charges[:3, j]
                        sigma_j = self.scf_driver.point_charges[4, j]
                        epsilon_j = self.scf_driver.point_charges[5, j]

                        vec_ij = xyz_j - xyz_i

                        # convert vector to nm
                        vec_ij *= bohr_in_angstrom() * 0.1

                        epsilon_ij = np.sqrt(epsilon_i * epsilon_j)
                        sigma_ij = 0.5 * (sigma_i + sigma_j)

                        if epsilon_ij == 0.0:
                            continue

                        rij = np.linalg.norm(vec_ij)
                        inv_r2 = 1.0 / rij**2
                        inv_r4 = inv_r2**2

                        sigma_r_6 = (sigma_ij / rij)**6
                        sigma_r_12 = sigma_r_6**2

                        coef_rr = (28.0 * sigma_r_12 - 8.0 * sigma_r_6) * inv_r4
                        coef_I = (2.0 * sigma_r_12 - sigma_r_6) * inv_r2

                        hess_ii = 24.0 * epsilon_ij * (
                            coef_rr * np.outer(vec_ij, vec_ij) -
                            coef_I * np.eye(3))

                        # convert hessian to atomic unit
                        hess_ii /= (hartree_in_kjpermol() *
                                    (10.0 / bohr_in_angstrom())**2)

                        hessian_point_charges[i, i, :, :] += hess_ii

            hessian_point_charges = self.comm.reduce(hessian_point_charges,
                                                     root=mpi_master())

        if self.rank == mpi_master():
            profiler.start_timer('Finalisation')
            hessian_nuclear_nuclear = self.hess_nuc_contrib(molecule)

            # Doing this post-hoc is much easier to implement, and the cost of the nuclear nuclear contribution is neglible
            if atom_pairs is not None:
                for i in range(natm):
                    for j in range(natm):
                        if (i, j) not in all_atom_pairs and (
                                j, i) not in all_atom_pairs:
                            hessian_nuclear_nuclear[i, j, :, :] = 0.0
                            if self._dft:
                                hessian_dft_xc[i * 3:i * 3 + 3,
                                               j * 3:j * 3 + 3] = 0.0

            self.hessian = (
                hessian_first_order_derivatives +
                hessian_2nd_order_derivatives.transpose(0, 2, 1, 3) +
                hessian_nuclear_nuclear.transpose(0, 2, 1, 3))

            if self.scf_driver.point_charges is not None:
                self.hessian += hessian_point_charges.transpose(0, 2, 1, 3)

            self.hessian = self.hessian.reshape(natm * 3, natm * 3)

            # add pe contr to hessian
            if self.scf_driver._pe:
                from .embedding import PolarizableEmbeddingHess
                embedding_drv = PolarizableEmbeddingHess(
                    molecule=molecule,
                    ao_basis=ao_basis,
                    options=self.scf_driver.embedding,
                    density=density * 2,
                    comm=self.comm)
                self.hessian += embedding_drv.compute_pe_energy_hess_contributions(
                    density_matrix=density * 2)

            if self._dft:
                self.hessian += hessian_dft_xc
            profiler.stop_timer('Finalisation')
            profiler.stop_timer('total')

        self.ostream.print_info('Second order derivative contributions' +
                                ' to the Hessian computed in' +
                                ' {:.2f} sec.'.format(tm.time() - t2))
        self.ostream.print_blank()
        self.ostream.flush()

        # Calculate the gradient of the dipole moment for IR intensities
        if self.do_dipole_gradient:
            self.compute_dipole_gradient(molecule, ao_basis, dist_cphf_ov)

    def determine_xc_hessian_grid_level(self, molecule, grid_level):
        """
        Determines XC Hessian grid level.

        :param molecule:
            The molecule.
        :param grid_level:
            The input grid_level.

        :return:
            The recommended grid_level for XC Hessian.
        """

        # make sure to use high grid level for DFT Hessian
        # see J. Chem. Phys. 119, 12763-12768 (2003)

        # determine grid level for XC Hessian based on default_grid_level
        # and max_elem_id

        default_grid_level = get_default_grid_level(self.scf_driver.xcfun)
        elem_ids = molecule.get_identifiers()
        max_elem_id = max(elem_ids)

        errmsg = 'Hessian calculation with '
        errmsg += self.scf_driver.xcfun.get_func_label().upper()
        errmsg += f' functional and max element id {max_elem_id} '
        errmsg += 'is not supported.'

        if default_grid_level <= 4:
            if max_elem_id <= 18:
                grid_level = max(6, grid_level)
            elif max_elem_id <= 36:
                grid_level = max(7, grid_level)
            else:
                assert_msg_critical(False, errmsg)

        elif default_grid_level in [5, 6]:
            if max_elem_id <= 10:
                grid_level = max(6, grid_level)
            elif max_elem_id <= 18:
                grid_level = max(7, grid_level)
            else:
                assert_msg_critical(False, errmsg)

        else:
            assert_msg_critical(False, errmsg)

        return grid_level

    def compute_dipole_gradient(self, molecule, ao_basis, dist_cphf_ov):
        """
        Computes the analytical gradient of the dipole moment.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        """

        # Number of atoms and atomic charges
        natm = molecule.number_of_atoms()
        nuclear_charges = molecule.get_element_ids()

        # Dipole integrals
        dipole_mats = compute_electric_dipole_integrals(molecule, ao_basis,
                                                        [0.0, 0.0, 0.0])
        # Note: compute_dipole_gradient uses r instead of mu for dipole operator
        dipole_ints = (
            -1.0 * dipole_mats[0],
            -1.0 * dipole_mats[1],
            -1.0 * dipole_mats[2],
        )

        if self.rank == mpi_master():
            # Initialize a local dipole gradient to zero
            dipole_gradient = np.zeros((3, natm, 3))

            # Put the nuclear contributions to the right place
            natm_zeros = np.zeros(natm)
            dipole_gradient[0, :, :] = np.vstack(
                (nuclear_charges, natm_zeros, natm_zeros)).T
            dipole_gradient[1, :, :] = np.vstack(
                (natm_zeros, nuclear_charges, natm_zeros)).T
            dipole_gradient[2, :, :] = np.vstack(
                (natm_zeros, natm_zeros, nuclear_charges)).T

            scf_tensors = self.scf_driver.scf_tensors
            density = scf_tensors['D_alpha']
            mo = scf_tensors['C_alpha']
            nocc = molecule.number_of_alpha_electrons()
            mo_occ = mo[:, :nocc]
            mo_vir = mo[:, nocc:]
            nvir = mo_vir.shape[1]

        ovlp_grad_drv = OverlapGeom100Driver()
        dip_grad_drv = ElectricDipoleMomentGeom100Driver()

        for a in range(natm):

            # overlap gradient

            gmats = ovlp_grad_drv.compute(molecule, ao_basis, a)

            ovlp_deriv_ao_a = []
            for x, label in enumerate(['X', 'Y', 'Z']):
                gmat = gmats.matrix_to_numpy(label)
                ovlp_deriv_ao_a.append(gmat + gmat.T)

            gmats = Matrices()

            # perturbed density

            perturbed_density_a = []

            for x in range(3):

                cphf_ov_ax = dist_cphf_ov[a * 3 + x].get_full_vector(0)

                if self.rank == mpi_master():

                    cphf_ov_ax = cphf_ov_ax.reshape(nocc, nvir)

                    perturbed_density_a.append(
                        # mj,xyij,ni->xymn
                        -np.linalg.multi_dot(
                            [density, ovlp_deriv_ao_a[x], density])
                        # ma,xyia,ni->xymn
                        + np.linalg.multi_dot([mo_vir, cphf_ov_ax.T, mo_occ.T])
                        # mi,xyia,na->xymn
                        + np.linalg.multi_dot([mo_occ, cphf_ov_ax, mo_vir.T]))

            # dipole gradient

            gmats_dip = dip_grad_drv.compute(molecule, ao_basis,
                                             [0.0, 0.0, 0.0], a)

            gmats_dip_components = [
                'X_X', 'X_Y', 'X_Z', 'Y_X', 'Y_Y', 'Y_Z', 'Z_X', 'Z_Y', 'Z_Z'
            ]

            comp_to_idx = {'X': 0, 'Y': 1, 'Z': 2}

            for i, label in enumerate(gmats_dip_components):
                gmat_dip = gmats_dip.matrix_to_numpy(label)

                icoord = comp_to_idx[label[0]]  # atom coordinate component
                icomp = comp_to_idx[label[-1]]  # dipole operator component

                if self.rank == mpi_master():

                    c, x = icomp, icoord

                    dipole_gradient[c, a, x] += -2.0 * (
                        np.sum(density * (gmat_dip + gmat_dip.T)) +
                        np.sum(perturbed_density_a[x] * dipole_ints[c]))

        if self.rank == mpi_master():
            self.dipole_gradient = dipole_gradient.reshape(3, 3 * natm)

    def compute_orbital_response(self, molecule, ao_basis):
        """
        TEMPORARY FUNCTION FOR PERFORMANCE TESTING
        Computes the CPHF orbital response.

        :param molecule:
            The molecule.
        :param ao_basis:
            Tha AO basis.
        """

        # get SCF tensors
        scf_tensors = self.scf_driver.scf_tensors

        # Set up a CPHF solver
        cphf_solver = HessianOrbitalResponse(self.comm, self.ostream)
        cphf_solver.update_settings(self.cphf_dict, self.method_dict)

        # Solve the CPHF equations
        cphf_solver.compute(molecule, ao_basis, scf_tensors, self.scf_driver)

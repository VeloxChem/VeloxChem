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
from copy import deepcopy
import numpy as np
import time
import math

from .veloxchemlib import AODensityMatrix, denmat
from .veloxchemlib import DenseMatrix
from .veloxchemlib import GradientScreeningData, GpuDevices
from .veloxchemlib import XCFunctional, MolecularGrid
from .veloxchemlib import compute_overlap_gradient_gpu
from .veloxchemlib import compute_kinetic_energy_gradient_gpu
from .veloxchemlib import compute_nuclear_potential_gradient_gpu
from .veloxchemlib import compute_fock_gradient_gpu, matmul_gpu
from .veloxchemlib import integrate_vxc_gradient_gpu
from .veloxchemlib import parse_xc_func
from .veloxchemlib import mpi_master
from .veloxchemlib import bohr_in_angstrom
from .molecularbasis import MolecularBasis
from .outputstream import OutputStream
from .gradientdriver import GradientDriver
from .errorhandler import assert_msg_critical
from .sanitychecks import (molecule_sanity_check, scf_results_sanity_check,
                           dft_sanity_check)


class ScfGradientDriver(GradientDriver):
    """
    Implements SCF gradient driver.

    :param scf_drv:
        The SCF driver.

    Instance variables
        - flag: The driver flag.
        - delta_h: The displacement for finite difference.
    """

    def __init__(self, scf_drv):
        """
        Initializes gradient driver.
        """

        super().__init__(scf_drv.comm, scf_drv.ostream)

        self.scf_driver = scf_drv
        self.flag = 'SCF Gradient Driver'

        self.eri_thresh = scf_drv.eri_thresh
        self.pair_thresh = scf_drv.pair_thresh
        self.density_thresh = scf_drv.density_thresh
        self.prelink_thresh = scf_drv.prelink_thresh

        self.timing = scf_drv.timing

        self.dispersion = scf_drv.dispersion

    def compute(self, molecule, basis, scf_results=None):
        """
        Performs calculation of gradient.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_results:
            The dictionary containing converged SCF results.
        """

        if scf_results is None:
            scf_results = self.scf_driver.scf_tensors

        start_time = time.time()
        self.print_header()

        if self.numerical:
            self.compute_numerical(molecule, basis, scf_results)

        else:
            # sanity checks
            molecule_sanity_check(molecule)
            scf_results_sanity_check(self, self.scf_driver.scf_tensors)
            dft_sanity_check(self, 'compute')

            if self.rank == mpi_master():
                scf_type = scf_results['scf_type']
            else:
                scf_type = None
            scf_type = self.comm.bcast(scf_type, root=mpi_master())

            if scf_type == 'restricted':
                self.compute_analytical_restricted(molecule, basis, scf_results)
            else:
                assert_msg_critical(
                    False,
                    'ScfGradientDriver: Not implemented for open-shell')

        # print gradient
        self.print_geometry(molecule)
        self.print_gradient(molecule)

        valstr = '*** Time spent in gradient calculation: '
        valstr += '{:.2f} sec ***'.format(time.time() - start_time)
        self.ostream.print_header(valstr)
        self.ostream.print_blank()
        self.ostream.flush()

    def compute_analytical_restricted(self, molecule, basis, scf_results):
        """
        Performs calculation of gradient for restricted SCF.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_results:
            The dictionary containing converged SCF results.
        """

        grad_timing = {
            'Screening': 0.0,
            'Overlap_grad': 0.0,
            'Kinetic_energy_grad': 0.0,
            'Nuclear_potential_grad': 0.0,
            'Point_charges_grad': 0.0,
            'Fock_grad': 0.0,
            'XC_grad': 0.0,
            'PE_grad': 0.0,
            'CPCM_grad': 0.0,
            'Classical': 0.0,
        }

        if self.rank == mpi_master():
            D = scf_results['D_alpha']
            nocc = molecule.number_of_alpha_electrons()
            ene_occ = scf_results['E_alpha'][:nocc]
            mo_occ = scf_results['C_alpha'][:, :nocc].copy()
            # multi_dot([mo_occ, np.diag(ene_occ), mo_occ.T])
            W = matmul_gpu(ene_occ * mo_occ, mo_occ.T)
            naos = D.shape[0]
        else:
            naos = None
        naos = self.comm.bcast(naos, root=mpi_master())

        if self.rank != mpi_master():
            D = np.zeros((naos, naos))
            W = np.zeros((naos, naos))

        self.comm.Bcast(D, root=mpi_master())
        self.comm.Bcast(W, root=mpi_master())

        dmat = AODensityMatrix([D], denmat.rest)
        wmat = DenseMatrix(W)

        natoms = molecule.number_of_atoms()

        self.gradient = np.zeros((natoms, 3))

        num_gpus_per_node = self.scf_driver._get_num_gpus_per_node()

        # screening for gradient

        t0 = time.time()

        grad_screener = GradientScreeningData(
            molecule, basis, dmat, wmat, num_gpus_per_node, self.pair_thresh,
            self.density_thresh, self.rank, self.nodes)

        grad_timing['Screening'] += time.time() - t0

        # kinetic energy contribution to gradient

        t0 = time.time()

        T_grad = compute_kinetic_energy_gradient_gpu(molecule, basis, grad_screener, self.rank, self.nodes)
        T_grad = T_grad.to_numpy()
        T_grad *= 2.0

        self.gradient += T_grad

        grad_timing['Kinetic_energy_grad'] += time.time() - t0

        # nuclear potential contribution to gradient

        t0 = time.time()

        V_grad = compute_nuclear_potential_gradient_gpu(molecule, basis, grad_screener, self.rank, self.nodes)
        V_grad = V_grad.to_numpy()
        V_grad *= 2.0

        self.gradient += V_grad

        grad_timing['Nuclear_potential_grad'] += time.time() - t0

        # TODO: point charges contribution

        t0 = time.time()

        if self.scf_driver._point_charges is not None:
            npoints = self.scf_driver._point_charges.shape[1]

            mm_coords = []
            mm_charges = []
            for p in range(npoints):
                xyz_p = self.scf_driver._point_charges[:3, p]
                chg_p = self.scf_driver._point_charges[3, p]
                mm_coords.append(xyz_p.copy())
                mm_charges.append(chg_p)

            assert_msg_critical(
                False,
                'ScfGradientDriver: Point charges not yet implemented')

            grad_timing['Point_charges_grad'] += time.time() - t0

        # orbital contribution to gradient
    
        t0 = time.time()

        S_grad = compute_overlap_gradient_gpu(molecule, basis, grad_screener, self.rank, self.nodes)
        S_grad = S_grad.to_numpy()
        S_grad *= -2.0

        self.gradient += S_grad
        
        grad_timing['Overlap_grad'] += time.time() - t0

        # ERI contribution to gradient

        if self.rank == mpi_master():
            use_dft = 'xcfun' in scf_results
        else:
            use_dft = False
        use_dft = self.comm.bcast(use_dft, root=mpi_master())

        # determine fock_type and exchange_scaling_factor
        if use_dft:
            if self.rank == mpi_master():
                xcfun = scf_results['xcfun']
            else:
                xcfun = None
            xcfun = self.comm.bcast(xcfun, root=mpi_master())
            xcfun = parse_xc_func(xcfun)

            if xcfun.is_hybrid():
                fock_type = '2jkx'
                exchange_scaling_factor = xcfun.get_frac_exact_exchange()
            else:
                fock_type = 'j'
                exchange_scaling_factor = 0.0
        else:
            fock_type = '2jk'
            exchange_scaling_factor = 1.0

        # further determine exchange_scaling_factor, erf_k_coef and omega
        need_omega = (use_dft and xcfun.is_range_separated())
        if need_omega:
            # range-separated hybrid
            full_k_coef = (xcfun.get_rs_alpha() + xcfun.get_rs_beta())
            erf_k_coef = -xcfun.get_rs_beta()
            omega = xcfun.get_rs_omega()
        else:
            full_k_coef = (xcfun.get_frac_exact_exchange() if use_dft else 1.0)
            erf_k_coef = None
            omega = None

        t0 = time.time()

        fock_grad = compute_fock_gradient_gpu(
            molecule, basis, dmat, 2.0, full_k_coef, 0.0, 'symm',
            self.eri_thresh, self.prelink_thresh, grad_screener,
            self.rank, self.nodes)

        self.gradient += fock_grad.to_numpy()

        # TODO: compute range-separted fock_grad in one shot
        if need_omega:
            fock_rs_grad = compute_fock_gradient_gpu(
                molecule, basis, dmat, 0.0, erf_k_coef, omega, 'symm',
                self.eri_thresh, self.prelink_thresh, grad_screener,
                self.rank, self.nodes)

            self.gradient += fock_rs_grad.to_numpy()

        grad_timing['Fock_grad'] += time.time() - t0

        # XC contribution to gradient

        t0 = time.time()

        if use_dft:
            if self.rank == mpi_master():
                xcfun_label = scf_results['xcfun']
            else:
                xcfun_label = None
            xcfun_label = self.comm.bcast(xcfun_label, root=mpi_master())

            xc_grad = integrate_vxc_gradient_gpu(molecule, basis, dmat, dmat, self.scf_driver._mol_grid, xcfun_label, num_gpus_per_node, self.rank, self.nodes)

            self.gradient += xc_grad.to_numpy()

            grad_timing['XC_grad'] += time.time() - t0

        else:
            xcfun_label = 'hf'

        # Embedding contribution to the gradient

        t0 = time.time()

        if self.scf_driver._pe:
            from .embedding import PolarizableEmbeddingGrad

            # pass along emb object from scf, or make a new one? -> for ind dipoles.
            self._embedding_drv = PolarizableEmbeddingGrad(
                molecule=molecule,
                ao_basis=basis,
                options=self.scf_driver.embedding,
                comm=self.comm)

            pe_grad = self._embedding_drv.compute_pe_contributions(
                density_matrix=2.0 * D)

            if self.rank == mpi_master():
                self.gradient += pe_grad

            grad_timing['PE_grad'] += time.time() - t0

        # nuclear contribution to gradient
        # and D4 dispersion correction if requested
        # (only added on master rank)

        t0 = time.time()

        if self.rank == mpi_master():
            self.gradient += self.grad_nuc_contrib(molecule)

            grad_timing['Classical'] += time.time() - t0

            t0 = time.time()

            if self.dispersion or (
                    self.scf_driver._dft and
                    'D4' in self.scf_driver.xcfun.get_func_label().upper()):
                assert_msg_critical(
                    False,
                    'Dispersion correction not yet implemented for the GPU version.')

            t0 = time.time()

        # nuclei-point charges contribution to gradient

        t0 = time.time()

        if self.scf_driver._point_charges is not None:

            assert_msg_critical(
                False,
                'ScfGradientDriver: Point charges not yet implemented')

            coords = molecule.get_coordinates_in_bohr()
            nuclear_charges = molecule.get_element_ids()
            npoints = self.scf_driver._point_charges.shape[1]

            for a in range(self.rank, natoms, self.nodes):
                z_a = nuclear_charges[a]
                r_a = coords[a]

                for p in range(npoints):
                    r_p = self.scf_driver._point_charges[:3, p]
                    q_p = self.scf_driver._point_charges[3, p]
                    r = np.linalg.norm(r_a - r_p)
                    f_ij = z_a * q_p * (r_p - r_a) / r**3

                    self.gradient[a] += f_ij

            if self.scf_driver.qm_vdw_params is not None:
                vdw_grad = np.zeros((natoms, 3))

                for a in range(self.rank, natoms, self.nodes):
                    xyz_i = coords[a]
                    sigma_i = self.scf_driver.qm_vdw_params[a, 0]
                    epsilon_i = self.scf_driver.qm_vdw_params[a, 1]

                    for p in range(npoints):
                        xyz_j = self.scf_driver._point_charges[:3, p]
                        sigma_j = self.scf_driver._point_charges[4, p]
                        epsilon_j = self.scf_driver._point_charges[5, p]

                        r_ij = xyz_j - xyz_i
                        distance_ij = np.linalg.norm(r_ij)
                        n_ij = r_ij / distance_ij

                        # bohr to nm
                        distance_ij *= bohr_in_angstrom() * 0.1

                        epsilon_ij = np.sqrt(epsilon_i * epsilon_j)
                        sigma_ij = 0.5 * (sigma_i + sigma_j)

                        sigma_r_6 = (sigma_ij / distance_ij)**6
                        sigma_r_12 = sigma_r_6**2

                        g = -24.0 * epsilon_ij * (2.0 * sigma_r_12 / distance_ij
                                                  - sigma_r_6 / distance_ij)

                        vdw_grad[a] += -g * n_ij

                # convert gradient to atomic unit
                vdw_grad /= (hartree_in_kjpermol() * 10.0 / bohr_in_angstrom())

                self.gradient += vdw_grad

        grad_timing['Classical'] += time.time() - t0

        # collect gradient

        self.gradient = self.comm.allreduce(self.gradient, op=MPI.SUM)

        if self.timing and self.rank == mpi_master():
            self.ostream.print_info('Gradient timing decomposition')
            for key, val in grad_timing.items():
                if val > 0.0:
                    self.ostream.print_info(f'    {key:<25}:  {val:.2f} sec')
            self.ostream.print_blank()

    def compute_energy(self, molecule, ao_basis, scf_results=None):
        """
        Computes the energy at current geometry.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param scf_results:
            The dictionary containing converged SCF results.

        :return:
            The energy.
        """

        if self.numerical:
            # disable restarting scf for numerical gradient
            self.scf_driver.restart = False
        else:
            # always try restarting scf for analytical gradient
            self.scf_driver.restart = True

        self.scf_driver.ostream.mute()
        new_scf_results = self.scf_driver.compute(molecule, ao_basis)
        self.scf_driver.ostream.unmute()

        assert_msg_critical(self.scf_driver.is_converged,
                            'ScfGradientDriver: SCF did not converge')

        if (self.rank == mpi_master()) and (scf_results is not None):
            scf_results.update(new_scf_results)

        return self.scf_driver.get_scf_energy()

    def __deepcopy__(self, memo):
        """
        Implements deepcopy.

        :param memo:
            The memo dictionary for deepcopy.

        :return:
            A deepcopy of self.
        """

        new_grad_drv = ScfGradientDriver(self.comm, self.ostream)

        for key, val in vars(self).items():
            if isinstance(val, (MPI.Intracomm, OutputStream)):
                pass
            elif isinstance(val, XCFunctional):
                new_grad_drv.key = XCFunctional(val)
            elif isinstance(val, MolecularGrid):
                new_grad_drv.key = MolecularGrid(val)
            else:
                new_grad_drv.key = deepcopy(val)

        return new_grad_drv

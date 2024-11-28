#
#                           VELOXCHEM 1.0-RC3
#         ----------------------------------------------------
#                     An Electronic Structure Code
#
#  Copyright © 2018-2022 by VeloxChem developers. All rights reserved.
#  Contact: https://veloxchem.org/contact
#
#  SPDX-License-Identifier: LGPL-3.0-or-later
#
#  This file is part of VeloxChem.
#
#  VeloxChem is free software: you can redistribute it and/or modify it under
#  the terms of the GNU Lesser General Public License as published by the Free
#  Software Foundation, either version 3 of the License, or (at your option)
#  any later version.
#
#  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
#  License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.

from mpi4py import MPI
import numpy as np
import math
import sys

from .veloxchemlib import bohr_in_angstrom, mpi_master
from .veloxchemlib import hartree_in_ev, bohr_in_angstroms
from .oneeints import compute_electric_dipole_integrals
from .subcommunicators import SubCommunicators
from .outputstream import OutputStream
from .errorhandler import assert_msg_critical

from .inputparser import (parse_input, print_keywords, print_attributes,
                          get_random_string_parallel)

class RixsDriver:
    """
    Implements RIXS driver.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables:
        - theta: Angle.
        - gamma_n: Broadening term.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes (CVS-)RIXS driver.
        """

        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(sys.stdout)
            else:
                ostream = OutputStream(None)

        # mpi information
        self.comm = comm
        self.rank = comm.Get_rank()
        self.nodes = comm.Get_size()

        # outputstream
        self.ostream = ostream

        # method settings
        self.nr_ve = None # nr_ce
        # CVS
        self.nr_ce = None # nr_ce
        self.nr_CO = None # nr_ce

        self.theta = 0
        self.gamma_n = .124/hartree_in_ev() # a.u.

        # input keywords
        self.input_keywords = {
            'rixs': {
                'angle': ('float', 'angle between incident polarization vector and propagation vector of outgoing'),
                'gamma': ('float', 'broadening term'),
            },
        }

    def update_settings(self, rixs_dict):
        """
        Updates settings in SCF driver.

        :param rixs_dict:
            The dictionary of rixs input.
        """

        #if method_dict is None:
        #    method_dict = {}

        rixs_keywords = {
            key: val[0] for key, val in self._input_keywords['rixs'].items()
        }

        parse_input(self, rixs_keywords, rixs_dict)

        if 'program_end_time' in rixs_dict:
            self.program_end_time = rixs_dict['program_end_time']
        if 'filename' in rixs_dict:
            self.filename = rixs_dict['filename']

        """
        if self.electric_field is not None:
            assert_msg_critical(
                len(self.electric_field) == 3,
                'SCF driver: Expecting 3 values in \'electric field\' input')
            assert_msg_critical(
                not self._pe,
                'SCF driver: \'electric field\' input is incompatible with ' +
                'polarizable embedding')
            # disable restart of calculation with static electric field since
            # checkpoint file does not contain information about the electric
            # field
            self.restart = False
        """
    
    def sts_tdm(self, molecule, scf_results, tda_res_val, tda_res_core):
        """
        The state-to-state transition-density-matrix

        :param molecule    : molecule object (could be removed w/ small modifications)
        :param ao_basis    : atomic orbital basis (could be removed)
        :param scf_results : results of the SCF calculation
        :tda_res_val       : linear-response results for the valence calculation
        :tda_res_core      : linear-response results for the core calculation

        returns the state-to-state transition-density-matrix in ao-basis
        """
        norb = scf_results['C_alpha'].shape[0]
        nocc = molecule.number_of_alpha_electrons()
        nvir = norb - nocc

        X_val_mo = tda_res_val['eigenvectors'].reshape(nocc,nvir,self.nr_ve)
        X_cor_mo = np.zeros((self.nocc, self.nvir, self.nr_ce))
        X_cor_mo[:self.nr_CO] = tda_res_core['eigenvectors'].reshape(self.nr_CO, self.nvir, self.nr_ce)

        # Electron- and hole density matrices
        gamma_ab = np.einsum('ipJ, iqI -> pqIJ', X_cor_mo, X_val_mo, optimize=True)
        gamma_ij = np.einsum('qaJ, paI -> pqIJ', X_cor_mo, X_val_mo, optimize=True)

        # Transform
        C_ab = scf_results['C_alpha'][:,nocc:]
        C_ij = scf_results['C_alpha'][:,:nocc]
        
        gamma_ab_ao = np.einsum('vp, pqIJ, wq -> vwIJ', C_ab, gamma_ab, C_ab, optimize=True)
        gamma_ij_ao = np.einsum('vp, pqIJ, wq -> vwIJ', C_ij, gamma_ij, C_ij, optimize=True)

        gamma_ao = gamma_ab_ao - gamma_ij_ao
        return gamma_ao

    def gts_tdm(self, molecule, scf_results, tda_res_core):
        """
        The GS-to-state transition-density-matrix (GS to core) --
        essentially transformed and reshaped eigenvector

        :param molecule    : molecule object (could be removed w/ small modifications)
        :param ao_basis    : atomic orbital basis (could be removed)
        :param scf_results : results of the SCF calculation (could be removed w/ small modifications)
        :tda_res_val       : linear-response results for the valence calculation
        :tda_res_core      : linear-response results for the core calculation

        returns the GS-to-state transition-density-matrix in ao-basis
        """
        norb = scf_results['C_alpha'].shape[0]
        nocc = molecule.number_of_alpha_electrons()
        nvir = norb - nocc

        X_cor_mo = np.zeros((self.nocc, self.nvir, self.nr_ce))
        X_cor_mo[:self.nr_CO] = tda_res_core['eigenvectors'].reshape(self.nr_CO, self.nvir, self.nr_ce)

        # Transform
        C_ab = scf_results['C_alpha'][:,nocc:]
        C_ij = scf_results['C_alpha'][:,:nocc]
        gamma_ng_ao = np.einsum('vp, pqJ, wq -> vwJ', C_ij, X_cor_mo, C_ab, optimize=True)

        return gamma_ng_ao


    def F_xy(self, w, f, gamma_ao,
         gamma_ng_ao, tda_res_core, tda_res_val,
         dipole_ints):
        """
        The RIXS scattering amptlitude (sum-over-states transition strength).
        
        :param w           : Energy of incoming photon 
        :param f           : Which, out of the nr_ve, final valence states to end
        :param gamma_factor: Broadening factor, or half-width of the core-excited state; 
                             in the resonant case corresponding to state with eigenvalue w
    
        returns the scattering amplitude/transition strength for all molecular directions and outgoing photon polarisation
        """
        
        eigenvalues = tda_res_core['eigenvalues'] # excitation energies
        e_n = (1 / (eigenvalues - w - self.gamma_n*1j)) # denominator
        scatt_amp = np.einsum('n, xij, ijn, yab, abn -> xy', e_n, -dipole_ints, gamma_ao[:,:,f], -dipole_ints, gamma_ng_ao, optimize='greedy')
        return scatt_amp
    
    def oscillator_strength(self, tdpm, tda_res_core, tda_res_val):
        d_E = [ce - tda_res_val['eigenvalues'] for ce in tda_res_core['eigenvalues']]
        unweighted_f = (2/3) * np.einsum('fnx->nf', tdpm ** 2)
        f = d_E * unweighted_f
        return f #np.array(d_E), f

    def rixs_xsection(self, w, f, gamma_ao,
         gamma_ng_ao, tda_res_core, tda_res_val,
         dipole_ints, w_p=None, theta=0, gamma_factor=.0124/hartree_in_ev()):
        """
        Calculate the RIXS cross-section, sigma.
    
        :param w            : Energy of the incident photon
        :param w_p          : Energy of the outgoing/scattered photon, if None: these are matched with valence energies
        :param f            : Index of the final/target excited state
        :param theta        : Scattering angle 
        :param gamma_factor : Broadening factor for Lorentzian profile
        :param tda_res(_val)  : The results tensor of the TDDFT/TDA ground to valence excited state
        
        returns the RIXS transition intensity for the given frequencies, w (and) and state f
        """
        if w_p == None:
            w_prime = w - tda_res_val['eigenvalues'][f] #w - w_f0 
        else:
            w_prime = w_p
            
        F = self.F_xy(w,f, gamma_ao,
         gamma_ng_ao, tda_res_core, tda_res_val,
         dipole_ints)
        
        sigma = w_prime/w * 1/15 * ((2 - (1/2) * np.sin(theta) ** 2) * np.sum(np.abs(F)**2) 
                   + ((3/4) * np.sin(theta) ** 2 - 1/2) * (np.sum(F * F.T.conj())
                                                   + np.trace(np.abs(F)**2)))
        return sigma

    #def rixs_map():
    def transition_dipole_mom(self, gamma_ao, dipole_integrals):
        T_fn = np.einsum('ijfn,xij->fnx', gamma_ao, -dipole_integrals)
        return T_fn

    def c_v_absenergy(self, tda_res_core, tda_res_val):
        d_E = [ce - tda_res_val['eigenvalues'] for ce in tda_res_core['eigenvalues']]
        return np.array(d_E)

    #def ene_loss()

    def compute(self, molecule, ao_basis, scf_results, tda_res_val, tda_res_core):
        self.norb = scf_results['C_alpha'].shape[0]
        self.nocc = molecule.number_of_alpha_electrons()
        self.nvir = self.norb - self.nocc
            
        gamma = self.sts_tdm(molecule, scf_results, tda_res_val, tda_res_core)
        gamma_gs = self.gts_tdm(molecule, scf_results,tda_res_core)

        dip_mats = compute_electric_dipole_integrals(molecule, ao_basis)
        dipole_ints = -1.0 * np.array([dip_mats[0],
                                       dip_mats[1],
                                       dip_mats[2]])

        omega_f = tda_res_val['eigenvalues']
        omega_n = tda_res_core['eigenvalues']
        photon_energy = omega_n # resonant
        #photon_energy = np.arange(np.min(omega_n)-.05,np.max(omega_n)+.05,.01)

        rixs_map = []
        # Calculate all elements of the full RIXS map, 
        # .real() is taken to handle eventual numerical inconsistencies
        for i, en_n in enumerate(photon_energy):
            for j in range(self.nr_ve):
                rixs_map.append(np.array([omega_f[j], en_n, self.rixs_xsection(en_n, j, gamma,
         gamma_gs, tda_res_core, tda_res_val,
         dipole_ints).real]))
        rixs_map = np.array(rixs_map)
        self.rixs_map = rixs_map
        
        T_fn = self.transition_dipole_mom(gamma, dipole_ints)
        f_f = self.oscillator_strength(T_fn, tda_res_core, tda_res_val)
        self.f_f = f_f
        emission_ene = self.c_v_absenergy(tda_res_core, tda_res_val)
        self.ene_loss_map = np.array([ce - emission_ene[i] for i,ce in enumerate(tda_res_core['eigenvalues'])])

    



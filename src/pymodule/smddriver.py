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
from io import StringIO
from contextlib import redirect_stderr
from pathlib import Path
import numpy as np
import sys
from itertools import combinations
import rdkit
import time


from .veloxchemlib import mpi_master, hartree_in_kjpermol, bohr_in_angstrom
from .outputstream import OutputStream
from .smdsolventproperties import get_smd_solvent_properties, get_sigma_properties, get_rzz_parameters
from .molecule import Molecule
from .molecularbasis import MolecularBasis
from .scfrestdriver import ScfRestrictedDriver
from .errorhandler import assert_msg_critical

## Note: Maybe adopt a keyword in scf_drv saying smd = True or similar, since then the only 
# thing extra to do is collecting the new parameters --> i.e., waking up smd VIA scf_drv, rather than
# the other way around 

class SMDDriver:
    """
    SMD (Solvation Model based on Density) calculations.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - 
    
    """

    def __init__(self, comm=None, ostream=None):

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

        # theory settings
        self.basis_set = '6-31G**'
        self.xcfun = 'B3LYP'
        self.ri_coulomb = False
        
        # solvent properties
        self.solute = None
        self.atom_radii = None
        self.cpcm_grid_per_sphere = 194
        self.smd_solvent_parameters = get_smd_solvent_properties() 
        self.sigma_param, self.sigma_water = get_sigma_properties()
        self.r_zz_dict = get_rzz_parameters()
        self.calculate_tanh = False

    def compute_solvation_energy(self, molecule, solvent='water'):
        """
        Perform SMD solvation.
        """

        self.solute = molecule
        self.solvent = solvent
        
        assert_msg_critical(
            self.solvent in self.smd_solvent_parameters,
            'Solvent {self.solvent} not available')
        
        self.smd_solvent_parameters = self.smd_solvent_parameters[solvent]
        self._get_intrinsic_coulomb_radii() # or execute this within the ENP function (if not used elsewhere)

        self.ostream.print_info("Computing the ENP Contribution...")
        self.ostream.flush()
        ENP_energy = self._get_ENP_contribution()
        
        self.ostream.print_info("Computing the CDS Contribution...")
        self.ostream.flush()
        CDS_energy = self._get_CDS_contribution()
        
        # DG_conc = DG_conc(gas-phase) + DG_conc(solv-phase) (=0 if 1 mol/L, =1.89 kcal/mol if gas-phase std state of 1 atm)

        # DG_solv = DG_ENP + DG_CDS (+ DG_conc)
        # N contribution in ENP: diff in total gas phase energy between gas-equil.geometry and liq-equil.geometry
        final_E = ENP_energy + CDS_energy 
        

        return final_E
    
    def _get_intrinsic_coulomb_radii(self):
        """
        Get the intrinsic Coulomb radii for the solute.

        """

        atom_labels = self.solute.get_labels() 
        bondi_radii = self.solute.vdw_radii_to_numpy()

        # Table 3 in SMD paper. 
        # If not in smd_elements, use Bondi radius (default in vlx); if not in Bondi, use 2.0 Å (also default in vlx)
        smd_elements = {
            'H':1.20,'C':1.85,'N':1.89,
            'O':1.52,'F':1.73,'Si':2.47,
            'P':2.12,'S':2.49,'Cl':2.38,
            'Br':3.06
        }

        alpha = self.smd_solvent_parameters['alpha']
        atom_radii = []

        for i, atom in enumerate(atom_labels):
            atom_radii.append(atom)
            
            if atom not in smd_elements.keys(): # Using Bondi Radius
                atom_radii.append(bondi_radii.tolist()[i])
            else:
                # re-parametrize O radius dep. on the solvent
                if atom == 'O':
                    atom_radii.append(1.52 if alpha < 0.43 else round((1.52 + 1.8*(0.43 - alpha)), 2))
                else:
                    atom_radii.append(smd_elements[atom])
                
        self.atom_radii = atom_radii

    def _get_ENP_contribution(self):
        """
        Get the ENP contribution to the solvation energy.
        """
        # Electrostatic contribution: DG_ENP = DG_EP if E_tot(gas-geometry in gas-phase) 
        # = E_tot(solv-geometry in gas-phase) (not the case in general for SMD model, 
        # but assumed for the paper)
        # --> DG_EP calculated fr MO self-consistent reaction field solver
        
        # TODO: check the cpcm-implemented definition of the cavity (with Erik)
        # TODO: Add N contribution? --> E-diff in gasphase E between gas-phase and solv-phase-optimized geometries 

        basis_set = self.basis_set
        basis = MolecularBasis.read(self.solute, basis_set)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.xcfun = self.xcfun
        scf_drv.ri_coulomb = self.ri_coulomb
        scf_drv.solvation_model = 'cpcm'
        scf_drv.cpcm_epsilon = self.smd_solvent_parameters['epsilon']
        scf_drv.cpcm_grid_per_sphere = self.cpcm_grid_per_sphere
        scf_drv.cpcm_custom_vdw_radii = self.atom_radii
        
        scf_drv.compute(self.solute, basis)
        ENP_energy = scf_drv.cpcm_epol * hartree_in_kjpermol()

        self.ostream.print_info(f"ENP Contribution: {ENP_energy:.4f} kJ/mol")
        self.ostream.flush()

        return ENP_energy

    #REMOVE:
    def _get_atomic_surface_tension_parameters(self):
        
        # Table 4:  Any possible surface tension parameter that is not in this table is set equal to zero in SMD. 
        # For example, there is no surface tension on P atoms in SMD. cal/molÅ^2
        
        param = self.smd_solvent_parameters
        n, alpha, beta = param['n'], param['alpha'], param['beta']
  
        # Create a sigma_i (tilde) 
        atoms = self.atom_radii[0::2]
        sigma_params = self.sigma_water if self.solvent == 'water' else self.sigma_param

        sigma_i = {}

        # sigma_k (single atom contribution)
        for atom in atoms:
            params = sigma_params.get((atom,))
            if params:
                if self.solvent == 'water':
                    sigma_i[atom] = params
                else:
                    s_n, s_alpha, s_beta = params
                    sigma_i[atom] = n * s_n + alpha * s_alpha + beta * s_beta
        
        # sigma_kk' (pair contributions) 
        for i, j in combinations(range(len(atoms)), 2):
            pair = (atoms[i], atoms[j])
            rev_pair = (atoms[j], atoms[i])
            
            if pair in sigma_params:
                params = sigma_params[pair]
            elif rev_pair in sigma_params:
                params = sigma_params[rev_pair]
            else:
                continue
            
            if self.solvent == 'water':
                sigma_i[pair] = params
            else:
                s_n, s_alpha, s_beta = params
                sigma_i[pair] = n * s_n + alpha * s_alpha + beta * s_beta
        
        self.sigma_i = sigma_i

    def _calculate_sigma_k(self):        
        # Table 4:  Any possible surface tension parameter that is not in this table is set equal to zero in SMD. 
        # For example, there is no surface tension on P atoms in SMD. cal/molÅ^2
        param = self.smd_solvent_parameters
        n, alpha, beta = param['n'], param['alpha'], param['beta']
        
        sigma_params = self.sigma_water if self.solvent == 'water' else self.sigma_param
        
        distance_matrix = self.solute.get_distance_matrix_in_angstrom()
        atoms = self.solute.get_labels()
        num_atoms = len(atoms)

        halogens = ['F', 'Si', 'S', 'Cl', 'Br'] #if Zk in [F, Si, S, Cl, Br], sigma_k = sigma_tilde_Z_k 
        all_smd_atoms = ['H', 'C', 'N', 'O'] + halogens #Z_k != H, C, N, O, F, Si, S, Cl, or Br: sigma_k = 0

        sigma_k = {}
        r_ZZ_param = self.r_zz_dict

        for k in range(num_atoms):

            Zk = atoms[k]
            sigma = 0.0

            if Zk not in all_smd_atoms:
                sigma_k[k] = 0.0
                continue

            # atomic surface tension
            if Zk in halogens:
                if self.solvent == 'water':
                    sigma += sigma_params.get((Zk,), 0.0)
                else:
                    s_n, s_alpha, s_beta = sigma_params.get((Zk,), (0.0, 0.0, 0.0))
                    sigma += n * s_n + alpha * s_alpha + beta * s_beta
                
                sigma_k[k] = sigma
                continue

            if self.solvent == 'water':
                sigma += sigma_params.get((Zk,), 0.0)
            else:
                s_n, s_alpha, s_beta = sigma_params.get((Zk,), (0.0, 0.0, 0.0))
                sigma += n * s_n + alpha * s_alpha + beta * s_beta

            # Initialize sums depending on atom type 
            if Zk == 'H':
                # Contributions from C and O neighbors
                for k2 in range(num_atoms):
                    if k == k2:
                        continue
                    Zk2 = atoms[k2]
                    if Zk2 not in ['C', 'O']:
                        continue
                    key = (Zk, Zk2)
                    if key not in r_ZZ_param or key not in sigma_params:
                        continue
                    R = distance_matrix[k, k2]
                    r_zz, dr_zz = r_ZZ_param[key]            
                    if R > r_zz + dr_zz:
                        tanh = np.exp(dr_zz / (R - dr_zz - r_zz))
                        if self.solvent == 'water':
                            sigma += sigma_params[key] * tanh
                        else:
                            s_n, s_alpha, s_beta = sigma_params[key]
                            sigma += (n * s_n + alpha * s_alpha + beta * s_beta) * tanh

            elif Zk == 'C':
                tanh_sum = 0.0
                for k2 in range(num_atoms):
                    if k == k2:
                        continue
                    Zk2 = atoms[k2]
                    key = (Zk, Zk2)
                    if key not in r_ZZ_param or key not in sigma_params:
                        continue
                    R = distance_matrix[k, k2]
                    r_zz, dr_zz = r_ZZ_param[key]
                    if R > r_zz + dr_zz:
                        tanh = np.exp(dr_zz / (R - dr_zz - r_zz))
                        if Zk2 == 'C':
                            if self.solvent == 'water':
                                sigma += sigma_params[key] * tanh
                            else:
                                s_n, s_alpha, s_beta = sigma_params[key]
                                sigma += (n * s_n + alpha * s_alpha + beta * s_beta) * tanh
                        else:
                            tanh_sum += tanh
                
                if self.solvent == 'water':
                    sigma += sigma_params.get(('C','N'), 0.0) * tanh_sum ** 2
                else:
                    s_n, s_alpha, s_beta = sigma_params.get(('C','N'), 0.0)
                    sigma += (n * s_n + alpha * s_alpha + beta * s_beta) * tanh_sum ** 2

            elif Zk == 'N':
                nested_cn_sum = 0.0
                nc3_sum = 0.0
                
                for k_prime in range(num_atoms):
                    if k_prime == k:
                        continue
                    Zk_prime = atoms[k_prime]
                    key_k_kp = ('N', Zk_prime)

                    # NC(3) term
                    if Zk_prime == 'C(3)' and key_k_kp in r_ZZ_param:
                        R_k_kp = distance_matrix[k, k_prime]
                        r_zz, dr_zz = r_ZZ_param[key_k_kp]
                        if R_k_kp > r_zz + dr_zz:
                            tanh_k_kp = np.exp(dr_zz / (R_k_kp - dr_zz - r_zz))
                            nc3_sum += tanh_k_kp
                    
                    # NC + C-X nested sum
                    elif Zk_prime == 'C':
                        R_k_kp = distance_matrix[k,k_prime]
                        r_nc, dr_nc = r_ZZ_param[('N', 'C')]
                        if R_k_kp <= r_nc + dr_nc:
                            continue

                        tanh_k_kp = np.exp(dr_nc / (R_k_kp - dr_nc - r_nc))

                        # inner sum over k'k"
                        inner_sum = 0.0
                        for k_pp in range(num_atoms):
                            if k_pp in (k, k_prime):
                                continue
                            Zk_pp = atoms[k_pp]
                            key_kp_kpp = ('C', Zk_pp)
                            if key_kp_kpp not in r_ZZ_param:
                                continue

                            R_kp_kpp = distance_matrix[k_prime, k_pp]
                            r_ckpp, dr_ckpp = r_ZZ_param[key_kp_kpp]
                            if R_kp_kpp > r_ckpp + dr_ckpp:
                                tanh_kp_kpp = np.exp(dr_ckpp / (R_kp_kpp - dr_ckpp - r_ckpp))
                                inner_sum += tanh_kp_kpp

                        nested_cn_sum += tanh_k_kp * (inner_sum ** 2)

                # Final sigma_k for nitrogen
                if self.solvent == 'water':
                    sigma += sigma_params.get(('N', 'C'), 0.0) * (nested_cn_sum ** 1.3)
                    sigma += sigma_params.get(('N', 'C(3)'), 0.0) * nc3_sum
                else:
                    n_1, alpha_1, beta_1 = sigma_params.get(('N', 'C'), 0.0)
                    sigma += (n * n_1 + alpha * alpha_1 + beta * beta_1) * (nested_cn_sum ** 1.3)
                    n_2, alpha_2, beta_2 = sigma_params.get(('N', 'C(3)'), 0.0)
                    sigma += (n * n_2 + alpha * alpha_2 + beta * beta_2) * nc3_sum

            elif Zk == 'O':
                for k2 in range(num_atoms):
                    if k == k2:
                        continue
                    Zk2 = atoms[k2]
                    if Zk2 not in ['C', 'N', 'O', 'P']:
                        continue
                    key = (Zk, Zk2)
                    if key not in r_ZZ_param or key not in sigma_params:
                        continue
                    R = distance_matrix[k, k2]
                    r_zz, dr_zz = r_ZZ_param[key]
                    if R > r_zz + dr_zz:
                        tanh = np.exp(dr_zz / (R - dr_zz - r_zz))
                        if self.solvent == 'water':
                            sigma += sigma_params[key] * tanh
                        else:
                            s_n, s_alpha, s_beta = sigma_params[key]
                            sigma += (n * s_n + alpha * s_alpha + beta * s_beta) * tanh

            sigma_k[k] = sigma

        self.ostream.print_info(f"Sigma_k: {sigma_k}")
        self.ostream.flush()
        
        return sigma_k
    
    def _get_SASA(self):
        """ 
        Computes the solvent accessible surface area (SASA) for each atom in the molecule

        """
        # Defined in: "The interpretation of protein structures: Estimation of static accessibility" 
        # https://www.sciencedirect.com/science/article/pii/002228367190324X?via%3Dihub
        # Using the Lee Richard algorithm available from FreeSASA via RDkit: 
        # OBS: The calculations for each atom are completely independent and can thus be parallelized over an arbitrary number of threads, whereas the calculation of adjacency lists has not been parallelized.
        from rdkit import Chem
        from rdkit.Chem import rdFreeSASA
    
        self.ostream.print_info("Computing Solvent Accessible Surface Area (SASA) per atom...")
        self.ostream.flush()
        
        vdw_radii = self.solute.vdw_radii_to_numpy() * bohr_in_angstrom()

        solute = Chem.MolFromXYZBlock(self.solute.get_xyz_string())
        r_s = 1.4 # Å (0.4 Å according to SMD paper (1.4 Å according to e.g., SM6, Lee&Richards etc.))
        
        opts = rdFreeSASA.SASAOpts(
            rdFreeSASA.SASAAlgorithm.LeeRichards,
            rdFreeSASA.SASAClassifier.Protor,  # optimized for small molecules
            r_s
        )

        rdFreeSASA.CalcSASA(solute, radii=vdw_radii.tolist(), confIdx=-1, opts=opts)

        atoms = solute.GetAtoms()
        SASA_list = [float(atoms[i].GetProp("SASA")) for i in range(len(atoms))]
        
        self.ostream.print_info(f"SASA list: {SASA_list}")
        self.ostream.print_info(f"Total SASA: {sum(SASA_list):.4f} Å^2")
        self.ostream.flush()
        
        return SASA_list

    def _get_molecular_surface_tension(self):
        # eq. 9 (OBS: sigma^[beta^2] = 0)
        # Calculate the full contribution here and add to the G_CDS expression
        
        # Table 5: solute-indep. param cal/molÅ^2
        sigma_gamma, sigma_phi2, sigma_psi2 = [0.35, -4.19, -6.68]
        param = self.smd_solvent_parameters
        
        self.sigma_M = (
            sigma_gamma * param['gamma'] 
            + sigma_phi2 * param['phi']**2 
            + sigma_psi2 * param['psi']**2
        )
 
    def _get_CDS_contribution(self):
        """
        Get the Cavity-Dispersion-Solvent-Structure (CDS) contribution to the solvation energy.
        """
        # Cavity-Dispersion-Solvent-Structure contribution: DG_CDS = see eq. 6 (simplified to only first term if water)
        
        CDS_energy = 0.0
        kJ_per_mol = 0.004184
        
        SASA_list = self._get_SASA()

        sigma_k = self._calculate_sigma_k()

        for k in range(len(self.solute.get_labels())):
            CDS_energy += sigma_k[k]*SASA_list[k]
        
        if self.solvent != 'water': #sigma[M] = 0 for water
            self._get_molecular_surface_tension()
            CDS_energy += self.sigma_M * sum(SASA_list)

        CDS_energy *= kJ_per_mol
        
        self.ostream.print_info(f"CDS Contribution: {CDS_energy:.4f} kJ/mol")
        self.ostream.flush()
        
        return CDS_energy
    
    
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
import sys

from .veloxchemlib import mpi_master, bohr_in_angstrom
from .outputstream import OutputStream
from .smdsolventproperties import get_smd_solvent_properties, get_sigma_properties, get_rzz_parameters
from .atomtypeidentifier import AtomTypeIdentifier
from .errorhandler import assert_msg_critical

class SmdDriver:
    """
    SMD (Solvation Model based on Density) calculations.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.
    Instance variables:
        - solvent (str): The solvent used for SMD calculations (default is 'water'). 
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes SMD driver.
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
        
        # solvent properties
        self.solvent = 'water'
        self.smd_solvent_parameters = get_smd_solvent_properties() 
        self.sigma_param, self.sigma_water = get_sigma_properties()
        self.r_zz_dict = get_rzz_parameters()

        self.solute = None

    
    def get_intrinsic_coulomb_radii(self):
        """
        Get the intrinsic Coulomb radii (Å) for the solute.
        """

        atom_labels = self.solute.get_labels() 
        vdw_radii = self.solute.vdw_radii_to_numpy() * bohr_in_angstrom()

        # Table 3 in SMD paper (Å). 
        # If not in smd_elements, use Bondi radius (default in vlx); if not in Bondi, use 2.0 Å (also default in vlx)
        smd_elements = {
            'H':1.20,'C':1.85,'N':1.89,
            'O':1.52,'F':1.73,'Si':2.47,
            'P':2.12,'S':2.49,'Cl':2.38,
            'Br':3.06
        }

        alpha = self.smd_solvent_parameters['alpha']
        atom_radii = []

        for atom, vdw_r in zip(atom_labels, vdw_radii):
            atom_radii.append(atom)
            
            if atom not in smd_elements:
                # Using Bondi Radius
                atom_radii.append(vdw_r)
            else:
                # re-parametrize O radius dep. on the solvent
                if atom == 'O':
                    atom_radii.append(1.52 if alpha >= 0.43 else round((1.52 + 1.8*(0.43 - alpha)), 4))
                else:
                    atom_radii.append(smd_elements[atom])


        return atom_radii

    def get_CDS_contribution(self):
        """
        Get the Cavity-Dispersion-Solvent-Structure (CDS) contribution to the solvation energy.
        """

        assert_msg_critical(
            self.solvent in self.smd_solvent_parameters,
            'Solvent {self.solvent} not available')
        
        self.smd_solvent_parameters = self.smd_solvent_parameters[self.solvent]
        self.epsilon = self.smd_solvent_parameters['epsilon']
        
        CDS_energy = 0.0
        cal_per_mol_to_hartree = 1/627509.5
        
        SASA_list = self._get_SASA()

        sigma_k = self._calculate_sigma_k()

        for k in range(len(self.solute.get_labels())):
            CDS_energy += sigma_k[k] * SASA_list[k]
        
        if self.solvent != 'water': # sigma[M] = 0 for water
            self._get_molecular_surface_tension()
            CDS_energy += self.sigma_M * sum(SASA_list)

        CDS_energy *= cal_per_mol_to_hartree
        
        return CDS_energy

    def _calculate_sigma_k(self):        
        # Table 4 (cal/molÅ^2):  Any possible surface tension parameter that is not in this table is set equal to zero in SMD. 
        param = self.smd_solvent_parameters
        
        # Extract all sigma_i for the solvent
        if self.solvent == 'water':
            sigma_params = self.sigma_water
        else:
            sigma_params = {
                key: vals[0] * param['n'] + vals[1] * param['alpha'] + vals[2] * param['beta']
                for key, vals in self.sigma_param.items()
            }
        
        distance_matrix = self.solute.get_distance_matrix_in_angstrom()
        atoms = self.solute.get_labels()
        num_atoms = len(atoms)
        
        atomtypeidentifier = AtomTypeIdentifier(self.comm)
        atomtypeidentifier.ostream.mute()
        atom_types = atomtypeidentifier.generate_gaff_atomtypes(self.solute)

        halogens = ['F', 'Si', 'S', 'Cl', 'Br'] #if Zk in halogens, sigma_k = sigma_tilde_Z_k 
        all_smd_atoms = ['H', 'C', 'N', 'O'] + halogens #if Z_k not in all_smd, sigma_k = 0

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
                sigma += sigma_params.get((Zk,), 0.0)
                sigma_k[k] = sigma
                continue

            sigma += sigma_params.get((Zk,), 0.0)

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
                    if R < r_zz + dr_zz:
                        tanh = np.exp(dr_zz / (R - dr_zz - r_zz))
                        sigma += sigma_params[key] * tanh

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
                    if R < r_zz + dr_zz:
                        tanh = np.exp(dr_zz / (R - dr_zz - r_zz))
                        if Zk2 == 'C':
                            sigma += sigma_params[key] * tanh
                        else:
                            tanh_sum += tanh
                
                sigma += sigma_params.get(('C','N'), 0.0) * (tanh_sum ** 2)

            elif Zk == 'N':

                nested_cn_sum = 0.0
                nc3_sum = 0.0
                
                for k_prime in range(num_atoms):
                    if k_prime == k:
                        continue
                    
                    Zk_prime = atoms[k_prime]
                    atom_type = atom_types[k_prime]

                    # NC(3) term - only in water
                    if atom_type == 'c3' and ('N', 'c3') in sigma_params:
                        R = distance_matrix[k, k_prime]
                        r_zz, dr_zz = r_ZZ_param[('N', 'c3')]
                        if R < r_zz + dr_zz:
                            tanh_k_kp = np.exp(dr_zz / (R - dr_zz - r_zz))
                            nc3_sum += tanh_k_kp
                    
                    # NC + C-X nested sum
                    elif Zk_prime == 'C' and atom_type != 'c3':
                        R = distance_matrix[k,k_prime]
                        r_nc, dr_nc = r_ZZ_param[('N', 'C')]
                        if R >= r_nc + dr_nc: 
                            continue

                        tanh_k_kp = np.exp(dr_nc / (R - dr_nc - r_nc))

                        # inner sum over k'k"
                        inner_sum = 0.0
                        for k_pp in range(num_atoms):
                            if k_pp in (k, k_prime):
                                continue
                            Zk_pp = atoms[k_pp]
                            key_kp_kpp = (Zk_prime, Zk_pp)
                            if key_kp_kpp not in r_ZZ_param:
                                continue

                            R = distance_matrix[k_prime, k_pp]
                            r_ckpp, dr_ckpp = r_ZZ_param[key_kp_kpp]
                            if R < r_ckpp + dr_ckpp:
                                tanh_kp_kpp = np.exp(dr_ckpp / (R - dr_ckpp - r_ckpp))
                                inner_sum += tanh_kp_kpp

                        nested_cn_sum += tanh_k_kp * (inner_sum ** 2)

                # Final sigma_k for nitrogen
                sigma += sigma_params.get(('N', 'C'), 0.0) * (nested_cn_sum ** 1.3)
                sigma += sigma_params.get(('N', 'c3'), 0.0) * nc3_sum

            elif Zk == 'O':
                for k2 in range(num_atoms):
                    if k == k2:
                        continue
                    Zk2 = atoms[k2]
                    if Zk2 not in ['C', 'N', 'O', 'P']:
                        continue
                    key = (Zk, Zk2)
                    if key not in sigma_params:
                        continue
                    R = distance_matrix[k, k2]
                    r_zz, dr_zz = r_ZZ_param[key]
                    if R < r_zz + dr_zz:
                        tanh = np.exp(dr_zz / (R - dr_zz - r_zz))
                        sigma += sigma_params[key] * tanh

            sigma_k[k] = sigma
        
        return sigma_k
    
    def _get_SASA(self):
        """ 
        Computes the solvent accessible surface area (SASA) for each atom in the molecule
        """

        # Defined in: "The interpretation of protein structures: Estimation of static accessibility" 
        # https://www.sciencedirect.com/science/article/pii/002228367190324X?via%3Dihub
        # Using the Lee Richard algorithm available from FreeSASA via RDkit 

        try:
            from rdkit import Chem
            from rdkit.Chem import rdFreeSASA
        except ImportError:
            raise ImportError('Unable to import rdkit.')
        
        solute = Chem.MolFromXYZBlock(self.solute.get_xyz_string())
        ptable = Chem.GetPeriodicTable()
        radii = [ptable.GetRvdw(atom.GetAtomicNum()) for atom in solute.GetAtoms()]
        r_s = 0.4 
        
        opts = rdFreeSASA.SASAOpts(
            rdFreeSASA.SASAAlgorithm.LeeRichards,
            rdFreeSASA.SASAClassifier.Protor,  
            r_s
        )

        rdFreeSASA.CalcSASA(solute, radii=radii, confIdx=-1, opts=opts)

        atoms = solute.GetAtoms()
        SASA_list = [float(atoms[i].GetProp("SASA")) for i in range(len(atoms))]
        
        return SASA_list

    def _get_molecular_surface_tension(self):
        # Eq. 9 (Note: sigma^[beta^2] = 0)
        
        # Table 5: solute-indep. param cal/molÅ^2
        sigma_gamma, sigma_phi2, sigma_psi2 = [0.35, -4.19, -6.68]
        param = self.smd_solvent_parameters
        
        self.sigma_M = (
            sigma_gamma * param['gamma'] 
            + sigma_phi2 * param['phi']**2 
            + sigma_psi2 * param['psi']**2
        )
 

    
    

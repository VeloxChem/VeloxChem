import os
import numpy as np
import veloxchem as vlx
import py3Dmol as p3d
from collections import defaultdict
from scipy import optimize


class MPC1:
    def __init__(self, molecule, deprotonated=False):
        self.molecule = molecule
        self.deprotonated = deprotonated 
        self.basis = None
        self.scf_results = None
        self.density_matrix = None
        self.surface_points = None
        self.atom_indices = None
        self.ea_values = None
        self.ie_values = None
        self.loprop_charges = None
        self.resp_charges = None
        self.max_esp_values = None
    
    
    def run_scf(self, functional='BLYP', solvation_model='CPCM'):
        if self.deprotonated == True:
            self.molecule.set_charge(-1)
        self.basis = vlx.MolecularBasis.read(self.molecule, "def2-SVP")
        self.scf_drv = vlx.ScfRestrictedDriver()
        self.scf_drv.xcfun = functional
        self.scf_drv.solvation_model = solvation_model
        # scf_drv.update_settings({'xcfun': 'BLYP'}, {'solvation_model': 'CPCM'})
        self.scf_drv.grid_level = 2
        self.scf_drv.conv_thresh = 1e-4
        self.scf_drv.cpcm_epsilon = 78.4
        self.scf_drv.ostream.mute()
        self.scf_drv.cpcm_grid_per_sphere = 194 # valid grid numbers are [50, 110, 194, 302, 434, 590, 770, 974, 2030]
        self.scf_drv.dispersion = True

        if functional.upper() =='B3LYP':
            self.scf_drv.ri_coulomb = False
        else:
            self.scf_drv.ri_coulomb = True
        self.scf_results = self.scf_drv.compute(self.molecule, self.basis)
    
    def optimize_geometry(self):

        # Run geometry optimization using BLYP functional
        self.run_scf(functional='BLYP')

        opt_drv = vlx.OptimizationDriver(self.scf_drv)
        opt_drv.ostream.mute()
        opt_drv.conv_energy = 1e-04
        opt_drv.conv_drms = 5e-03
        opt_drv.conv_dmax = 1e-02
        opt_drv.conv_grms = 1e-03
        opt_drv.conv_gmax = 3e-03
        optimize_res = opt_drv.compute(self.molecule, self.basis, self.scf_results)
        self.molecule = vlx.Molecule.read_xyz_string(optimize_res['final_geometry'])
        if self.deprotonated == True:
            self.molecule.set_charge(-1)
        # Use B3LYP functional for further calculations
        self.run_scf(functional='B3LYP', solvation_model='CPCM')


    def get_xyz_string(self):
        return self.molecule.get_xyz_string()
    
    def compute_electron_density(self):

        grid_drv = vlx.GridDriver()
        grid_drv.set_level(4)
        molgrid = grid_drv.generate(self.molecule)
        chi_g = vlx.XCIntegrator().compute_gto_values(self.molecule, self.basis, molgrid)

        D = self.scf_results["D_alpha"] + self.scf_results["D_beta"]
        G = np.einsum("ab,bg->ag", D, chi_g)
        n_g = np.einsum("ag,ag->g", chi_g, G)
        
        # Extract surface points based on density range
        indices = np.where((n_g > 0.0095) & (n_g < 0.0105))
        self.surface_points = np.column_stack((molgrid.x_to_numpy()[indices],
                                               molgrid.y_to_numpy()[indices],
                                               molgrid.z_to_numpy()[indices]))
        
        
    def _get_surface_point_atom_indices(self):

        atom_positions = self.molecule.get_coordinates_in_bohr()
        atom_vdw_radii = self.molecule.vdw_radii_to_numpy()  # Get van der Waals radii for each atom in bohr
        atom_indices = []

        for point in self.surface_points:

            # Calculate distances from this point to all atoms
            distances = np.linalg.norm(atom_positions - point, axis=1)

            # Adjust distances by subtracting van der Waals radii
            adjusted_distances = distances - atom_vdw_radii

            # Find the closest atom considering van der Waals radii
            closest_atom_index = np.argmin(adjusted_distances)
            atom_indices.append(closest_atom_index)
        
        self.atom_indices = atom_indices

    def compute_molecular_orbitals(self):

        grid_drv = vlx.GridDriver() # Handles the numerical grid used to evaluate the electron density
        molgrid = grid_drv.generate(self.molecule)
        xc_drv = vlx.XCIntegrator()
        
        # generate AOs on the grid points
        chi_g = xc_drv.compute_gto_values(self.molecule, self.basis, molgrid) # Computes the values of atomic orbitals (GTOs) at grid points


        mol_orb_coeff = self.scf_results["C_alpha"]
        occ_mo_val, unocc_mo_val = [], []
        
        occ_orb = sum(self.scf_results['occ_alpha'])
        tot_orb = len(self.scf_results['occ_alpha'])
        print(f"Number of occupied orbitals: {occ_orb}")
        print(f"Total number of orbitals: {tot_orb}")
        
        for i in range(int(occ_orb)):
            mo_val = np.dot(chi_g.T, mol_orb_coeff[:, i])
            occ_mo_val.append(mo_val)
        
        for i in range(int(occ_orb), int(tot_orb)):
            mo_val = np.dot(chi_g.T, mol_orb_coeff[:, i])
            unocc_mo_val.append(mo_val)
        
        return occ_mo_val, unocc_mo_val
    
    def compute_local_properties(self):
        print(f'The total amount of surface points is {len(self.surface_points)}')

        occ_mo_val, unocc_mo_val = self.compute_molecular_orbitals()
        energy_occ = self.scf_results['E_alpha'][:len(occ_mo_val)]
        energy_unocc = self.scf_results['E_alpha'][len(occ_mo_val):]
        # den_val = self.compute_electron_density() # Electron density values at surface points
              

        EA, IE = [], []
        for idx, point in enumerate(self.surface_points):

            ea_point_val, ie_point_val = 0, 0 # Initialize EA and IE values for this point

            for j in range(len(unocc_mo_val)):
              
                unocc_mo_val_squared = unocc_mo_val[j][idx]**2
                ea_point_val += unocc_mo_val_squared * energy_unocc[j] # Sum over all unoccupied MOs

            ea_point_val /= sum(unocc_mo_val[j][idx]**2 for j in range(len(unocc_mo_val))) * (-1)


            for j in range(len(occ_mo_val)):
                
                occ_mo_val_squared = occ_mo_val[j][idx]**2
                ie_point_val += occ_mo_val_squared * abs(energy_occ[j]) # Sum over all occupied MOs

            ie_point_val /= sum(occ_mo_val[j][idx]**2 for j in range(len(occ_mo_val)))

            EA.append(ea_point_val)
            IE.append(ie_point_val)

        
        self.ea_values = EA
        self.ie_values = IE
    
        data_bank = [{'atom': int(self.atom_indices[idx]), 'IE': IE[idx], 'EA': EA[idx]} for idx in range(len(IE))]
        # print(f"Data bank: {data_bank}")

	    # Organize data into a dictionary where keys are atoms and values are lists of (EA, index) tuples
        self.EA_data = defaultdict(lambda: {"EA": [], "indices": []})
        self.IE_data = defaultdict(lambda: {"IE": [], "indices": []})

        for idx, entry in enumerate(data_bank):
            atom = entry["atom"]
            self.EA_data[atom]["EA"].append(entry["EA"])
            self.EA_data[atom]["indices"].append(idx)  # Store the index of the point
            self.IE_data[atom]["IE"].append(entry["IE"])
            self.IE_data[atom]["indices"].append(idx)  # Store the index of the point
        
    def _local_property_EA(self):
        self.EA_results = {}
        self.min_EA_points = []  # List to store indices of min EA points

        for atom, values in self.EA_data.items():
            # Assuming 'values' contains 'EA' (list of EA values) and 'indices' (corresponding surface point indices)
            if values["EA"]:  # Ensure there's at least one EA value
                min_EA_value = min(values["EA"])  # Get min EA value
                min_EA_index = values["indices"][values["EA"].index(min_EA_value)]  # Find corresponding surface point index
                max_EA_value = max(values["EA"])  # Get max EA value
                mean_EA_value = np.mean(values["EA"])  # Calculate mean EA value
                median_EA_value = np.median(values["EA"])  # Calculate median EA value

                self.EA_results[atom] = {
                    "min_EA": min_EA_value,
                    "max_EA": max_EA_value,
                    "mean_EA": mean_EA_value,
                    "median_EA": median_EA_value,
                    "min_EA_index": min_EA_index
                }

                self.min_EA_points.append(min_EA_index)  # Store min EA point index

        return self.EA_results, self.min_EA_points

    def _local_property_IE(self):
        self.IE_results = {}
        self.min_IE_points = []  # List to store indices of min IE points

        for atom, values in self.IE_data.items():
            # Assuming 'values' contains 'IE' (list of IE values) and 'indices' (corresponding surface point indices)
            if values["IE"]:  # Ensure there's at least one IE value
                min_IE_value = min(values["IE"])  # Get min IE value
                min_IE_index = values["indices"][values["IE"].index(min_IE_value)]  # Find corresponding surface point index
                max_IE_value = max(values["IE"])  # Get max IE value
                mean_IE_value = np.mean(values["IE"])  # Calculate mean IE value
                median_IE_value = np.median(values["IE"])  # Calculate median IE value

                self.IE_results[atom] = {
                    "min_IE": min_IE_value,
                    "max_IE": max_IE_value,
                    "mean_IE": mean_IE_value,
                    "median_IE": median_IE_value,
                    "min_IE_index": min_IE_index
                }

                self.min_IE_points.append(min_IE_index)  # Store min IE point index

        return self.IE_results, self.min_IE_points

    def _get_max_min_EA_IE(self):
        # EA_data and IE_data should be structured as dictionaries where the atom index maps to a list of EA/IE values
        self.result_EA, self.min_EA_points = self._local_property_EA()
        self.result_IE, self.min_IE_points = self._local_property_IE()

        return self.result_EA, self.result_IE, self.min_EA_points, self.min_IE_points


    def _compute_resp_loprop(self):
        resp_drv = vlx.RespChargesDriver()
        resp_drv.ostream.mute()
        self.resp_charges = resp_drv.compute(self.molecule, self.basis, self.scf_results, 'resp')
        
        if 'solvation_model' in self.scf_results:
            del self.scf_results['solvation_model']

        loprop_drv = vlx.LoPropDriver()
        loprop_drv.ostream.mute()
        self.loprop_charges = loprop_drv.compute(self.molecule, self.basis, self.scf_results)

        localized_polarizabilities = self.loprop_charges['localized_polarizabilities']
        self.local_polarizabilities = []

        # Compute the polarizability vector for each atom
        for polarizability_tensor in localized_polarizabilities:
            # Reshape the tensor components into a 3x3 matrix
            polarizability_matrix = np.array([
                [polarizability_tensor[0], polarizability_tensor[1], polarizability_tensor[2]],
                [polarizability_tensor[1], polarizability_tensor[3], polarizability_tensor[4]],
                [polarizability_tensor[2], polarizability_tensor[4], polarizability_tensor[5]]
            ])

            # Compute the eigenvalues and eigenvectors of the polarizability matrix
            eigenvalues, eigenvectors = np.linalg.eigh(polarizability_matrix)

            # Find the principal polarizability vector (corresponding to the largest eigenvalue)
            max_index = np.argmax(eigenvalues)
            principal_vector = eigenvectors[:, max_index]

            # Store the principal polarizability vector for the atom
            self.local_polarizabilities.append(principal_vector)

    
    def compute_esp(self):

        respdrv = vlx.RespChargesDriver()
        respdrv.ostream.mute()
        self.max_esp_values, self.min_esp_values = respdrv.get_electrostatic_potential_with_maximum_and_minimum_values(
            self.surface_points, self.molecule, self.basis, self.scf_results, self.atom_indices
        )
    
    def generate_data_matrix(self):
        num_atoms = len(self.resp_charges)
        data_matrix = np.zeros((num_atoms, 18))

        for i, atom_idx in enumerate(range(num_atoms)):
            data_matrix[i, 0] = atom_idx + 1
            data_matrix[i, 1] = self.molecule.get_element_ids()[atom_idx]
            data_matrix[i, 2] = self.result_EA[atom_idx]['min_EA']
            data_matrix[i, 3] = self.result_EA[atom_idx]['max_EA']
            data_matrix[i, 4] = self.result_EA[atom_idx]['mean_EA']
            data_matrix[i, 5] = self.result_EA[atom_idx]['median_EA']
            data_matrix[i, 6] = self.result_IE[atom_idx]['min_IE']
            data_matrix[i, 7] = self.result_IE[atom_idx]['max_IE']
            data_matrix[i, 8] = self.result_IE[atom_idx]['mean_IE']
            data_matrix[i, 9] = self.result_IE[atom_idx]['median_IE']
            data_matrix[i, 10] = self.loprop_charges['localized_charges'][i]
            data_matrix[i, 11] = self.resp_charges[i]
            data_matrix[i, 12] = self.max_esp_values.get(atom_idx, np.nan)
            data_matrix[i, 13] = self.min_esp_values.get(atom_idx, np.nan)
            data_matrix[i, 14:17] = self.local_polarizabilities[i]  # Principal polarizability vector
            data_matrix[i,17] = self.scf_results['scf_energy'] # SCF energy

        return data_matrix

    def run_all_calculations(self):
        self.optimize_geometry()
        print('Successfully optimized the geometry')
        self.compute_electron_density()
        print('Successfully computed the electron density')
        self._get_surface_point_atom_indices()
        print('Successfully computed the surface point atom indices')
        self.compute_molecular_orbitals()
        print('Successfully computed the molecular orbitals')
        self.compute_local_properties()
        print('Successfully computed the local properties')
        self._local_property_EA()
        print('Successfully computed the local property EA')
        self._local_property_IE()
        print('Successfully computed the local property IE')
        self._get_max_min_EA_IE()
        print('Successfully computed the max and min EA and IE')
        self._compute_resp_loprop()
        print('Successfully computed the RESP and LoProp charges')
        self.compute_esp()
        print('Successfully computed the ESP')
        return self.generate_data_matrix()
    print('Successfully generated the data matrix')
    

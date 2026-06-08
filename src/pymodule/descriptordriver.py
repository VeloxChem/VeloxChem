from collections import defaultdict
import numpy as np
import veloxchem as vlx


class DescriptorDriver:
    def __init__(self):
        self.molecule = None       #molecule for descriptor calculation
        self.basis = None          
        self.scf_drv = None        
        self.scf_results = None    
        
        self.molgrid = None                     
        self.surface_points = None            
        self.point_atom_indices = None        

        self.occ_mo_points_amplitude = []      #MO amplitude for occupied orbitals in every point
        
        self.unocc_mo_points_amplitude = []    #MO amplitude for unoccupied orbitals in every point


        self.ea_values = []
        self.ie_values = []

        self.atomtypes = []                    #List containing atomtypes of all atoms according to amber atomtypes in GAFF (atom indices 1 - n)
        self.SASA = None                       #Solvent-accessible surface area (Å^2) 
        self.log_p = None

        self.resp_charges = None
        self.atom_indices = None
        

        self.ea_ie_results = None

        self.results_dict = {}





    def compute_respcharges(self):
                # RESP computation
        resp_drv = vlx.RespChargesDriver()
        resp_basis = vlx.MolecularBasis.read(self.molecule, "6-31G*")
        self.resp_charges = resp_drv.compute(self.molecule, resp_basis)
        
        
    
   
    

    #Work-around for vlx.XCIntegrator in order to apply it to cpcm grid. 
    def compute_gto_values_at_custom_points(self, custom_points: np.ndarray):
        """
        custom_points: numpy array of shape (N, 3) in bohr
        """
       
        weights = np.ones(len(custom_points))

        # 2. Stack them into a (4, N) array
        # Use np.vstack to stack the x, y, z, and weights vertically
        molgrid_points = np.vstack((
        custom_points[:, 0],  # row 0: x
        custom_points[:, 1],  # row 1: y
        custom_points[:, 2],  # row 2: z
        weights               # row 3: weights
        ))
        from veloxchem.veloxchemlib import DenseMatrix, MolecularGrid

       


        coords_mat = DenseMatrix(molgrid_points)
        mol_grid = MolecularGrid(coords_mat)


        mol_grid.partition_grid_points()
        mol_grid.distribute_counts_and_displacements(0, 1) 

        xc_drv = vlx.XCIntegrator()
        chi_g = xc_drv.compute_gto_values(self.molecule, self.basis, mol_grid)
        
        return chi_g
    



    def create_surface_points(self):
        
        cpcm_drv = vlx.CpcmDriver()
        
        raw_cpcm_points = cpcm_drv.generate_cpcm_grid(self.molecule)
        self.molgrid = raw_cpcm_points[0][:, :3]
        chi_g = self.compute_gto_values_at_custom_points(self.molgrid)   #finds which basis functions (GTOs) in what grid point?

        D = self.scf_results['D_alpha'] + self.scf_results['D_beta'] #total electron density matrix (alpha and beta spin)
        G = np.einsum("ab,bg->ag", D, chi_g)       #Here we evaluate each grid-point with weighted basis functions. N_basis X N_gridpoints
        n_g = np.einsum("ag,ag->g", chi_g, G)      #1-D array of electorn densities in the different points 
    
        self.surface_points = np.column_stack((self.molgrid[:, 0], self.molgrid[:, 1], self.molgrid[:, 2])) 
        print(f'No. of surface points created: {len(self.surface_points)}')
        





    def assign_points_to_atoms(self):
        atom_positions = self.molecule.get_coordinates_in_bohr()
        atom_vdw_radii = self.molecule.vdw_radii_to_numpy()  # Get van der Waals radii for each atom in bohr

        point_atom_indices = []
        for point in self.surface_points:
            distances = np.linalg.norm(atom_positions - point, axis=1)     # Calculate distances from this point to all atoms
            adjusted_distances = distances - atom_vdw_radii      # Adjust distances by subtracting van der Waals radii
            closest_atom_index = np.argmin(adjusted_distances)       # Find the index of the closest atom considering van der Waals radii
            point_atom_indices.append(closest_atom_index)

        self.point_atom_indices = point_atom_indices
        
        #check that amount of surface points equals the no. of point-to-atom indices 
        assert len(self.surface_points) == len(self.point_atom_indices), \
        f"Mismatch! Created {len(self.surface_points)} surface points " \
        f"but generated {len(self.point_atom_indices)} point-to-atom indices"



    def molecular_orbital_amplitude_per_point(self):
        
        chi_g = self.compute_gto_values_at_custom_points(self.molgrid)
        mol_orb_coeff = self.scf_results["C_alpha"]

        occ_orb = sum(self.scf_results['occ_alpha'])   
        tot_orb = len(self.scf_results['occ_alpha'])

        for i in range(int(occ_orb)):                              #for occupied orbitals
            mo_val = np.dot(chi_g.T, mol_orb_coeff[:, i])   # Atomic orbitals amplitude (per point in gid) * Coeffients for linnearily combining AO into MO = MO amplitude in each point
            self.occ_mo_points_amplitude.append(mo_val)

        for i in range(int(occ_orb), int(tot_orb)):                     #same
            mo_val = np.dot(chi_g.T, mol_orb_coeff[:, i])
            self.unocc_mo_points_amplitude.append(mo_val)
                

    def compute_IE_and_EA(self):
        energy_occ = self.scf_results['E_alpha'][:int(sum(self.scf_results['occ_alpha']))]    # Energies (Hartree) of occupied alpha MOs as calculated by scfdrver    
        energy_unocc = self.scf_results['E_alpha'][int(sum(self.scf_results['occ_alpha'])):] # Energies (Hartree) of virtual (empty) alpha MOs as calculated by scfdrver
         
        for idx, point in enumerate(self.surface_points):     #loop over all the surface points, keeping idx to find the correct surface point for every MO in "self.unocc_mo_points_amplitude"
            
            ea_point_numerator, ea_point_denominator  = 0, 0 # Initialize 

            for j in range(len(self.unocc_mo_points_amplitude)):  #Loop over  all unoccupied MOs; LUMO, LUMO+1, etc and their amplitude in every point [[MO1: p1, p2 ... pn], [MO2 p1, p2 ... pn]]

                unocc_mo_ampl_squared = (self.unocc_mo_points_amplitude[j][idx])**2         #squared cumulative amplitude for LUMO1 point1, LUMO2 point 1, LUMO3 point1 until LUMOn, point n

                ea_point_numerator += (unocc_mo_ampl_squared * energy_unocc[j])        #sum of squared amplitudes x virtual orbital energy
                ea_point_denominator += unocc_mo_ampl_squared                          #sum of squared amplitudes

            ea = (-ea_point_numerator / ea_point_denominator )* 27.211407953              #Hartree to eV
            self.ea_values.append(ea)  


            ie_point_numerator, ie_point_denominator  = 0, 0 # Initialize 

            for j in range(len(self.occ_mo_points_amplitude)):  #Loop over  all occupied MOs and their amplitude in every point [[MO1: p1, p2 ... pn], [MO2 p1, p2 ... pn]]

                occ_mo_ampl_squared = (self.occ_mo_points_amplitude[j][idx])**2         #squared cumulative amplitude 

                ie_point_numerator += (occ_mo_ampl_squared * abs(energy_occ[j]))        #sum of squared amplitudes x orbital energy
                ie_point_denominator += occ_mo_ampl_squared                          #sum of squared amplitudes

                 
            ie = (ie_point_numerator/ie_point_denominator)* 27.211407953            #Hartree to eV
            self.ie_values.append(ie)

        print(f'{len(self.ea_values)}, {len(self.ie_values)} values of electron affinity and ionization energy calculated for {len(self.surface_points)} surface points ')




    def compute_atom_statistics(self):

        # Create point-level data
        point_data = []
        for i in range(len(self.surface_points)):
            point_data.append({
                'point_index': i + 1,
                'atom_label': self.molecule.get_label(self.point_atom_indices[i]),
                'atom_index': int(self.point_atom_indices[i]),
                'IE': self.ie_values[i],
                'EA': self.ea_values[i]
            })
        
        # Group by atom for statistics
        ea_data = defaultdict(list)
        ie_data = defaultdict(list)
        
        for point in point_data:
            atom_idx = point['atom_index']
            ea_data[atom_idx].append(point['EA'])
            ie_data[atom_idx].append(point['IE'])
        
        # Calculate per-atom statistics
        atom_statistics = {}
        
        for atom in range(len(self.molecule.get_labels())):
            print(f'Atom : {atom}')
            if not ea_data[atom] or not ie_data[atom]:
                atom_statistics[atom] = {
                    'min_EA': np.nan,
                    'max_EA': np.nan,
                    'mean_EA': np.nan,
                    'std_EA': np.nan,
                    'n_points': 0,
                    'min_IE': np.nan,
                    'max_IE': np.nan,
                    'mean_IE': np.nan,
                    'std_IE': np.nan
                }
                continue
            atom_statistics[atom] = {
                'min_EA': min(ea_data[atom]),
                'max_EA': max(ea_data[atom]),
                'mean_EA': np.mean(ea_data[atom]),
                'std_EA': np.std(ea_data[atom]),
                'n_points': len(ea_data[atom]),
                'min_IE': min(ie_data[atom]),
                'max_IE': max(ie_data[atom]),
                'mean_IE': np.mean(ie_data[atom]),
                'std_IE': np.std(ie_data[atom])
            }
        
        # Store in self.results
        self.ea_ie_results = atom_statistics
        
        

        
    def define_atom_types(self):                           #Uses AtomTypeIdentifier to categorize each atom according to GAFF
        atomtypeidentifier = vlx.AtomTypeIdentifier()
        self.atomtypes = atomtypeidentifier.generate_gaff_atomtypes(self.molecule)




    def compute_log_p_and_SASA(self): 
    

        #SASA calculation
        smd = vlx.SmdDriver()
        smd.solute = self.molecule
        sasa_list = smd._get_SASA()
        self.sasa = sum(sasa_list) # Å^2

        #Calc scf energy for water solvation
        basis = vlx.MolecularBasis.read(self.molecule, 'def2-svp')

        scf_drv = vlx.ScfRestrictedDriver()
        scf_drv.xcfun = 'b3lyp'
        scf_drv.solvation_model = 'smd'
        scf_drv.smd_solvent = 'water'

        scf_results_water = scf_drv.compute(self.molecule, basis)

        #Calc scf energy for octanol solvation
        basis = vlx.MolecularBasis.read(self.molecule, 'def2-svp')

        scf_drv = vlx.ScfRestrictedDriver()
        scf_drv.xcfun = 'b3lyp'
        scf_drv.solvation_model = 'smd'
        scf_drv.smd_solvent = '1-octanol'

        scf_results_octanol = scf_drv.compute(self.molecule, basis)
        

        #calculate log p
        hartree_to_j_mol = 2625500.2
        R = 8.3144626  # J / (mol * K)
        T = 298.15     # Kelvin

 
        ddg_solv_hartree  = scf_results_octanol['scf_energy'] - scf_results_water['scf_energy']
        ddg_solv_j_mol = ddg_solv_hartree * hartree_to_j_mol

        self.log_p = -ddg_solv_j_mol / (np.log(10) * R * T)



    def compile_results(self):
        
        self.results_dict['IE_EA'] = self.ea_ie_results
        self.results_dict['atomtypes'] = self.atomtypes
        self.results_dict['log_p'] = self.log_p
        self.results_dict['sasa'] = self.sasa
        self.results_dict['RESP_charges'] = self.resp_charges
        self.results_dict['ie_surface_average'] = np.mean(self.ie_values)
        self.results_dict['ea_surface_average'] = np.mean(self.ea_values)
        self.results_dict['scf_results'] = self.scf_results


        









    def compute_descriptors(self, molecule, basis=None, scf_drv=None, scf_results = None):
        self.molecule = molecule


        # Create basis if not provided
        if basis is None:
            basis = vlx.MolecularBasis.read(molecule, 'def2-svpd')
        self.basis = basis
    
        # instantiate SCF driver if not provided  
        if scf_drv is None:
            scf_drv = vlx.ScfRestrictedDriver()
            scf_drv.xcfun = 'b3lyp'
        self.scf_drv = scf_drv
        
        if scf_results is None:
            scf_results = scf_drv.compute(self.molecule, self.basis)  #electronic structure needed to find surface grid
        self.scf_results = scf_results
 
      
        # Compute RESP charges for all geometries
        print(f"computing RESP-charges...")
        self.compute_respcharges()    

        #compute surface points for molecule
        print(f"Creating surface grid...")
        self.create_surface_points()
        
        #assign each point to closest atom
        print(f"Assigning grid points to atoms...")
        self.assign_points_to_atoms()

        #calculate molecular orbital amplitude at each point
        print(f"Computing atomic IE and EA enegies...")
        self.molecular_orbital_amplitude_per_point()

        self.compute_IE_and_EA()    #eV

        # Compute atom statistics and store in self.results
        
        self.compute_atom_statistics()



        # Define atomtypes in molecule
        self.define_atom_types()


       # self.compute_log_p_and_SASA()
        self.compute_log_p_and_SASA() 
        self.compile_results()


        return self.results_dict

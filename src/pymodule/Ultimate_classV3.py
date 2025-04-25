import os
import time
import h5py
import shutil
import logging
import numpy as np
import pandas as pd
import veloxchem as vlx

from rdkit import Chem
from rdkit.Chem import rdDetermineBonds

logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] %(message)s')

class MoleculeProcessor2:
    def __init__(self, folder_paths, output_folder='results', deprotonate = False, protonate = False):
        self.folder_paths = folder_paths
        self.output_folder = output_folder
        self.deprotonate = deprotonate
        self.protonate = protonate
        self.iodine_files = []
        self.non_iodine_files = []
        self.failed_deprot = []
        self.failed_mapping = []
        self.h5_deprotonated_dir = 'h5_deprotonated_molecules'
        self.deprotonated_dir = 'deprotonated_molecules'
        self.protonated_dir = 'protonated_molecules'
        self.h5_protonated_dir = 'h5_protonated_molecules'

        os.makedirs(self.output_folder, exist_ok=True)
        os.makedirs(self.h5_deprotonated_dir, exist_ok=True)
        os.makedirs(self.deprotonated_dir, exist_ok=True)
        os.makedirs(self.h5_protonated_dir, exist_ok=True)
        os.makedirs(self.protonated_dir, exist_ok=True)

        print(f"Created output directory: {self.output_folder}")
        print(f"Created H5 deprotonated directory: {self.h5_deprotonated_dir}")
        print(f"Created deprotonated molecules directory: {self.deprotonated_dir}")
        print(f"Created H5 protonated directory: {self.h5_protonated_dir}")
        print(f"Created protonated molecules directory: {self.protonated_dir}")

    def classify_molecules(self):
        for folder_path in self.folder_paths:
            xyz_files = [f for f in os.listdir(folder_path) if f.endswith(".xyz")]
            for file in xyz_files:
                with open(os.path.join(folder_path, file), "r") as f:
                    content = f.read()
                    if " I " in content or "I\n" in content:
                        self.iodine_files.append(os.path.join(folder_path, file))
                    else:
                        self.non_iodine_files.append(os.path.join(folder_path, file))
        logging.info(f"Found {len(self.iodine_files)} iodine and {len(self.non_iodine_files)} non-iodine files.")

    def calculate_properties(self):
        for xyz_file in self.non_iodine_files:
            folder, filename = os.path.split(xyz_file)

            h5_file_path = os.path.join(self.output_folder, f'{filename}.h5')
            if os.path.exists(h5_file_path):
                logging.info(f"Skipping {filename} as its h5 file already exists.")
                continue

            try:
                molecule = vlx.Molecule.read_xyz_file(xyz_file)
            
                calc = vlx.MPC1(molecule)

                print('Starting calculations for', filename)

                data_matrix = calc.run_all_calculations()

                print('Yippie! calculations done')

                df = pd.DataFrame(data_matrix, columns=[
                    'Atom idx', 'Atom numbers', 'min_EA', 'max_EA', 'mean_EA', 'median_EA', 'min_IE',
                    'max_IE', 'mean_IE', 'median_IE', 'min_AE', 'max_AE', 'mean_AE', 'median_AE',
                    'localized_charge', 'resp_charge', 'max_esp_values', 'min_esp_values',
                    'local_polarizability X', 'local_polarizability Y', 'local_polarizability Z', 'Total energy'
                ])
                with h5py.File(os.path.join(self.output_folder, f'{filename}.h5'), 'w') as f:
                    for key, value in df.items():
                        f.create_dataset(key, data=value)
                    xyz_coords = np.array(calc.get_xyz_string(), dtype=h5py.string_dtype(encoding='utf-8'))
                    f.create_dataset('xyz coordinates', data=xyz_coords)

                    print('Yippie! h5 file created')

            except Exception as e:
                logging.error(f"Error processing {filename}: {e}")
                
                print('Failed to process', filename)

    def deprotonate_and_process(self):
        if not self.deprotonate:
            return
        self.deprot_info = [] # This is just a list to later see if a molecule containing nitrogen was deprotonated
        logging.info("Deprotonating molecules...")
        
        # Sort the non iodine files
        self.non_iodine_files = sorted(self.non_iodine_files)

        for file in self.non_iodine_files:
            folder, filename = os.path.split(file)
            mol = Chem.MolFromXYZFile(file)
            try:
                rdDetermineBonds.DetermineBonds(mol, 0)
                deprot = vlx.OxygenDeprotonation(mol, filename, output_folder=self.deprotonated_dir)
                self.succes = deprot.Deprotonate() # Succes will be True if the molecule had a deprotonatable hydrogen (OH)
                
                print('The return from Deprotonation.py is', self.succes)
                if self.succes: 
                    self.deprot_info.append(True)
                else:
                    self.deprot_info.append(False)

            except Exception as e:
                logging.error(f"Failed to deprotonate {filename}: {e}")
                self.failed_deprot.append(filename)

        logging.info("Calculating properties for deprotonated molecules...")

        for deprot_file in os.listdir(self.deprotonated_dir):
            if not deprot_file.endswith('.xyz'):
                continue

            prefix = deprot_file[:12] if len(deprot_file) >= 12 else deprot_file[:8]
            deprot_path = os.path.join(self.deprotonated_dir, deprot_file)

            for folder_path in self.folder_paths:
                neutral_file = next((f for f in os.listdir(folder_path) if f.startswith(prefix)), None)
                if neutral_file is None:
                    continue

                neutral_path = os.path.join(folder_path, neutral_file)
                neutral_mol = Chem.MolFromXYZFile(neutral_path)
                deprot_mol = Chem.MolFromXYZFile(deprot_path)

                if neutral_mol is None or deprot_mol is None:
                    logging.warning(f"Could not load molecules: {neutral_file} or {deprot_file}")
                    continue

                try:
                    mapping = vlx.MCS(neutral_mol, deprot_mol, deprotonated = True).get_atom_mapping()
                    if mapping is None:
                        raise ValueError("Atom mapping failed.")
                except Exception as e:
                    logging.error(f"Mapping failed for {deprot_file}: {e}")
                    self.failed_mapping.append(deprot_file)
                    continue

                try:
                    mol = vlx.Molecule.read_xyz_file(deprot_path)
                    calc = vlx.MPC1(mol, deprotonated=True)
                    data_matrix = calc.run_all_calculations()
                    xyz = np.array(calc.get_xyz_string(), dtype=h5py.string_dtype(encoding='utf-8'))
                    df = pd.DataFrame(data_matrix, columns=[
                        'Atom idx', 'Atom numbers', 'deprot_min_EA', 'deprot_max_EA', 'deprot_mean_EA', 'deprot_median_EA',
                        'deprot_min_IE', 'deprot_max_IE', 'deprot_mean_IE', 'deprot_median_IE', 'deprot_min_AE', 'deprot_max_AE',
                        'deprot_mean_AE', 'deprot_median_AE', 'deprot_localized_charge', 'deprot_resp_charge',
                        'deprot_max_esp_values', 'deprot_min_esp_values', 'deprot_local_polarizability X',
                        'deprot_local_polarizability Y', 'deprot_local_polarizability Z', 'deprot_Total energy'
                    ])

                    df = df.iloc[:, 2:]
                    idx_map = {d: p for p, d in mapping}
                    reordered = df.iloc[list(idx_map.keys())]
                    reordered.index = [idx_map[i] for i in idx_map]
                    reordered = reordered.sort_index()

                    neutral_h5_file = next((f for f in os.listdir(self.output_folder) if f.startswith(prefix) and f.endswith('.h5')), None)
                    if neutral_h5_file:
                        deprot_h5_filename = deprot_file.replace(".xyz", ".h5")
                        h5_path = os.path.join(self.h5_deprotonated_dir, deprot_h5_filename)
                        shutil.copy2(os.path.join(self.output_folder, neutral_h5_file), h5_path)
                        with h5py.File(h5_path, 'a') as f:
                            for key, value in reordered.items():
                                f.create_dataset(key, data=value)
                            f.create_dataset("xyz coordinates deprot", data=xyz)
                        logging.info(f"Added deprotonated data to {h5_path}")
                except Exception as e:
                    logging.error(f"Failed to process deprotonated {deprot_file}: {e}")

    def protonate_and_process(self):
        if not self.protonate:
            return

        logging.info("Protonating molecules...")
        number_of_versions=[]

        # Sort the non iodine files
        self.non_iodine_files = sorted(self.non_iodine_files)

        for file in self.non_iodine_files:
            number = 0
            folder, filename = os.path.split(file)
            mol = Chem.MolFromXYZFile(file)
            try:
                rdDetermineBonds.DetermineBonds(mol, 0)
                prot_drv = vlx.NitrogenProtonation(mol, filename, output_folder=self.protonated_dir)
                number = prot_drv.Protonate()
                number_of_versions.append(number)

            except Exception as e:
                logging.error(f"Failed to protonate {filename}: {e}")
                self.failed_deprot.append(filename)

        zero_indices = [i for i, v in enumerate(number_of_versions) if v == 0]
        non_zero_values = [v for v in number_of_versions if v != 0]

        print("Non-zero values:", non_zero_values)
        print("Indices of zeros:", zero_indices)

        logging.info("Calculating properties for protonated molecules...")

        # Similar processing as in deprotonate_and_process...
        idx4info = 0
        
        for prot_file in sorted(os.listdir(self.protonated_dir)):

            if not prot_file.endswith('.xyz'):
                continue
            
                        
            prefix = prot_file[:12] if len(prot_file) >= 12 else prot_file[:8]
            prot_path = os.path.join(self.protonated_dir, prot_file)


            try:
                mol = vlx.Molecule.read_xyz_file(prot_path)
                calc = vlx.MPC1(mol, protonated=True)
                data_matrix = calc.run_all_calculations()
                xyz = np.array(calc.get_xyz_string(), dtype=h5py.string_dtype(encoding='utf-8'))
                df = pd.DataFrame(data_matrix, columns=[
                    'Atom idx', 'Atom numbers', 'prot_min_EA', 'prot_max_EA', 'prot_mean_EA', 'prot_median_EA',
                    'prot_min_IE', 'prot_max_IE', 'prot_mean_IE', 'prot_median_IE',
                    'prot_min_AE', 'prot_max_AE', 'prot_mean_AE', 'prot_median_AE',
                    'prot_localized_charge', 'prot_resp_charge', 'prot_max_esp_values', 'prot_min_esp_values',
                    'prot_local_polarizability X', 'prot_local_polarizability Y',
                    'prot_local_polarizability Z', 'prot_Total energy'
                    ])

                print('THIS IS THE self.deprot_info', self.deprot_info)
                print('This is how many elements that are in the list:', len(self.deprot_info))

                    
                filtered_list = [v for i, v in enumerate(self.deprot_info) if i not in zero_indices]

                print("Filtered list:", filtered_list)

                # Update the filtered list to account for different protonation sites of the molecule
                updated_filtered_list = [v for count,v in zip(non_zero_values, filtered_list) for _ in range(count)]

                print('This is the updated filtered list:', updated_filtered_list)

                print('idx4info is:', idx4info)
                #for i in range(non_zero_values[idx4info]):
                #    print('Entering the first loop to loop through each index (True/False) this many times (should be 2 in this example):', non_zero_values[idx4info])
                #    print('The element of the filtered list is:', filtered_list[idx4info])

                if updated_filtered_list[idx4info] is False:
                    print('The index was False')
                             
                    # Look in h5 folder for the neutral molecule for a match and add the data
                    print('The nitrogen containing molecule did NOT have an OH-group')
                    neutral_h5_file = next((f for f in os.listdir(self.output_folder) if f.startswith(prefix) and f.endswith('.h5')), None)
                    if neutral_h5_file:
                        print('Found matching h5 file in the /results folder', prot_path)
                        prot_h5_filename = prot_file.replace(".xyz", ".h5")
                        h5_path = os.path.join(self.h5_protonated_dir, prot_h5_filename)
                        shutil.copy2(os.path.join(self.output_folder, neutral_h5_file), h5_path)
                        with h5py.File(h5_path, 'a') as f:
                            print('Adding info to the neutral h5 file')
                            # for key, value in reordered.items():
                            for key, value in df.items():
                                f.create_dataset(key, data=value)
                            f.create_dataset("xyz coordinates deprot", data=xyz)
                        logging.info(f"Added deprotonated data to {h5_path}")

                    else:
                            print('Could not find a matching neutral molecule (to edit the h5 file) for molecule:', prot_file)

                else:
                    # Look in deprotonated folder for the neutral molecule for a match and add the data
                    print('The nitrogen containing compound DID also contain an OH-group')
                    neutral_h5_file = next((f for f in os.listdir(self.h5_deprotonated_dir) if f.startswith(prefix) and f.endswith('.h5')), None)
                    if neutral_h5_file:
                        print('Found matching h5 file in the /h5_deprot folder', prot_path)
                        prot_h5_filename = prot_file.replace(".xyz", ".h5")
                        h5_path = os.path.join(self.h5_protonated_dir, prot_h5_filename)
                        shutil.copy2(os.path.join(self.h5_deprotonated_dir, neutral_h5_file), h5_path)
                        with h5py.File(h5_path, 'a') as f:
                            print('Adding info to the deprotonated h5 file')
                            # for key, value in reordered.items():
                            for key, value in df.items():
                                f.create_dataset(key, data=value)
                            f.create_dataset("xyz coordinates deprot", data=xyz)
                        logging.info(f"Added deprotonated data to {h5_path}")

                print('All protonation saker done for The Molle, döda mig')
            
            except Exception as e:
                logging.error(f"Failed to process deprotonated {prot_file}: {e}")
    
            finally:
                idx4info += 1

    def copy_h5_files(self, final_output_folder):
        # Create the final output folder if it doesn't exist
        os.makedirs(final_output_folder, exist_ok=True)
        
        # Define the directories to copy from in the specified order
        directories_to_copy = [
            self.h5_protonated_dir,
            self.h5_deprotonated_dir,
            self.output_folder
        ]
        
        # Iterate over each directory and copy the h5 files
        for directory in directories_to_copy:
            for file_name in os.listdir(directory):
                if file_name.endswith('.h5'):
                    source_path = os.path.join(directory, file_name)
                    dest_path = os.path.join(final_output_folder, file_name)

                    # Check if the file already exists in the final folder
                    if not os.path.exists(dest_path):
                        try:
                            shutil.copy2(source_path, dest_path)
                            logging.info(f"Copied {file_name} from {directory} to {final_output_folder}.")
                        except Exception as e:
                            logging.error(f"Failed to copy {file_name} from {directory} to {final_output_folder}: {e}")
                    else:
                        logging.info(f"Skipping {file_name}, already exists in {final_output_folder}.")

                    
    def run(self):
        start = time.time()
        self.classify_molecules()
        self.calculate_properties()
        self.deprotonate_and_process()
        self.protonate_and_process()
        logging.info(f"Finished processing in {time.time() - start:.2f} seconds")

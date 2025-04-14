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
# from veloxchem import Molecule
# from Atom_mapping import MCS
# from Class import MolecularPropertyCalculator
# from Deprotonation import OxygenDeprotonation
# from Protonation import NitrogenProtonation

logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] %(message)s')

class MoleculeProcessor:
    def __init__(self, folder_paths, output_folder='results', deprotonate=False, protonate = False):
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
            
                calc = vlx.MPC1(molecule, xyz_filename=filename)

                print('Starting calculations for', filename)

                data_matrix = calc.run_all_calculations()

                print('Yippie! calculations done')

                df = pd.DataFrame(data_matrix, columns=[
                    'Atom idx', 'Atom numbers', 'min_EA', 'max_EA', 'mean_EA', 'median_EA', 'min_IE',
                    'max_IE', 'mean_IE', 'median_IE', 'localized_charge', 'resp_charge', 'max_esp_values',
                    'min_esp_values', 'local_polarizability X', 'local_polarizability Y',
                    'local_polarizability Z', 'Total energy'
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

        logging.info("Deprotonating molecules...")

        for file in self.non_iodine_files:
            folder, filename = os.path.split(file)
            mol = Chem.MolFromXYZFile(file)
            try:
                rdDetermineBonds.DetermineBonds(mol, 0)
                deprot = OxygenDeprotonation(mol, filename, output_folder=self.deprotonated_dir)
                deprot.Deprotonate()
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
                    mapping = MCS(neutral_mol, deprot_mol).get_atom_mapping()
                    if mapping is None:
                        raise ValueError("Atom mapping failed.")
                except Exception as e:
                    logging.error(f"Mapping failed for {deprot_file}: {e}")
                    self.failed_mapping.append(deprot_file)
                    continue

                try:
                    mol = Molecule.read_xyz_file(deprot_path)
                    calc = MolecularPropertyCalculator(mol, deprotonated=True)
                    data_matrix = calc.run_all_calculations()
                    xyz = np.array(calc.get_xyz_string(), dtype=h5py.string_dtype(encoding='utf-8'))
                    df = pd.DataFrame(data_matrix, columns=[
                        'Atom idx', 'Atom numbers', 'deprot_min_EA', 'deprot_max_EA', 'deprot_mean_EA', 'deprot_median_EA',
                        'deprot_min_IE', 'deprot_max_IE', 'deprot_mean_IE', 'deprot_median_IE', 'deprot_localized_charge',
                        'deprot_resp_charge', 'deprot_max_esp_values', 'deprot_min_esp_values',
                        'deprot_local_polarizability X', 'deprot_local_polarizability Y',
                        'deprot_local_polarizability Z', 'deprot_Total energy'
                    ])

                    df = df.iloc[:, 2:]
                    idx_map = {d: p for p, d in mapping}
                    reordered = df.iloc[list(idx_map.keys())]
                    reordered.index = [idx_map[i] for i in idx_map]
                    reordered = reordered.sort_index()

                    neutral_h5_file = next((f for f in os.listdir(self.output_folder) if f.startswith(prefix) and f.endswith('.h5')), None)
                    if neutral_h5_file:
                        h5_path = os.path.join(self.h5_deprotonated_dir, neutral_h5_file)
                        shutil.copy2(os.path.join(self.h5_protonated_dir, neutral_h5_file), h5_path)
                        with h5py.File(h5_path, 'a') as f:
                            for key, value in reordered.items():
                                f.create_dataset(f"deprotonated_{key}", data=value)
                            f.create_dataset("xyz coordinates deprot", data=xyz)
                        logging.info(f"Added deprotonated data to {h5_path}")
                except Exception as e:
                    logging.error(f"Failed to process deprotonated {deprot_file}: {e}")

    def protonate_and_process(self):
        if not self.protonate:
            return

        logging.info("Protonating molecules...")

        for file in self.non_iodine_files:
            folder, filename = os.path.split(file)
            mol = Chem.MolFromXYZFile(file)
            try:
                rdDetermineBonds.DetermineBonds(mol, 0)
                prot_drv = NitrogenProtonation(mol, filename, output_folder=self.protonated_dir)
                prot_drv.Protonate()
            except Exception as e:
                logging.error(f"Failed to protonate {filename}: {e}")
                self.failed_deprot.append(filename)

        logging.info("Calculating properties for protonated molecules...")

        # Similar processing as in deprotonate_and_process...   
        for prot_file in os.listdir(self.protonated_dir):
            if not prot_file.endswith('.xyz'):
                continue

            prefix = prot_file[:12] if len(prot_file) >= 12 else prot_file[:8]
            prot_path = os.path.join(self.protonated_dir, prot_file)

            for folder_path in self.folder_paths:
                neutral_file = next((f for f in os.listdir(folder_path) if f.startswith(prefix)), None)
                if neutral_file is None:
                    continue

                neutral_path = os.path.join(folder_path, neutral_file)
                neutral_mol = Chem.MolFromXYZFile(neutral_path)
                prot_mol = Chem.MolFromXYZFile(prot_path)

                if prot_mol is None or neutral_mol is None:
                    logging.warning(f"Could not load molecules: {prot_file} or {prot_file}")
                    continue

                try:
                    mapping = MCS(prot_mol, neutral_mol).get_atom_mapping()
                    if mapping is None:
                        raise ValueError("Atom mapping failed.")
                except Exception as e:
                    logging.error(f"Mapping failed for {prot_file}: {e}")
                    self.failed_mapping.append(prot_file)
                    continue

                try:
                    mol = Molecule.read_xyz_file(prot_path)
                    calc = MolecularPropertyCalculator(mol, protonated=True)
                    data_matrix = calc.run_all_calculations()
                    xyz = np.array(calc.get_xyz_string(), dtype=h5py.string_dtype(encoding='utf-8'))
                    df = pd.DataFrame(data_matrix, columns=[
                        'Atom idx', 'Atom numbers', 'prot_min_EA', 'prot_max_EA', 'prot_mean_EA', 'prot_median_EA',
                        'prot_min_IE', 'prot_max_IE', 'prot_mean_IE', 'prot_median_IE', 'prot_localized_charge',
                        'prot_resp_charge', 'prot_max_esp_values', 'prot_min_esp_values',
                        'prot_local_polarizability X', 'prot_local_polarizability Y',
                        'prot_local_polarizability Z', 'Total energy'
                    ])

                    df = df.iloc[:, 2:]
                    idx_map = {d: p for p, d in mapping}
                    reordered = df.iloc[list(idx_map.keys())]
                    reordered.index = [idx_map[i] for i in idx_map]
                    reordered = reordered.sort_index()

                    neutral_h5_file = next((f for f in os.listdir(self.h5_protonated_dir) if f.startswith(prefix) and f.endswith('.h5')), None)
                    if neutral_h5_file:
                        h5_path = os.path.join(self.h5_protonated_dir, neutral_h5_file)
                        shutil.copy2(os.path.join(self.output_folder, neutral_h5_file), h5_path)
                        with h5py.File(h5_path, 'a') as f:
                            for key, value in reordered.items():
                                f.create_dataset(f"deprotonated_{key}", data=value)
                            f.create_dataset("xyz coordinates deprot", data=xyz)
                        logging.info(f"Added deprotonated data to {h5_path}")
                except Exception as e:
                    logging.error(f"Failed to process deprotonated {neutral_file}: {e}")

                    
    def run(self):
        start = time.time()
        self.classify_molecules()
        self.calculate_properties()
        self.deprotonate_and_process()
        self.protonate_and_process()
        logging.info(f"Finished processing in {time.time() - start:.2f} seconds")

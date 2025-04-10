import os
import time
import h5py
import shutil
import logging
import numpy as np
import pandas as pd
import veloxchem as vlx
import sys

from rdkit import Chem
from rdkit.Chem import rdDetermineBonds


logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] %(message)s', stream = sys.stdout)

class MoleculeProcessor:
    def __init__(self, folder_paths, output_folder='Results4protonated_molecules', deprotonate=True):
        self.folder_paths = folder_paths
        self.output_folder = output_folder
        self.deprotonate = deprotonate
        self.iodine_files = []
        self.non_iodine_files = []
        self.failed_deprot = []
        self.failed_mapping = []
        self.h5_deprotonated_dir = 'h5_deprotonated_molecules'
        self.deprotonated_dir = 'deprotonated_molecules'

        os.makedirs(self.output_folder, exist_ok=True)
        os.makedirs(self.h5_deprotonated_dir, exist_ok=True)
        os.makedirs(self.deprotonated_dir, exist_ok=True)

        print(f"Created output directory: {self.output_folder}")
        print(f"Created H5 deprotonated directory: {self.h5_deprotonated_dir}")
        print(f"Created deprotonated molecules directory: {self.deprotonated_dir}")

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
                # calc = vlx.MolecularPropertyCalculator(folder, xyz_filename=filename)
                calc = vlx.MPI1(molecule, deprotonated=False)

                print('Starting calculations for', filename)

                data_matrix = calc.run_all_calculations()

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

            except Exception as e:
                logging.error(f"Error processing {filename}: {e}")
                
                print('Failed to process', filename)

    def deprotonate_and_process(self):
        if not self.deprotonate:
            return

        logging.info("Deprotonating molecules...")


        existing_h5_files = {f[:12] for f in os.listdir(self.h5_deprotonated_dir) if f.endswith('.h5')}
        self.non_iodine_files = [
            file for file in self.non_iodine_files
            if os.path.split(file)[1][:12] not in existing_h5_files
        ]

        for file in self.non_iodine_files:
            folder, filename = os.path.split(file)
            mol = Chem.MolFromXYZFile(file)
            try:
                rdDetermineBonds.DetermineBonds(mol, 0)
                deprot = vlx.OxygenDeprotonation(mol, filename, output_folder=self.deprotonated_dir)
                deprot.Deprotonate()
            except Exception as e:
                logging.error(f"Failed to deprotonate {filename}: {e}")
                self.failed_deprot.append(filename)

        logging.info("Calculating properties for deprotonated molecules...")

        for deprot_file in os.listdir(self.deprotonated_dir):
            if not deprot_file.endswith('.xyz'):
                continue

            prefix = deprot_file[:12]
            deprot_path = os.path.join(self.deprotonated_dir, deprot_file)

            for folder_path in self.folder_paths:
                prot_file = next((f for f in os.listdir(folder_path) if f.startswith(prefix)), None)
                if prot_file is None:
                    continue

                prot_path = os.path.join(folder_path, prot_file)
                prot_mol = Chem.MolFromXYZFile(prot_path)
                deprot_mol = Chem.MolFromXYZFile(deprot_path)

                if prot_mol is None or deprot_mol is None:
                    logging.warning(f"Could not load molecules: {prot_file} or {deprot_file}")
                    continue

                try:
                    mapping = vlx.MCS(prot_mol, deprot_mol).get_atom_mapping()
                    if mapping is None:
                        raise ValueError("Atom mapping failed.")
                except Exception as e:
                    logging.error(f"Mapping failed for {deprot_file}: {e}")
                    self.failed_mapping.append(deprot_file)
                    continue

                try:
                    mol = vlx.Molecule.read_xyz_file(deprot_path)
                    calc = vlx.MolecularPropertyCalculator(mol, deprotonated=True)
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

                    prot_h5_file = next((f for f in os.listdir(self.h5_protonated_dir) if f.startswith(prefix) and f.endswith('.h5')), None)
                    if prot_h5_file:
                        h5_path = os.path.join(self.h5_deprotonated_dir, prot_h5_file)
                        shutil.copy2(os.path.join(self.h5_protonated_dir, prot_h5_file), h5_path)
                        with h5py.File(h5_path, 'a') as f:
                            for key, value in reordered.items():
                                f.create_dataset(key, data=value)
                            f.create_dataset("xyz coordinates deprot", data=xyz)
                        logging.info(f"Added deprotonated data to {h5_path}")
                except Exception as e:
                    logging.error(f"Failed to process deprotonated {deprot_file}: {e}")

    def run(self):
        start = time.time()
        self.classify_molecules()
        self.calculate_properties()
        self.deprotonate_and_process()
        logging.info(f"Finished processing in {time.time() - start:.2f} seconds")

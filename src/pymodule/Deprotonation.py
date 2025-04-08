import rdkit
import os
from rdkit import Chem
from rdkit.Chem import Draw
import veloxchem as vlx
import numpy as np
from rdkit.Chem import AllChem
from io import BytesIO
from PIL import Image

class OxygenDeprotonation:

    """
    Class for deprotonation of oxygen atoms,
    and saving the new molecules in an xyz file.
    """
    def __init__(self, molecule, molecule_name, output_folder="deprotonated_molecules"):
        
        # self.folder_path = folder_path
        # self.xyz_filename = xyz_filename
        self.molecule = None
        self.mol = molecule
        self.basis = None
        self.deprotonated_molecules = None
        self.molecule_name = molecule_name
        self.output_folder = output_folder

        

    # def load_molecule(self):
    #     """
    #     Load the molecule from the xyz file.
    #     """
    #     file_path = os.path.join(self.folder_path, self.xyz_filename)
    #     if os.path.exists(file_path):
    #         self.molecule = vlx.Molecule.read_xyz_file(file_path)
    #         # self.mol = Chem.MolFromXYZFile(file_path)
    #         # self.molecule.show(atom_indices=True)
    #         self.basis = vlx.MolecularBasis.read(self.molecule, "def2-SVP")
    #     else:
    #         raise FileNotFoundError(f"Error: The file {file_path} does not exist.")
    
        
    def _find_hydroxyl_hydrogens(self):
        """
        Find the hydroxyl hydrogens in the molecule.
        """
        hydroxyl_hydrogens = [] # List of indices of hydroxyl hydrogens
        for atom in self.mol.GetAtoms():
            if atom.GetSymbol() == "O":
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetSymbol() == "H":
                        hydroxyl_hydrogens.append(neighbor.GetIdx())   
        
        # Convert the list into an array
        oxygen_hydrogens = np.array(hydroxyl_hydrogens)

        print(f"Hydrogen(s) indice(s) bonded to oxygen: {oxygen_hydrogens} with atom indices starting from 0")
        # print(f"Hydrogen(s) indice(s) bonded to oxygen: {oxygen_hydrogens + 1} with atom indices starting from 1")
        return oxygen_hydrogens
    
    def remove_hydroxyl_hydrogen(self):
        """
        Remove the hydroxyl hydrogen from the molecule.
        """
        

        self.deprotonated_molecules = [] # List to save the deprotonated molecules
        oxygen_hydrogen_index = self._find_hydroxyl_hydrogens()

        for idx in sorted(oxygen_hydrogen_index, reverse=True):
            # Convert into an editable molecule
            new_molecule = Chem.RWMol(self.mol)
            atom_index_to_remove = int(idx)
            new_molecule.RemoveAtom(atom_index_to_remove)
            self.deprotonated_molecules.append(new_molecule)

    def save_deprotonated_molecules_as_xyz(self, filename):
        """
        Save the deprotonated molecules as an xyz file.
        """
        for deprotonated_molecule in self.deprotonated_molecules:

            AllChem.EmbedMolecule(deprotonated_molecule, AllChem.ETKDG()) # Generate 3D coordinates
            conf = deprotonated_molecule.GetConformer() # Get the conformer

            # Write the file
            with open(filename, "w") as f:
                f.write(f"{deprotonated_molecule.GetNumAtoms()}\n")
                f.write("\n")
                for atom in deprotonated_molecule.GetAtoms():
                    pos = conf.GetAtomPosition(atom.GetIdx())
                    f.write(f"{atom.GetSymbol()} {pos.x} {pos.y} {pos.z}\n")

    def Deprotonate(self):
        """
        Deprotonate the oxygen atom and save the new molecules.
        """
        # self.load_molecule()
        self.remove_hydroxyl_hydrogen()

        # saved_files = []

        for i, file_name in enumerate(self.deprotonated_molecules):
            # filapath = os.path.join(self.folder_path, f"/Users/simonisaksson/Examensarbete/veloxchem/Atom mapping/deprotonated_molecules/ deprotonated_{molecule}.xyz")
            self.save_deprotonated_molecules_as_xyz(f"{self.output_folder}/{str(self.molecule_name.strip())}_{i}.xyz")
            # saved_files.append(file_name)
            print(f"Deprotonated molecule {i} saved as {self.molecule_name.strip()}.xyz")

        print('All Molecules with OH-groups have been deprotonated, molecules without this group is not saved as an xyz file.')
        print('Remember to set the charge to -1 for the molecule when loaded with VeloxChem')
        # return saved_files
        
        
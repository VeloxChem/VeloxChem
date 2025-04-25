import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
import os

class NitrogenProtonation:
    """
    Class for protonation of nitrogen atoms,
    and saving the new molecules in an xyz file.
    """

    def __init__(self, molecule, molecule_name, output_folder="protonated_molecules"):
        self.mol = molecule
        self.molecule_name = molecule_name
        self.output_folder = output_folder
        self.protonated_molecules = []

    def _find_nitrogens(self):
        """
        Find the nitrogen atoms in the molecule.
        """
        nitrogen_indices = [atom.GetIdx() for atom in self.mol.GetAtoms() if atom.GetSymbol() == "N"]
        print(f"Nitrogen atom indices: {nitrogen_indices}")
        return nitrogen_indices
    
    def protonate_nitrogens(self):
        nitrogen_indices = self._find_nitrogens()

        for idx in nitrogen_indices:
            new_mol = None

            # Create a deep copy of the molecule to avoid modifying the original one
            mol_copy = Chem.Mol(self.mol)  # Create a new molecule object
            mol_copy = Chem.Mol(mol_copy)  # Ensure a deep copy is made

            # Add hydrogens to the molecule
            mol_with_hs = Chem.AddHs(mol_copy)

            # Make it editable for modifications
            editable_mol = Chem.RWMol(mol_with_hs)

            # Add a new hydrogen atom to the nitrogen
            h_idx = editable_mol.AddAtom(Chem.Atom("H"))

            # Bond it to the nitrogen
            editable_mol.AddBond(idx, h_idx, order=Chem.rdchem.BondType.SINGLE)

            # Set the nitrogen to have a +1 formal charge
            editable_mol.GetAtomWithIdx(idx).SetFormalCharge(1)

            # Finalize and sanitize
            new_mol = editable_mol.GetMol()
            Chem.SanitizeMol(new_mol)

            # Generate 3D coordinates
            AllChem.EmbedMolecule(new_mol, AllChem.ETKDG())
            # new_mol = Chem.MolFromSmiles(Chem.MolToSmiles(new_mol))
            # Store the new protonated molecule
            self.protonated_molecules.append(new_mol)

            print(f"Protonated nitrogen at index {idx}.")
            print(f"Total protonated molecules: {len(self.protonated_molecules)}")


    def save_mol_as_xyz(self, mol, filename, directory=""):

        # Ensure hydrogens are added and coordinates exist
        mol = Chem.AddHs(mol)
        if mol.GetNumConformers() == 0:
            AllChem.EmbedMolecule(mol, AllChem.ETKDG())

        conf = mol.GetConformer()
        num_atoms = mol.GetNumAtoms()

        # Create the full path
        full_path = os.path.join(directory, filename)

        with open(full_path, "w") as f:
            f.write(f"{num_atoms}\n\n")  # XYZ header
            for i in range(num_atoms):
                atom = mol.GetAtomWithIdx(i)
                pos = conf.GetAtomPosition(i)
                f.write(f"{atom.GetSymbol()} {pos.x:.6f} {pos.y:.6f} {pos.z:.6f}\n")

    def Protonate(self):
        """
        Protonate and save new molecules.
        """
        self.protonate_nitrogens()

        for i, _ in enumerate(self.protonated_molecules):
            protonated_molecules = self.protonated_molecules
            
            base_name = os.path.splitext(self.molecule_name.strip())[0]
            filename = f"{base_name}_protonated_{i}.xyz"

            self.save_mol_as_xyz(protonated_molecules[i], filename, self.output_folder)


            #self.save_mol_as_xyz(protonated_molecules[i],
            #                      self.molecule_name.strip() + "_protonated_" + str(i) + '.xyz', self.output_folder)
               
            print(f"Protonated molecule {i} saved as:", filename)

        print("All nitrogen atoms have been individually protonated. Remember to adjust the charge if needed.")     
        
        return len(self.protonated_molecules)

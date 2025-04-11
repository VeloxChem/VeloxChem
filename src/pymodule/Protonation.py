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
        mol_with_hs = Chem.AddHs(self.mol)
        nitrogen_indices = self._find_nitrogens()

        for idx in nitrogen_indices:
            editable_mol = Chem.RWMol(mol_with_hs)

            # Add a new hydrogen atom
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

            # Store the new protonated molecule
            self.protonated_molecules.append(new_mol)


    def save_protonated_molecules_as_xyz(self, filename):
        """
        Save the protonated molecules as an xyz file.
        """
        # Ensure the directory exists
        directory = os.path.dirname(filename)
        if not os.path.exists(directory):
            os.makedirs(directory)

        for mol in self.protonated_molecules:
            AllChem.EmbedMolecule(mol, AllChem.ETKDG())
            conf = mol.GetConformer()

            with open(filename, "w") as f:
                f.write(f"{mol.GetNumAtoms()}\n\n")
                for atom in mol.GetAtoms():
                    pos = conf.GetAtomPosition(atom.GetIdx())
                    f.write(f"{atom.GetSymbol()} {pos.x:.6f} {pos.y:.6f} {pos.z:.6f}\n")

    def Protonate(self):
        """
        Protonate and save new molecules.
        """
        self.protonate_nitrogens()

        for i, _ in enumerate(self.protonated_molecules):
            self.save_protonated_molecules_as_xyz(
                f"{self.output_folder}/{self.molecule_name.strip()}_protonated_{i}.xyz"
            )
            print(f"Protonated molecule {i} saved as {self.molecule_name.strip()}_protonated_{i}.xyz")

        print("All nitrogen atoms have been individually protonated. Remember to adjust the charge if needed.")

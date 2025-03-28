from veloxchem.conformergenerator import ConformerGenerator
from veloxchem.molecule import Molecule
conf = ConformerGenerator()
conf.molecule = Molecule.read_smiles("CC1([C@@H](N2[C@H](S1)[C@@H](C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C")
conf.number_of_conformers_to_select = 1
conf.top_file_name = "MOL.top" #default top file name
conf.save_xyz = True #if True then save the xyz file in filtered folder
conf.em_tolerance_value = 1 #default is 10
conf.generate()
conf.show_global_minimum()
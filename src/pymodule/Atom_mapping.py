
import rdkit
from rdkit import Chem
from rdkit.Chem import Draw
import veloxchem as vlx
import numpy as np
from rdkit.Chem.Draw import MolDraw2DCairo
from io import BytesIO
from PIL import Image
# import cv2
from rdkit.Chem import rdFMCS
from rdkit.Chem import rdDetermineBonds

class MCS:
    def __init__(self, Protonated_molecule, Deprotonated_molecule, deprotonated = False, protonated = False):
        self.molecule1 = Protonated_molecule
        self.molecule2 = Deprotonated_molecule
        self.deprotonated = deprotonated
        self.protonated = protonated

    def give_the_molecules_personality(self):
        if self.deprotonated:
            self.mol_1 = self.molecule1
            rdDetermineBonds.DetermineBonds(self.mol_1, charge=0) # Determine bonds 
            self.mol_1.SetProp('_Name', 'Molecule 1')

            self.mol_2 = self.molecule2
            rdDetermineBonds.DetermineBonds(self.mol_2, charge=-1)
            self.mol_2.SetProp('_Name', 'Molecule 2')

            molecules = [self.mol_1, self.mol_2]
            return molecules
        
        if self.protonated:
            self.mol_1 = self.molecule1
            rdDetermineBonds.DetermineBonds(self.mol_1, charge=1) # Determine bonds 
            self.mol_1.SetProp('_Name', 'Molecule 1')

            self.mol_2 = self.molecule2
            rdDetermineBonds.DetermineBonds(self.mol_2, charge=0)
            self.mol_2.SetProp('_Name', 'Molecule 2')

            molecules = [self.mol_1, self.mol_2]
            return molecules


    def parameters(self):
        
        # Make sure that the molecules are sanitized
        Chem.SanitizeMol(self.mol_1)
        Chem.SanitizeMol(self.mol_2)

        params = rdFMCS.MCSParameters()
        params.AtomTyper = rdFMCS.AtomCompare.CompareElements
        params.BondTyper = rdFMCS.BondCompare.CompareOrder

        if self.molecule1.GetRingInfo().NumRings()>0: # if molecule1 has rings
            if self.molecule2.GetRingInfo().NumRings()>0: # if molecule2 has rings
                params.RingMatchesRingOnly = True
                params.CompleteRingsOnly = True
        else:
            params.RingMatchesRingOnly = False
            params.CompleteRingsOnly = False

        return params
    
    def show_mol(self, d2d,mol,legend='',highlightAtoms=[]):
        d2d.DrawMolecule(mol,legend=legend, highlightAtoms=highlightAtoms)
        d2d.FinishDrawing()
        bio = BytesIO(d2d.GetDrawingText())
        return Image.open(bio)

    # Function to Combine Images Horizontally
    def show_images(self, imgs,buffer=5):
        height = 0
        width = 0
        for img in imgs:
            height = max(height,img.height)
            width += img.width
        width += buffer*(len(imgs)-1)
        res = Image.new("RGBA",(width,height))
        x = 0
        for img in imgs:
            res.paste(img,(x,0))
            x += img.width + buffer
        return res

    def get_res(self): # Get the result
        molecules = self.give_the_molecules_personality()
        self.res = rdFMCS.FindMCS(molecules, self.parameters())


    def final_info_mcs(self):

        imgs = []

        d2d = MolDraw2DCairo(350,300)
        d2d.drawOptions().addAtomIndices = True

        dopts = d2d.drawOptions()
        dopts.setBackgroundColour((0,.9,.9,.3))

        highlightAtoms_mol_1 = self.mol_1.GetSubstructMatch(self.res.queryMol)

        imgs.append(self.show_mol(d2d, self.mol_1, self.mol_1.GetProp('_Name'), highlightAtoms_mol_1))

        # imgs.append(self.show_mol(d2d,self.mol_1,legend = self.mol_1.GetProp('_Name'), highlightAtoms = highlightAtoms_mol_1 ))

        d2d = MolDraw2DCairo(350,300)
        d2d.drawOptions().addAtomIndices = True

        dopts = d2d.drawOptions()
        dopts.setBackgroundColour((0,.9,.9,.3))

        highlightAtoms_product = self.mol_2.GetSubstructMatch(self.res.queryMol)

        imgs.append(self.show_mol(d2d,self.mol_2,legend = self.mol_2.GetProp('_Name'), highlightAtoms = highlightAtoms_product ))


        # alternative 
        self.show_images(imgs).save('output.png')

        # cv2.namedWindow("Output")
        # display(self.show_images(imgs))

        # cv2.imwrite('output.png' , cv2.cvtColor(np.array(self.show_images(imgs)), cv2.COLOR_BGR2RGB))


        print('\n\n')
        print('res.queryMol : ' , self.res.queryMol)
        print('res.queryMol.GetAtoms : ' , [ i.GetIdx() for i in self.res.queryMol.GetAtoms()])

        print('\n\n')
        print('Protonated Molecule : ' , highlightAtoms_mol_1)

        print('\n\n')
        print('Deprotonated Molecule : ' , highlightAtoms_product)

        print('\n\n')
        list_of_tuples = list(zip(highlightAtoms_mol_1, highlightAtoms_product)) 
        print(list_of_tuples)
        if len(list_of_tuples) != self.mol_2.GetNumAtoms():
            print('Could not find a match for all atoms in the molecule') 
            return None
        else:
            return list_of_tuples

    def get_atom_mapping(self):
        self.give_the_molecules_personality()
        self.parameters()
        self.get_res()
        return self.final_info_mcs()
    

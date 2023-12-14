import math
import numpy as np
# from veloxchem import mathconst_pi, hartree_in_kcalpermol, bohr_in_angstrom, Molecule
from .veloxchemlib import mathconst_pi
from .veloxchemlib import hartree_in_kcalpermol
from .veloxchemlib import bohr_in_angstrom
from .molecule import Molecule
from .seminario import Seminario
from .atomtypeidentifier import AtomTypeIdentifier

class Topology():

    """
    example:
    >>> top = Topology()
    >>> top.lead_forcefield(molecule,"path_to_gaff.dat")
    >>> top.load_charge(charge_array)
    >>> top.reparameterise(molecule, hessian)
    >>> top.write_itp(filename,output_folder)
    >>> top.write_top(filename,output_folder)
    """
    def __init__(self):
        ""

    def load_forcefield(self,molecule,ff_file_path):
        """
        Populates the topology based on the force-field parameters corresponding to a given molecule
        """
        ao_identifier = AtomTypeIdentifier()
        ao_identifier.generate_gaff_atomtypes(molecule)
        ff_data = ao_identifier.generate_force_field_dict(ff_file_path)
        self.atoms = ff_data['atomtypes']
        self.bonds = ff_data['bonds']
        self.angles = ff_data['angles']
        self.dihedrals = ff_data['dihedrals']
        self.impropers = ff_data['impropers']
        self.pairs = ff_data['pairs']
        return ao_identifier

    def update_charge(self,charges):
        """
        Takes list list of charges, and assigns those to the atoms 
        Ordering of the coordinates is assumed to be the same as the atoms
        """
        for atom_id,charge in zip(self.atoms,charges):
            self.atoms[atom_id]['charge']=float(charge)

    def get_itp_string(self,RES="Mol"):
        """
        generate an itp string based on this topology.

        Res: residue name used in the topology, defaults to "MOL"
        """
        # Initialize the itp file
        itp_string = '; Created by VeloxChem\n\n'

        # Print the atomtypes section
        itp_string += '\n[ atomtypes ]\n; name  bond_type  mass  charge  ptype  sigma  epsilon\n'
        for atom_id, atom_data in self.atoms.items():
            itp_string += f'{atom_data["type"]:<6} {atom_data["type"]:<10} {atom_data["mass"]:<7.2f} 0.0000  A {atom_data["sigma"]:11.4e} {atom_data["epsilon"]:11.4e}\n'

        # Print the moleculetype section
        itp_string += f"[moleculetype]\n; name  nrexcl\n{RES}  3\n\n"
        
        # Print the atoms section
        itp_string += '[ atoms ]\n'
        itp_string += '; nr    type  resnr  residue  atom  cgnr  charge  mass type2 charge2; comment - comment2\n'
        for atom_id, atom_data in self.atoms.items():
            itp_string += f'{atom_id:<6} {atom_data["type"]:>4}  1     {RES}     {atom_data["type"]:>4}  1     {0.00000:9.5f} {atom_data["mass"]:9.5f} ; {atom_data["comment"]}\n' #todo fix charge

        # Print the bonds section
        itp_string += '\n[ bonds ]\n; ai  aj  funct  r0 (nm)  fc (kJ/(mol nm2)) r02 fc2; comment - comment2\n'
        for bond_key, bond_data in self.bonds.items():
            itp_string += f'{bond_key[0]:<4} {bond_key[1]:<4} 1 {bond_data["eq"]:>6.6f} {bond_data["fc"]:>8.3f}; {bond_data["comment"]}\n'

        # Print the angles section
        itp_string += '\n[ angles ]\n; ai  aj  ak  funct  theta0 (degr)  fc (kJ/(mol rad2)) theta02 fc2; comment - comment2\n'
        for angle_key, angle_data in self.angles.items():
            itp_string += f'{angle_key[0]:<4} {angle_key[1]:<4} {angle_key[2]:<4} 1 {angle_data["eq"]:>6.3f} {angle_data["fc"]:>8.3f}; {angle_data["comment"]}\n'

        # Print the dihedrals section
        itp_string += '\n[ dihedrals ]\n; ai  aj  ak  al  funct  theta  k  mult; comment\n'#todo deal with other types of dihedrals
        for dihedral_key, dihedral_data in self.dihedrals.items():
            itp_string += f'{dihedral_key[0]:<4} {dihedral_key[1]:<4} {dihedral_key[2]:<4} {dihedral_key[3]:<4} 1 {dihedral_data["eq"]:>8.3f} {dihedral_data["fc"]:>8.3f}  1; {dihedral_data["comment"]}\n'

        # Print the impropers section
        itp_string += '\n[ dihedrals ] ; Improper dihedral section\n; ai  aj  ak  al  type  phi0  fc  n; comment\n'
        for improper_key, improper_data in self.impropers.items():
            itp_string += f'{improper_key[0]:<4} {improper_key[1]:<4} {improper_key[2]:<4} {improper_key[3]:<4} 4 {improper_data["eq"]:>8.3f} {improper_data["fc"]:>8.3f} {improper_data["periodicity"]:>2}; {improper_data["comment"]}\n'

        # Print the pairs section
        itp_string += '\n[ pairs ]\n; ai  aj  funct\n'
        for pair_key in self.pairs:
            itp_string += f'{pair_key[0]:<4} {pair_key[1]:<4} 1\n'
        return itp_string

    def write_itp(self,filename,output_folder,RES="MOL"):
        """
        Generate a GROMACS itp file based on this topology. Filename and RES will be set as instance variables

        Res: residue name used in the topology, defaults to "MOL". Should match with the residue name in the .top file.
        """
        if hasattr(self,filename):
            if self.filename != filename:
                print(f"WARNING: {self.filename} was set earlier as filename (probably in write_top), while {filename} was passed to write_itp")
        if hasattr(self,RES):
            if self.RES != RES:
                print(f"WARNING: {self.RES} was set earlier as residue name (probably in write_top), while {RES} was passed to write_itp")
        
        self.filename = filename
        self.RES=RES
        itp_string = self.get_itp_string(RES)
        with open(f"{output_folder}/{filename}.itp", "w") as f:
            f.write(itp_string)

    def get_top_string(self,filename,RES="MOL"):
        """
        Generate a GROMACS top file based on this topology. Filename and RES will be set as instance variables

        Res: residue name used in the topology, defaults to "MOL". Should match with the residue name in the .top file.
        """
        # Construct the topol.top file string
        # Initialize the top file
        top_string = '; Created by VeloxChem\n\n'

        # Print the include section for including the itp file
        top_string += f"""#include "{filename}.itp"\n"""

        # Print the defaults section
        top_string += '\n[ defaults ]\n; nbfunc  comb-rule  gen-pairs  fudgeLJ  fudgeQQ\n1  3  yes  0.5  0.8333\n'

        # Print the system directive
        top_string += f"\n[ system ]\n; name\n{RES}\n"

        # Print the molecules section
        top_string += f"\n[ molecules ]\n; name  number\n{RES}  1\n"
        return top_string

    def write_top(self,filename,output_folder,RES="MOL"):
        """
        Generate a GROMACS top file based on this topology. Filename and RES will be set as instance variables, and should match with the .itp file

        Res: residue name used in the topology, defaults to "MOL".
        """
        if hasattr(self,filename):
            if self.filename != filename:
                print(f"WARNING: {self.filename} was set earlier as filename (probably in write_itp), while {filename} was passed to write_top")
        if hasattr(self,RES):
            if self.RES != RES:
                print(f"WARNING: {self.RES} was set earlier as residue name (probably in write_itp), while {RES} was passed to write_top")

        top_string = self.get_top_string(filename,RES)

        with open(f"{output_folder}/{filename}.top", "w") as f:
            f.write(top_string)
    
    def reparameterise(self,mol,hes,keys=None,element=None,print_all=False,only_eq=False):
        """
        Reparameterises all unknown parameters with the seminario method using the given hessian

        mol: veloxchem molecule, used for geometry information
        hes: cartesian hessian matrix as 2D array
        keys: one or multiple tuples corresponding to bonds angles or dihedrals that need to be reparameterised. If passed, the unknown parameters will be ignored
        print_all: if true, all parameters will be printed, and not just the adjusted ones
        only_eq: only adjust the eq values, and not the fc values. Useful for transition states
        """
        print("\nVeloxchem force-field reparameterisation\n\n----------------------------------------")
        A_to_nm=0.1 #1 angstrom is 0.1 nm
        Bohr_to_nm = bohr_in_angstrom()*A_to_nm
        cal_to_joule = 4.184 #1 calorie is 4.184 joule
        Hartree_to_kJmol = hartree_in_kcalpermol()*cal_to_joule

        if (type(keys) == tuple):
            keys = [keys]

        coords = mol.get_coordinates_in_bohr() #Coordinates in Bohr
        sem = Seminario(hes,coords)

        def process_parameter(parameters, label):
            for ids in parameters.keys():
                types = parameters[ids]['types']
                eq = parameters[ids]["eq"]
                fc = parameters[ids]["fc"]

                if label == "bond":
                    neweq = AtomTypeIdentifier.measure_length(coords[ids[0]-1],coords[ids[1]-1])
                    neweq *= Bohr_to_nm #Convert to nm

                    newfc = sem.bond_fc(ids[0]-1,ids[1]-1) #fc in H/Bohr^2
                    newfc *= Hartree_to_kJmol/(Bohr_to_nm**2) #Convert to kJ/mol nm^2
                elif label == "angle":
                    neweq = AtomTypeIdentifier.measure_angle(coords[ids[0]-1],coords[ids[1]-1],coords[ids[2]-1])
                    newfc = sem.angle_fc(ids[0]-1,ids[1]-1,ids[2]-1) #Returns in H/rad^2 i think
                    newfc *= Hartree_to_kJmol #Convert to kJmol^-1/rad^2
                elif label == "dihedral":
                    phase = self.dihedrals[ids]["phase"]
                    perio = self.dihedrals[ids]["periodicity"]
                elif label == "improper":
                    perio = self.impropers[ids]["periodicity"]
                    neweq = AtomTypeIdentifier.measure_dihedral(coords[ids[0]-1],coords[ids[1]-1],coords[ids[2]-1],coords[ids[3]-1])

                    newfc = sem.dihed_fc(ids[0]-1,ids[1]-1,ids[2]-1,ids[3]-1) #Returns in H/rad^2 i think
                    newfc *= Hartree_to_kJmol #Convert to kJmol^-1/rad^2

                comment = parameters[ids]["comment"]
                
                repar = False
                if keys is not None:
                    if ids in keys:
                        repar = True
                elif keys is None and element is not None:
                    for id in ids:
                        if element == mol.get_labels()[id-1]:
                            repar = True
                        
                elif keys is None and element is None and comment =="unknown":
                    repar = True
                

                if repar:
                    if not label == "dihedral":
                        parameters[ids]["eq"] = neweq
                        if not only_eq:
                            parameters[ids]["fc"] = newfc
                            new_comment = f"reparameterised from (eq,fc)=({eq:>.1f},{fc:.1f}) to (eq,fc)=({neweq:.1f},{newfc:.1f})"
                        else:
                            new_comment = f"reparameterised from eq={eq:>.1f} to eq={neweq:.1f}"
                        parameters[ids]["comment"] = new_comment
                        comment = new_comment
                    else:
                        comment+=", reparameterisation of proper dihedrals should be done with other method"

                if repar or print_all:
                    if label == "bond":
                        print(f"{f'{ids}{types}:':<{20}}\t{eq:>6f}\t{fc:>6f}\t\t{neweq:>6f}\t{newfc:>6f}\t\t{comment}")
                    if label == "angle":
                        print(f"{f'{ids}{types}:':<{30}}\t{eq:>6f}\t{fc:>6f}\t\t{neweq:>6f}\t{newfc:>6f}\t\t{comment}")
                    if label == "dihedral":
                        print(f"{f'{ids}{types}: ':<{40}}{eq:>6f}\t{fc:>6f}\t\t{phase:>6f}\t{perio}\t\t{comment}")
                    if label == "improper":
                        print(f"{f'{ids}{types}: ':<{40}}{eq:>6f}\t{fc:>6f}\t\t{perio}\t\t{neweq:.<3f}\t{newfc:.<3f}\t\t{perio}\t\t{comment}")

        print("Bond \t\t\tr0(nm) \t\tfc(kJmol^-1/nm^2) \tNew r0(nm) \tNew fc(kJmol^-1/nm^2) \tComment")
        process_parameter(self.bonds, "bond")

        print("\nAngle \t\t\t\ttheta0(deg)  \tfc(kJmol^-1/rad^2) \tNew theta0(deg)\tNew fc(kJmol^-1/rad^2) \tComment")
        process_parameter(self.angles, "angle")

        print("\nDihedrals\t\t\t\teq(deg)  \tfc(kJmol^-1/rad^2?) \tPhase(deg) \tPeriodicity \tComment")
        process_parameter(self.dihedrals, "dihedral")

        print("\nImpropers\t\t\t\teq(deg)  \tfc(kJmol^-1/rad^2?) \tperiodicity \tNew eq(deg) \tNew fc(kJmol^-1/rad^2?)\tNew periodicity\tComment")
        process_parameter(self.impropers, "improper")

    @staticmethod
    def write_FEP_itp(filename,output_folder,topA,topB,molA,molB,RES="MOL"):
        """
        Writes a gromacs FEP incorporating all bonds in both topA and topB
        """
        topA.filename = filename
        topA.RES=RES
        
        # Initialize the itp file
        itp_AB = '; Created by VeloxChem\n\n'
        itp_A = itp_AB
        itp_B = itp_AB

        coordsA = molA.get_coordinates()
        coordsB = molB.get_coordinates()

        # Print the atomtypes section
        atomtypes_sec = '\n[ atomtypes ]\n; name  bond_type  mass  charge  ptype  sigma  epsilon\n'
        #Atoms in both topologies should be the same #todo implement checks
        for atom_id, atom_data in topA.atoms.items():
            atomtypes_sec += f'{atom_data["type"]:<6} {atom_data["type"]:<10} {atom_data["mass"]:<7.2f} 0.0000  A {atom_data["sigma"]:11.4e} {atom_data["epsilon"]:11.4e}\n'

        # Print the moleculetype section
        moltype_sec = f"[moleculetype]\n; name  nrexcl\{RES}  3\n\n"

        itp_AB += atomtypes_sec + moltype_sec 
        itp_A  += atomtypes_sec + moltype_sec 
        itp_B  += atomtypes_sec + moltype_sec 

        # Print the atoms section
        atoms_sec = '[ atoms ]\n ; nr    type  resnr  residue  atom  cgnr  charge  mass type2 charge2; comment - comment2\n'
        atoms_secA = atoms_sec
        atoms_secB = atoms_sec
        #Atoms in both topologies should be the same #todo implement checks
        for atom_id, atom_data in topA.atoms.items():
            atomA = f'{atom_id:<6} {atom_data["type"]:>4}  1     {RES}     {atom_data["type"]:>4}  1     {0.00000:9.5f} {atom_data["mass"]:9.5f}' #todo fix charge
            atomB = f'{atom_id:<6} {atom_data["type"]:>4}  1     {RES}     {atom_data["type"]:>4}  1     {0.00000:9.5f} {atom_data["mass"]:9.5f}' #todo fix charge
            atom = atomA+f' {topB.atoms[atom_id]["type"]:>4} {0.00000:9.5f}'
            
            commentA =f' {topA.atoms[atom_id]["comment"]}'
            commentB =f' {topB.atoms[atom_id]["comment"]}'
            
            atoms_secA += atomA+";"+commentA+'\n'
            atoms_secB += atomB+";"+commentB+'\n'
            atoms_sec += atom+";"+commentA+" -"+commentB+'\n'
        
        itp_AB += atoms_sec
        itp_A  += atoms_secA
        itp_B  += atoms_secB

        #todo refactor this code
        # Print the bonds section
        sec = '\n[ bonds ]\n; ai  aj  funct  r0 (nm)  fc (kJ/(mol nm2)) r02 fc2; comment - comment2\n'
        secA = sec
        secB = sec
        keys = topA.bonds.keys() | topB.bonds.keys()
        for key in keys:
            id = f'{key[0]:<4} {key[1]:<4} 1'
            
            Aflag = key in topA.bonds.keys()
            Bflag = key in topB.bonds.keys()
            
            if Aflag: 
                parA = f' {topA.bonds[key]["eq"]:>6.6f} {topA.bonds[key]["fc"]:>8.3f}'
                commentA =f' {topA.bonds[key]["comment"]}'
            else: 
                bondlength = AtomTypeIdentifier.measure_length(coordsA[key[0]-1],coordsA[key[1]-1])
                bondlength *= bohr_in_angstrom()*0.1
                parA = f' {bondlength:>6.6f} 0'
                commentA = ' Nonbonded'

            if Bflag: 
                parB = f' {topB.bonds[key]["eq"]:>6.6f} {topB.bonds[key]["fc"]:>8.3f}'
                commentB =f' {topB.bonds[key]["comment"]}'
            else: 
                bondlength = AtomTypeIdentifier.measure_length(coordsB[key[0]-1],coordsB[key[1]-1])
                bondlength *= bohr_in_angstrom()*0.1
                parB = f' {bondlength:>6.6f} 0'
                commentB = ' Nonbonded'
            
            secA += id+parA+";"+commentA+'\n'
            secB += id+parB+";"+commentB+'\n'
            sec += id+parA+parB+";"+commentA+" -"+commentB+'\n'
        
        itp_AB += sec
        itp_A  += secA
        itp_B  += secB

        # Print the angles section
        sec = '\n[ angles ]\n; ai  aj  ak  funct  theta0 (degr)  fc (kJ/(mol rad2)) theta02 fc2; comment - comment2\n'
        secA = sec
        secB = sec
        keys = topA.angles.keys() | topB.angles.keys()
        for key in keys:
            Aflag = key in topA.angles.keys()
            Bflag = key in topB.angles.keys()

            id = f'{key[0]:<4} {key[1]:<4} {key[2]:<4} 1 '
            if Aflag: 
                parA = f' {topA.angles[key]["eq"]:>6.6f} {topA.angles[key]["fc"]:>8.3f}'
                commentA = f' {topA.angles[key]["comment"]}'
            else:
                angle = AtomTypeIdentifier.measure_angle(coordsA[key[0]-1],coordsA[key[1]-1],coordsA[key[2]-1])
                parA = f' {angle:>6.6f} 0'
                commentA = ' Nonbonded'

            if Bflag:
                parB = f' {topB.angles[key]["eq"]:>6.6f} {topB.angles[key]["fc"]:>8.3f}'
                commentB = f' {topB.angles[key]["comment"]}'
            else:
                angle = AtomTypeIdentifier.measure_angle(coordsB[key[0]-1],coordsB[key[1]-1],coordsB[key[2]-1])
                parB = f' {topA.angles[key]["eq"]:>6.6f} 0'
                commentB = ' Nonbonded'

            secA += id+parA+';'+commentA+'\n'
            secB += id+parB+';'+commentB+'\n'
            sec += id+parA+parB+';'+commentA+" -"+commentB+'\n'

        itp_AB += sec
        itp_A+=secA
        itp_B+=secB
        
        # Print the dihedrals section
        sec = '\n[ dihedrals ]\n; ai  aj  ak  al  funct  theta  k  mult theta2 k2 mult2; comment - comment2 ; WARNING: theta2, k2 and mult2 might be wrong in the wrong order or might be missing parameters\n'
        secA = sec
        secB = sec

        keys = topA.dihedrals.keys() | topB.dihedrals.keys()
        for key in  keys:
            Aflag = key in topA.dihedrals.keys()
            Bflag = key in topB.dihedrals.keys()

            id = f'{key[0]:<4} {key[1]:<4} {key[2]:<4} {key[3]:<4} 1'

            if Aflag:
                parA = f' {topA.dihedrals[key]["eq"]:>8.3f} {topA.dihedrals[key]["fc"]:>8.3f} 1'
                commentA = f' {topA.dihedrals[key]["comment"]}' 
            else: 
                angle = AtomTypeIdentifier.measure_dihedral(coordsA[key[0]-1],coordsA[key[1]-1],coordsA[key[2]-1],coordsA[key[3]-1])
                parA = f' {angle:>8.3f} 0 1'
                commentA = ' Nonbonded'
            if Bflag: 
                parB = f' {topB.dihedrals[key]["eq"]:>8.3f} {topB.dihedrals[key]["fc"]:>8.3f} 1'
                commentB = f' {topB.dihedrals[key]["comment"]}'
            else: 
                angle = AtomTypeIdentifier.measure_dihedral(coordsB[key[0]-1],coordsB[key[1]-1],coordsB[key[2]-1],coordsB[key[3]-1])
                parB = f' {angle:>8.3f} 0 1'
                commentB = ' Nonbonded'
        
            secA += id+parA+';'+commentA+'\n'
            secB += id+parB+';'+commentB+'\n'
            sec += id+parA+parB+';'+commentA+" -"+commentB+'\n'

        itp_AB += sec
        itp_A+=secA
        itp_B+=secB

        # Print the impropers section
        sec = '\n[ dihedrals ] ; Improper dihedral section\n; ai  aj  ak  al  type  phi0  fc  n phi02 fc2 n2; comment - comment2 ; WARNING: phi02, fc2 and n2 might be wrong in the wrong order or might be missing parameters\n'
        secA = sec
        secB = sec
        keys = topA.impropers.keys() | topB.impropers.keys()
        for key in keys:
            Aflag = key in topA.impropers.keys()
            Bflag = key in topB.impropers.keys()

            id = f'{key[0]:<4} {key[1]:<4} {key[2]:<4} {key[3]:<4} 4'
            
            if Aflag:
                parA = f' {topA.impropers[key]["eq"]:>8.3f} {topA.impropers[key]["fc"]:>8.3f} {topA.impropers[key]["periodicity"]:>2}'
                commentA = f' {topA.impropers[key]["comment"]}'
            else:
                angle = AtomTypeIdentifier.measure_dihedral(coordsA[key[0]-1],coordsA[key[1]-1],coordsA[key[2]-1],coordsA[key[3]-1])
                parA = f' {angle:>8.3f} 0 1' 
                commentA = ' Nonbonded'
            if Bflag: 
                parB = f' {topB.impropers[key]["eq"]:>8.3f} {topB.impropers[key]["fc"]:>8.3f} {topB.impropers[key]["periodicity"]:>2}'
                commentA = f' {topB.impropers[key]["comment"]}'
            else:
                angle = AtomTypeIdentifier.measure_dihedral(coordsB[key[0]-1],coordsB[key[1]-1],coordsB[key[2]-1],coordsB[key[3]-1])
                parB = f' {angle} 0 1' 
                commentB = ' Nonbonded'

            secA += id+parA+';'+commentA+'\n'
            secB += id+parB+';'+commentB+'\n'
            sec += id+parA+parB+';'+commentA+" -"+commentB+'\n'

        itp_AB += sec
        itp_A+=secA
        itp_B+=secB

        # Print the pairs section
        itp_AB += '\n[ pairs ]\n; ai  aj  funct\n'
        for pair_key in topA.pairs:
            itp_AB += f'{pair_key[0]:<4} {pair_key[1]:<4} 1\n'
            #todo add pairs for topB, whattodo with the pairs at all?
            
        with open(f"{output_folder}/{filename}AB.itp", "w") as f:
            f.write(itp_AB)
        with open(f"{output_folder}/{filename}A.itp", "w") as f:
            f.write(itp_A)
        with open(f"{output_folder}/{filename}B.itp", "w") as f:
            f.write(itp_B)
    #Utility functions for measuring parameters
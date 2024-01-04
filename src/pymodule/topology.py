#
#                           VELOXCHEM 1.0-RC3
#         ----------------------------------------------------
#                     An Electronic Structure Code
#
#  Copyright © 2018-2022 by VeloxChem developers. All rights reserved.
#  Contact: https://veloxchem.org/contact
#
#  SPDX-License-Identifier: LGPL-3.0-or-later
#
#  This file is part of VeloxChem.
#
#  VeloxChem is free software: you can redistribute it and/or modify it under
#  the terms of the GNU Lesser General Public License as published by the Free
#  Software Foundation, either version 3 of the License, or (at your option)
#  any later version.
#
#  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
#  License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.

import numpy as np
from .veloxchemlib import mathconst_pi
from .veloxchemlib import hartree_in_kcalpermol
from .veloxchemlib import bohr_in_angstrom
from .molecule import Molecule
from .seminario import Seminario
from .atomtypeidentifier import AtomTypeIdentifier

class Topology:

    """
    Class for creating a topology based on a molecule and a force-field file

    example:
    >>> top = Topology()
    >>> top.identify_atomtypes(molecule,"path_to_gaff.dat")
    >>> top.create_topology()
    >>> top.write_gromacs_files("mol","output_folder")
    """

    def __init__(self):
        """
        Initialize the topology
        """
        self.atoms = {}
        self.bonds = {}
        self.angles = {}
        self.dihedrals = {}
        self.impropers = {}
        self.pairs = []
        self.coordinates = []
        self.filename = ""
        self.RES = ""   
    
    # TODO: Add the ffld parser function
    def parse_ffld(self, molecule, itp, ffld):
        ffld_data = ffld_parser.parse(itp,ffld,molecule)
        masses = molecule.masses_to_numpy()
        for i,id in enumerate(ffld_data['atoms']): #todo do i need to add masses manually?
            ffld_data['atoms'][id]['mass'] =masses[i]
        self.atoms = ffld_data['atoms']
        self.bonds = ffld_data['bonds']
        self.angles = ffld_data['angles']
        self.dihedrals = ffld_data['dihedrals']
        self.impropers = ffld_data['impropers']
        self.pairs = ffld_data['pairs']

    def identify_atomtypes(self,molecule,ff_file_path):
        """
        Uses the atomtype identifier to identify the atomtypes of a molecule
        This method can be used to correct missidentifications in the force-field file

        Example:
        >>> top = Topology()
        >>> top.identify_atomtypes(molecule,"path_to_gaff.dat")
        >>> top.atomtypes #list of atomtypes
        If the atom type in the 4th atom is wrong, you can correct it with:
        >>> top.atomtypes[3] = "ca"

        molecule: Molecule object from veloxchem
        ff_file_path: path to the force-field file

        Returns: a list of atomtypes       
        """

        self.ao_identifier = AtomTypeIdentifier()
        self.atomtypes = self.ao_identifier.generate_gaff_atomtypes(molecule)
        self.ao_identifier.check_for_bad_assignations(ff_file_path)

        self.molecule = molecule
        self.ff_file_path = ff_file_path

    def create_topology(self):
        """
        Creates a topology based on a molecule and a force-field file
        The atomtypes should be identified first with identify_atomtypes
        """

        # Print an error if the atomtypes are not identified
        if self.ao_identifier is None:
            raise Exception("Atomtypes are not identified yet. Use the identify_atomtypes method first")
        
        # Generate the force-field dictionary
        # Method from atomtypeidentifier.py to analize the topology of the molecule 
        ff_data = self.ao_identifier.generate_force_field_dict(self.ff_file_path)
        masses = self.molecule.masses_to_numpy()

        for i,id in enumerate(ff_data['atomtypes']):
            ff_data['atomtypes'][id]['mass'] =masses[i]
        self.atoms = ff_data['atomtypes']
        self.bonds = ff_data['bonds']
        self.angles = ff_data['angles']
        self.dihedrals = ff_data['dihedrals']
        self.impropers = ff_data['impropers']
        self.pairs = ff_data['pairs']
        self.coordinates = self.molecule.get_coordinates_in_angstrom() * 0.1
    
    # Methods for updating the topology
    def update_charge(self,charges):
        """
        Takes list list of charges, and assigns those to the atoms 
        Ordering of the coordinates is assumed to be the same as the atoms
        """
        for atom_id,charge in zip(self.atoms,charges):
            self.atoms[atom_id]['charge']=float(charge)

    # Methods for generating GROMACS the topology
            
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
            itp_string += f'{atom_data["type"]:<6} {atom_data["type"]:<10} {atom_data["mass"]:<7.2f} {atom_data["charge"]:9.5f}  A {atom_data["sigma"]:11.4e} {atom_data["epsilon"]:11.4e}\n'

        # Print the moleculetype section
        itp_string += f"[moleculetype]\n; name  nrexcl\n{RES}  3\n\n"
        
        # Print the atoms section
        itp_string += '[ atoms ]\n'
        itp_string += '; nr    type  resnr  residue  atom  cgnr  charge  mass; comment\n'
        for atom_id, atom_data in self.atoms.items():
            itp_string += f'{atom_id:<6} {atom_data["type"]:>4}  1     {RES}     {atom_data["type"]:>4}  1     {atom_data["charge"]:9.5f} {atom_data["mass"]:9.5f} ; {atom_data["comment"]}\n'

        # Print the bonds section
        itp_string += '\n[ bonds ]\n; ai  aj  funct  r0 (nm)  fc (kJ/(mol nm2)); comment\n'
        for bond_key, bond_data in self.bonds.items():
            itp_string += f'{bond_key[0]:<4} {bond_key[1]:<4} 1 {bond_data["eq"]:>6.6f} {bond_data["fc"]:>8.3f}; {bond_data["comment"]}\n'

        # Print the angles section
        itp_string += '\n[ angles ]\n; ai  aj  ak  funct  theta0 (degr)  fc (kJ/(mol rad2)); comment\n'
        for angle_key, angle_data in self.angles.items():
            itp_string += f'{angle_key[0]:<4} {angle_key[1]:<4} {angle_key[2]:<4} 1 {angle_data["eq"]:>6.3f} {angle_data["fc"]:>8.3f}; {angle_data["comment"]}\n'

        # Print the dihedrals section
        itp_string += '\n[ dihedrals ]\n; ai  aj  ak  al  funct  theta  k  mult; comment\n'
        for dihedral_key, dihedral_data in self.dihedrals.items():
            if dihedral_data["type"]==9:
                itp_string += f'{dihedral_key[0]:<4} {dihedral_key[1]:<4} {dihedral_key[2]:<4} {dihedral_key[3]:<4} {dihedral_data["type"]} {dihedral_data["eq"]:>8.3f} {dihedral_data["fc"]:>8.3f}  {abs(dihedral_data["periodicity"])}; {dihedral_data["comment"]}\n'
            elif dihedral_data["type"]==5:
                itp_string += f'{dihedral_key[0]:<4} {dihedral_key[1]:<4} {dihedral_key[2]:<4} {dihedral_key[3]:<4} {dihedral_data["type"]} {dihedral_data["c1"]} {dihedral_data["c2"]} {dihedral_data["c3"]} {dihedral_data["c4"]} {dihedral_data["c5"]} ; {dihedral_data["comment"]}\n'
            else:
                print(f"ERROR: dihedral type {dihedral_data['type']} not supported") #todo raise exception?

        # Print the impropers section
        itp_string += '\n[ dihedrals ] ; Improper dihedral section\n; ai  aj  ak  al  type  phi0  fc  n; comment\n'
        for improper_key, improper_data in self.impropers.items():
            if improper_data['type']==4:
                itp_string += f'{improper_key[0]:<4} {improper_key[1]:<4} {improper_key[2]:<4} {improper_key[3]:<4} {improper_data["type"]} {improper_data["eq"]:>8.3f} {improper_data["fc"]:>8.3f} {improper_data["periodicity"]:>2}; {improper_data["comment"]}\n'
            elif improper_data['type']==2:
                itp_string += f'{improper_key[0]:<4} {improper_key[1]:<4} {improper_key[2]:<4} {improper_key[3]:<4} {improper_data["type"]} {improper_data["eq"]:>8.3f} {improper_data["fc"]:>8.3f}; {improper_data["comment"]}\n'
            else:
                print(f"ERROR: dihedral type {improper_data['type']} not supported") #todo raise exception?

        # Print the pairs section
        itp_string += '\n[ pairs ]\n; ai  aj  funct\n'
        for pair_key in self.pairs:
            itp_string += f'{pair_key[0]:<4} {pair_key[1]:<4} 1\n'
        return itp_string

    def get_top_string(self,filename,RES="MOL"):
        """
        Generate a GROMACS top file based on this topology. Filename and RES will be set as instance variables

        Res: residue name used in the topology, defaults to "MOL". Should match with the residue name in the .top file.
        """
        # Construct the topol.top file string
        # Initialize the top file
        top_string = '; Created by VeloxChem\n\n'

        # Print the include section for including the itp file
        top_string += " "
        top_string += f"""#include "{filename}.itp"\n"""

        # Print the defaults section all aligned
        top_string += '[ defaults ]'
        top_string += ';nbfunc comb-rule gen-pairs fudgeLJ fudgeQQ'
        # Values aligned with the comment
        top_string += '\n1        2        yes       0.5     0.8333\n'

        # Print the system directive
        top_string += f"\n[ system ]\n; name\n{RES}\n"

        # Print the molecules section
        top_string += f"\n[ molecules ]\n; name  number\n{RES}  1\n"
        return top_string
    
    # Methods for writing the GROMACS files
    
    def write_itp(self, filename, output_folder, RES="MOL"):
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

    def write_top(self, filename, output_folder, RES="MOL"):
        """
        Generate a GROMACS top file based on this topology. Filename and RES will be set as instance variables

        Res: residue name used in the topology, defaults to "MOL". Should match with the residue name in the .itp file.
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

    def write_gro(self, filename, output_folder, RES='MOL'):
        """
        Writes a GROMACS gro file based on the coordinates of the molecule
        and the atom types in the topology.

        :param filename: Name of the gro file
        :param output_folder: Output directory for the gro file
        :param RES: Name of the residue
        """
        # Check if coordinates are available
        if not hasattr(self, 'coordinates'):
            raise ValueError("Coordinates are not set for the topology.")

        # Construct the gro file content
        gro_string = f'Generated by VeloxChem\n'
        gro_string += f'{len(self.atoms)}\n'

        # Atom lines formatting
        for atom_id, atom_data in enumerate(self.atoms.values(), start=1):
            x, y, z = self.coordinates[atom_id - 1]
            gro_string += f'{1 % 100000:>5d}{RES:<5}{atom_data["type"]:>5}{atom_id % 100000:>5d}{x:8.3f}{y:8.3f}{z:8.3f}\n'
        
        # TODO: Replace with actual box dimensions if available
        gro_string += '   2.00000   2.00000   2.00000\n'

        # Write to file
        with open(f"{output_folder}/{filename}.gro", "w") as f:
            f.write(gro_string)

    # This method is a wrapper for the other write methods

    def write_gromacs_files(self, filename, output_folder, RES="MOL"):
        """
        Generate a GROMACS top and itp file based on this topology. Filename and RES will be set as instance variables, and should match with the .itp file

        Res: residue name used in the topology, defaults to "MOL".
        """
        self.write_itp(filename,output_folder,RES)
        self.write_top(filename,output_folder,RES)
        self.write_gro(filename,output_folder,RES)

    # Methods for reparameterizing the force-field
    
    def reparameterize(self, molecule, hessian, origin='ORCA', keys=None, element=None, print_all=False, only_eq=False, no_repar=False, repar_imp=False):
        """
        Reparameterizes all unknown parameters with the seminario method using the given hessian

        molecule: Molecule object from veloxchem
        hessian: Hessian matrix in a 2D numpy array or path to the hessian file
        origin: origin of the hessian, either 'ORCA' or 'VeloxChem'
            If the origin is ORCA, the hessian will be parsed from the file
            If the origin is VeloxChem, the hessian should be already a 2D numpy array
        keys: list of tuples of atom ids, for which the parameters should be reparameterised
        element: element for which all parameters should be reparameterised
        print_all: print all parameters, not just the ones that are reparameterised
        only_eq: only reparameterise the equilibrium values, not the force constants
        no_repar: do not reparameterise anything, just print the parameters
        repar_imp: reparameterise impropers as well (experimental)

        """

        print(
        """
        VeloxChem Force-Field Reparameterization
        Based on the Seminario method

        ------------------------------------------------------------------------------------------
        References:
        J. M. Seminario (1996) Calculation of intramolecular force fields from second-derivative tensors. 
        Internat. J. Quant. Chem. 60:1271-1277
        """
        )

        A_to_nm=0.1 #1 angstrom is 0.1 nm
        Bohr_to_nm = bohr_in_angstrom()*A_to_nm
        cal_to_joule = 4.184 #1 calorie is 4.184 joule
        Hartree_to_kJmol = hartree_in_kcalpermol()*cal_to_joule

        if (type(keys) == tuple):
            keys = [keys]

        if origin == 'ORCA':
            hessian = Seminario.parse_orca_hessian(hessian)
        elif origin == 'VeloxChem':
            hessian = np.array(hessian)
        else:
            raise Exception("Please specify the origin of the hessian as 'ORCA' or 'VeloxChem'")
        
        coords = molecule.get_coordinates_in_bohr() #Coordinates in Bohr
        sem = Seminario(hessian,coords)

        def process_parameter(parameters, label):
            for ids in parameters.keys():
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
                    perio = self.dihedrals[ids]["periodicity"]
                elif label == "improper":
                    perio = self.impropers[ids]["periodicity"]
                    if repar_imp:
                        neweq = AtomTypeIdentifier.measure_dihedral(coords[ids[0]-1],coords[ids[1]-1],coords[ids[2]-1],coords[ids[3]-1])
                        newfc = sem.dihed_fc(ids[0]-1,ids[1]-1,ids[2]-1,ids[3]-1,avg_improper=True) #Returns in H/rad^2 i think
                        newfc *= Hartree_to_kJmol #Convert to kJmol^-1/rad^2
                    else:
                        neweq = 0
                        newfc = 0

                comment = parameters[ids]["comment"] #GAFF2 or unknown
                
                repar = False
                if no_repar:
                    ""
                elif keys is not None:
                    if ids in keys:
                        repar = True
                elif keys is None and element is not None:
                    for id in ids:
                        if element == molecule.get_labels()[id-1]:
                            repar = True
                elif keys is None and element is None and "unknown" in comment:
                    repar = True

                if repar:
                    if label =="dihedral" and not repar_imp:
                        comment += ", not reparameterizing impropers"
                    elif not label == "dihedral":
                        parameters[ids]["eq"] = neweq
                        if not only_eq:
                            parameters[ids]["fc"] = newfc
                            comment += f", reparameterized from (eq,fc)=({eq:>.1f},{fc:.1f}) to (eq,fc)=({neweq:.1f},{newfc:.1f})"
                        else:
                            comment += f", reparameterized from eq={eq:>.1f} to eq={neweq:.1f}"
                        parameters[ids]["comment"] = comment

                        if label =="improper":
                            parameters[ids]["type"] =2    
                            parameters[ids]["periodicity"] = ""                        
                    else:
                        comment+=", reparameterization of proper dihedrals should be done with other method"

                if repar or print_all:
                    if label == "bond":
                        print(f"{f'{ids}:':<{10}}\t{eq:>6f}\t{fc:>6f}\t\t{neweq:>6f}\t{newfc:>6f}\t\t{comment}")
                    if label == "angle":
                        print(f"{f'{ids}:':<{15}}\t{eq:>6f}\t{fc:>6f}\t\t{neweq:>6f}\t{newfc:>6f}\t\t{comment}")
                    if label == "dihedral":
                        print(f"{f'{ids}: ':<{20}}{eq:>6f}\t{fc:>6f}\t\t{abs(perio)}\t\t{comment}")
                    if label == "improper":
                        print(f"{f'{ids}: ':<{20}}{eq:>6f}\t{fc:>6f}\t\t{perio}\t\t{neweq:.<3f}\t{newfc:.<3f}\t\t{perio}\t\t{comment}")

        print("Bond \t\tr0(nm) \t\tfc(kJmol^-1/nm^2) \tNew r0(nm) \tNew fc(kJmol^-1/nm^2) \tComment")
        process_parameter(self.bonds, "bond")

        print("\nAngle \t\ttheta0(deg)  \tfc(kJmol^-1/rad^2) \tNew theta0(deg)\tNew fc(kJmol^-1/rad^2) \tComment")
        process_parameter(self.angles, "angle")

        print("\nDihedrals\t\teq(deg)  \tfc(kJmol^-1) \tPeriodicity \tComment")
        process_parameter(self.dihedrals, "dihedral")

        print("\nImpropers\t\teq(deg)  \tfc(kJmol^-1/rad^2?) \tperiodicity \tNew eq(deg) \tNew fc(kJmol^-1/rad^2?)\tNew periodicity\tComment")
        process_parameter(self.impropers, "improper")

    # TODO Check this method because the specific libraries to be imported are 
    # not clear.
        
    def test_force_field(self,filename, output_folder, save_trajectory=False, show_output=False):
        """
        This method will perfom a short MD simulation and measure the RMSD
        of the molecule during the simulation compared to the original molecule.

        This method requires the GROMACS force-field files to be present in the current directory.
        """
        # TODO: USE VELOXCHEM INSTEAD OF OPENMM

        from openmm.app import GromacsGroFile, GromacsTopFile, PDBReporter, StateDataReporter, Simulation, PME, HBonds
        from openmm import LangevinMiddleIntegrator
        from openmm.unit import kelvin, picosecond, nanometer

        from sys import stdout
        import mdtraj as md
        import numpy as np
        import os
    
        # Load GROMACS files
        # Check if the GROMACS files are present
        if not os.path.isfile(f"{filename}.gro"):
            raise ValueError(f"{filename}.gro is not present in the current directory")
        if not os.path.isfile(f"{filename}.top"):
            raise ValueError(f"{filename}.top is not present in the current directory")
        
        gro = GromacsGroFile(f"{filename}.gro")
        top = GromacsTopFile(f"{filename}.top", periodicBoxVectors=gro.getPeriodicBoxVectors())

        # Create the simulation
        system = top.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer,
                constraints=HBonds)
        integrator = LangevinMiddleIntegrator(298*kelvin, 1/picosecond, 0.004*picosecond)
        simulation = Simulation(top.topology, system, integrator)
        simulation.context.setPositions(gro.positions)

        # Run the simulation
        simulation.minimizeEnergy()
        simulation.reporters.append(PDBReporter('trajectory.pdb', 10))
        if show_output:
            simulation.reporters.append(StateDataReporter(stdout, 10, step=True,
                    potentialEnergy=True, temperature=True))
        simulation.step(1000)

        # Load the trajectory
        traj = md.load('trajectory.pdb', top=f"{filename}.gro")
        traj.superpose(traj, 0)
        rmsd = md.rmsd(traj, traj, 0)
        
        # Print a table with the RMSD, Avg, Max, Min and StdDev
        print("Runnning a short MD simulation to test the force-field")
        print("RMSD of the molecule during the simulation compared to the original molecule")
        print("RMSD report:")
        print("Avg\tMax\tMin\tStdDev")
        print(f"{np.mean(rmsd):.3f}\t{np.max(rmsd):.3f}\t{np.min(rmsd):.3f}\t{np.std(rmsd):.3f}")

        # Print warnings if the RMSD is too high
        if np.mean(rmsd) > 0.1:
            print("WARNING: RMSD is too high, check the force-field parameters")
        else:
            print("RMSD is under 0.1 nm, the force-field parameters seem to be stable")
        
        # Save the trajectory if requested
        if save_trajectory == False:
            os.remove("trajectory.pdb")
        else:
            print("Trajectory saved as trajectory.pdb")




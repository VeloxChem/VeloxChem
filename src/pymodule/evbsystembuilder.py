import sys
import typing
import copy
import math

import numpy as np
import openmm as mm
import openmm.app as mmapp
import openmm.unit as mmunit

from mpi4py import MPI

from .veloxchemlib import hartree_in_kcalpermol, bohr_in_angstrom
from .veloxchemlib import mpi_master
from .outputstream import OutputStream
from .forcefieldgenerator import ForceFieldGenerator
from .atomtypeidentifier import AtomTypeIdentifier
from .systembuilder import SystemBuilder
from .molecule import Molecule


class EvbSystemBuilder():

    def __init__(self, comm=None, ostream=None):
        '''
        Initialize the EVB driver class.
        '''

        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(sys.stdout)
            else:
                ostream = OutputStream(None)
        
        # output stream
        self.ostream = ostream

        # MPI information
        self.comm = comm
        self.rank = comm.Get_rank()
        self.size = comm.Get_size()

        self.temperature: float = 300
        self.Lambda: list[float]

        self.soft_core_coulomb_ref = True
        self.soft_core_coulomb_run = False
        self.soft_core_lj_ref = True
        self.soft_core_lj_run = False

        self.sc_alpha_lj: float = 0.85
        self.sc_alpha_q: float = 0.3
        self.sc_sigma_q: float = 1.0
        self.sc_power: float = 1 / 6
        self.morse_D_default: float = 500  # kj/mol, default dissociation energy if none is given
        self.morse_couple: float = 1  # kj/mol, scaling for the morse potential to emulate a coupling between two overlapping bonded states
        self.restraint_k: float = 1000  # kj/mol nm^2, force constant for the position restraints
        self.restraint_r_default: float = 0.5  # nm, default position restraint distance if none is given
        self.restraint_r_offset: float = 0.1  # nm, distance added to the measured distance in a structure to set the position restraint distance
        self.coul14_scale: float = 0.833
        self.lj14_scale: float = 0.5

        self.verbose = False
        # self.write_xml = True

        self.constraints: list[dict] = []

        self.k = 4.184 * hartree_in_kcalpermol() * 0.1 * bohr_in_angstrom()

        self.deg_to_rad: float = np.pi / 180

        self.data_folder: str | None = None
        self.run_folder: str | None = None

    def build_systems(
        self,
        reactant: ForceFieldGenerator,
        product: ForceFieldGenerator,
        Lambda: list,
        configuration: dict,
        constraints: list = []
    ):

        self.temperature=configuration.get("temperature", self.temperature)
        NPT=configuration.get("NPT", False)
        pressure=configuration.get("pressure", 1)
        solvent=configuration.get("solvent", None)
        padding=configuration.get("padding", 1.2)
        CNT=configuration.get("CNT", False)
        Graphene=configuration.get("graphene", False)
        M=configuration.get("M", 5)
        N=configuration.get("N", 9)
        ion_count=configuration.get("ion_count", 0)
        no_reactant=configuration.get("no_reactant", False)
        E_field=configuration.get("E_field", [0, 0, 0])

        self.constraints = constraints
        
        system = mm.System()
        topology = mmapp.Topology()
        reaction_chain = topology.addChain()
        reaction_residue = topology.addResidue(name="REA", chain=reaction_chain)
        reaction_atoms = []
        nb_force = mm.NonbondedForce()
        nb_force.setName("General nonbonded force")

        if not no_reactant:
            # add atoms of the solute to the topology and the system
            elements = reactant.molecule.get_labels()
            for i, atom in enumerate(reactant.atoms.values()):
                mm_element = mmapp.Element.getBySymbol(elements[i])
                name = f"{elements[i]}{i}"
                reaction_atom = topology.addAtom(name, mm_element, reaction_residue)
                reaction_atoms.append(reaction_atom)
                system.addParticle(mm_element.mass)
                nb_force.addParticle(0, 1, 0)  #Placeholder values, actual values depend on lambda and will be set later

        if not no_reactant:
            system_mol = Molecule(reactant.molecule)
            positions = list(system_mol.get_coordinates_in_angstrom())  #A
        else:
            system_mol = Molecule()
            positions = []
        X: float = 0
        Y: float = 0

        if CNT or Graphene:
            if CNT and Graphene:
                raise ValueError("CNT and Graphene cannot be used simultaneously, pick one please")
            cc_eq = 0.14080860183951827  # nm
            cc_eq = 0.1418
            x_disp = 3 * cc_eq * 10  # A
            y_disp = math.sqrt(3) * cc_eq * 10  # A
            X = M * x_disp  # A
            Y = N * y_disp  # A

            # system_mol.set_coordinates_in_angstrom(positions)
            R = Y / (2 * np.pi)
            if Graphene:
                self.ostream.print_info(f"Building graphene sheet with X: {X/10:.3f}nm, Y: {Y/10:.3f}nm")
            if CNT:
                self.ostream.print_info(f"Building CNT with X: {X/10:.3f}nm, R: {R/10.:3f}nm")

            graphene_cell = Molecule.read_xyz_file("input_files/Graphene_cell.xyz")

            C_chain = topology.addChain()
            C_residue = topology.addResidue("CCC", chain=reaction_chain)
            CC_atoms = []

            # https://pubs.acs.org/doi/10.1021/jp011344u
            graphene_bond: dict = {
                'type': 'morse',
                'force_constant': 4579.87 * 100,
                'equilibrium': cc_eq,
                'comment': 'Graphene bond, cg-cg GAFF parameters'
            }
            graphene_angle: dict = {
                'type': 'harmonic',
                'force_constant': 504.75776,
                'equilibrium': 120,
                'comment': 'Graphene angle, ce-ce-ce GAFF parameters'
            }

            graphene_dih_0: dict = {
                "type": "Fourier",
                "barrier": 25,
                "phase": 0.0,
                "periodicity": 2,
                "comment": 'Graphene torsion, X -X -c -o GAFF parameters'
            }

            graphene_dih_180: dict = {
                "type": "Fourier",
                "barrier": 25,
                "phase": 180.0,
                "periodicity": 2,
                "comment": 'Graphene torsion, X -X -c -o GAFF parameters'
            }

            sigma = 0.3851
            epsilon = 0.4396

            bonds = {}
            angles = {}
            dihedrals = {}
            impropers = {}
            carbon_atoms = []

            def cell(_m, _n):
                return 4 * (_m % M * N + _n % N)

            index = 0
            for m in range(0, M):
                for n in range(0, N):
                    for i, (element_label, coord) in enumerate(
                            zip(graphene_cell.get_labels(), graphene_cell.get_coordinates_in_angstrom())):

                        #todo make sure the molecule and the added CNT do not clip
                        if CNT:
                            x_offset = X / 5
                            y_offset = X / 5
                            z_offset = X / 5
                        if Graphene:
                            x_offset = X / 2
                            y_offset = Y / 2
                            z_offset = X / 5
                        x = coord[0] + m * x_disp + x_offset
                        y = coord[1] + n * y_disp + y_offset
                        z = coord[2] + 5 + z_offset
                        phi = 2 * np.pi * y / Y
                        new_coord_graphene = np.array([x, y, z])  # A
                        new_coord_CNT = np.array(
                            [x, R * math.sin(phi) + X / 2, R * math.cos(phi) + X / 2]
                        )  # A, X/2 factor is to center the CNT in the middle of the box. The X length is used to get the proper box size
                        if CNT:
                            new_coord = new_coord_CNT
                        else:
                            new_coord = new_coord_graphene

                        system_mol.add_atom(element_label, *(new_coord / bohr_in_angstrom()))  # Bohr
                        positions.append(new_coord)
                        mm_element = mmapp.Element.getBySymbol(element_label)
                        name = index
                        index += 1
                        atom = topology.addAtom(f"{name}", mm_element, C_residue)
                        carbon_atoms.append(atom)
                        CC_atoms.append(atom)
                        system.addParticle(mm_element.mass)
                        #todo placeholder values from gaff ca

                        nb_force.addParticle(0, sigma, epsilon)

                    # Obtain the starting indices for the current unit cell and all adjecent unit cells
                    center = cell(m, n)
                    left = cell(m - 1, n)
                    right = cell(m + 1, n)
                    up = cell(m, n + 1)
                    upright = cell(m + 1, n + 1)
                    down = cell(m, n - 1)

                    # Add bonds and angles
                    bonds.update({(center + 0, center + 1): graphene_bond})
                    bonds.update({(center + 1, center + 2): graphene_bond})
                    bonds.update({(center + 2, center + 3): graphene_bond})
                    bonds.update({(center + 3, right + 0): graphene_bond})
                    bonds.update({(center + 0, up + 1): graphene_bond})
                    bonds.update({(center + 3, up + 2): graphene_bond})

                    angles.update({(center + 0, center + 1, center + 2): graphene_angle})
                    angles.update({(center + 1, center + 2, center + 3): graphene_angle})
                    angles.update({(up + 1, center + 0, center + 1): graphene_angle})
                    angles.update({(center + 2, center + 3, up + 2): graphene_angle})
                    angles.update({(left + 3, center + 0, up + 1): graphene_angle})
                    angles.update({(left + 3, center + 0, center + 1): graphene_angle})
                    angles.update({(center + 0, center + 1, down + 0): graphene_angle})
                    angles.update({(down + 0, center + 1, center + 2): graphene_angle})
                    angles.update({(center + 1, center + 2, down + 3): graphene_angle})
                    angles.update({(center + 3, center + 2, down + 3): graphene_angle})
                    angles.update({(center + 2, center + 3, right + 0): graphene_angle})
                    angles.update({(right + 0, center + 3, up + 2): graphene_angle})

                    dihedrals.update({(left + 3, center + 0, center + 1, down + 0): graphene_dih_0})
                    dihedrals.update({(up + 1, center + 0, center + 1, down + 0): graphene_dih_180})
                    dihedrals.update({(left + 3, center + 0, center + 1, center + 2): graphene_dih_180})
                    dihedrals.update({(left + 3, center + 0, center + 1, center + 2): graphene_dih_0})

                    dihedrals.update({(center + 0, center + 1, center + 2, center + 3): graphene_dih_0})
                    dihedrals.update({(down + 0, center + 1, center + 2, center + 3): graphene_dih_180})
                    dihedrals.update({(center + 0, center + 1, center + 2, down + 3): graphene_dih_180})
                    dihedrals.update({(down + 0, center + 1, center + 2, down + 3): graphene_dih_0})

                    dihedrals.update({(center + 1, center + 2, center + 3, up + 2): graphene_dih_0})
                    dihedrals.update({(down + 3, center + 2, center + 3, up + 2): graphene_dih_180})
                    dihedrals.update({(center + 1, center + 2, center + 3, right + 0): graphene_dih_180})
                    dihedrals.update({(down + 3, center + 2, center + 3, right + 0): graphene_dih_0})

                    dihedrals.update({(center + 2, center + 3, right + 0, right + 1): graphene_dih_0})
                    dihedrals.update({(up + 2, center + 3, right + 0, right + 1): graphene_dih_180})
                    dihedrals.update({(center + 2, center + 3, right + 0, upright + 1): graphene_dih_180})
                    dihedrals.update({(up + 2, center + 3, right + 0, upright + 1): graphene_dih_0})

                    dihedrals.update({(left + 3, center + 0, up + 1, up + 0): graphene_dih_0})
                    dihedrals.update({(center + 1, center + 0, up + 1, up + 0): graphene_dih_180})
                    dihedrals.update({(left + 3, center + 0, up + 1, up + 2): graphene_dih_180})
                    dihedrals.update({(center + 1, center + 0, up + 1, up + 2): graphene_dih_0})

                    dihedrals.update({(center + 2, center + 3, up + 2, up + 1): graphene_dih_0})
                    dihedrals.update({(center + 2, center + 3, up + 2, up + 3): graphene_dih_180})
                    dihedrals.update({(right + 0, center + 3, up + 2, up + 1): graphene_dih_180})
                    dihedrals.update({(right + 0, center + 3, up + 2, up + 3): graphene_dih_0})

                    impropers.update({(left + 3, center + 0, center + 1, up + 2): graphene_dih_180})
                    impropers.update({(center + 0, center + 1, down + 0, center + 2): graphene_dih_180})
                    impropers.update({(center + 1, center + 2, down + 3, center + 3): graphene_dih_180})
                    impropers.update({(center + 2, center + 3, right + 0, up + 2): graphene_dih_180})

            # carbon_harmonic_bond_force = mm.HarmonicBondForce()
            # carbon_harmonic_bond_force.setName("Carbon harmonic bond")
            # carbon_harmonic_bond_force.setUsesPeriodicBoundaryConditions(True)

            carbon_harmonic_bond_force = mm.CustomBondForce("0.5*graphene_fc_factor*k*(r-r0)^2")
            carbon_harmonic_bond_force.setName("Carbon harmonic bond")
            carbon_harmonic_bond_force.addPerBondParameter("r0")
            carbon_harmonic_bond_force.addPerBondParameter("k")
            carbon_harmonic_bond_force.addGlobalParameter("graphene_fc_factor", 1)
            carbon_harmonic_bond_force.setUsesPeriodicBoundaryConditions(True)

            system.addForce(carbon_harmonic_bond_force)
            for key, bond in bonds.items():
                carbon_harmonic_bond_force.addBond(
                    CC_atoms[key[0]].index,
                    CC_atoms[key[1]].index,
                    [bond["equilibrium"], bond["force_constant"]],
                )

            carbon_harmonic_angle_force = mm.HarmonicAngleForce()
            carbon_harmonic_angle_force.setName("Carbon harmonic angle")
            carbon_harmonic_angle_force.setUsesPeriodicBoundaryConditions(True)
            system.addForce(carbon_harmonic_angle_force)
            for key, angle in angles.items():
                carbon_harmonic_angle_force.addAngle(
                    CC_atoms[key[0]].index,
                    CC_atoms[key[1]].index,
                    CC_atoms[key[2]].index,
                    angle["equilibrium"] * np.pi / 180,
                    angle["force_constant"],
                )

            carbon_fourier_dihedral_force = mm.PeriodicTorsionForce()
            carbon_fourier_dihedral_force.setName("Carbon improper fourier torsion")
            carbon_fourier_dihedral_force.setUsesPeriodicBoundaryConditions(True)

            for key, dihedral in dihedrals.items():
                carbon_fourier_dihedral_force.addTorsion(
                    CC_atoms[key[0]].index,
                    CC_atoms[key[1]].index,
                    CC_atoms[key[2]].index,
                    CC_atoms[key[3]].index,
                    dihedral["periodicity"],
                    dihedral["phase"],
                    dihedral["barrier"],
                )
            for key, improper in impropers.items():
                carbon_fourier_dihedral_force.addTorsion(
                    CC_atoms[key[0]].index,
                    CC_atoms[key[1]].index,
                    CC_atoms[key[2]].index,
                    CC_atoms[key[3]].index,
                    improper["periodicity"],
                    improper["phase"],
                    improper["barrier"],
                )

            for i in range(len(carbon_atoms)):
                for j in range(len(carbon_atoms)):
                    if j > i:
                        nb_force.addException(carbon_atoms[i].index, carbon_atoms[j].index, 0, 1, 0)

        if CNT:
            box_x: float = X * 0.1
            box_y: float = X * 0.1
            box_z: float = X * 0.1
        elif Graphene:
            box_x = X * 0.1
            box_y = Y * 0.1
            box_z = X * 0.1
        else:
            box_x = 2 * padding
            box_y = 2 * padding
            box_z = 2 * padding
        box: list[float] = [box_x, box_y, box_z]  # nm
        self.ostream.print_info(f"Building system in box with dimensions {box_x:.3f} x {box_y:.3f} x {box_z:.3f} nm")

        if solvent is None:
            self.positions = np.array(positions) * 0.1
            num_solvent_molecules = 0
            num_solvent_atoms_per_molecule = 0

        else:
            vlxsysbuilder = SystemBuilder()
            # vlxsysbuilder.solvate(reactant.molecule, solvent, padding)

            box_volume: float = box_x * box_y * box_z  # nm^3
            box_A: list[float] = [10 * elem for elem in box]  # A

            na_mol = Molecule.read_smiles("[Na+]")
            cl_mol = Molecule.read_smiles("[Cl-]")

            charge = reactant.molecule.get_charge()
            charge_quantities = [0, 0]
            self.ostream.print_info(f"Charge of the system: {charge}")
            if charge > 0:
                charge_quantities = [ion_count, ion_count + charge]
            elif charge < 0:
                charge_quantities = [ion_count - charge, ion_count]

            solute_volume: float = vlxsysbuilder._get_volume(system_mol) * 0.001  #nm^3
            solvent_volume: float = box_volume - solute_volume
            solvent_mols_per_nm3, solvent_density, solvent_smiles_code = vlxsysbuilder._solvent_properties(solvent)
            solvent_molecule = Molecule.read_smiles(solvent_smiles_code)
            Na_volume = (4.0 / 3.0) * np.pi * 0.332840**3
            Na_fac = Na_volume / solvent_volume
            Cl_volume = (4.0 / 3.0) * np.pi * 0.440104**3
            Cl_fac = Cl_volume / solvent_volume
            self.ostream.print_info(f"Na+ factor: {Na_fac}, Cl- factor: {Cl_fac}")
            quantities = [
                int(solvent_volume * solvent_mols_per_nm3 - charge_quantities[0] * Na_fac -
                    charge_quantities[1] * Cl_fac), charge_quantities[0], charge_quantities[1]
            ]
            quantities = [
                int(solvent_volume * solvent_mols_per_nm3 - charge_quantities[0] * Na_fac -
                    charge_quantities[1] * Cl_fac), charge_quantities[0], charge_quantities[1]
            ]
            solvents: list[Molecule] = [solvent_molecule, na_mol, cl_mol]
            self.ostream.print_info(
                f"Solute volume: {solute_volume:.3f} nm^3, box valume: {box_volume:.3f} nm^3, solvent volume: {solvent_volume:.3f} nm^3"
            )

            self.ostream.print_info(
                f"Adding {quantities[0]} molecules of {solvent} with density {solvent_mols_per_nm3} mol/nm^3 to the system"
            )
            if charge != 0:
                self.ostream.print_info(f"Adding {quantities[1]} Na+ and {quantities[2]} Cl- ions to the system to neutralize the charge")

            vlxsysbuilder.custom_solvate(
                system_mol,
                solvents,
                quantities,
                box_A,
            )
            self.positions = vlxsysbuilder.system_molecule.get_coordinates_in_angstrom() * 0.1

            # box: list[float] = [elem * 0.1 for elem in vlxsysbuilder.box]
            for i, (num_solvent_molecules,
                    vlx_solvent_molecule) in enumerate(zip(vlxsysbuilder.added_solvent_counts, vlxsysbuilder.solvents)):

                mol_index = -1
                if vlx_solvent_molecule == solvent_molecule:
                    mol_index = 0
                elif vlx_solvent_molecule == na_mol:
                    mol_index = 1
                elif vlx_solvent_molecule == cl_mol:
                    mol_index = 2
                else:
                    self.ostream.print_warning("Could not find molecule")
                # num_solvent_molecules = vlxsysbuilder.added_solvent_counts[0]
                # num_solvent_atoms_per_molecule = vlxsysbuilder.solvents[0].number_of_atoms()
                # solvent_molecule = vlxsysbuilder.solvents[0]
                num_solvent_atoms_per_molecule = vlx_solvent_molecule.number_of_atoms()
                elements = vlx_solvent_molecule.get_labels()

                self.ostream.print_info("Generating solvent forcefield")
                solvent_ff = ForceFieldGenerator()
                if mol_index == 0:
                    if solvent == 'spce':
                        #todo maybe don't hardcode this
                        self.ostream.print_info('Using SPCE from amber03 FF')
                        solvent_ff.bonds = {
                            (0, 1): {
                                'type': 'harmonic',
                                'force_constant': 462750.4,
                                'equilibrium': 0.1,
                                'comment': 'SPCE water'
                            },
                            (0, 2): {
                                'type': 'harmonic',
                                'force_constant': 462750.4,
                                'equilibrium': 0.1,
                                'comment': 'SPCE water'
                            }
                        }
                        solvent_ff.angles = {
                            (1, 0, 2): {
                                'type': 'harmonic',
                                'force_constant': 836.8,
                                'equilibrium': 1.91061193216 / self.deg_to_rad,
                                'comment': 'SPCE water'
                            }
                        }
                        solvent_ff.atoms = {
                            0: {
                                'type': 'ow',
                                'name': 'O',
                                'mass': 15.994915,
                                'charge': -0.8476,
                                'sigma': 0.31657195050398818,
                                'epsilon': 0.6497752,
                                'equivalent_atom': 'SPCE water'
                            },
                            1: {
                                'type': 'hw',
                                'name': 'H1',
                                'mass': 1.007825,
                                'charge': 0.4238,
                                'sigma': 1.0,
                                'epsilon': 0.0,
                                'equivalent_atom': 'SPCE water'
                            },
                            2: {
                                'type': 'hw',
                                'name': 'H2',
                                'mass': 1.007825,
                                'charge': 0.4238,
                                'sigma': 1.0,
                                'epsilon': 0.0,
                                'equivalent_atom': 'SPCE water'
                            }
                        }
                        solvent_ff.dihedrals = {}
                        solvent_ff.impropers = {}
                    elif solvent == 'tip3p':
                        self.ostream.print_info('Using TIP3 from amber03 FF')
                        solvent_ff.bonds = {
                            (0, 1): {
                                'type': 'harmonic',
                                'force_constant': 462750.4,
                                'equilibrium': 0.09572,
                                'comment': 'TIP-3P water'
                            },
                            (0, 2): {
                                'type': 'harmonic',
                                'force_constant': 462750.4,
                                'equilibrium': 0.09572,
                                'comment': 'TIP-3P water'
                            }
                        }
                        solvent_ff.angles = {
                            (1, 0, 2): {
                                'type': 'harmonic',
                                'force_constant': 836.8000000000001,
                                'equilibrium': 1.82421813418 / self.deg_to_rad,
                                'comment': 'TIP-3P water'
                            }
                        }
                        solvent_ff.atoms = {
                            0: {
                                'type': 'ow',
                                'name': 'O',
                                'mass': 15.994915,
                                'charge': -0.834,
                                'sigma': 0.31507524065751241,
                                'epsilon': 0.635968,
                                'equivalent_atom': 'TIP-3P water'
                            },
                            1: {
                                'type': 'hw',
                                'name': 'H1',
                                'mass': 1.007825,
                                'charge': 0.417,
                                'sigma': 1.0,
                                'epsilon': 0.0,
                                'equivalent_atom': 'TIP-3P water'
                            },
                            2: {
                                'type': 'hw',
                                'name': 'H2',
                                'mass': 1.007825,
                                'charge': 0.417,
                                'sigma': 1.0,
                                'epsilon': 0.0,
                                'equivalent_atom': 'TIP-3P water'
                            }
                        }
                        solvent_ff.dihedrals = {}
                        solvent_ff.impropers = {}
                    else:
                        solvent_ff = ForceFieldGenerator()
                        solvent_ff.create_topology(vlx_solvent_molecule)
                elif mol_index == 1:
                    # https://github.com/gromacs/gromacs/blob/main/share/top/amber03.ff/ffnonbonded.itp
                    # solvent_ff.partial_charges = [1]
                    solvent_ff.atoms = {
                        0: {
                            'type': 'Na',
                            'name': 'Na+',
                            'mass': 22.99,
                            'charge': 1,
                            'sigma': 0.332840,
                            'epsilon': 0.015897,
                        }
                    }
                    solvent_ff.bonds = {}
                    solvent_ff.angles = {}
                    solvent_ff.dihedrals = {}
                    solvent_ff.impropers = {}
                elif mol_index == 2:
                    # https://github.com/gromacs/gromacs/blob/main/share/top/amber03.ff/ffnonbonded.itp
                    # solvent_ff.partial_charges = [1]
                    solvent_ff.atoms = {
                        0: {
                            'type': 'Cl',
                            'name': 'Cl-',
                            'mass': 35.45,
                            'charge': -1,
                            'sigma': 0.440104,
                            'epsilon': 0.418400,
                        }
                    }
                    solvent_ff.bonds = {}
                    solvent_ff.angles = {}
                    solvent_ff.dihedrals = {}
                    solvent_ff.impropers = {}

                harmonic_bond_force = mm.HarmonicBondForce()
                harmonic_bond_force.setName(f"Solvent {i} harmonic bond")
                harmonic_angle_force = mm.HarmonicAngleForce()
                harmonic_angle_force.setName(f"Solvent {i} harmonic angle")
                fourier_force = mm.PeriodicTorsionForce()
                fourier_force.setName(f"Solvent {i} proper fourier torsion")
                RB_force = mm.RBTorsionForce()
                RB_force.setName(f"Solvent {i} proper RB torsion")
                fourier_imp_force = mm.PeriodicTorsionForce()
                fourier_imp_force.setName(f"Solvent {i} improper fourier torsion")
                system.addForce(harmonic_bond_force)
                system.addForce(harmonic_angle_force)
                system.addForce(fourier_force)
                system.addForce(RB_force)
                system.addForce(fourier_imp_force)

                # Loop over all solvent molecules
                self.ostream.print_info(f"Adding {num_solvent_molecules} molecules of {solvent} to the system")
                solvent_nb_atom_count = 0
                solvent_system_atom_count = 0
                for i in range(num_solvent_molecules):

                    # Create a new residue and add it to the topology
                    resname = "SOL"
                    if mol_index == 0:
                        if solvent == 'spce' or solvent == 'tip3p':
                            resname = "HOH"
                        else:
                            resname = "SOL"
                    else:
                        resname = "ION"
                    solvent_residue = topology.addResidue(name=resname, chain=reaction_chain)
                    solvent_atoms = []
                    # Loop over all atoms in the solvent molecule
                    for j in range(num_solvent_atoms_per_molecule):
                        # Figure out the element of the atom
                        mm_element = mmapp.Element.getBySymbol(elements[j])
                        name = f"{elements[j]}{j}"
                        # add the atom to the topology
                        solvent_atom = topology.addAtom(name, mm_element, solvent_residue)
                        solvent_atoms.append(solvent_atom)
                        # Add the atom as a particle to the system
                        system.addParticle(mm_element.mass)
                        solvent_system_atom_count += 1
                        # Add this atom to the nonbonded force

                    # add all bonded interactions in this molecule
                    #todo can this be merged with _add_reaction_forces?
                    for key, bond in solvent_ff.bonds.items():
                        harmonic_bond_force.addBond(solvent_atoms[key[0]].index, solvent_atoms[key[1]].index,
                                                    bond["equilibrium"], bond["force_constant"])
                    for key, angle in solvent_ff.angles.items():
                        harmonic_angle_force.addAngle(solvent_atoms[key[0]].index, solvent_atoms[key[1]].index,
                                                      solvent_atoms[key[2]].index,
                                                      angle["equilibrium"] * self.deg_to_rad, angle["force_constant"])
                    for key, dihedral in solvent_ff.dihedrals.items():
                        if dihedral["type"] == "RB":
                            RB_force.addTorsion(
                                solvent_atoms[key[0]].index,
                                solvent_atoms[key[1]].index,
                                solvent_atoms[key[2]].index,
                                solvent_atoms[key[3]].index,
                                dihedral["RB_coefficients"][0],
                                dihedral["RB_coefficients"][1],
                                dihedral["RB_coefficients"][2],
                                dihedral["RB_coefficients"][3],
                                dihedral["RB_coefficients"][4],
                                dihedral["RB_coefficients"][5],
                            )
                        elif dihedral["type"] == "Fourier":
                            fourier_force.addTorsion(
                                solvent_atoms[key[0]].index,
                                solvent_atoms[key[1]].index,
                                solvent_atoms[key[2]].index,
                                solvent_atoms[key[3]].index,
                                dihedral["periodicity"],
                                dihedral["phase"] * self.deg_to_rad,
                                dihedral["barrier"],
                            )
                        else:
                            assert False, "Unknown dihedral type"
                    for key, dihedral in solvent_ff.impropers.items():
                        if dihedral["type"] == "Fourier":
                            fourier_imp_force.addTorsion(
                                solvent_atoms[key[0]].index,
                                solvent_atoms[key[1]].index,
                                solvent_atoms[key[2]].index,
                                solvent_atoms[key[3]].index,
                                dihedral["periodicity"],
                                dihedral["phase"] * self.deg_to_rad,
                                dihedral["barrier"],
                            )
                        else:
                            assert False, "Unknown dihedral type"
                    exceptions = self._create_exceptions_from_bonds(solvent_ff)
                    for key, atom in solvent_ff.atoms.items():

                        sigma = atom["sigma"]
                        if sigma == 0:
                            if atom["epsilon"] == 0:
                                sigma = 1
                            else:
                                raise ValueError("Sigma is 0 while epsilon is not, which will cause division by 0")
                        solvent_nb_atom_count += 1
                        nb_force.addParticle(atom["charge"], atom["sigma"], atom["epsilon"])
                    # #Loop over all atoms, and check if their id's are part of any exceptions
                    for i in solvent_ff.atoms.keys():
                        for j in solvent_ff.atoms.keys():
                            if i < j:
                                key = (i, j)
                                if key in exceptions.keys():
                                    epsilon = exceptions[key]["epsilon"]
                                    sigma = exceptions[key]["sigma"]
                                    qq = exceptions[key]["qq"]

                                    nb_force.addException(
                                        solvent_atoms[i].index,
                                        solvent_atoms[j].index,
                                        qq,
                                        sigma,
                                        epsilon,
                                    )
                self.ostream.print_info(
                    f"Added {solvent_nb_atom_count} atoms to the nonbonded force and {solvent_system_atom_count} atoms to the system"
                )

        cmm_remover = mm.CMMotionRemover()
        cmm_remover.setName("CMM remover")
        system.addForce(cmm_remover)

        nb_force.setNonbondedMethod(mm.NonbondedForce.PME)
        # solvent_nb_force.setNonbondedMethod(mm.NonbondedForce.CutoffNonPeriodic)
        cutoff = min(1.0,
                     min(min(box_x, box_y), box_z) *
                     0.4)  # nm, 0.4 factor to accomodate for shrinkage of the box in NPT simulations
        nb_force.setCutoffDistance(cutoff)
        nb_force.setUseSwitchingFunction(True)
        self.ostream.print_info(f"Setting nonbonded cutoff to {cutoff:.3f} nm")
        nb_force.setSwitchingDistance(0.9 * cutoff)
        system.addForce(nb_force)

        if NPT:
            system.addForce(
                mm.MonteCarloFlexibleBarostat(
                    pressure * mmunit.bar,  # type: ignore
                    self.temperature * mmunit.kelvin,  # type: ignore
                ))

        E_field_force = mm.CustomExternalForce("-q*(Ex*x+Ey*y+Ez*z)")
        E_field_force.addGlobalParameter("Ex", E_field[0])
        E_field_force.addGlobalParameter("Ey", E_field[1])
        E_field_force.addGlobalParameter("Ez", E_field[2])
        E_field_force.addPerParticleParameter("q")
        E_field_force.setName("Electric field")
        for i in range(system.getNumParticles()):
            E_field_force.addParticle(i, [0])  # Actual charges are lambda dependent, and are set later
        system.addForce(E_field_force)

        #Add the reactant to the nonbonded force
        if not no_reactant:
            for i, atom in enumerate(reactant.atoms.values()):
                #Make sure the solute does not interact with itself through, as there will be another nonbonded force to take care of this
                for j in range(len(reactant.atoms.values())):
                    if j > i:
                        nb_force.addException(
                            reaction_atoms[i].index,  #todo rewrite loop to reaction_atoms 
                            reaction_atoms[j].index,
                            0.0,
                            1.0,
                            0.0,
                        )

        vector_box: list[mm.Vec3] = [mm.Vec3(box[0], 0, 0), mm.Vec3(0, box[1], 0), mm.Vec3(0, 0, box[2])]
        expanded_box = [[box[0], 0, 0], [0, box[1], 0], [0, 0, box[2]]]
        system.setDefaultPeriodicBoxVectors(*vector_box)
        topology.setPeriodicBoxVectors(expanded_box)

        self.systems: typing.Dict = {}
        for lam in Lambda:
            total_charge = 0
            for i, (reactant_atom, product_atom) in enumerate(zip(reactant.atoms.values(), product.atoms.values())):
                charge = (1 - lam) * reactant_atom["charge"] + lam * product_atom["charge"]
                total_charge += charge
                sigma = (1 - lam) * reactant_atom["sigma"] + lam * product_atom["sigma"]
                epsilon = (1 - lam) * reactant_atom["epsilon"] + lam * product_atom["epsilon"]
                if sigma == 0:
                    if epsilon == 0:
                        sigma = 1
                    else:
                        raise ValueError("Sigma is 0 while epsilon is not, which will cause division by 0")
                nb_force.setParticleParameters(i, charge, sigma, epsilon)

            for i in range(system.getNumParticles()):
                charge = nb_force.getParticleParameters(i)[0]
                E_field_force.setParticleParameters(i, i, [charge])

            if not round(total_charge, 5).is_integer():
                self.ostream.print_warning(f"Warning: total charge for lambda {lam} is {total_charge} and is not a whole number")

            self.systems[lam] = copy.deepcopy(system)
            if lam == 0.0:
                self.systems["reactant"] = copy.deepcopy(system)
            if lam == 1.0:
                self.systems["product"] = copy.deepcopy(system)

        self.reactant = reactant
        self.product = product

        self.reaction_atoms = reaction_atoms  #Used in _add_reaction_forces
        self.topology: mmapp.Topology = topology
        # self.Lambda = Lambda
        # self.temperature = temperature
        self.system_mol = system_mol

        if not no_reactant:

            for lam in Lambda:
                self._add_reaction_forces(self.systems[lam], lam, reference_state=False, lj_soft_core=self.soft_core_lj_run, coul_soft_core=self.soft_core_coulomb_run)
            self._add_reaction_forces(self.systems["reactant"], 0, reference_state=True, lj_soft_core=self.soft_core_lj_ref, coul_soft_core=self.soft_core_coulomb_ref)
            self._add_reaction_forces(self.systems["product"], 1, reference_state=True, lj_soft_core=self.soft_core_lj_ref, coul_soft_core=self.soft_core_coulomb_ref)
        #Give all forces a unique forcegroup so that they can be recalculated separately later on
        for system in self.systems.values():
            for i, force in enumerate(system.getForces()):
                force.setForceGroup(i)

        return self.systems, self.topology, self.positions
        # for i, force in enumerate(self.systems["reactant"].getForces()):
        #     if "CNT" in force.getName():
        #         self.systems["reactant"].removeForce(i)

        # for i, force in enumerate(self.systems["product"].getForces()):
        #     if "CNT" in force.getName():
        #         self.systems["product"].removeForce(i)

    def _add_reaction_forces(self, system, lam, reference_state=False, lj_soft_core=False, coul_soft_core=False) -> mm.System:
        reference_state = reference_state

        forces: list[mm.Force] = []
        forces = forces + self._create_bond_forces(lam, reference_state)
        forces = forces + self._create_angle_forces(lam)
        forces = forces + self._create_proper_torsion_forces(lam)
        forces = forces + self._create_improper_torsion_forces(lam)
        forces = forces + self._create_nonbonded_forces(lam, lj_soft_core, coul_soft_core)
        forces = forces + self._create_constraint_forces(lam, reference_state)

        for i, force in enumerate(forces):
            #     force.setForceGroup(i + 1)  # The 0th force group are the solvent forces
            system.addForce(force)
        return system

    def _create_bond_forces(self, lam, reference_state) -> list[mm.Force]:

        harmonic_force = mm.HarmonicBondForce()
        harmonic_force.setName("Reaction harmonic bond")

        morse_expr = "D*(1-exp(-a*(r-re)))^2;"
        morse_force = mm.CustomBondForce(morse_expr)
        morse_force.setName("Reaction morse bond")
        morse_force.addPerBondParameter("D")
        morse_force.addPerBondParameter("a")
        morse_force.addPerBondParameter("re")

        max_dist_expr = "k*step(r-rmax)*(r-rmax)^2"
        max_distance = mm.CustomBondForce(max_dist_expr)
        max_distance.setName("Reaction distance restraint")
        max_distance.addPerBondParameter("rmax")
        max_distance.addPerBondParameter("k")

        bond_keys = list(set(self.reactant.bonds) | set(self.product.bonds))
        static_bond_keys = list(set(self.reactant.bonds) & set(self.product.bonds))
        broken_bond_keys = list(set(self.reactant.bonds) - set(self.product.bonds))
        formed_bond_keys = list(set(self.product.bonds) - set(self.reactant.bonds))

        for key in bond_keys:
            # mm_top.addBond(mm_atoms[key[0]], mm_atoms[key[1]])
            # if in a and b: interpolate harmonic

            if key in static_bond_keys:

                bondA = self.reactant.bonds[key]
                bondB = self.product.bonds[key]

                if round(bondA["equilibrium"], 3) == round(bondB["equilibrium"], 3):
                    harmonic_force.addBond(
                        self._reaction_to_total_atomid(key[0]),
                        self._reaction_to_total_atomid(key[1]),
                        bondA["equilibrium"],
                        (1 - lam) * bondA["force_constant"] + lam * bondB["force_constant"],
                    )
                else:
                    harmonic_force.addBond(
                        self._reaction_to_total_atomid(key[0]),
                        self._reaction_to_total_atomid(key[1]),
                        bondA["equilibrium"],
                        (1 - lam) * bondA["force_constant"],
                    )
                    harmonic_force.addBond(
                        self._reaction_to_total_atomid(key[0]),
                        self._reaction_to_total_atomid(key[1]),
                        bondB["equilibrium"],
                        lam * bondB["force_constant"],
                    )
            # if in a or b:
            else:
                if key in broken_bond_keys:
                    scale = 1 - lam
                    bond = self.reactant.bonds[key]
                    if hasattr(self.product, "molecule"):
                        coords = self.product.molecule.get_coordinates_in_angstrom()
                        r = (AtomTypeIdentifier.measure_length(coords[key[0]], coords[key[1]]) * 0.1)
                    else:
                        r = self.restraint_r_default
                        if self.verbose:
                            self.ostream.print_info(f"No product geometry given, defaulting position restraint to {r} nm")
                elif key in formed_bond_keys:
                    scale = lam
                    bond = self.product.bonds[key]
                    coords = self.reactant.molecule.get_coordinates_in_angstrom()
                    r = (AtomTypeIdentifier.measure_length(coords[key[0]], coords[key[1]]) * 0.1)
                else:

                    assert (False), "A bond can either be static, or dynamic, in  which case it can be broken or formed"
                # if scale > 0:
                #     # scale morse
                if "D" not in bond.keys():
                    D = self.morse_D_default
                    if self.verbose:
                        self.ostream.print_info(
                            f"No D value associated with bond {key[0]} {key[1]}. Setting to default value {self.morse_D_default}"
                        )
                else:
                    D = bond["D"]
                a = math.sqrt(bond["force_constant"] / (2 * D))
                re = bond["equilibrium"]
                morse_force.addBond(
                    self._reaction_to_total_atomid(key[0]),
                    self._reaction_to_total_atomid(key[1]),
                    [scale * D, a, re],
                )
                    
                # if scale > 0:
                #     # scale morse
                if "D" not in bond.keys():
                    D = self.morse_D_default
                    if self.verbose:
                        self.ostream.print_info(
                            f"No D value associated with bond {key[0]} {key[1]}. Setting to default value {self.morse_D_default}"
                        )
                else:
                    D = bond["D"]
                a = math.sqrt(bond["force_constant"] / (2 * D))
                re = bond["equilibrium"]
                morse_force.addBond(
                    self._reaction_to_total_atomid(key[0]),
                    self._reaction_to_total_atomid(key[1]),
                    [scale * D, a, re],
                )
                    

                if not reference_state:
                    k = self.restraint_k
                    rmax = r + self.restraint_r_offset
                    max_distance.addBond(
                        self._reaction_to_total_atomid(key[0]),
                        self._reaction_to_total_atomid(key[1]),
                        [rmax, (1 - scale) * k],
                    )
                    if self.verbose:
                        self.ostream.print_info(
                            f"Adding maximum distance {rmax} with k {(1-scale)*k} to atoms {key[0]} and {key[1]} of for lambda {lam}"
                        )
        return [harmonic_force, morse_force, max_distance]

    def _create_angle_forces(self, lam: float) -> typing.List[mm.Force]:
        harmonic_force = mm.HarmonicAngleForce()
        harmonic_force.setName("Reaction harmonic angle")
        angle_keys = list(set(self.reactant.angles) | set(self.product.angles))
        for key in angle_keys:
            if key in self.reactant.angles.keys() and key in self.product.angles.keys():
                angleA = self.reactant.angles[key]
                angleB = self.product.angles[key]
                if round(angleA["equilibrium"], 3) == round(angleB["equilibrium"], 3):
                    harmonic_force.addAngle(
                        self._reaction_to_total_atomid(key[0]),
                        self._reaction_to_total_atomid(key[1]),
                        self._reaction_to_total_atomid(key[2]),
                        angleA["equilibrium"] * self.deg_to_rad,
                        (1 - lam) * angleA["force_constant"] + lam * angleB["force_constant"],
                    )
                else:
                    harmonic_force.addAngle(
                        self._reaction_to_total_atomid(key[0]),
                        self._reaction_to_total_atomid(key[1]),
                        self._reaction_to_total_atomid(key[2]),
                        angleA["equilibrium"] * self.deg_to_rad,
                        (1 - lam) * angleA["force_constant"],
                    )
                    harmonic_force.addAngle(
                        self._reaction_to_total_atomid(key[0]),
                        self._reaction_to_total_atomid(key[1]),
                        self._reaction_to_total_atomid(key[2]),
                        angleB["equilibrium"] * self.deg_to_rad,
                        lam * angleB["force_constant"],
                    )
            else:
                if key in self.reactant.angles.keys():
                    scale = 1 - lam
                    angle = self.reactant.angles[key]
                else:
                    scale = lam
                    angle = self.product.angles[key]
                if scale > 0:
                    harmonic_force.addAngle(
                        self._reaction_to_total_atomid(key[0]),
                        self._reaction_to_total_atomid(key[1]),
                        self._reaction_to_total_atomid(key[2]),
                        angle["equilibrium"] * self.deg_to_rad,
                        scale * angle["force_constant"],
                    )
        return [harmonic_force]

    def _create_proper_torsion_forces(self, lam) -> typing.List[mm.Force]:
        fourier_force = mm.PeriodicTorsionForce()
        fourier_force.setName("Reaction proper fourier torsion")
        RB_force = mm.RBTorsionForce()
        RB_force.setName("Reaction proper RB torsion")

        dihedral_keys = list(set(self.reactant.dihedrals) | set(self.product.dihedrals))
        for key in dihedral_keys:
            total_atom_id = [self._reaction_to_total_atomid(reaction_atomid) for reaction_atomid in key]
            if (key in self.reactant.dihedrals.keys() and key in self.product.dihedrals.keys()):
                dihedA = self.reactant.dihedrals[key]
                dihedB = self.product.dihedrals[key]
                self._add_torsion(fourier_force, dihedA, total_atom_id, 1 - lam)
                self._add_torsion(fourier_force, dihedB, total_atom_id, lam)
            else:
                if key in self.reactant.dihedrals.keys():
                    scale = 1 - lam
                    dihed = self.reactant.dihedrals[key]
                else:
                    scale = lam
                    dihed = self.product.dihedrals[key]
                if scale > 0:
                    self._add_torsion(fourier_force, dihed, total_atom_id, scale)
        return [fourier_force, RB_force]

    def _create_improper_torsion_forces(self, lam) -> typing.List[mm.Force]:

        fourier_force = mm.PeriodicTorsionForce()
        fourier_force.setName("Reaction improper fourier torsion")

        dihedral_keys = list(set(self.reactant.impropers) | set(self.product.impropers))
        for key in dihedral_keys:
            total_atom_id = [self._reaction_to_total_atomid(reaction_atomid) for reaction_atomid in key]
            if (key in self.reactant.impropers.keys() and key in self.product.impropers.keys()):
                dihedA = self.reactant.impropers[key]
                self._add_torsion(fourier_force, dihedA, total_atom_id, 1 - lam)
                dihedB = self.product.impropers[key]
                self._add_torsion(fourier_force, dihedB, total_atom_id, lam)
            else:
                if key in self.reactant.impropers.keys():
                    scale = 1 - lam
                    dihed = self.reactant.impropers[key]
                else:
                    scale = lam
                    dihed = self.product.impropers[key]
                if scale > 0:
                    self._add_torsion(fourier_force, dihed, total_atom_id, scale)
        return [fourier_force]

    def _add_torsion(self, fourier_force, torsion_dict, total_atom_id, barrier_scaling):
        if torsion_dict["type"] == "Fourier":
            if torsion_dict.get("multiple", False):
                for periodicity, phase, barrier in zip(torsion_dict["periodicity"], torsion_dict["phase"], torsion_dict["barrier"]):
                    fourier_force.addTorsion(
                        total_atom_id[0],
                        total_atom_id[1],
                        total_atom_id[2],
                        total_atom_id[3],
                        abs(periodicity),
                        phase * self.deg_to_rad,
                        barrier_scaling * barrier,
                    )
            else:
                fourier_force.addTorsion(
                    total_atom_id[0],
                    total_atom_id[1],
                    total_atom_id[2],
                    total_atom_id[3],
                    torsion_dict["periodicity"],
                    torsion_dict["phase"] * self.deg_to_rad,
                    barrier_scaling * torsion_dict["barrier"],
                )
        else:
            assert False, "Unknown dihedral type"

    def _get_couloumb_expression(self, soft_core: bool):
        soft_core_expression = (
            " (1-l) * CoultotA + l  * CoultotB; "
            "CoultotA   = step(r-rqA) * CoulA + step(rqA-r) * CoullinA;"
            "CoultotB   = step(r-rqB) * CoulB + step(rqB-r) * CoullinB;"
            "CoullinA   = k * ( ( qqA / rqAdiv^3 ) * r^2 - 3 * ( qqA / rqAdiv^2 ) * r + 3 * ( qqA / rqAdiv ) );"
            "CoullinB   = k * ( ( qqB / rqBdiv^3 ) * r^2 - 3 * ( qqB / rqBdiv^2 ) * r + 3 * ( qqB / rqBdiv ) );"
            "rqAdiv     = select(rqA,rqA,1);"
            "rqBdiv     = select(rqB,rqB,1);"
            "rqA        = alphaq * ( 1 + sigmaq * qqA )  * 1 ^ pow;"
            "rqB        = alphaq * ( 1 + sigmaq * qqB )  * 1 ^ pow;"
            "CoulA      = k*qqA/r; "
            "CoulB      = k*qqB/r; "
            ""
            f"k         = {self.k};"
            ""
            f"alphalj   = {self.sc_alpha_lj};"
            f"alphaq    = {self.sc_alpha_q};"
            f"sigmaq    = {self.sc_sigma_q};"
            f"pow       = {self.sc_power};")

        hard_core_expression = (
            " (1-l) * k*qqA/r + l * k*qqB/r; "
            ""
            f"k         = {self.k};")
            
        if soft_core:
            return soft_core_expression
        else:
            return hard_core_expression

    def _get_lj_expression(self, soft_core: bool):
        soft_core_expression = (
            " (1-l) * LjtotA + l * LjtotB; "
            ""
            "LjtotA     = (step(r - rljA) * LjA + step(rljA - r) * LjlinA);"
            "LjtotB     = (step(r - rljB) * LjB + step(rljB - r) * LjlinB);"
            "LjlinA     = ( (78*A12) / (rljAdiv^14) - (21*A6) / (rljAdiv^8) )*r^2 - ( (168*A12) / (rljAdiv^13) - (48*A6) / (rljAdiv^7) )*r + ( (91*A12) / (rljAdiv^12) - (28*A6) / (rljAdiv^6) );"
            "LjlinB     = ( (78*B12) / (rljBdiv^14) - (21*B6) / (rljBdiv^8) )*r^2 - ( (168*B12) / (rljBdiv^13) - (48*B6) / (rljBdiv^7) )*r + ( (91*B12) / (rljBdiv^12) - (28*B6) / (rljBdiv^6) );"
            # if rljA = 0, returns 1, otherwise returns rljA. Prevents division by 0 while the step factor is already 0
            "rljAdiv    = select(rljA,rljA,1);"
            "rljBdiv    = select(rljB,rljB,1);"
            "rljA       = alphalj * ( (26/7 ) * A6  * 1 ) ^ pow;"
            "rljB       = alphalj * ( (26/7 ) * B6  * 1 ) ^ pow;"
            "LjA        = A12 / r^12 - A6 / r^6;"
            "LjB        = B12 / r^12 - B6 / r^6;"
            "A12        = 4 * epsilonA * sigmaA ^ 12; "
            "B12        = 4 * epsilonB * sigmaB ^ 12; "
            "A6         = 4 * epsilonA * sigmaA ^ 6; "
            "B6         = 4 * epsilonB * sigmaB ^ 6; "
            ""
            f"alphalj   = {self.sc_alpha_lj};"
            f"alphaq    = {self.sc_alpha_q};"
            f"sigmaq    = {self.sc_sigma_q};"
            f"pow       = {self.sc_power};")

        hard_core_expression = (
            " (1-l) * LjtotA "
            "  + l  * LjtotB; "
            ""
            "LjtotA     = A12 / r^12 - A6 / r^6;"
            "LjtotB     = B12 / r^12 - B6 / r^6;"
            "A12        = 4 * epsilonA * sigmaA ^ 12; "
            "B12        = 4 * epsilonB * sigmaB ^ 12; "
            "A6         = 4 * epsilonA * sigmaA ^ 6; "
            "B6         = 4 * epsilonB * sigmaB ^ 6; ")

        if soft_core:
            return soft_core_expression
        else:
            return hard_core_expression

    def _create_nonbonded_forces(self, lam, lj_soft_core=False, coul_soft_core=False) -> typing.List[mm.Force]:
        coulomb_force = mm.CustomBondForce(self._get_couloumb_expression(coul_soft_core))
        if coul_soft_core:
            coulomb_force.setName("Reaction internal coulomb")
        else:
            coulomb_force.setName("Reaction internal coulomb")
        coulomb_force.addPerBondParameter("qqA")
        coulomb_force.addPerBondParameter("qqB")
        coulomb_force.addGlobalParameter("l", lam)

        lj_force = mm.CustomBondForce(self._get_lj_expression(lj_soft_core))
        if lj_soft_core:
            lj_force.setName("Reaction internal lennard-jones soft-core")
        else:
            lj_force.setName("Reaction internal lennard-jones")
        lj_force.addPerBondParameter("sigmaA")
        lj_force.addPerBondParameter("sigmaB")
        lj_force.addPerBondParameter("epsilonA")
        lj_force.addPerBondParameter("epsilonB")
        lj_force.addGlobalParameter("l", lam)

        

        reactant_exceptions = self._create_exceptions_from_bonds(self.reactant)
        product_exceptions = self._create_exceptions_from_bonds(self.product)

        #Loop over all atoms, and check if their id's are part of any exceptions
        for i in self.reactant.atoms.keys():
            for j in self.reactant.atoms.keys():
                if i < j:
                    key = (i, j)
                    # Remove any exception from the nonbondedforce
                    # and add it instead to the exception bond force
                    if key in reactant_exceptions.keys() :
                        epsilonA = reactant_exceptions[key]["epsilon"]
                        sigmaA = reactant_exceptions[key]["sigma"]
                        qqA = reactant_exceptions[key]["qq"]
                    
                    else:
                        atomA1 = self.reactant.atoms[key[0]]
                        atomA2 = self.reactant.atoms[key[1]]
                        epsilonA = math.sqrt(atomA1["epsilon"] * atomA2["epsilon"])
                        sigmaA = 0.5 * (atomA1["sigma"] + atomA2["sigma"])
                        qqA = atomA1["charge"] * atomA2["charge"]

                    if key in product_exceptions.keys():
                        epsilonB = product_exceptions[key]["epsilon"]
                        sigmaB = product_exceptions[key]["sigma"]
                        qqB = product_exceptions[key]["qq"]
                    else:
                        atomB1 = self.product.atoms[key[0]]
                        atomB2 = self.product.atoms[key[1]]
                        epsilonB = math.sqrt(atomB1["epsilon"] * atomB2["epsilon"])
                        sigmaB = 0.5 * (atomB1["sigma"] + atomB2["sigma"])
                        qqB = atomB1["charge"] * atomB2["charge"]

                    if sigmaA == 1.0:
                        sigmaA = sigmaB
                    elif sigmaB == 1.0:
                        sigmaB = sigmaA
                    if not (qqA == 0.0 and qqB == 0.0):
                        coulomb_force.addBond(
                            self._reaction_to_total_atomid(key[0]),
                            self._reaction_to_total_atomid(key[1]),
                            [qqA, qqB],
                        )
                    if not (epsilonA == 0.0 and epsilonB == 0.0):
                        lj_force.addBond(
                            self._reaction_to_total_atomid(key[0]),
                            self._reaction_to_total_atomid(key[1]),
                            [sigmaA, sigmaB, epsilonA, epsilonB],
                        )

        return [lj_force,coulomb_force]

    def _create_exceptions_from_bonds(self, molecule: ForceFieldGenerator) -> dict[tuple[int, int], dict[str, float]]:
        particles = molecule.atoms

        exclusions = [set() for _ in range(len(particles))]
        bonded12 = [set() for _ in range(len(particles))]
        exceptions = {}

        # Populate bonded12 with bonds.
        for bond in molecule.bonds.keys():
            bonded12[bond[0]].add(bond[1])
            bonded12[bond[1]].add(bond[0])

        # Find particles separated by 1, 2, or 3 bonds.
        for i in range(len(particles)):
            self._add_exclusions_to_set(bonded12, exclusions[i], i, i, 2)

        # Find particles separated by 1 or 2 bonds and create exceptions.
        for i in range(len(exclusions)):
            bonded13 = set()
            self._add_exclusions_to_set(bonded12, bonded13, i, i, 1)
            for j in exclusions[i]:
                if j < i:
                    if j not in bonded13:
                        # This is a 1-4 interaction.
                        particle1 = particles[j]
                        particle2 = particles[i]
                        charge_prod = (self.coul14_scale * particle1["charge"] * particle2["charge"])
                        sigma = 0.5 * (particle1["sigma"] + particle2["sigma"])
                        epsilon = (self.lj14_scale * (particle1["epsilon"] * particle2["epsilon"])**0.5)
                        exceptions[tuple(sorted((i, j)))] = {
                            "qq": charge_prod,
                            "sigma": sigma,
                            "epsilon": epsilon,
                        }
                    else:
                        # This interaction should be completely excluded.
                        exceptions[tuple(sorted((i, j)))] = {
                            "qq": 0.0,
                            "sigma": 1.0,
                            "epsilon": 0.0,
                        }
        return exceptions

    @staticmethod
    def _add_exclusions_to_set(bonded12, exclusions, base_particle, from_particle, current_level):
        for i in bonded12[from_particle]:
            if i != base_particle:
                exclusions.add(i)
            if current_level > 0:
                EvbSystemBuilder._add_exclusions_to_set(bonded12, exclusions, base_particle, i, current_level - 1)

    def _create_constraint_forces(self, lam, reference_state: bool) -> typing.List[mm.Force]:
        bond_constraint = mm.HarmonicBondForce()
        bond_constraint.setName("Bond constraint")
        angle_constraint = mm.HarmonicAngleForce()
        angle_constraint.setName("Angle constraint")
        torsion_constraint = mm.CustomTorsionForce("0.5*k*(theta-theta0)^2")
        torsion_constraint.setName("Harmonic torsion constraint")
        torsion_constraint.addPerTorsionParameter("theta0")
        torsion_constraint.addPerTorsionParameter("k")
        if len(self.constraints) == 0:
            self.ostream.print_info(f"No constraints found")
        else:
            self.ostream.print_info(f"Adding constraints: {self.constraints}")
        if (
                not reference_state
        ):  # Return the reference state forces empty so that we have the same number of forces in the reference state and run state
            for constraint in self.constraints:
                key = list(constraint.keys())[0]
                if "lambda_couple" in constraint[key].keys():
                    if constraint[key]["lambda_couple"] == 1:
                        scale = lam
                    elif constraint[key]["lambda_couple"] == -1:
                        scale = 1 - lam
                    else:
                        scale = 1
                else:
                    scale = 1

                if len(key) == 2:
                    bond_constraint.addBond(
                        self._reaction_to_total_atomid(key[0]),
                        self._reaction_to_total_atomid(key[1]),
                        constraint[key]["equilibrium"],
                        constraint[key]["force_constant"] * scale,
                    )
                if len(key) == 3:
                    angle_constraint.addAngle(
                        self._reaction_to_total_atomid(key[0]),
                        self._reaction_to_total_atomid(key[1]),
                        self._reaction_to_total_atomid(key[2]),
                        constraint[key]["equilibrium"] * self.deg_to_rad,
                        constraint[key]["force_constant"] * scale,
                    )
                if len(key) == 4:
                    torsion_constraint.addTorsion(
                        self._reaction_to_total_atomid(key[0]),
                        self._reaction_to_total_atomid(key[1]),
                        self._reaction_to_total_atomid(key[2]),
                        self._reaction_to_total_atomid(key[3]),
                        [
                            constraint[key]["equilibrium"] * self.deg_to_rad,
                            constraint[key]["force_constant"] * scale,
                        ],
                    )
        return [bond_constraint, angle_constraint, torsion_constraint]

    def _reaction_to_total_atomid(self, reaction_atom_id: int) -> int:
        return self.reaction_atoms[reaction_atom_id].index

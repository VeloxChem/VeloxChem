#
#                                   VELOXCHEM
#              ----------------------------------------------------
#                          An Electronic Structure Code
#
#  SPDX-License-Identifier: BSD-3-Clause
#
#  Copyright 2018-2025 VeloxChem developers
#
#  Redistribution and use in source and binary forms, with or without modification,
#  are permitted provided that the following conditions are met:
#
#  1. Redistributions of source code must retain the above copyright notice, this
#     list of conditions and the following disclaimer.
#  2. Redistributions in binary form must reproduce the above copyright notice,
#     this list of conditions and the following disclaimer in the documentation
#     and/or other materials provided with the distribution.
#  3. Neither the name of the copyright holder nor the names of its contributors
#     may be used to endorse or promote products derived from this software without
#     specific prior written permission.
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
#  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
#  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
#  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
#  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
#  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
#  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
#  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
#  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

from mpi4py import MPI
import numpy as np
import typing
import copy
import math
import sys
from enum import Enum, auto

from .veloxchemlib import hartree_in_kcalpermol, bohr_in_angstrom, Point
from .veloxchemlib import mpi_master
from .outputstream import OutputStream
from .mmforcefieldgenerator import MMForceFieldGenerator
from .atomtypeidentifier import AtomTypeIdentifier
from .solvationbuilder import SolvationBuilder
from .molecule import Molecule
from .errorhandler import assert_msg_critical

try:
    import openmm as mm
    import openmm.app as mmapp
    import openmm.unit as mmunit
except ImportError:
    pass


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

        self.sc_alpha_lj: float = 0.7
        self.sc_alpha_q: float = 0.3
        self.sc_sigma_q: float = 1.0
        self.sc_power: float = 1 / 6
        self.morse_D_default: float = 10000  # kj/mol, default dissociation energy if none is given
        self.morse_couple: float = 1  # kj/mol, scaling for the morse potential to emulate a coupling between two overlapping bonded states
        self.restraint_k: float = 1000  # kj/mol nm^2, force constant for the position restraints
        self.restraint_r_default: float = 0.5  # nm, default position restraint distance if none is given
        self.restraint_r_offset: float = 0.1  # nm, distance added to the measured distance in a structure to set the position restraint distance
        self.coul14_scale: float = 0.833
        self.lj14_scale: float = 0.5
        self.minimal_nb_cutoff: float = 1  # nm, minimal cutoff for the nonbonded force

        self.soft_core_coulomb_pes = True
        self.soft_core_lj_pes = True

        self.soft_core_coulomb_int = False
        self.soft_core_lj_int = False
        self.pressure: float = -1.
        self.solvent: str = None  #type: ignore
        self.padding: float = 1.
        self.no_reactant: bool = False
        self.E_field: list[float] = [0, 0, 0]
        self.neutralize: bool = False

        self.bonded_integration: bool = True  # If the integration potential should use bonded (harmonic/morse) forces for forming/breaking bonds, instead of replacing them with nonbonded potentials
        self.bonded_integration_fac: float = 0.25  # Scaling factor for the bonded integration forces.

        self.verbose = False

        self.constraints: list[dict] = []

        self.k = 4.184 * hartree_in_kcalpermol() * 0.1 * bohr_in_angstrom(
        )  # Coulombic pre-factor

        self.deg_to_rad: float = np.pi / 180

        self.data_folder: str | None = None
        self.run_folder: str | None = None

        self.keywords = {
            "temperature": {
                "type": float
            },
            "minimal_nb_cutoff": {
                "type": float
            },
            "bonded_integration": {
                "type": bool
            },
            "soft_core_coulomb_pes": {
                "type": bool
            },
            "soft_core_lj_pes": {
                "type": bool
            },
            "soft_core_coulomb_int": {
                "type": bool
            },
            "soft_core_lj_int": {
                "type": bool
            },
            "bonded_integration_fac": {
                "type": float
            },
            "pressure": {
                "type": float
            },
            "solvent": {
                "type": str
            },
            "padding": {
                "type": float
            },
            "no_reactant": {
                "type": bool
            },
            "E_field": {
                "type": list
            },
            "neutralize": {
                "type": bool
            },
            "morse_D_default": {
                "type": float
            },
            "morse_couple": {
                "type": float
            },
            "restraint_k": {
                "type": float
            },
            "restraint_r_default": {
                "type": float
            },
            "restraint_r_offset": {
                "type": float
            },
            "coul14_scale": {
                "type": float
            },
            "lj14_scale": {
                "type": float
            },
        }

    def build_systems(
        self,
        reactant: MMForceFieldGenerator,
        product: MMForceFieldGenerator,
        Lambda: list,
        configuration: dict,
        constraints: list = [],
    ):

        assert_msg_critical('openmm' in sys.modules,
                            'openmm is required for EvbSystemBuilder.')

        for keyword, value in self.keywords.items():
            if keyword in configuration:
                if (not isinstance(configuration[keyword], value["type"])
                        and not (isinstance(configuration[keyword], int)
                                 and value["type"] == float)):
                    raise ValueError(
                        f"Configuration option {keyword} should be of type {value['type']}"
                    )
                else:
                    setattr(self, keyword, configuration[keyword])
                    self.ostream.print_info(
                        f"{keyword}: {getattr(self, keyword)}")

            else:
                self.ostream.print_info(
                    f"{keyword}: {getattr(self, keyword)} (default)")

        # CNT = configuration.get("CNT", False)
        # CNT = False  # todo fix the exploding CNT
        # Graphene = configuration.get("graphene", False)
        # graphene_size = configuration.get("graphene_size", 2)
        # CNT_radius = configuration.get("CNT_radius", 0.5)
        # ion_count = configuration.get("ion_count", 0)

        self.constraints = constraints

        system = mm.System()
        topology = mmapp.Topology()
        reaction_chain = topology.addChain()
        reaction_residue = topology.addResidue(name="REA", chain=reaction_chain)
        reaction_atoms = []
        nb_force = mm.NonbondedForce()
        nb_force.setName("General nonbonded force")

        if not self.no_reactant:
            # add atoms of the solute to the topology and the system
            elements = reactant.molecule.get_labels()
            for i, atom in enumerate(reactant.atoms.values()):
                mm_element = mmapp.Element.getBySymbol(elements[i])
                name = f"{elements[i]}{i}"
                reaction_atom = topology.addAtom(name, mm_element,
                                                 reaction_residue)
                reaction_atoms.append(reaction_atom)
                system.addParticle(mm_element.mass)
                nb_force.addParticle(
                    0, 1, 0
                )  #Placeholder values, actual values depend on lambda and will be set later

        if not self.no_reactant:
            system_mol = Molecule(reactant.molecule)
            positions = system_mol.get_coordinates_in_angstrom()
            x_size = 0.1 * (max(positions[:, 0]) - min(positions[:, 0]))
            y_size = 0.1 * (max(positions[:, 1]) - min(positions[:, 1]))
            z_size = 0.1 * (max(positions[:, 2]) - min(positions[:, 2]))
            box = [
                2 * self.padding + x_size, 2 * self.padding + y_size,
                2 * self.padding + z_size
            ]
            self.ostream.print_info(
                f"Size of the molecule: {x_size:.3f} x {y_size:.3f} x {z_size:.3f}, padding: {self.padding:.3f} nm."
            )
        else:
            box = [1, 1, 1]
            system_mol = Molecule()
            positions = np.array()

        # if CNT or Graphene:
        #     assert False, "Rethink CNT/graphene input"
        #     box = self._add_CNT_graphene(system, CNT, Graphene, nb_force, topology, system_mol, positions, box,
        #                                  graphene_size, CNT_radius)

        # self.ostream.print_info(f"Building system in box with dimensions {box[0]:.3f} x {box[1]:.3f} x {box[2]:.3f} nm")

        if self.solvent:
            box = self._add_solvent(system, system_mol, self.solvent, topology,
                                    nb_force, self.neutralize, self.padding)

        else:
            self.positions = np.array(positions) * 0.1

        cmm_remover = mm.CMMotionRemover()
        cmm_remover.setName("CMM remover")
        cmm_remover.setForceGroup(EvbForceGroup.CMM_REMOVER.value)
        system.addForce(cmm_remover)

        nb_force.setNonbondedMethod(mm.NonbondedForce.PME)

        cutoff = min(
            self.minimal_nb_cutoff,
            min(min(box[0], box[1]), box[2]) * 0.4
        )  # nm, 0.4 factor to accomodate for shrinkage of the box in NPT simulations
        nb_force.setCutoffDistance(cutoff)
        nb_force.setUseSwitchingFunction(True)
        self.ostream.print_info(f"Setting nonbonded cutoff to {cutoff:.3f} nm")
        nb_force.setSwitchingDistance(0.9 * cutoff)
        nb_force.setForceGroup(EvbForceGroup.NB_FORCE.value)
        system.addForce(nb_force)

        if self.pressure > 0:
            barostat = mm.MonteCarloBarostat(
                self.pressure * mmunit.bar,  # type: ignore
                self.temperature * mmunit.kelvin,  # type: ignore
            )
            barostat.setForceGroup(EvbForceGroup.BAROSTAT.value)
            system.addForce(barostat)

        if np.any(np.array(self.E_field) > 0.001):
            E_field_force = self._create_E_field(system, self.E_field)
            E_field_force.setForceGroup(EvbForceGroup.E_FIELD.value)
            system.addForce(E_field_force)

        #Add the reactant to the nonbonded force
        if not self.no_reactant:
            for i, atom in enumerate(reactant.atoms.values()):
                #Make sure the solute does not interact with itself through, as there will be another nonbonded force to take care of this
                for j in range(len(reactant.atoms.values())):
                    if j > i:
                        nb_force.addException(
                            reaction_atoms[i].index,
                            reaction_atoms[j].index,
                            0.0,
                            1.0,
                            0.0,
                        )

        vector_box: list[mm.Vec3] = [
            mm.Vec3(box[0], 0, 0),
            mm.Vec3(0, box[1], 0),
            mm.Vec3(0, 0, box[2])
        ]
        expanded_box = [[box[0], 0, 0], [0, box[1], 0], [0, 0, box[2]]]
        system.setDefaultPeriodicBoxVectors(*vector_box)
        topology.setPeriodicBoxVectors(expanded_box)

        self.systems: typing.Dict = {}
        for lam in Lambda:
            total_charge = 0
            for i, (reactant_atom, product_atom) in enumerate(
                    zip(reactant.atoms.values(), product.atoms.values())):
                charge = (1 - lam) * reactant_atom[
                    "charge"] + lam * product_atom["charge"]
                total_charge += charge
                sigma = (
                    1 -
                    lam) * reactant_atom["sigma"] + lam * product_atom["sigma"]
                epsilon = (1 - lam) * reactant_atom[
                    "epsilon"] + lam * product_atom["epsilon"]
                if sigma == 0:
                    if epsilon == 0:
                        sigma = 1
                    else:
                        raise ValueError(
                            "Sigma is 0 while epsilon is not, which will cause division by 0"
                        )
                nb_force.setParticleParameters(i, charge, sigma, epsilon)

            for i in range(system.getNumParticles()):
                charge = nb_force.getParticleParameters(i)[0]
                if np.any(np.array(self.E_field) > 0.001):
                    E_field_force.setParticleParameters(i, i, [charge])

            if not round(total_charge, 5).is_integer():
                self.ostream.print_warning(
                    f"Warning: total charge for lambda {lam} is {total_charge} and is not a whole number"
                )

            self.systems[lam] = copy.deepcopy(system)
            # if lam == 0.0:
            #     self.systems["reactant"] = copy.deepcopy(system)
            # if lam == 1.0:
            #     self.systems["product"] = copy.deepcopy(system)

        self.reactant = reactant
        self.product = product

        self.reaction_atoms = reaction_atoms  #Used in _add_reaction_forces
        self.topology: mmapp.Topology = topology
        self.system_mol = system_mol

        if not self.no_reactant:
            for lam in Lambda:
                self._add_reaction_forces(self.systems[lam], lam)

        self.ostream.flush()
        return self.systems, self.topology, self.positions

    def _add_CNT_graphene(self, system, CNT: bool, Graphene: bool, nb_force,
                          topology, system_mol, positions, box, graphene_size,
                          CNT_radius):

        assert_msg_critical('openmm' in sys.modules,
                            'openmm is required for EvbSystemBuilder.')

        if CNT and Graphene:
            raise ValueError(
                "CNT and Graphene cannot be used simultaneously, pick one please"
            )

        # cc_eq = 0.14080860183951827  # nm
        cc_eq = 0.1418  # nm, distance between two carbon atoms in graphene
        x_disp = 3 * cc_eq * 10  # A, size of the graphene unit cell in the x direction
        y_disp = math.sqrt(
            3
        ) * cc_eq * 10  # A, size of the graphene unit cell in the y direction
        if Graphene:
            x_minim = max(graphene_size, box[0]) * 10
            y_minim = max(graphene_size, box[1]) * 10
            M = math.ceil(x_minim / x_disp)
            N = math.ceil(y_minim / y_disp)
            X = M * x_disp  # A
            Y = N * y_disp  # A
            self.ostream.print_info(
                f"Box size x: {box[0]:.3f} y: {box[1]:.3f} nm, minimum graphene size: {graphene_size:.3f} nm. Using largest to determine the size of the graphene sheet and size of the box"
            )
            self.ostream.print_info(
                f"Building graphene sheet with X: {X/10:.3f} nm, Y: {Y/10:.3f} nm"
            )
            box[0] = X / 10
            box[1] = Y / 10
        else:
            x_minim = box[0] * 10
            y_minim = CNT_radius * 2 * np.pi * 10
            M = math.ceil(x_minim / x_disp)
            N = math.ceil(y_minim / y_disp)
            X = M * x_disp  # A
            Y = N * y_disp  # A
            R = Y / (2 * np.pi)  # A
            self.ostream.print_info(
                f"Box size x: {box[0]:.3f}, minimum CNT radius: {CNT_radius:.3f} nm."
            )
            self.ostream.print_info(
                f"Building CNT with X: {X/10:.3f} nm, R: {R/10.:3f} nm")

        graphene_xyz = """4

        C        0.00000        1.21944        0.00000
        C        0.70404        0.00000        0.00000
        C        2.11212        0.00000        0.00000
        C        2.81616        1.21944        0.00000"""
        graphene_cell = Molecule.read_xyz_string(graphene_xyz)

        C_chain = topology.addChain()
        C_residue = topology.addResidue("CCC", chain=C_chain)
        CC_atoms = []

        # https://pubs.acs.org/doi/10.1021/jp011344u
        graphene_bond: dict = {
            'type': 'harmonic',
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

        graphene_dihedral: dict = {
            "type": "Fourier",
            "barrier": 25,
            "phase": 180,
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
        mol_positions = system_mol.get_coordinates_in_angstrom()
        middle_x = (max(mol_positions[:, 0]) + min(mol_positions[:, 0])) / 2
        middle_y = (max(mol_positions[:, 1]) + min(mol_positions[:, 1])) / 2
        min_z = min(mol_positions[:, 2])
        min_y = min(mol_positions[:, 1])
        if Graphene:
            x_offset = middle_x - X / 2
            y_offset = middle_y - Y / 2
            z_offset = min_z - 2
        else:
            x_offset = middle_x - X / 2
            y_offset = min_y - 0.7 * (R + 1)
            z_offset = min_z - 0.7 * (R + 1)

        for m in range(0, M):
            for n in range(0, N):
                for i, (atomic_number, coord) in enumerate(
                        zip(graphene_cell.get_element_ids(),
                            graphene_cell.get_coordinates_in_angstrom())):

                    coord[0] += m * x_disp
                    coord[1] += n * y_disp
                    if CNT:
                        phi = 2 * np.pi * coord[1] / Y
                        coord = np.array(
                            [coord[0], R * math.sin(phi), R * math.cos(phi)]
                        )  # A, X/2 factor is to center the CNT in the middle of the box. The X length is used to get the proper box size
                    coord += np.array([x_offset, y_offset, z_offset])

                    mm_element = mmapp.Element.getByAtomicNumber(atomic_number)
                    name = index
                    index += 1
                    atom = topology.addAtom(f"{name}", mm_element, C_residue)
                    carbon_atoms.append(atom)
                    CC_atoms.append(atom)

                    system_mol.add_atom(int(atomic_number), Point(coord),
                                        'angstrom')  # Bohr
                    positions = np.vstack([positions, coord])
                    system.addParticle(mm_element.mass)
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

                angles.update({
                    (center + 0, center + 1, center + 2):
                    graphene_angle
                })
                angles.update({
                    (center + 1, center + 2, center + 3):
                    graphene_angle
                })
                angles.update({
                    (up + 1, center + 0, center + 1): graphene_angle
                })
                angles.update({
                    (center + 2, center + 3, up + 2): graphene_angle
                })
                angles.update({(left + 3, center + 0, up + 1): graphene_angle})
                angles.update({
                    (left + 3, center + 0, center + 1): graphene_angle
                })
                angles.update({
                    (center + 0, center + 1, down + 0): graphene_angle
                })
                angles.update({
                    (down + 0, center + 1, center + 2): graphene_angle
                })
                angles.update({
                    (center + 1, center + 2, down + 3): graphene_angle
                })
                angles.update({
                    (center + 3, center + 2, down + 3): graphene_angle
                })
                angles.update({
                    (center + 2, center + 3, right + 0):
                    graphene_angle
                })
                angles.update({(right + 0, center + 3, up + 2): graphene_angle})

                dihedrals.update({
                    (left + 3, center + 0, center + 1, down + 0):
                    graphene_dihedral
                })
                dihedrals.update({
                    (up + 1, center + 0, center + 1, down + 0):
                    graphene_dihedral
                })
                dihedrals.update({
                    (left + 3, center + 0, center + 1, center + 2):
                    graphene_dihedral
                })
                dihedrals.update({
                    (left + 3, center + 0, center + 1, center + 2):
                    graphene_dihedral
                })

                dihedrals.update({
                    (center + 0, center + 1, center + 2, center + 3):
                    graphene_dihedral
                })
                dihedrals.update({
                    (down + 0, center + 1, center + 2, center + 3):
                    graphene_dihedral
                })
                dihedrals.update({
                    (center + 0, center + 1, center + 2, down + 3):
                    graphene_dihedral
                })
                dihedrals.update({
                    (down + 0, center + 1, center + 2, down + 3):
                    graphene_dihedral
                })

                dihedrals.update({
                    (center + 1, center + 2, center + 3, up + 2):
                    graphene_dihedral
                })
                dihedrals.update({
                    (down + 3, center + 2, center + 3, up + 2):
                    graphene_dihedral
                })
                dihedrals.update({
                    (center + 1, center + 2, center + 3, right + 0):
                    graphene_dihedral
                })
                dihedrals.update({
                    (down + 3, center + 2, center + 3, right + 0):
                    graphene_dihedral
                })

                dihedrals.update({
                    (center + 2, center + 3, right + 0, right + 1):
                    graphene_dihedral
                })
                dihedrals.update({
                    (up + 2, center + 3, right + 0, right + 1):
                    graphene_dihedral
                })
                dihedrals.update({
                    (center + 2, center + 3, right + 0, upright + 1):
                    graphene_dihedral
                })
                dihedrals.update({
                    (up + 2, center + 3, right + 0, upright + 1):
                    graphene_dihedral
                })

                dihedrals.update({
                    (left + 3, center + 0, up + 1, up + 0):
                    graphene_dihedral
                })
                dihedrals.update({
                    (center + 1, center + 0, up + 1, up + 0):
                    graphene_dihedral
                })
                dihedrals.update({
                    (left + 3, center + 0, up + 1, up + 2):
                    graphene_dihedral
                })
                dihedrals.update({
                    (center + 1, center + 0, up + 1, up + 2):
                    graphene_dihedral
                })

                dihedrals.update({
                    (center + 2, center + 3, up + 2, up + 1):
                    graphene_dihedral
                })
                dihedrals.update({
                    (center + 2, center + 3, up + 2, up + 3):
                    graphene_dihedral
                })
                dihedrals.update({
                    (right + 0, center + 3, up + 2, up + 1):
                    graphene_dihedral
                })
                dihedrals.update({
                    (right + 0, center + 3, up + 2, up + 3):
                    graphene_dihedral
                })

                impropers.update({
                    (left + 3, center + 0, center + 1, up + 2):
                    graphene_dihedral
                })
                impropers.update({
                    (center + 0, center + 1, down + 0, center + 2):
                    graphene_dihedral
                })
                impropers.update({
                    (center + 1, center + 2, down + 3, center + 3):
                    graphene_dihedral
                })
                impropers.update({
                    (center + 2, center + 3, right + 0, up + 2):
                    graphene_dihedral
                })

        carbon_harmonic_bond_force = mm.HarmonicBondForce()
        carbon_harmonic_bond_force.setName("Carbon harmonic bond")
        carbon_harmonic_bond_force.setUsesPeriodicBoundaryConditions(True)

        carbon_harmonic_angle_force = mm.HarmonicAngleForce()
        carbon_harmonic_angle_force.setName("Carbon harmonic angle")
        carbon_harmonic_angle_force.setUsesPeriodicBoundaryConditions(True)

        carbon_fourier_dihedral_force = mm.PeriodicTorsionForce()
        carbon_fourier_dihedral_force.setName("Carbon improper fourier torsion")
        carbon_fourier_dihedral_force.setUsesPeriodicBoundaryConditions(True)

        for key, bond in bonds.items():
            atom_ids = self._key_to_id(key, CC_atoms)
            self._add_bond(
                carbon_harmonic_bond_force,
                atom_ids,
                bond['equilibrium'],
                bond['force_constant'],
            )

        for key, angle in angles.items():
            atom_ids = self._key_to_id(key, CC_atoms)
            self._add_angle(
                carbon_harmonic_angle_force,
                atom_ids,
                angle['equilibrium'],
                angle['force_constant'],
            )

        for key, dihedral in dihedrals.items():
            atom_ids = self._key_to_id(key, CC_atoms)
            self._add_torsion(carbon_fourier_dihedral_force, dihedral, atom_ids)

        for key, improper in impropers.items():
            atom_ids = self._key_to_id(key, CC_atoms)
            self._add_torsion(carbon_fourier_dihedral_force, improper, atom_ids)

        for i in range(len(carbon_atoms)):
            for j in range(len(carbon_atoms)):
                if j > i:
                    nb_force.addException(carbon_atoms[i].index,
                                          carbon_atoms[j].index, 0, 1, 0)

        carbon_harmonic_bond_force.setForceGroup(EvbForceGroup.CARBON.value)
        carbon_harmonic_angle_force.setForceGroup(EvbForceGroup.CARBON.value)
        carbon_fourier_dihedral_force.setForceGroup(EvbForceGroup.CARBON.value)
        system.addForce(carbon_harmonic_bond_force)
        system.addForce(carbon_harmonic_angle_force)
        system.addForce(carbon_fourier_dihedral_force)

        return box

    def _add_solvent(self, system, system_mol, solvent, topology, nb_force,
                     neutralize, padding):

        assert_msg_critical('openmm' in sys.modules,
                            'openmm is required for EvbSystemBuilder.')

        vlxsysbuilder = SolvationBuilder()
        mols_per_nm3, density, smiles_code = vlxsysbuilder._solvent_properties(
            solvent)
        vlxsysbuilder.solvate(
            solute=system_mol,
            solvent=solvent,
            padding=padding,
            neutralize=neutralize,
            target_density=density * 0.95,
        )

        self.positions = vlxsysbuilder.system_molecule.get_coordinates_in_angstrom(
        ) * 0.1
        box = [side * 0.1 for side in vlxsysbuilder.box]

        solvents = vlxsysbuilder.solvents
        solvent_counts = vlxsysbuilder.added_solvent_counts
        resnames = ["SOL"] * len(vlxsysbuilder.added_solvent_counts)
        if hasattr(vlxsysbuilder, "added_counterions"):
            if vlxsysbuilder.added_counterions > 0:
                solvents = solvents + [vlxsysbuilder._counterion_molecules()
                                       ] * vlxsysbuilder.added_counterions
                solvent_counts = solvent_counts + [
                    vlxsysbuilder.added_counterions
                ]
                resnames = resnames + ["ION"]

        for i, (solvent_count, vlx_solvent_molecule,
                resname) in enumerate(zip(solvent_counts, solvents, resnames)):

            num_solvent_atoms_per_molecule = vlx_solvent_molecule.number_of_atoms(
            )
            elements = vlx_solvent_molecule.get_labels()

            harmonic_bond_force = mm.HarmonicBondForce()
            harmonic_bond_force.setName(f"Solvent {i} harmonic bond")
            harmonic_angle_force = mm.HarmonicAngleForce()
            harmonic_angle_force.setName(f"Solvent {i} harmonic angle")
            fourier_force = mm.PeriodicTorsionForce()
            fourier_force.setName(f"Solvent {i} proper fourier torsion")
            fourier_imp_force = mm.PeriodicTorsionForce()
            fourier_imp_force.setName(f"Solvent {i} improper fourier torsion")

            # Loop over all solvent molecules
            solvent_nb_atom_count = 0
            solvent_system_atom_count = 0
            solvent_chain = topology.addChain()
            solvent_ff = MMForceFieldGenerator()
            solvent_ff.create_topology(vlx_solvent_molecule)

            for atom in solvent_ff.atoms.values():
                if atom['type'] == 'ow':
                    sigma = 1.8200 * 2**(-1 / 6) * 2 / 10
                    epsilon = 0.0930 * 4.184
                    atom['type'] = 'oh'
                    atom['sigma'] = sigma
                    atom['epsilon'] = epsilon
                    atom['comment'] = "Reaction-water oxygen"
                elif atom['type'] == 'hw':
                    sigma = 0.3019 * 2**(-1 / 6) * 2 / 10
                    epsilon = 0.0047 * 4.184
                    atom['type'] = 'ho'
                    atom['sigma'] = sigma
                    atom['epsilon'] = epsilon
                    atom['comment'] = "Reaction-water hydrogen"

            for i in range(solvent_count):

                solvent_residue = topology.addResidue(name=resname,
                                                      chain=solvent_chain)
                solvent_atoms = []
                # Loop over all atoms in the solvent molecule
                for j in range(num_solvent_atoms_per_molecule):
                    # Figure out the element of the atom
                    mm_element = mmapp.Element.getBySymbol(elements[j])
                    name = f"{elements[j]}{j}"
                    # add the atom to the topology
                    solvent_atom = topology.addAtom(name, mm_element,
                                                    solvent_residue)
                    solvent_atoms.append(solvent_atom)
                    # Add the atom as a particle to the system
                    system.addParticle(mm_element.mass)
                    solvent_system_atom_count += 1
                    # Add this atom to the nonbonded force

                # add all bonded interactions in this molecule
                for key, bond in solvent_ff.bonds.items():
                    atom_ids = self._key_to_id(key, solvent_atoms)
                    self._add_bond(
                        harmonic_bond_force,
                        atom_ids,
                        bond['equilibrium'],
                        bond['force_constant'],
                    )

                for key, angle in solvent_ff.angles.items():
                    atom_ids = self._key_to_id(key, solvent_atoms)
                    self._add_angle(
                        harmonic_angle_force,
                        atom_ids,
                        angle['equilibrium'],
                        angle['force_constant'],
                    )

                for key, dihedral in solvent_ff.dihedrals.items():
                    atom_ids = self._key_to_id(key, solvent_atoms)
                    self._add_torsion(fourier_force, dihedral, atom_ids)
                for key, dihedral in solvent_ff.impropers.items():
                    atom_ids = self._key_to_id(key, solvent_atoms)
                    self._add_torsion(fourier_imp_force, dihedral, atom_ids)

                exceptions = self._create_exceptions_from_bonds(
                    solvent_ff.atoms, solvent_ff.bonds)

                for key, atom in solvent_ff.atoms.items():

                    sigma = atom["sigma"]
                    if sigma == 0:
                        if atom["epsilon"] == 0:
                            sigma = 1
                        else:
                            raise ValueError(
                                "Sigma is 0 while epsilon is not, which will cause division by 0"
                            )
                    solvent_nb_atom_count += 1
                    nb_force.addParticle(atom["charge"], atom["sigma"],
                                         atom["epsilon"])
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

            if harmonic_bond_force.getNumBonds() > 0:
                harmonic_bond_force.setForceGroup(EvbForceGroup.SOLVENT.value)
                system.addForce(harmonic_bond_force)
            if harmonic_angle_force.getNumAngles() > 0:
                harmonic_angle_force.setForceGroup(EvbForceGroup.SOLVENT.value)
                system.addForce(harmonic_angle_force)
            if fourier_force.getNumTorsions() > 0:
                fourier_force.setForceGroup(EvbForceGroup.SOLVENT.value)
                system.addForce(fourier_force)
            if fourier_imp_force.getNumTorsions() > 0:
                fourier_imp_force.setForceGroup(EvbForceGroup.SOLVENT.value)
                system.addForce(fourier_imp_force)
            self.ostream.print_info(
                f"Added {solvent_nb_atom_count} atoms to the nonbonded force and {solvent_system_atom_count} atoms to the system"
            )
        return box

    def _add_reaction_forces(self, system, lam):

        assert_msg_critical('openmm' in sys.modules,
                            'openmm is required for EvbSystemBuilder.')
        bonded_harmonic, bonded_integration, morse_force, max_distance, = self._create_bond_forces(
            lam)
        angle, angle_integration = self._create_angle_forces(lam)
        torsion = self._create_proper_torsion_forces(lam)
        improper = self._create_improper_torsion_forces(lam)

        # The bonded_integration flag causes the nonbonded exceptions to be created as if the bonds of the reactant and product are all present all the time, and thus to exclude these nonbonded interactions over the entire lambda vector.
        # This is compensated through the extra bonded interactions (labeled with integration)
        # This only affects the integration potential as the hard core interactions are used for integration
        intlj, intcoul = self._create_nonbonded_forces(
            lam,
            bonded=self.bonded_integration,
            lj_soft_core=self.soft_core_lj_int,
            coul_soft_core=self.soft_core_coulomb_int,
        )
        # The soft core forces are used for the PES calculations, and thus always have the 'correct' exceptions
        peslj, pescoul = self._create_nonbonded_forces(
            lam,
            bonded=False,
            lj_soft_core=self.soft_core_lj_pes,
            coul_soft_core=self.soft_core_coulomb_pes,
        )
        # The hard
        bond_constraint, constant_force, angle_constraint, torsion_constraint = self._create_constraint_forces(
            lam)

        bonded_harmonic.setForceGroup(EvbForceGroup.REACTION_BONDED.value)
        angle.setForceGroup(EvbForceGroup.REACTION_BONDED.value)
        torsion.setForceGroup(EvbForceGroup.REACTION_BONDED.value)
        improper.setForceGroup(EvbForceGroup.REACTION_BONDED.value)

        intlj.setForceGroup(EvbForceGroup.INTLJ.value)
        intcoul.setForceGroup(EvbForceGroup.INTCOUL.value)
        peslj.setForceGroup(EvbForceGroup.PESLJ.value)
        pescoul.setForceGroup(EvbForceGroup.PESCOUL.value)
        bond_constraint.setForceGroup(EvbForceGroup.CONSTRAINT.value)
        constant_force.setForceGroup(EvbForceGroup.CONSTRAINT.value)
        angle_constraint.setForceGroup(EvbForceGroup.CONSTRAINT.value)
        torsion_constraint.setForceGroup(EvbForceGroup.CONSTRAINT.value)

        system.addForce(bonded_harmonic)
        system.addForce(angle)
        system.addForce(torsion)
        system.addForce(improper)
        system.addForce(intlj)
        system.addForce(intcoul)
        system.addForce(peslj)
        system.addForce(pescoul)
        system.addForce(bond_constraint)
        system.addForce(constant_force)
        system.addForce(angle_constraint)
        system.addForce(torsion_constraint)

        if self.bonded_integration:
            bonded_integration.setForceGroup(EvbForceGroup.INTEGRATION.value)
            angle_integration.setForceGroup(EvbForceGroup.INTEGRATION.value)
            morse_force.setForceGroup(EvbForceGroup.PES.value)
            system.addForce(bonded_integration)
            system.addForce(angle_integration)
            system.addForce(morse_force)
        else:
            morse_force.setForceGroup(EvbForceGroup.REACTION_BONDED.value)
            max_distance.setForceGroup(EvbForceGroup.RESTRAINT.value)
            system.addForce(morse_force)
            system.addForce(max_distance)
        return system

    def _create_E_field(self, system, E_field):

        assert_msg_critical('openmm' in sys.modules,
                            'openmm is required for EvbSystemBuilder.')

        E_field_force = mm.CustomExternalForce("-q*(Ex*x+Ey*y+Ez*z)")
        E_field_force.addGlobalParameter("Ex", E_field[0])
        E_field_force.addGlobalParameter("Ey", E_field[1])
        E_field_force.addGlobalParameter("Ez", E_field[2])
        E_field_force.addPerParticleParameter("q")
        E_field_force.setName("Electric field")
        for i in range(system.getNumParticles()):
            E_field_force.addParticle(
                i,
                [0])  # Actual charges are lambda dependent, and are set later
        return E_field_force

    def _create_bond_forces(self, lam):

        assert_msg_critical('openmm' in sys.modules,
                            'openmm is required for EvbSystemBuilder.')

        harmonic_force = mm.HarmonicBondForce()
        harmonic_force.setName("Reaction harmonic bond")

        integration_force = mm.HarmonicBondForce()
        integration_force.setName("Integration reaction harmonic bond")

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
        static_bond_keys = list(
            set(self.reactant.bonds) & set(self.product.bonds))
        broken_bond_keys = list(
            set(self.reactant.bonds) - set(self.product.bonds))
        formed_bond_keys = list(
            set(self.product.bonds) - set(self.reactant.bonds))

        for key in bond_keys:
            atom_ids = self._key_to_id(key, self.reaction_atoms)
            if key in static_bond_keys:
                bondA = self.reactant.bonds[key]
                bondB = self.product.bonds[key]
                self._add_bond(
                    harmonic_force,
                    atom_ids,
                    bondA['equilibrium'],
                    bondA['force_constant'] * (1 - lam),
                )
                self._add_bond(
                    harmonic_force,
                    atom_ids,
                    bondB['equilibrium'],
                    bondB['force_constant'] * lam,
                )
            else:
                if key in broken_bond_keys:
                    scale = 1 - lam
                    bond = self.reactant.bonds[key]
                    broken_coords = self.product.molecule.get_coordinates_in_angstrom(
                    )
                elif key in formed_bond_keys:
                    scale = lam
                    bond = self.product.bonds[key]
                    broken_coords = self.reactant.molecule.get_coordinates_in_angstrom(
                    )
                else:
                    assert (
                        False
                    ), "A bond can either be static, or dynamic, in  which case it can be broken or formed"
                broken_length = (self.measure_length(
                    broken_coords[key[0]], broken_coords[key[1]]) * 0.1)

                # self._add_harm_bond(harmonic_dynamic_force, bond, atom_ids, scale)
                self._add_morse_bond(morse_force, bond, atom_ids, scale)
                self._add_distance_restraint(max_distance, atom_ids,
                                             broken_length, 1 - scale)

                # This the force constant is not scaled down (fully) with lambda, because if bonded_integration is turned on, the morse force is not included in the integration potential
                # The bonded integration factor is scaled up based on lambda, to make sure that there is a full bonded interaction when the bond is present
                self._add_bond(
                    integration_force,
                    atom_ids,
                    bond['equilibrium'] * scale + broken_length * (1 - scale),
                    bond['force_constant'] *
                    (1 - scale * self.bonded_integration_fac),
                )
        return harmonic_force, integration_force, morse_force, max_distance

    def _create_angle_forces(self, lam: float):

        assert_msg_critical('openmm' in sys.modules,
                            'openmm is required for EvbSystemBuilder.')

        harmonic_force = mm.HarmonicAngleForce()
        harmonic_force.setName("Reaction angle")
        integration_force = mm.HarmonicAngleForce()
        integration_force.setName("Integration reaction angle")

        angle_keys = list(set(self.reactant.angles) | set(self.product.angles))
        for key in angle_keys:
            atom_ids = self._key_to_id(key, self.reaction_atoms)
            if key in self.reactant.angles.keys(
            ) and key in self.product.angles.keys():
                angleA = self.reactant.angles[key]
                angleB = self.product.angles[key]
                self._add_angle(harmonic_force, atom_ids, angleA['equilibrium'],
                                angleA['force_constant'] * (1 - lam))
                self._add_angle(harmonic_force, atom_ids, angleB['equilibrium'],
                                angleB['force_constant'] * (lam))
            else:
                if key in self.reactant.angles.keys():
                    scale = 1 - lam
                    angle = self.reactant.angles[key]
                    broken_coords = self.product.molecule.get_coordinates_in_angstrom(
                    )
                else:
                    scale = lam
                    angle = self.product.angles[key]
                    broken_coords = self.reactant.molecule.get_coordinates_in_angstrom(
                    )
                broken_equil = self.measure_angle(broken_coords[key[0]],
                                                  broken_coords[key[1]],
                                                  broken_coords[key[2]])

                self._add_angle(harmonic_force, atom_ids, angle['equilibrium'],
                                angle['force_constant'] * scale)
                self._add_angle(integration_force, atom_ids, broken_equil,
                                angle['force_constant'] * (1 - scale)*self.bonded_integration_fac)

        return harmonic_force, integration_force

    def _create_proper_torsion_forces(self, lam):

        assert_msg_critical('openmm' in sys.modules,
                            'openmm is required for EvbSystemBuilder.')

        fourier_force = mm.PeriodicTorsionForce()
        fourier_force.setName("Reaction proper fourier torsion")

        dihedral_keys = list(
            set(self.reactant.dihedrals) | set(self.product.dihedrals))
        for key in dihedral_keys:
            atom_ids = self._key_to_id(key, self.reaction_atoms)
            if (key in self.reactant.dihedrals.keys()
                    and key in self.product.dihedrals.keys()):
                dihedA = self.reactant.dihedrals[key]
                dihedB = self.product.dihedrals[key]
                self._add_torsion(fourier_force, dihedA, atom_ids, 1 - lam)
                self._add_torsion(fourier_force, dihedB, atom_ids, lam)
            else:
                if key in self.reactant.dihedrals.keys():
                    scale = 1 - lam
                    dihed = self.reactant.dihedrals[key]
                else:
                    scale = lam
                    dihed = self.product.dihedrals[key]
                if scale > 0:
                    self._add_torsion(fourier_force, dihed, atom_ids, scale)
        return fourier_force

    def _create_improper_torsion_forces(self, lam):

        assert_msg_critical('openmm' in sys.modules,
                            'openmm is required for EvbSystemBuilder.')

        fourier_force = mm.PeriodicTorsionForce()
        fourier_force.setName("Reaction improper fourier torsion")

        dihedral_keys = list(
            set(self.reactant.impropers) | set(self.product.impropers))
        for key in dihedral_keys:
            atom_ids = self._key_to_id(key, self.reaction_atoms)
            if (key in self.reactant.impropers.keys()
                    and key in self.product.impropers.keys()):
                dihedA = self.reactant.impropers[key]
                self._add_torsion(fourier_force, dihedA, atom_ids, 1 - lam)
                dihedB = self.product.impropers[key]
                self._add_torsion(fourier_force, dihedB, atom_ids, lam)
            else:
                if key in self.reactant.impropers.keys():
                    scale = 1 - lam
                    dihed = self.reactant.impropers[key]
                else:
                    scale = lam
                    dihed = self.product.impropers[key]
                if scale > 0:
                    self._add_torsion(fourier_force, dihed, atom_ids, scale)
        return fourier_force

    def _add_bond(self, bond_force, atom_id, equil, fc):
        if fc > 0:
            bond_force.addBond(
                atom_id[0],
                atom_id[1],
                equil,
                fc,
            )

    def _add_morse_bond(self,
                        bond_force,
                        bond_dict,
                        atom_id,
                        barrier_scaling=1.):
        if "D" not in bond_dict.keys():
            D = self.morse_D_default
            if self.verbose:
                self.ostream.print_info(
                    f"No D value associated with bond {atom_id[0]} {atom_id[1]}. Setting to default value {self.morse_D_default}"
                )
        else:
            D = bond_dict["D"]

        assert "force_constant" in bond_dict.keys(
        ), "No force constant found for bond"
        assert "equilibrium" in bond_dict.keys(
        ), "No equilibrium distance found for bond"

        a = math.sqrt(bond_dict["force_constant"] / (2 * D))
        re = bond_dict["equilibrium"]

        if barrier_scaling * D > 0:
            bond_force.addBond(
                atom_id[0],
                atom_id[1],
                [barrier_scaling * D, a, re],
            )

    def _add_distance_restraint(self,
                                restraint_force,
                                atom_id,
                                broken_length,
                                barrier_scaling=1.):
        rmax = broken_length

        rmax += self.restraint_r_offset
        k = self.restraint_k * barrier_scaling

        if k > 0:
            restraint_force.addBond(
                atom_id[0],
                atom_id[1],
                [rmax, k],
            )

    def _add_angle(self, angle_force, atom_id, equil, fc):
        if fc > 0:
            angle_force.addAngle(
                atom_id[0],
                atom_id[1],
                atom_id[2],
                equil * self.deg_to_rad,
                fc,
            )

    def _add_torsion(self,
                     fourier_force,
                     torsion_dict,
                     atom_id,
                     barrier_scaling=1.):
        assert torsion_dict["type"] == "Fourier", "Unknown dihedral type"
        if barrier_scaling > 0:
            if torsion_dict.get("multiple", False):
                for periodicity, phase, barrier in zip(
                        torsion_dict["periodicity"], torsion_dict["phase"],
                        torsion_dict["barrier"]):
                    if barrier > 0:
                        fourier_force.addTorsion(
                            atom_id[0],
                            atom_id[1],
                            atom_id[2],
                            atom_id[3],
                            abs(periodicity),
                            phase * self.deg_to_rad,
                            barrier_scaling * barrier,
                        )
            else:
                if torsion_dict["barrier"] > 0:
                    fourier_force.addTorsion(
                        atom_id[0],
                        atom_id[1],
                        atom_id[2],
                        atom_id[3],
                        torsion_dict["periodicity"],
                        torsion_dict["phase"] * self.deg_to_rad,
                        barrier_scaling * torsion_dict["barrier"],
                    )

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

        hard_core_expression = (" (1-l) * k*qqA/r + l * k*qqB/r; "
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

        hard_core_expression = (" (1-l) * LjtotA "
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

    def _create_nonbonded_forces(self,
                                 lam,
                                 bonded=False,
                                 lj_soft_core=False,
                                 coul_soft_core=False):

        assert_msg_critical('openmm' in sys.modules,
                            'openmm is required for EvbSystemBuilder.')

        coulomb_force = mm.CustomBondForce(
            self._get_couloumb_expression(coul_soft_core))
        name = "Reaction internal coulomb"
        if coul_soft_core:
            name += " soft-core"

        coulomb_force.setName(name)
        coulomb_force.addPerBondParameter("qqA")
        coulomb_force.addPerBondParameter("qqB")
        coulomb_force.addGlobalParameter("l", lam)

        lj_force = mm.CustomBondForce(self._get_lj_expression(lj_soft_core))
        name = "Reaction internal lennard-jones"
        if lj_soft_core:
            name += " soft-core"
        lj_force.setName(name)
        lj_force.addPerBondParameter("sigmaA")
        lj_force.addPerBondParameter("sigmaB")
        lj_force.addPerBondParameter("epsilonA")
        lj_force.addPerBondParameter("epsilonB")
        lj_force.addGlobalParameter("l", lam)

        if not bonded:
            reactant_bonds = self.reactant.bonds
            product_bonds = self.product.bonds
        else:
            bonds = list(set(self.reactant.bonds) | set(self.product.bonds))
            reactant_bonds = bonds
            product_bonds = bonds
        reactant_exceptions = self._create_exceptions_from_bonds(
            self.reactant.atoms, reactant_bonds)
        product_exceptions = self._create_exceptions_from_bonds(
            self.product.atoms, product_bonds)

        #Loop over all atoms, and check if their id's are part of any exceptions
        for i in self.reactant.atoms.keys():
            for j in self.reactant.atoms.keys():
                if i < j:
                    key = (i, j)
                    # Remove any exception from the nonbondedforce
                    # and add it instead to the exception bond force
                    if key in reactant_exceptions.keys():
                        epsilonA = reactant_exceptions[key]["epsilon"]
                        sigmaA = reactant_exceptions[key]["sigma"]
                        qqA = reactant_exceptions[key]["qq"]

                    else:
                        atomA1 = self.reactant.atoms[key[0]]
                        atomA2 = self.reactant.atoms[key[1]]
                        epsilonA = math.sqrt(atomA1["epsilon"] *
                                             atomA2["epsilon"])
                        sigmaA = 0.5 * (atomA1["sigma"] + atomA2["sigma"])
                        qqA = atomA1["charge"] * atomA2["charge"]

                    if key in product_exceptions.keys():
                        epsilonB = product_exceptions[key]["epsilon"]
                        sigmaB = product_exceptions[key]["sigma"]
                        qqB = product_exceptions[key]["qq"]
                    else:
                        atomB1 = self.product.atoms[key[0]]
                        atomB2 = self.product.atoms[key[1]]
                        epsilonB = math.sqrt(atomB1["epsilon"] *
                                             atomB2["epsilon"])
                        sigmaB = 0.5 * (atomB1["sigma"] + atomB2["sigma"])
                        qqB = atomB1["charge"] * atomB2["charge"]

                    if sigmaA == 1.0:
                        sigmaA = sigmaB
                    elif sigmaB == 1.0:
                        sigmaB = sigmaA
                    atom_ids = self._key_to_id(key, self.reaction_atoms)
                    if not (qqA == 0.0 and qqB == 0.0):
                        coulomb_force.addBond(
                            atom_ids[0],
                            atom_ids[1],
                            [qqA, qqB],
                        )
                    if not (epsilonA == 0.0 and epsilonB == 0.0):
                        lj_force.addBond(
                            atom_ids[0],
                            atom_ids[1],
                            [sigmaA, sigmaB, epsilonA, epsilonB],
                        )

        return lj_force, coulomb_force

    def _create_exceptions_from_bonds(
            self, particles, bonds) -> dict[tuple[int, int], dict[str, float]]:

        exclusions = [set() for _ in range(len(particles))]
        bonded12 = [set() for _ in range(len(particles))]
        exceptions = {}

        # Populate bonded12 with bonds.
        for bond in bonds:
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
                        charge_prod = (self.coul14_scale * particle1["charge"] *
                                       particle2["charge"])
                        sigma = 0.5 * (particle1["sigma"] + particle2["sigma"])
                        epsilon = (
                            self.lj14_scale *
                            (particle1["epsilon"] * particle2["epsilon"])**0.5)
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
    def _add_exclusions_to_set(bonded12, exclusions, base_particle,
                               from_particle, current_level):
        for i in bonded12[from_particle]:
            if i != base_particle:
                exclusions.add(i)
            if current_level > 0:
                EvbSystemBuilder._add_exclusions_to_set(bonded12, exclusions,
                                                        base_particle, i,
                                                        current_level - 1)

    def _create_constraint_forces(self, lam):

        assert_msg_critical('openmm' in sys.modules,
                            'openmm is required for EvbSystemBuilder.')

        bond_constraint = mm.HarmonicBondForce()
        bond_constraint.setName("Bond constraint")

        constant_force = mm.CustomBondForce("-k*r")
        constant_force.setName("Constant bond force constraint")
        constant_force.addPerBondParameter("k")

        angle_constraint = mm.HarmonicAngleForce()
        angle_constraint.setName("Angle constraint")
        torsion_constraint = mm.CustomTorsionForce("0.5*k*(theta-theta0)^2")
        torsion_constraint.setName("Harmonic torsion constraint")
        torsion_constraint.addPerTorsionParameter("theta0")
        torsion_constraint.addPerTorsionParameter("k")
        if len(self.constraints) > 0:
            if self.verbose:
                self.ostream.print_info(
                    f"Adding constraints: {self.constraints}")

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

            atomids = self._key_to_id(key, self.reaction_atoms)
            if len(key) == 2:
                if constraint[key]['type'] == 'harmonic':
                    self._add_bond(
                        bond_constraint,
                        atomids,
                        constraint[key]["equilibrium"],
                        constraint[key]["force_constant"] * scale,
                    )
                elif constraint[key]['type'] == 'linear':
                    constant_force.addBond(
                        atomids[0],
                        atomids[1],
                        [constraint[key]["force_constant"] * scale],
                    )
                else:
                    raise ValueError(
                        f"Unknown constraint bond type {constraint[key]['type']}"
                    )
            if len(key) == 3:
                self._add_angle(
                    angle_constraint,
                    atomids,
                    constraint[key]["equilibrium"],
                    constraint[key]["force_constant"] * scale,
                )

            if len(key) == 4:
                torsion_constraint.addTorsion(
                    atomids[0],
                    atomids[1],
                    atomids[2],
                    atomids[3],
                    [
                        constraint[key]["equilibrium"] * self.deg_to_rad,
                        constraint[key]["force_constant"] * scale,
                    ],
                )
        return bond_constraint, constant_force, angle_constraint, torsion_constraint

    def _key_to_id(self, key: tuple[int, ...], atom_list: list) -> list[int]:
        return [atom_list[key[i]].index for i in range(len(key))]

    @staticmethod
    def measure_length(v1, v2):
        """
        Calculates the distance between v1 and v2
        """

        return np.linalg.norm(np.array(v1) - np.array(v2))

    @staticmethod
    def measure_angle(v1, v2, v3):
        """
        Calculates the angle between v1 and v2 and v3
        """

        vector_u = (v1 - v2) / np.linalg.norm(v1 - v2)
        vector_v = (v3 - v2) / np.linalg.norm(v3 - v2)

        # Calculate the angle using the dot product
        dot_product = np.dot(vector_u, vector_v)
        angle = np.arccos(dot_product)
        angle = angle * 180 / np.pi  # Convert to degrees

        return angle


# The EVB procedure uses two different potentials, one for the integration of the EOMs to explore phase space, and another one for the calculation of the PES
# The integration potential is optimized to explore all relevant areas of phase space in an efficient manner. This includes constraints and distance restraints.
# The PES potential is optimized to calculate the potential energy surface accurate. Constraints and restraints are ommitted for this.
# Soft core long range potentials provide faster convergence for averages and are thus included in the PES potential. They cause unstable integration though.
class EvbForceGroup(Enum):
    DEFAULT = auto(
    )  # Default force group, included for both integration and energy calculations
    # THERMOSTAT = auto()  # Thermostat
    INTLJ = auto()  # Integration lennard-jones potential
    INTCOUL = auto()  # Integration coulombic potential
    PESLJ = auto()  # Lennard-jones potential
    PESCOUL = auto()  # Coulombic potential
    CONSTRAINT = auto()

    CMM_REMOVER = auto()  # Center of mass motion remover
    NB_FORCE = auto()  # Solvent-solvent and solvent-solute nb force
    BAROSTAT = auto()  # Barostat
    E_FIELD = auto()  # Electric field force
    REACTION_BONDED = auto()  # Bonded forces for the reaction atoms

    # Constraints that also should be included in the PES calculations. Currently only used for the linear bond constraint
    PES_CONSTRAINT = auto()
    RESTRAINT = auto()
    SOLVENT = auto(
    )  # All solvent-solvent interactions. Does not include the solute-solvent long range interaction
    CARBON = auto()  # Graphene and CNTs
    INTEGRATION = auto(
    )  # All leftover forces that should only be used for integration
    PES = auto(
    )  # All leftover forces that should only be used for the calculation of the PES
    NONE = auto(
    )  # Forces that should not be included in the integration or PES calculations
    DEBUG1PES = auto()  # Debugging force group 1
    DEBUG2PES = auto()  # Debugging force group 2
    DEBUG1INT = auto()  # Debugging force group 1
    DEBUG2INT = auto()  # Debugging force group 2
    DEBUG1 = auto()  # Debugging force group 1
    DEBUG2 = auto()  # Debugging force group 2

    # Both methods return classes because integrator.setIntegrationForceGroups() takes a set as argument
    @classmethod
    def integration_force_groups(cls):
        return set([
            cls.DEFAULT.value,
            cls.CMM_REMOVER.value,
            cls.NB_FORCE.value,
            cls.BAROSTAT.value,
            cls.E_FIELD.value,
            cls.REACTION_BONDED.value,
            cls.INTLJ.value,
            cls.INTCOUL.value,
            cls.CONSTRAINT.value,
            cls.PES_CONSTRAINT.value,
            cls.RESTRAINT.value,
            cls.SOLVENT.value,
            cls.CARBON.value,
            cls.INTEGRATION.value,
            # cls.DEBUG1INT.value,
            # cls.DEBUG2INT.value,
            # cls.DEBUG1.value,
            # cls.DEBUG2.value,
        ])

    @classmethod
    def NVT_integration_force_groups(cls):
        int_fg = cls.integration_force_groups()
        int_fg.remove(cls.BAROSTAT.value)
        return int_fg

    @classmethod
    def pes_force_groups(cls):
        return set([
            cls.DEFAULT.value,
            cls.CMM_REMOVER.value,
            cls.NB_FORCE.value,
            cls.BAROSTAT.value,
            cls.E_FIELD.value,
            cls.REACTION_BONDED.value,
            cls.PES_CONSTRAINT.value,
            cls.PESLJ.value,
            cls.PESCOUL.value,
            cls.SOLVENT.value,
            cls.CARBON.value,
            cls.PES.value,
            # cls.DEBUG1PES.value,
            # cls.DEBUG2PES.value,
            # cls.DEBUG1.value,
            # cls.DEBUG2.value,
        ])

    @classmethod
    def get_header(cls):
        header = ""
        integration_forcegroups = cls.integration_force_groups()
        pes_forcegroups = cls.pes_force_groups()
        for fg in cls:
            in_int = fg.value in integration_forcegroups
            in_pes = fg.value in pes_forcegroups
            fg_cat = 'b'
            if in_int and not in_pes:
                fg_cat = 'i'
            if not in_int and in_pes:
                fg_cat = 'p'
            header += f"{fg.name}({fg.value}-{fg_cat}), "
        header = header[:-2] + '\n'
        return header

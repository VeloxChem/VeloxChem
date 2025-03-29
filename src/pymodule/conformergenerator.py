from mpi4py import MPI
from pathlib import Path
from copy import deepcopy
import numpy as np
import itertools
import time
import re

from .veloxchemlib import bohr_in_angstrom, mpi_master
from .molecule import Molecule
from .mmforcefieldgenerator import MMForceFieldGenerator

from openmm.app import GromacsTopFile, NoCutoff, Simulation
from openmm import LangevinIntegrator
from openmm import Platform  # Platform added
from openmm.unit import nanometer, md_unit_system, kelvin, picoseconds, picosecond


class ConformerGenerator:

    def __init__(self, comm=None):

        if comm is None:
            comm = MPI.COMM_WORLD

        self._comm = comm
        self._rank = comm.Get_rank()
        self._size = comm.Get_size()

        self.molecule = None
        self.number_of_conformers_to_select = 50

        self.top_file_name = None

        self.save_xyz_files = False
        self.save_path = None

        self.em_tolerance = 1.0

    def _convert_molecule_xyz_string(self, labels, coords, comment=""):

        if len(labels) != len(coords):
            raise ValueError("The length of labels and coords should be the same")
        xyz_string = ""
        xyz_string += str(len(labels)) + "\n"
        xyz_string += comment + "\n"  # add comment line to second line
        for i in range(len(labels)):
            xyz_string += (
                labels[i]
                + "  "
                + str(coords[i][0])
                + "  "
                + str(coords[i][1])
                + "  "
                + str(coords[i][2])
                + "\n"
            )
        return xyz_string

    def _write_molecule_xyz_file(self, path, labels, coords, comment=""):

        xyz_string = self._convert_molecule_xyz_string(labels, coords, comment)
        with open(path, "w") as f:
            f.write(xyz_string)

    # TODO: use equivalent atoms
    def _check_methyl_in_dihedrals(self, dihedral_indices, dihedrals_dict, atom_info_dict):

        dih_comment = dihedrals_dict[dihedral_indices]["comment"]
        if isinstance(dih_comment, list):  # if there are multiple comments as a list
            dih_comment = dih_comment[0]

        dih_comment = re.sub(r" ", "", dih_comment)  # comment out
        dih_comment = dih_comment.split("-")  # comment out
        if "c3" not in dih_comment:
            return False
        else:
            if (
                dih_comment[1] == "c3" or dih_comment[2] == "c3"
            ):  # check if starts from methyl group
                # if c3 connected to 3H
                # we check twice because the order of comment could be wrong
                # round1
                c3_index = dihedral_indices[1] + 1  # convert to 1 based index
                c3_bonded_atoms = atom_info_dict[c3_index]["ConnectedAtoms"]
                # if 3H in the connected atoms
                if len(c3_bonded_atoms) == 4 and c3_bonded_atoms.count("H") == 3:
                    return True
                # round2
                c3_index = dihedral_indices[2] + 1
                c3_bonded_atoms = atom_info_dict[c3_index]["ConnectedAtoms"]
                if len(c3_bonded_atoms) == 4 and c3_bonded_atoms.count("H") == 3:
                    return True
        return False

    def _get_dihedral_candidates(self, molecule, top_file_name):

        mmff_gen = MMForceFieldGenerator(self._comm)
        mmff_gen.ostream.mute()
        # TODO: double check partial charge
        mmff_gen.partial_charges = molecule.get_partial_charges(molecule.get_charge())
        mmff_gen.create_topology(molecule)
        mmff_gen.write_gromacs_files(filename=top_file_name)
        self._comm.barrier()

        atom_info_dict = deepcopy(mmff_gen.atom_info_dict)
        rotatable_bonds = deepcopy(mmff_gen.rotatable_bonds)
        dihedrals_dict = deepcopy(mmff_gen.dihedrals)

        rotatable_bonds = [[i[0] - 1, i[1] - 1] for i in rotatable_bonds]
        rotatable_dihedrals_dict = {}

        # only pick one dihedral for each rotatable bond
        for key, value in dihedrals_dict.items():
            bond = [key[1], key[2]]
            sorted_bond = sorted(bond)

            if [sorted_bond[0], sorted_bond[1]] in rotatable_bonds:
                if (sorted_bond[0], sorted_bond[1]) not in rotatable_dihedrals_dict:
                    rotatable_dihedrals_dict[(sorted_bond[0], sorted_bond[1])] = value
                    rotatable_dihedrals_dict[(sorted_bond[0], sorted_bond[1])][
                        "dihedral_indices"
                    ] = key
                else:
                    continue

        dihedrals_candidates = []

        def fetch_periodicity_value(periodicity):
            if isinstance(periodicity, list):
                return max([abs(p) for p in periodicity])
            else:
                return periodicity

        for k, v in rotatable_dihedrals_dict.items():
            periodicity_value = fetch_periodicity_value(v["periodicity"])
            if periodicity_value == 3:
                dih_angle = [60, 180, 300]
            elif periodicity_value == 2:
                dih_angle = [0, 180]
            elif periodicity_value == 4:
                dih_angle = [0, 90, 180, 270]
            else:
                continue

            dih_index = v["dihedral_indices"]
            if self._check_methyl_in_dihedrals(dih_index, dihedrals_dict, atom_info_dict):
                continue

            dihedrals_candidates.append((dih_index, dih_angle))

        return dihedrals_candidates, atom_info_dict, dihedrals_dict

    def _get_dihedral_combinations(self, dihedrals_candidates):

        # assemble all possible combinations of dihedrals, each line pick a value from a list [180,0,300] for example

        iter_stime = time.time()

        test = [i[1] for i in dihedrals_candidates]
        dihedrals_combinations = list(itertools.product(*test))
        dihedral_list = [i[0] for i in dihedrals_candidates]

        print(f"Combination of dihedrals took {time.time() - iter_stime:.2f} sec")

        return dihedrals_combinations, dihedral_list

    def _get_mol_comb(self, molecule, top_file_name, dihedrals_candidates, atom_info_dict, dihedrals_dict):

        dihedrals_combinations, dihedral_list = self._get_dihedral_combinations(dihedrals_candidates)

        # assemble the dihedral and the angle to a dict
        conformation_dih_dict = []
        for i in range(len(dihedrals_combinations)):
            combo = np.array(dihedrals_combinations[i]).reshape(-1, 1)
            conformation_dih_dict.append(np.hstack((dihedral_list, combo)))

        # make an array to store the dihedral conformation_dih_dict for broadcast
        dih_comb_array = np.array(conformation_dih_dict)

        # should be aware that dihedral_dict count atom index from 0, but molecule to set dihedral count from 1
        return dih_comb_array

    def _init_openmm_system(self, topology_file):

        top = GromacsTopFile(topology_file)
        system = top.createSystem(NoCutoff)

        # platform settings for small molecule
        platform = Platform.getPlatformByName("CPU")
        platform.setPropertyDefaultValue("Threads", "1")

        integrator = LangevinIntegrator(300 * kelvin, 1 / picosecond, 0.001 * picoseconds)
        simulation = Simulation(top.topology, system, integrator, platform)

        return simulation

    def _minimize_energy(self, molecule, simulation, em_tolerance):

        coords_nm = molecule.get_coordinates_in_angstrom() * 0.1
        simulation.context.setPositions(coords_nm * nanometer)

        simulation.minimizeEnergy(tolerance=em_tolerance)

        state = simulation.context.getState(
            getPositions=True,
            getEnergy=True,
            getForces=True,
        )

        energy = state.getPotentialEnergy().value_in_unit_system(md_unit_system)
        # convert to angstrom
        optimized_coords = state.getPositions(asNumpy=True).value_in_unit_system(md_unit_system) * 10

        return energy, optimized_coords

    def _preoptimize_molecule(self, molecule, top_file_name, em_tolerance):

        mmff_gen = MMForceFieldGenerator(self._comm)
        mmff_gen.ostream.mute()
        # TODO: double check partial charge
        mmff_gen.partial_charges = molecule.get_partial_charges(molecule.get_charge())
        mmff_gen.create_topology(molecule)
        mmff_gen.write_gromacs_files(filename=self.top_file_name)
        self._comm.barrier()

        if self._rank == mpi_master():
            simulation = self._init_openmm_system(top_file_name)
            energy, opt_coords = self._minimize_energy(molecule, simulation, em_tolerance)
        else:
            opt_coords = None
        opt_coords = self._comm.bcast(opt_coords, root=mpi_master())

        new_molecule = Molecule(molecule)
        for i in range(len(opt_coords)):
            new_molecule.set_atom_coordinates(i, opt_coords[i] / bohr_in_angstrom())

        return new_molecule

    def generate(self, molecule, filename=None):

        if filename is None:
            self.top_file_name = "MOL.top"
        else:
            self.top_file_name = f"{filename}.top"

        self.molecule = molecule

        molecule = self._preoptimize_molecule(self.molecule, self.top_file_name, self.em_tolerance)

        comm = self._comm
        rank = self._comm.Get_rank()
        size = self._comm.Get_size()

        dihedrals_candidates, atom_info_dict, dihedrals_dict = (
            self._get_dihedral_candidates(molecule, self.top_file_name))

        if rank == mpi_master():
            conformation_dih_arr = self._get_mol_comb(
                molecule, self.top_file_name, dihedrals_candidates,
                atom_info_dict, dihedrals_dict)
        else:
            conformation_dih_arr = None
        conformation_dih_arr = comm.bcast(conformation_dih_arr, root=mpi_master())

        dih_comb_arr_rank = conformation_dih_arr[rank::size, :, :].copy()

        # generate conformers based on the dihedral combinations and based on
        # previous conformer to save time

        # each rank will generate the conformers based on the assigned dihedral
        # combinations

        conf_start_time = time.time()

        conformations = []

        for i in range(dih_comb_arr_rank.shape[0]):
            if i > 0:
                new_molecule = Molecule(conformations[-1])
                old_dih_settings = dih_comb_arr_rank[i - 1][:, 4]
                new_dih_settings = dih_comb_arr_rank[i][:, 4]
                # compare the difference between the two sets and only update the dihedrals that are different
                diff_dih_ind = np.where(old_dih_settings != new_dih_settings)[0]
            else:
                new_molecule = Molecule(molecule)
                diff_dih_ind = np.arange(dih_comb_arr_rank[i].shape[0])

            value_atom_index = dih_comb_arr_rank[i, :, 0:4] + 1

            for j in diff_dih_ind[::-1]:
                new_molecule.set_dihedral_in_degrees(value_atom_index[j], dih_comb_arr_rank[i, j, 4])

            conformations.append(Molecule(new_molecule))

        print(f"Rank {rank} generated {len(dih_comb_arr_rank)} conformations in {time.time() - conf_start_time:.2f} sec")

        # optimize energy and coordinates for each conformation

        opt_start_time = time.time()

        simulation = self._init_openmm_system(self.top_file_name)

        energy_coords = []

        for mol_i in range(len(conformations)):
            energy, opt_coords = self._minimize_energy(
                conformations[mol_i], simulation, self.em_tolerance)
            energy_coords.append([energy, opt_coords])

        print(f"Rank {rank} finished minimization of {len(conformations)} conformations in {time.time() - opt_start_time:.2f} sec")

        # sort and select energy_coords
        sorted_energy_coords = sorted(energy_coords)[:self.number_of_conformers_to_select]

        # gather energy and opt_coords
        gathered_energy_coords = comm.gather(sorted_energy_coords, root=mpi_master())

        if rank == mpi_master():
            # now we have all optimized conformer and energy, so we can analyze
            # the energy and get the lowest energy conformer reshape the
            # all_energy_coords to a list
            all_sel_energy_coords = [
                ene_coord for local_energy_coords in gathered_energy_coords
                for ene_coord in local_energy_coords]

            # sort and select all_energy_coords
            all_sorted_energy_coords = sorted(all_sel_energy_coords)[:self.number_of_conformers_to_select]

            # get the lowest energy conformer
            global_minimum_energy, global_minimum_conformer  = all_sorted_energy_coords[0]
            print(f"Global minimum energy is {global_minimum_energy}")

            # return conformers and coordinates
            selected_conformers = []
            for i in range(len(all_sorted_energy_coords)):
                conformer = {}
                conformer["energy"] = all_sorted_energy_coords[i][0]
                conformer["labels"] = molecule.get_labels()
                conformer["coordinates"] = all_sorted_energy_coords[i][1]
                selected_conformers.append(conformer)

            # save the selected conformers to file
            if self.save_xyz_files:
                if self.save_path is None:
                    save_path = Path("selected_conformers")
                else:
                    save_path = Path(self.save_path)
                save_path.mkdir(parents=True, exist_ok=True)

                for i, conf in enumerate(selected_conformers):
                    xyz_path = str(save_path / f"conformer_{i}.xyz")
                    # assemble coordinates and atom labels to xyz file
                    self._write_molecule_xyz_file(
                        xyz_path,
                        conf["labels"],
                        conf["coordinates"],
                        comment="Energy:" + str(conf["energy"]),
                    )

            self.global_minimum_conformer = global_minimum_conformer
            self.global_minimum_energy = global_minimum_energy
            self.selected_conformers = selected_conformers

    def show_global_minimum(self, atom_indices=False, atom_labels=False):

        if self._rank == mpi_master():
            min_conformer_xyz_string = self._convert_molecule_xyz_string(
                self.molecule.get_labels(),
                self.global_minimum_conformer,
                comment="Energy:" + str(self.global_minimum_energy),
            )
            min_conformer = Molecule.read_xyz_string(min_conformer_xyz_string)
            print("Global minimum energy:", self.global_minimum_energy)
            min_conformer.show(atom_indices=atom_indices, atom_labels=atom_labels)

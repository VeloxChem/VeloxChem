from .molecule import Molecule
from .mmforcefieldgenerator import MMForceFieldGenerator
from pathlib import Path
import re
import numpy as np
from mpi4py import MPI
import time
import itertools
from openmm.app import GromacsTopFile, NoCutoff, Simulation
from openmm import LangevinIntegrator
from openmm import Platform  # Platform added
from openmm.unit import nanometer, md_unit_system, kelvin, picoseconds, picosecond

class ConformerGenerator:
    def __init__(self,comm=None):
        if comm is None:
            comm = MPI.COMM_WORLD
            self._comm = comm
            self._rank = comm.Get_rank()
            self._size = comm.Get_size()
            self.molecule = None
            self.number_of_conformers_to_select = None
            self.top_file_name = None
            self.save_xyz = None
            self.em_tolerance_value = None
    

    # use atom label and coordinates to create a xyz string
    def _convert_molecule_xyz_string(self,labels, coords, comment=""):
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

    def _write_molecule_xyz_file(self,path,labels, coords, comment=""):
        xyz_string = self._convert_molecule_xyz_string(labels, coords, comment)
        with open(path, "w") as f:
            f.write(xyz_string)

    # use the ring identifier in the molecule set dihedral to filter out the dihedrals that are not in the ring and add periodicity information
    def _check_methyl_in_dihedrals(self,dihedral_indices, dihedrals_dict, atom_info_dict):
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


    def _rotate_check(self,molecule, zero_based_dih_indices):
        i = zero_based_dih_indices[1]
        j = zero_based_dih_indices[2]
        connectivity_matrix = molecule.get_connectivity_matrix()
        connectivity_matrix[i, j] = 0
        connectivity_matrix[j, i] = 0
        atoms_connected_to_j = molecule._find_connected_atoms(j, connectivity_matrix)

        if i in atoms_connected_to_j:
            return False
        return True

    def _append_peridocity_to_real_rotatablebonds(self,molecule, top_file_name="MOL"):
        mmff_gen = MMForceFieldGenerator()  # ADD MM
        mmff_gen.partial_charges = molecule.get_partial_charges(molecule.get_charge())
        mmff_gen.create_topology(molecule)  ##REMOVE resp=False
        mmff_gen.write_gromacs_files(filename=top_file_name)
        atom_info_dict = mmff_gen.atom_info_dict
        rotatable_bonds = mmff_gen.rotatable_bonds
        dihedrals_dict = mmff_gen.dihedrals
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
            # check if periodicity is a list or a single value
            if isinstance(periodicity, list):
                periodicity = periodicity
            else:
                periodicity = [periodicity]

            if len(periodicity) == 1:
                return abs(int(periodicity[0]))
            else:
                # return max and abs value
                max_value = max([abs(int(i)) for i in periodicity])
                return max_value

        for k, v in rotatable_dihedrals_dict.items():
            if v["periodicity"] is None:
                continue
            periodicity_value = fetch_periodicity_value(v["periodicity"])
            if periodicity_value == 3:
                dih_angle = [60, 180, 300]
                if not self._rotate_check(molecule, v["dihedral_indices"]):
                    continue

            elif periodicity_value == 2:
                dih_angle = [0, 180]
                if not self._rotate_check(molecule, v["dihedral_indices"]):
                    continue
            elif periodicity_value == 4:
                dih_angle = [0, 90, 180, 270]
                if not self._rotate_check(molecule, v["dihedral_indices"]):
                    continue
            else:
                continue

            dih_index = v["dihedral_indices"]
            if self._check_methyl_in_dihedrals(dih_index, dihedrals_dict, atom_info_dict):
                continue
            dihedrals_candidates.append((dih_index, dih_angle))

        return dihedrals_candidates, atom_info_dict, dihedrals_dict


    # assemble all possible combinations of dihedrals, each line pick a value from a list [180,0,300] for example
    def _get_dihedral_combinations(self,dihedrals_candidates):
        iter_stime = time.time()
        test = [i[1] for i in dihedrals_candidates]
        dihedrals_combinations = list(itertools.product(*test))
        dihedral_list = [i[0] for i in dihedrals_candidates]
        print(f"iteration for combinations cost time: {time.time() - iter_stime} seconds")
        return dihedrals_combinations, dihedral_list


    def _get_mol_comb(self,molecule, top_file_name): # this should in rank 0. and broadcast the combination to all ranks
        dihedrals_candidates,atom_info_dict, dihedrals_dict = (
            self._append_peridocity_to_real_rotatablebonds(
                molecule, top_file_name=top_file_name
            )
        )
        dihedrals_combinations, dihedral_list = self._get_dihedral_combinations(
            dihedrals_candidates
        )
        #assemble the dihedral and the angle to a dict
        conformation_dih_dict = {}
        if len(dihedrals_combinations) == 1:
            dihedral_list.append(dihedrals_combinations[0])
            dihedral_combinations = np.array(dihedral_list).reshape(-1, 1)
            conformation_dih_dict[0] = dihedral_combinations
        else:
            for i in range(len(dihedrals_combinations)):
                combo = np.array(dihedrals_combinations[i]).reshape(-1, 1)
                dihedral_combinations = np.hstack([dihedral_list, combo])
                conformation_dih_dict[i] = dihedral_combinations
        #make an array to store the dihedral conformation_dih_dict for broadcast
        dih_comb_array = np.array([conformation_dih_dict[i] for i in conformation_dih_dict.keys()])
        # should be aware that dihedral_dict count atom index from 0, but molecule to set dihedral count from 1
        return dih_comb_array

    def _init_openmm_system(self,topology_file):
        top = GromacsTopFile(topology_file)
        system = top.createSystem(NoCutoff)
        platform = Platform.getPlatformByName(
            "CPU"
        )  # CPU is faster than OpenCL for small systems
        platform.setPropertyDefaultValue(
            "Threads", "1"
        )  # set the number of threads to 1
        integrator = LangevinIntegrator(
            300 * kelvin, 1 / picosecond, 0.001 * picoseconds
        )
        simulation = Simulation(
            top.topology, system, integrator, platform
        )  # platform added
        return simulation

    def _get_molecule_energy_coords(self,molecule, simulation,em_tolerance_value):
        coords_nm = molecule.get_coordinates_in_angstrom() * 0.1
        simulation.context.setPositions(coords_nm * nanometer)
        simulation.minimizeEnergy(tolerance = em_tolerance_value) 
        state = simulation.context.getState(
            getPositions=True,
            getEnergy=True,
            getForces=True,
        )

        energy = state.getPotentialEnergy().value_in_unit_system(md_unit_system)
        optimized_coords = (
            state.getPositions(asNumpy=True).value_in_unit_system(md_unit_system) * 10
        )  # convert to angstrom
        return energy, optimized_coords


        
    def generate(self):#gen(self,molecule, top_file_name, number_of_conformers_to_select, save_xyz):
        if self.molecule is None:
            raise ValueError("Molecule is not set")
        else:
            molecule = self.molecule
        if self.number_of_conformers_to_select is None:
            number_of_conformers_to_select = 10
        else:
            number_of_conformers_to_select = self.number_of_conformers_to_select
        if self.top_file_name is None:
            top_file_name = "MOL.top"
        else:
            top_file_name = self.top_file_name
        if self.save_xyz is None:
            save_xyz = True
        else:
            save_xyz = self.save_xyz

        if self.em_tolerance_value is None:
            em_tolerance_value = 10
        else:
            if isinstance(self.em_tolerance_value, int):
                em_tolerance_value = self.em_tolerance_value
            elif isinstance(self.em_tolerance_value, float):
                em_tolerance_value = self.em_tolerance_value
            else:
                raise ValueError("em_tolerance_value should be an integer or float")

        comm = self._comm
        rank = self._comm.Get_rank()
        size = self._comm.Get_size()
        if rank == 0:
            conformation_dih_arr=self._get_mol_comb(molecule, top_file_name)
            # broadcast the arr to all ranks
            # decide the index to all ranks based on the size
            avg,res = divmod(conformation_dih_arr.shape[0],size)
            count = [avg + (1 if i < res else 0) for i in range(size)] # remove the 0 because split will start from 0
            displs = [sum(count[:i]) for i in range(1,size)] 
            #generate indices for each rank
            index = np.arange(conformation_dih_arr.shape[0])
            assigned_indices = np.split(index,displs)

            # for each rank, send out the subarray of dih and topfilename based on the assigned index
            for i in range(1,size):
                if len(assigned_indices[i]) == 0:
                    arr_shape = (0,0,5) # empty array
                    comm.send(arr_shape, dest=i, tag=i) 
                    dih_comb_arr_rank = np.empty((0,0,5))
                    comm.send(dih_comb_arr_rank, dest=i, tag=i) 
                    continue
                dih_comb_arr_rank = conformation_dih_arr[assigned_indices[i]]
                arr_shape = dih_comb_arr_rank.shape
                comm.send(arr_shape, dest=i, tag=i) 
                comm.send(dih_comb_arr_rank, dest=i, tag=i) 
                #print(f"rank 0 sent array shape {arr_shape} to rank {i}")\
                
            # make the dih_comb_arr_rank for rank 0
            dih_comb_arr_rank = conformation_dih_arr[assigned_indices[0]]
            #print(f"rank 0 get array shape {dih_comb_arr_rank.shape} to rank 0")

        else:
            arr_shape = None
            dih_comb_arr_rank = np.empty((0,0,0))
            arr_shape = comm.recv(source=0, tag=rank)
            dih_comb_arr_rank = np.empty(arr_shape)
            dih_comb_arr_rank = comm.recv(source=0, tag=rank)
            #print(f"rank {rank} received array shape {arr_shape}")
            #print(f"rank {rank} received array {dih_comb_arr_rank,dih_comb_arr_rank.shape}")
        

        # will start to generate conformers based on the dihedral combinations and based on previous conformer to save time
        conformations = [] # each rank will generate the conformers based on the assigned dihedral combinations
        start_time = time.time()
        if dih_comb_arr_rank.shape[0] == 0:
            print(f"rank {rank} received empty array, skipping conformer generation")
            energy_coords = []
            selected_energy_coords = []
            #MPI.COMM_WORLD.Barrier()  # Synchronize all ranks
            #exit()  # Exit the process for this rank
        else:
            for i in range(dih_comb_arr_rank.shape[0]):

                if i > 1:
                    old_dih_settings = dih_comb_arr_rank[i - 1][:, 4]
                    new_dih_settings = dih_comb_arr_rank[i][:, 4]
                    # compare the difference between the two sets and only update the dihedrals that are different
                    diff_dih_ind = np.where(old_dih_settings != new_dih_settings)[0]
                else:
                    diff_dih_ind = np.arange(dih_comb_arr_rank[i].shape[0])
                value_atom_index = dih_comb_arr_rank[i][:, 0:4] + 1

                # only loop through the dihedrals that are different from the previous set
                for j in diff_dih_ind:
                    molecule.set_dihedral_in_degrees(value_atom_index[j], dih_comb_arr_rank[i][j,4])

                conformations.append(molecule)
            print(f"rank {rank} generated {len(conformations)} conformations, with time {time.time() - start_time} seconds")

        # each rank will start to minimize the conformations and get the energy and opt_coords
        #initialize the openmm simulation
        simulation = self._init_openmm_system(top_file_name)
        # get the energy and opt_coords for each conformation
        energy_coords = []
        min_stime = time.time()
        for mol_i in range(len(conformations)):
            energy, opt_coords = self._get_molecule_energy_coords(conformations[mol_i], simulation,em_tolerance_value)
            energy_coords.append([energy, opt_coords])
        print(f"rank {rank} finished minimization for {len(conformations)} conformation with time: {time.time() - min_stime} seconds")

        # sort the energy_coords based on the energy 
        sorted_energy_coords = sorted(energy_coords, key=lambda x: x[0])
        if number_of_conformers_to_select < len(sorted_energy_coords):
            selected_energy_coords = sorted_energy_coords[:number_of_conformers_to_select]
        else:
            selected_energy_coords = sorted_energy_coords

        # gather all the energy and opt_coords of all 
        all_selected_energy_coords = comm.gather(selected_energy_coords, root=0)
        if rank == 0:
            # now we have all optimized conformer and energy, so we can analyze the energy and get the lowest energy conformer
            #reshape the all_energy_coords to a list
            all_sel_energy_coords =[]
            for i in range(size):
                if len(all_selected_energy_coords[i]) ==0:
                    print(f"rank {i} has no conformer")
                    continue
                print(f"rank {i} has {len(all_selected_energy_coords[i])} conformers")
                all_sel_energy_coords.extend(all_selected_energy_coords[i])
            #sort the all_energy_coords based on the energy again 
            all_sorted_energy_coords = sorted(all_sel_energy_coords, key=lambda x: x[0])
            #truncate the conformers to the selected number
            if number_of_conformers_to_select < len(all_sorted_energy_coords):
                all_sorted_energy_coords = all_sorted_energy_coords[:number_of_conformers_to_select]
            else:
                all_sorted_energy_coords = all_sorted_energy_coords
            #get the lowest energy conformer
            global_minimum_energy = all_sorted_energy_coords[0][0]
            global_minimum_conformer = all_sorted_energy_coords[0][1]
            print(f"Global minimum energy is {global_minimum_energy}")
            #save them all
            if save_xyz:
                molecule_label = molecule.get_labels() # get the atom labels for all 
                for i in range(len(all_sorted_energy_coords)):
                    xyz_path = str(Path("filtered", "/conformers" + str(i) + ".xyz"))
                    # assemble coordinates and atom labels to xyz file
                    self._write_molecule_xyz_file(
                        xyz_path,
                        molecule_label,
                        all_sorted_energy_coords[i][1],
                        comment="Energy:" + str(all_sorted_energy_coords[i][0]),
                    )
                
                self.global_minimum_conformer = global_minimum_conformer
                self.global_minimum_energy = global_minimum_energy
            
    def show_global_minimum(self,atom_indices=False, atom_labels=False, energy_print=True):
        rank = self._rank
        if rank == 0:
            min_conformer_xyz_string = self._convert_molecule_xyz_string(
                self.molecule.get_labels(),
                self.global_minimum_conformer,
                comment="Energy:" + str(self.global_minimum_energy),
            )
            min_conformer = Molecule.read_xyz_string(min_conformer_xyz_string)
            if energy_print:
                print("Global minimum energy:", self.global_minimum_energy)
            min_conformer.show(atom_indices=atom_indices, atom_labels=atom_labels)
            

        

if __name__ == "__main__":
    conf = ConformerGenerator()
    conf.molecule = Molecule.read_smiles("CC1([C@@H](N2[C@H](S1)[C@@H](C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C")
    conf.number_of_conformers_to_select = 1
    conf.top_file_name = "MOL.top" #default top file name
    conf.save_xyz = True #if True then save the xyz file in filtered folder
    conf.em_tolerance_value = 1 #default is 10
    conf.generate()
    conf.show_global_minimum()

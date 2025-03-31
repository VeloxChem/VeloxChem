from mpi4py import MPI
from pathlib import Path
import numpy as np
import itertools
import time
import sys
import re
from copy import deepcopy

from .veloxchemlib import bohr_in_angstrom, mpi_master
from .outputstream import OutputStream
from .molecule import Molecule
from .atomtypeidentifier import AtomTypeIdentifier
from .mmforcefieldgenerator import MMForceFieldGenerator
from .errorhandler import assert_msg_critical

try:
    from openmm import LangevinIntegrator, Platform
    from openmm.app import GromacsTopFile, NoCutoff, Simulation
    from openmm.unit import (nanometer, md_unit_system, kelvin, picoseconds,
                             picosecond)
except ImportError:
    pass


class ConformerGenerator:

    def __init__(self, comm=None, ostream=None):

        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(sys.stdout)
            else:
                ostream = OutputStream(None)

        self._comm = comm
        self._rank = comm.Get_rank()
        self._size = comm.Get_size()

        self.ostream = ostream

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
                + "   "
                + "%10.5f" % round(coords[i][0], 5)  # format the coordinates to 3 decimal places
                + "   "
                + "%10.5f" % round(coords[i][1], 5)
                + "   "
                + "%10.5f" % round(coords[i][2], 5)
                + "\n"
            )
        return xyz_string



    def _write_molecule_xyz_file(self, path, labels, coords, comment=""):

        xyz_string = self._convert_molecule_xyz_string(labels, coords, comment)
        with open(path, "w") as f:
            f.write(xyz_string)


    def analyze_equiv(self,molecule):
        def get_equiv(mol):
            idtf = AtomTypeIdentifier()
            idtf.ostream.mute()
            atom_type = idtf.generate_gaff_atomtypes(mol)
            idtf.identify_equivalences()
            equivalent_charges = idtf.equivalent_charges
            return atom_type, equivalent_charges

        atom_type, equiv = get_equiv(molecule)
        equiv_atoms_groups = []
        # Split the string by comma
        substrings = equiv.split(",")
        # Split each substring by "="
        for substr in substrings:
            unit = substr.split("=")
            equiv_atoms_groups.append(unit)
        # map str to int
        one_based_equiv_atoms_groups = [list(map(int, x)) for x in equiv_atoms_groups]
        return  atom_type, one_based_equiv_atoms_groups
        
    # use equivalent atoms
    def _check_equivside_in_dihedrals(self, dihedral_indices,atom_info_dict,one_based_equiv_atoms_groups):
        
        # according to the connected_atom and equivalent list, if one side is connected to a equivalent group "like methyl" then no need to sample and return False

        #i,j,k,l we check side_j and side_k
        # we check twice because any side works
        # round1
        side_j_index = dihedral_indices[1] + 1  # convert to 1 based index
        one_based_connected_atom_numbers = atom_info_dict[side_j_index]["ConnectedAtomsNumbers"]
        connected_set = set(one_based_connected_atom_numbers)
        connected_set = connected_set - {dihedral_indices[2]+1}  # remove the dihedral atom_k
        for equiv_g in one_based_equiv_atoms_groups:
            if connected_set.issubset(set(equiv_g)):
                return True
        # round2
        side_k_index = dihedral_indices[2] + 1  # convert to 1 based index
        one_based_connected_atom_numbers = atom_info_dict[side_k_index]["ConnectedAtomsNumbers"]
        connected_set = set(one_based_connected_atom_numbers)
        connected_set = connected_set - {dihedral_indices[1]+1}  # remove the dihedral atom_j
        for equiv_g in one_based_equiv_atoms_groups:
            if connected_set.issubset(set(equiv_g)):
                return True
        # if not found, return False
        return False
    

    def _useMMFF_generator(self, comm, mol, top_file_name="MOL"):
        # use MMFF generator to generate the topology and dihedrals
        mmff_gen = MMForceFieldGenerator(comm=comm)
        mmff_gen.partial_charges = mol.get_partial_charges(mol.get_charge())
        mmff_gen.create_topology(mol)
        mmff_gen.write_gromacs_files(filename=top_file_name)
        atom_info_dict = mmff_gen.atom_info_dict
        rotatable_bonds = mmff_gen.rotatable_bonds
        dihedrals_dict = mmff_gen.dihedrals
        return atom_info_dict, rotatable_bonds, dihedrals_dict

    def _get_dihedral_candidates(self, molecule, top_file_name,atom_info_dict, rotatable_bonds, dihedrals_dict):
        _,one_based_equiv_atoms_groups = self.analyze_equiv(molecule)
        self._comm.barrier()

        rotatable_bonds_zero_based = [(i - 1, j - 1) for (i, j) in rotatable_bonds]
        rotatable_dihedrals_dict = {}

        def get_max_periodicity(periodicity):
            if isinstance(periodicity, list):
                return max([abs(p) for p in periodicity])
            else:
                return periodicity

        # only pick one dihedral for each rotatable bond # dihedral with max periodicity absolute value
        for (i, j, k, l), dih in dihedrals_dict.items():

            sorted_bond = tuple(sorted([j, k]))
            max_periodicity = get_max_periodicity(dih["periodicity"])

            if sorted_bond in rotatable_bonds_zero_based:
                if sorted_bond not in rotatable_dihedrals_dict:
                    rotatable_dihedrals_dict[sorted_bond] = deepcopy(dih)
                    rotatable_dihedrals_dict[sorted_bond]["dihedral_indices"] = (i, j, k, l)
                    rotatable_dihedrals_dict[sorted_bond]["max_periodicity"] = max_periodicity
                else:
                    curr_periodicity = rotatable_dihedrals_dict[sorted_bond]["max_periodicity"]
                    rotatable_dihedrals_dict[sorted_bond]["max_periodicity"] = max(
                        curr_periodicity, max_periodicity)

        dihedrals_candidates = []

        for k, v in rotatable_dihedrals_dict.items():
            max_periodicity = v["max_periodicity"]
            if max_periodicity == 3:
                dih_angle = [60, 180, 300]
            elif max_periodicity == 2:
                dih_angle = [0, 180]
            elif max_periodicity == 4:
                dih_angle = [0, 90, 180, 270]
            else:
                continue

            dih_index = v["dihedral_indices"]
            if self._check_equivside_in_dihedrals(dih_index, atom_info_dict, one_based_equiv_atoms_groups):
                continue

            dihedrals_candidates.append((dih_index, dih_angle))

        return dihedrals_candidates

    def _get_dihedral_combinations(self, dihedrals_candidates):

        # assemble all possible combinations of dihedrals
        test = [i[1] for i in dihedrals_candidates]
        dihedrals_combinations = list(itertools.product(*test))
        dihedral_list = [i[0] for i in dihedrals_candidates]

        self.ostream.print_info(f"{len(dihedrals_combinations)} conformers will be generated.")
        self.ostream.flush()

        return dihedrals_combinations, dihedral_list
    

    
    def _get_mol_comb(self, dihedrals_candidates):

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

        assert_msg_critical('openmm' in sys.modules,
                            'OpenMM is required for ConformerGenerator.')

        top = GromacsTopFile(topology_file)
        system = top.createSystem(NoCutoff)

        # platform settings for small molecule
        platform = Platform.getPlatformByName("CPU")
        platform.setPropertyDefaultValue("Threads", "1")

        integrator = LangevinIntegrator(300 * kelvin, 1 / picosecond, 0.001 * picoseconds)
        simulation = Simulation(top.topology, system, integrator, platform)

        return simulation

    def _minimize_energy(self, molecule, simulation, em_tolerance):

        assert_msg_critical('openmm' in sys.modules,
                            'OpenMM is required for ConformerGenerator.')

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
    
    def _preoptimize_molecule(self, mol, top_file_name, em_tolerance_value):
        # use MMFF generator to generate the topology and dihedrals in all ranks
        _, _, _ = self._useMMFF_generator(self._comm, mol, top_file_name)
        self._comm.barrier()

        if self._rank == mpi_master():
            simulation = self._init_openmm_system(top_file_name)
            energy, opt_coords = self._minimize_energy(mol, simulation, em_tolerance_value)
        else:
            energy = None
            opt_coords = None
        energy = self._comm.bcast(energy, root=mpi_master())
        opt_coords = self._comm.bcast(opt_coords, root=mpi_master())
        # update the coordinates of the molecule to new_mol
        new_mol = Molecule(mol)
        for i in range(len(opt_coords)):
            new_mol.set_atom_coordinates(i, opt_coords[i] / bohr_in_angstrom())

        return energy,new_mol
    


    def generate(self, molecule):

        conf_gen_t0 = time.time()

        top_file_name = self.top_file_name
        if top_file_name is None:
            top_file_name = "MOL.top"

        self.molecule = molecule

        pre_molecule_energy,pre_molecule = self._preoptimize_molecule(self.molecule, top_file_name, self.em_tolerance)

        comm = self._comm
        rank = self._comm.Get_rank()
        size = self._comm.Get_size()


        atom_info_dict, rotatable_bonds, dihedrals_dict = self._useMMFF_generator(
            comm, pre_molecule, top_file_name
        )
        
        if rank == mpi_master():
            dihedrals_candidates = self._get_dihedral_candidates(pre_molecule, top_file_name,atom_info_dict, rotatable_bonds, dihedrals_dict)
            num_dih_candidates = len(dihedrals_candidates)

        else:
            num_dih_candidates = None

        num_dih_candidates = comm.bcast(num_dih_candidates, root=mpi_master())
        if num_dih_candidates == 0:
            if rank == mpi_master():
                self.ostream.print_info("No rotatable bonds found, no new conformers will be generated.")
                self.ostream.flush()
                pre_molecule_conformer = {}
                pre_molecule_conformer["energy"] = pre_molecule_energy
                pre_molecule_conformer["labels"] = pre_molecule.get_labels()
                pre_molecule_conformer["coordinates"] = pre_molecule.get_coordinates_in_angstrom()
                self.global_minimum_conformer = pre_molecule_conformer["coordinates"]
                self.global_minimum_energy = pre_molecule_energy
                self.selected_conformers = [pre_molecule_conformer]
            return


        if rank == mpi_master():
            conformation_dih_arr = self._get_mol_comb(
                dihedrals_candidates)
        else:
            conformation_dih_arr = None
        conformation_dih_arr = comm.bcast(conformation_dih_arr, root=mpi_master())

        ave, rem = divmod(len(conformation_dih_arr), size)
        counts = [ave + 1 if p < rem else ave for p in range(size)]
        displs = [sum(counts[:p]) for p in range(size)]

        num_total_conformers = len(conformation_dih_arr)

        dih_comb_arr_rank = conformation_dih_arr[
            displs[rank]:displs[rank] + counts[rank]].copy()

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

            for j in diff_dih_ind:
                new_molecule.set_dihedral_in_degrees(value_atom_index[j], dih_comb_arr_rank[i, j, 4])

            conformations.append(Molecule(new_molecule))

        conf_dt = time.time() - conf_start_time

        info = f"{num_total_conformers} conformers generated in {conf_dt:.2f} sec"
        if comm.Get_size() > 1:
            dt_list = comm.gather(conf_dt, root=mpi_master())
            if comm.Get_rank() == mpi_master():
                load_imb = 1.0 - sum(dt_list) / (len(dt_list) * max(dt_list))
                info += f' (load imb.: {load_imb * 100:.1f}%)'
        self.ostream.print_info(info)
        self.ostream.flush()

        # optimize energy and coordinates for each conformation

        opt_start_time = time.time()

        simulation = self._init_openmm_system(top_file_name)

        energy_coords = []

        for mol_i in range(len(conformations)):
            energy, opt_coords = self._minimize_energy(
                conformations[mol_i], simulation, self.em_tolerance)
            energy_coords.append([energy, opt_coords])

        opt_dt = time.time() - opt_start_time

        info = f"Energy minimization of {num_total_conformers} conformers took {opt_dt:.2f} sec"
        if comm.Get_size() > 1:
            dt_list = comm.gather(opt_dt, root=mpi_master())
            if comm.Get_rank() == mpi_master():
                load_imb = 1.0 - sum(dt_list) / (len(dt_list) * max(dt_list))
                info += f' (load imb.: {load_imb * 100:.1f}%)'
        self.ostream.print_info(info)
        self.ostream.flush()

        # sort and select energy_coords
        sorted_energy_coords = sorted(
            energy_coords, key=lambda x: x[0])[:self.number_of_conformers_to_select]

        # gather energy and opt_coords
        gathered_energy_coords = comm.gather(sorted_energy_coords, root=mpi_master())

        if self.save_xyz_files:
            if self.save_path is None:
                save_path = Path("selected_conformers")
            else:
                save_path = Path(self.save_path)
            save_path.mkdir(parents=True, exist_ok=True)

        if rank == mpi_master():
            # now we have all optimized conformer and energy, so we can analyze
            # the energy and get the lowest energy conformer reshape the
            # all_energy_coords to a list
            all_sel_energy_coords = [
                ene_coord for local_energy_coords in gathered_energy_coords
                for ene_coord in local_energy_coords]

            # sort and select all_energy_coords
            all_sorted_energy_coords = sorted(
                all_sel_energy_coords, key=lambda x: x[0])[:self.number_of_conformers_to_select]

            # get the lowest energy conformer
            global_minimum_energy, global_minimum_conformer = all_sorted_energy_coords[0]
            self.ostream.print_info(f"Global minimum energy: {global_minimum_energy:.3f} kJ/mol")

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
                for i, conf in enumerate(selected_conformers):
                    xyz_path = str(save_path / f"conformer_{i + 1}.xyz")
                    # assemble coordinates and atom labels to xyz file
                    self._write_molecule_xyz_file(
                        xyz_path,
                        conf["labels"],
                        conf["coordinates"],
                        comment=f"Energy: {conf['energy']:.3f} kJ/mol"
                    )

            self.global_minimum_conformer = global_minimum_conformer
            self.global_minimum_energy = global_minimum_energy
            self.selected_conformers = selected_conformers

            if self.save_xyz_files:
                self.ostream.print_info(
                    f"{self.number_of_conformers_to_select} conformers with " +
                    f"the lowest energies are saved in folder {str(save_path)}")

            self.ostream.print_info(
                "Total time spent in generating conformers: " +
                f"{time.time() - conf_gen_t0:.2f} sec")

    def show_global_minimum(self, atom_indices=False, atom_labels=False):

        if self._rank == mpi_master():
            min_conformer_xyz_string = self._convert_molecule_xyz_string(
                self.molecule.get_labels(),
                self.global_minimum_conformer,
                comment="Energy:" + str(self.global_minimum_energy),
            )
            min_conformer = Molecule.read_xyz_string(min_conformer_xyz_string)
            min_conformer.show(atom_indices=atom_indices, atom_labels=atom_labels)

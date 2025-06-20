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
from pathlib import Path
from copy import deepcopy
import numpy as np
import itertools
import time
import sys

from .veloxchemlib import bohr_in_angstrom, mpi_master
from .outputstream import OutputStream
from .molecule import Molecule
from .atomtypeidentifier import AtomTypeIdentifier
from .mmforcefieldgenerator import MMForceFieldGenerator
from .errorhandler import assert_msg_critical
from .mofutils import svd_superimpose
from .molecularbasis import MolecularBasis
from .respchargesdriver import RespChargesDriver

try:
    import openmm
    from openmm import LangevinIntegrator, Platform
    from openmm.app import NoCutoff, Simulation, PDBFile, ForceField, GromacsTopFile
    from openmm.unit import nanometer, md_unit_system, kelvin, picoseconds, picosecond
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
        self.number_of_conformers_to_select = None

        self.top_file_name = None
        self.partial_charges = None
        self.resp_charges = True
        self.resp_charges_driver = RespChargesDriver()

        self.save_xyz_files = False
        self.save_path = None

        # thresholds for energy minimization
        self.em_tolerance = 1.0

        # thresholds for removing duplicate conformers
        # rmsd threshold in Angstrom, energy threshold in kJ/mol
        # TODO: double check the thresholds
        self.rmsd_threshold = 1.2
        self.energy_threshold = 1.2

        self.implicit_solvent_model = None
        self.solute_dielectric = 1.0
        self.solvent_dielectric = 78.39

        self.use_gromacs_files = False

    def _analyze_equiv(self, molecule):

        idtf = AtomTypeIdentifier()
        idtf.ostream.mute()

        atom_type = idtf.generate_gaff_atomtypes(molecule)

        idtf.identify_equivalences()
        equivalent_charges = idtf.equivalent_charges
        if len(equivalent_charges) == 0:
            return []

        equiv_atoms_groups = []
        for substr in equivalent_charges.split(","):
            equiv_atoms_groups.append(substr.split("="))

        one_based_equiv_atoms_groups = [
            list(map(int, x)) for x in equiv_atoms_groups
        ]

        return one_based_equiv_atoms_groups

    def _check_equivside_in_dihedrals(self, dihedral_indices, atom_info_dict,
                                      one_based_equiv_atoms_groups):

        # according to the connected_atom and equivalent list, if one side is
        # connected to a equivalent group "like methyl" then no need to sample

        # for i,j,k,l we check side_j and side_k
        # we check twice because any side works

        side_j_index = dihedral_indices[1] + 1  # convert to 1 based index
        side_k_index = dihedral_indices[2] + 1  # convert to 1 based index

        side_equiv = False
        max_equiv_atoms = 0

        # check j and k side of dihedral
        for a, b in [(side_j_index, side_k_index),
                     (side_k_index, side_j_index)]:

            # check equiv_atoms on one dihedral atom
            if atom_info_dict[a]["AtomicSymbol"] == "C":
                one_based_connected_atom_numbers = atom_info_dict[a][
                    "ConnectedAtomsNumbers"]
                connected_set = set(one_based_connected_atom_numbers)
                connected_set = connected_set - {
                    b
                }  # remove the other dihedral atom

                for equiv_g in one_based_equiv_atoms_groups:
                    if connected_set.issubset(set(equiv_g)):
                        side_equiv = True
                        max_equiv_atoms = max(max_equiv_atoms,
                                              len(connected_set))
                        break

        return side_equiv, max_equiv_atoms

    def _check_methyl_group(self, dihedral_indices, atom_info_dict):

        side_j_index = dihedral_indices[1] + 1  # convert to 1 based index
        side_k_index = dihedral_indices[2] + 1  # convert to 1 based index

        # check j and k side of dihedral
        for a, b in [(side_j_index, side_k_index),
                     (side_k_index, side_j_index)]:

            if atom_info_dict[a]["AtomicSymbol"] == "C":
                one_based_connected_atom_numbers = atom_info_dict[a][
                    "ConnectedAtomsNumbers"]
                connected_set = set(one_based_connected_atom_numbers)
                connected_set = connected_set - {
                    b
                }  # remove the other dihedral atom
                connected_elements = [
                    atom_info_dict[idx]["AtomicSymbol"]
                    for idx in connected_set
                ]
                if tuple(connected_elements) == ("H", "H", "H"):
                    return True

        return False

    def _get_dihedral_candidates(self, molecule, top_file_name,
                                 partial_charges):

        mmff_gen = MMForceFieldGenerator(self._comm)
        mmff_gen.ostream.mute()
        if partial_charges is None:
            warn_text = "ConformerGenerator: Partial charges not provided. "
            warn_text += "Will use a quick (and likely inaccurate) estimation of partial charges."
            self.ostream.print_warning(warn_text)
            mmff_gen.partial_charges = molecule.get_partial_charges(
                molecule.get_charge())
        else:
            mmff_gen.partial_charges = partial_charges
        mmff_gen.create_topology(molecule)
        if self.use_gromacs_files:
            mmff_gen.write_gromacs_files(filename=top_file_name)
        else:
            mmff_gen.write_openmm_files(filename=top_file_name)
        # make sure to sync the ranks after the top file is written
        self._comm.barrier()

        atom_info_dict = deepcopy(mmff_gen.atom_info_dict)
        rotatable_bonds = deepcopy(mmff_gen.rotatable_bonds)
        dihedrals_dict = deepcopy(mmff_gen.dihedrals)

        rotatable_bonds_zero_based = [(i - 1, j - 1)
                                      for (i, j) in rotatable_bonds]
        rotatable_dihedrals_dict = {}

        def get_max_periodicity(periodicity):
            if isinstance(periodicity, list):
                return max([abs(p) for p in periodicity])
            else:
                return periodicity

        # only pick one dihedral for each rotatable bond
        for (i, j, k, l), dih in dihedrals_dict.items():

            sorted_bond = tuple(sorted([j, k]))
            max_periodicity = get_max_periodicity(dih["periodicity"])

            if sorted_bond in rotatable_bonds_zero_based:
                if sorted_bond not in rotatable_dihedrals_dict:
                    rotatable_dihedrals_dict[sorted_bond] = deepcopy(dih)
                    rotatable_dihedrals_dict[sorted_bond][
                        "dihedral_indices"] = (i, j, k, l)
                    rotatable_dihedrals_dict[sorted_bond][
                        "max_periodicity"] = max_periodicity
                else:
                    curr_periodicity = rotatable_dihedrals_dict[sorted_bond][
                        "max_periodicity"]
                    rotatable_dihedrals_dict[sorted_bond][
                        "max_periodicity"] = max(curr_periodicity,
                                                 max_periodicity)

        dihedrals_candidates = []

        # needed by _check_equivside_in_dihedrals
        one_based_equiv_atoms_groups = self._analyze_equiv(molecule)

        for k, v in rotatable_dihedrals_dict.items():
            max_periodicity = v["max_periodicity"]
            if max_periodicity == 2:
                dih_angle = [0, 180]
            elif max_periodicity == 3:
                dih_angle = [60, 180, 300]
            elif max_periodicity == 4:
                dih_angle = [0, 90, 180, 270]
            elif max_periodicity == 5:
                dih_angle = [36, 54, 84, 144, 324]
            else:
                continue

            dih_index = v["dihedral_indices"]

            # skip dihedral angle involving methyl group
            if self._check_methyl_group(dih_index, atom_info_dict):
                continue

            # look for equiv_atoms, and skip if number of equiv_atoms equals max_periodicity
            side_equiv, max_equiv_atoms = self._check_equivside_in_dihedrals(
                dih_index, atom_info_dict, one_based_equiv_atoms_groups)
            if side_equiv and max_equiv_atoms == max_periodicity:
                continue

            dihedrals_candidates.append((dih_index, dih_angle))

        return dihedrals_candidates, atom_info_dict, dihedrals_dict

    def _get_dihedral_combinations(self, dihedrals_candidates):

        # assemble all possible combinations of dihedrals

        dih_angles = [i[1] for i in dihedrals_candidates]
        dihedrals_combinations = list(itertools.product(*dih_angles))
        dihedral_list = [i[0] for i in dihedrals_candidates]

        self.ostream.print_info(
            f"{len(dihedrals_combinations)} conformers will be generated.")
        self.ostream.flush()

        return dihedrals_combinations, dihedral_list

    def _get_mol_comb(self, molecule, top_file_name, dihedrals_candidates):

        dihedrals_combinations, dihedral_list = self._get_dihedral_combinations(
            dihedrals_candidates)

        # assemble the dihedral and the angle to a dict
        conformation_dih_dict = []
        for i in range(len(dihedrals_combinations)):
            combo = np.array(dihedrals_combinations[i]).reshape(-1, 1)
            conformation_dih_dict.append(np.hstack((dihedral_list, combo)))

        # make an array to store the dihedral conformation_dih_dict for broadcast
        dih_comb_array = np.array(conformation_dih_dict)

        # should be aware that dihedral_dict count atom index from 0, but molecule to set dihedral count from 1
        return dih_comb_array

    def _init_openmm_system(self, topology_file, implicit_solvent_model):

        assert_msg_critical('openmm' in sys.modules,
                            'OpenMM is required for ConformerGenerator.')

        assert_msg_critical(
            not (self.use_gromacs_files
                 and implicit_solvent_model is not None),
            'ConformerGenerator: Please set to use_gromacs_files to False ' +
            'when using implicit solvent model.')

        if self.use_gromacs_files:
            if topology_file.endswith(".top"):
                topology_file = str(Path(topology_file))
            else:
                topology_file = str(Path(topology_file).with_suffix(".top"))
            top = GromacsTopFile(topology_file)
            system = top.createSystem(NoCutoff)
        else:
            pdb_file = str(Path(topology_file).with_suffix(".pdb"))
            xml_file = str(Path(topology_file).with_suffix(".xml"))
            pdb = PDBFile(pdb_file)
            if implicit_solvent_model is None:
                forcefield = ForceField(xml_file)
                system = forcefield.createSystem(pdb.topology,
                                                 nonbondedMethod=NoCutoff)
            else:
                implicit_fpath = Path(
                    openmm.__file__).parent / "app" / "data" / "implicit"
                implicit_model_fname = str(implicit_fpath /
                                           f"{implicit_solvent_model}.xml")
                assert_msg_critical(
                    Path(implicit_model_fname).is_file(),
                    f"ConformerGenerator: Could not find file {implicit_model_fname}"
                )
                forcefield = ForceField(xml_file, implicit_model_fname)
                system = forcefield.createSystem(
                    pdb.topology,
                    nonbondedMethod=NoCutoff,
                    soluteDielectric=self.solute_dielectric,
                    solventDielectric=self.solvent_dielectric)
                self.ostream.print_info(
                    f"Using implicit solvent model {implicit_solvent_model}")
                self.ostream.flush()

        # platform settings for small molecule
        platform = Platform.getPlatformByName("CPU")
        platform.setPropertyDefaultValue("Threads", "1")

        integrator = LangevinIntegrator(300 * kelvin, 1 / picosecond,
                                        0.001 * picoseconds)
        if self.use_gromacs_files:
            simulation = Simulation(top.topology, system, integrator, platform)
        else:
            simulation = Simulation(pdb.topology, system, integrator, platform)

        return simulation

    def show_available_implicit_solvent_models(self):

        if self._rank == mpi_master():
            implicit_folder_path = Path(
                openmm.__file__).parent / "app" / "data" / "implicit"
            implicit_solvent_files = [
                f for f in implicit_folder_path.iterdir()
                if f.is_file() and f.suffix == ".xml"
            ]
            print("Available implicit solvent files:")
            for f in implicit_solvent_files:
                print(f.name)

    def _minimize_energy(self, molecule, simulation, em_tolerance):

        assert_msg_critical('openmm' in sys.modules,
                            'OpenMM is required for ConformerGenerator.')

        coords_nm = molecule.get_coordinates_in_angstrom() * 0.1
        simulation.context.setPositions(coords_nm * nanometer)

        simulation.minimizeEnergy(tolerance=em_tolerance)

        state = simulation.context.getState(
            getPositions=True,
            getEnergy=True,
            # getForces=True,
        )

        energy = state.getPotentialEnergy().value_in_unit_system(
            md_unit_system)
        # convert to angstrom
        optimized_coords = state.getPositions(
            asNumpy=True).value_in_unit_system(md_unit_system) * 10

        return energy, optimized_coords

    def _optimize_molecule(self, molecule, top_file_name, em_tolerance,
                           partial_charges, implicit_solvent_model):

        mmff_gen = MMForceFieldGenerator(self._comm)
        mmff_gen.ostream.mute()
        if partial_charges is None:
            warn_text = "ConformerGenerator: Partial charges not provided. "
            warn_text += "Will use a quick (and likely inaccurate) estimation of partial charges."
            self.ostream.print_warning(warn_text)
            mmff_gen.partial_charges = molecule.get_partial_charges(
                molecule.get_charge())
        else:
            mmff_gen.partial_charges = partial_charges
        mmff_gen.create_topology(molecule)
        if self.use_gromacs_files:
            mmff_gen.write_gromacs_files(filename=top_file_name)
        else:
            mmff_gen.write_openmm_files(filename=top_file_name)
        # make sure to sync the ranks after the top file is written
        self._comm.barrier()

        if self._rank == mpi_master():
            simulation = self._init_openmm_system(top_file_name,
                                                  implicit_solvent_model)
            energy, opt_coords = self._minimize_energy(molecule, simulation,
                                                       em_tolerance)
        else:
            energy, opt_coords = None, None
        energy, opt_coords = self._comm.bcast((energy, opt_coords),
                                              root=mpi_master())

        new_molecule = Molecule(molecule)
        for i in range(len(opt_coords)):
            new_molecule.set_atom_coordinates(
                i, opt_coords[i] / bohr_in_angstrom())

        return energy, new_molecule

    def generate(self, molecule):

        # sanity check
        if self.implicit_solvent_model is not None:
            self.use_gromacs_files = False

        conf_gen_t0 = time.time()

        top_file_name = self.top_file_name
        if top_file_name is None:
            top_file_name = "MOL"

        self.molecule = molecule

        comm = self._comm
        rank = self._comm.Get_rank()
        size = self._comm.Get_size()
        #default partial charges are RESP charges: self.resp_charges is True and partial_charges is None as default 
        #if the user provides partial charges, then skip RESP charges calculation
        if self.resp_charges and (self.partial_charges is None):
            basis = MolecularBasis.read(molecule, "6-31g*")
            self.partial_charges =self.resp_charges_driver.compute(
                molecule, basis, 'resp')

        dihedrals_candidates, atom_info_dict, dihedrals_dict = (
            self._get_dihedral_candidates(molecule, top_file_name,
                                          self.partial_charges))

        # exit early if there is no candidate dihedral to rotate
        if not dihedrals_candidates:
            self.ostream.print_info(
                "No rotatable bond found, no new conformer will be generated.")
            self.ostream.flush()

            energy, molecule = self._optimize_molecule(
                self.molecule, top_file_name, self.em_tolerance,
                self.partial_charges, self.implicit_solvent_model)

            if rank == mpi_master():
                self.global_minimum_conformer = molecule
                self.global_minimum_energy = energy
                return {
                    'energies': [energy],
                    'molecules': [Molecule(molecule)],
                    'geometries': [
                        molecule.get_xyz_string(
                            comment=f"Energy: {energy:.3f} kJ/mol")
                    ],
                }
            else:
                return None

        if rank == mpi_master():
            conformation_dih_arr = self._get_mol_comb(molecule, top_file_name,
                                                      dihedrals_candidates)
        else:
            conformation_dih_arr = None
        conformation_dih_arr = comm.bcast(conformation_dih_arr,
                                          root=mpi_master())

        ave, rem = divmod(len(conformation_dih_arr), size)
        counts = [ave + 1 if p < rem else ave for p in range(size)]
        displs = [sum(counts[:p]) for p in range(size)]

        num_total_conformers = len(conformation_dih_arr)

        dih_comb_arr_rank = conformation_dih_arr[displs[rank]:displs[rank] +
                                                 counts[rank]].copy()

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
                diff_dih_ind = np.where(
                    old_dih_settings != new_dih_settings)[0]
            else:
                new_molecule = Molecule(molecule)
                diff_dih_ind = np.arange(dih_comb_arr_rank[i].shape[0])

            value_atom_index = dih_comb_arr_rank[i, :, 0:4] + 1

            for j in diff_dih_ind:
                new_molecule.set_dihedral_in_degrees(
                    value_atom_index[j], dih_comb_arr_rank[i, j, 4])

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

        simulation = self._init_openmm_system(top_file_name,
                                              self.implicit_solvent_model)

        energy_coords = []

        for mol_i in range(len(conformations)):
            energy, opt_coords = self._minimize_energy(conformations[mol_i],
                                                       simulation,
                                                       self.em_tolerance)
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
        if self.number_of_conformers_to_select is None:
            self.number_of_conformers_to_select = num_total_conformers
        sorted_energy_coords = sorted(
            energy_coords,
            key=lambda x: x[0])[:self.number_of_conformers_to_select]

        # gather energy and opt_coords
        gathered_energy_coords = comm.gather(sorted_energy_coords,
                                             root=mpi_master())

        if rank == mpi_master():
            # now we have all optimized conformer and energy, so we can analyze
            # the energy and get the lowest energy conformer reshape the
            # all_energy_coords to a list
            all_sel_energy_coords = [
                ene_coord for local_energy_coords in gathered_energy_coords
                for ene_coord in local_energy_coords
            ]

            # sort all_energy_coords
            all_sorted_energy_coords = sorted(all_sel_energy_coords,
                                              key=lambda x: x[0])
            # get the lowest energy conformer
            min_energy, min_coords_angstrom = all_sorted_energy_coords[0]
            min_mol = Molecule(molecule)
            for iatom in range(min_mol.number_of_atoms()):
                min_mol.set_atom_coordinates(
                    iatom, min_coords_angstrom[iatom] / bohr_in_angstrom())
            self.ostream.print_info(
                f"Global minimum energy: {min_energy:.3f} kJ/mol")
            self.ostream.flush()

            self.global_minimum_conformer = min_mol
            self.global_minimum_energy = min_energy

            # return conformers info
            conformers_dict = {
                'energies': [],
                'molecules': [],
                'geometries': [],
            }

            for conf_energy, conf_coords_angstrom in all_sorted_energy_coords:
                conformers_dict["energies"].append(conf_energy)
                mol_copy = Molecule(molecule)
                for iatom in range(mol_copy.number_of_atoms()):
                    mol_copy.set_atom_coordinates(
                        iatom,
                        conf_coords_angstrom[iatom] / bohr_in_angstrom())
                conformers_dict["molecules"].append(mol_copy)
                conformers_dict["geometries"].append(
                    mol_copy.get_xyz_string(
                        comment=f"Energy: {conf_energy:.3f} kJ/mol"))

            num_conformers = len(conformers_dict["energies"])

            equiv_conformer_pairs = []

            for i in range(num_conformers):
                xyz_i = conformers_dict["molecules"][
                    i].get_coordinates_in_angstrom()
                ene_i = conformers_dict["energies"][i]

                for j in range(i + 1, num_conformers):
                    xyz_j = conformers_dict["molecules"][
                        j].get_coordinates_in_angstrom()
                    ene_j = conformers_dict["energies"][j]

                    if abs(ene_i - ene_j) < self.energy_threshold:
                        rmsd, rot, trans = svd_superimpose(xyz_j, xyz_i)
                        if rmsd < self.rmsd_threshold:
                            equiv_conformer_pairs.append((i, j))

            duplicate_conformers = [j for i, j in equiv_conformer_pairs]
            duplicate_conformers = sorted(list(set(duplicate_conformers)))

            filtered_energies = [
                e for i, e in enumerate(conformers_dict["energies"])
                if i not in duplicate_conformers
            ]
            filtered_molecules = [
                Molecule(m) for i, m in enumerate(conformers_dict["molecules"])
                if i not in duplicate_conformers
            ]
            filtered_geometries = [
                g for i, g in enumerate(conformers_dict["geometries"])
                if i not in duplicate_conformers
            ]

            conformers_dict = {
                "energies":
                filtered_energies[:self.number_of_conformers_to_select],
                "molecules":
                filtered_molecules[:self.number_of_conformers_to_select],
                "geometries":
                filtered_geometries[:self.number_of_conformers_to_select],
            }

            # save the selected conformers to file
            if self.save_xyz_files:
                if self.save_path is None:
                    save_path = Path("selected_conformers")
                else:
                    save_path = Path(self.save_path)
                save_path.mkdir(parents=True, exist_ok=True)

                for i in range(len(conformers_dict["energies"])):
                    conf_energy = conformers_dict["energies"][i]
                    conf_mol = conformers_dict["molecules"][i]
                    xyz_path = save_path / f"conformer_{i + 1}.xyz"
                    with xyz_path.open("w") as fh:
                        fh.write(
                            conf_mol.get_xyz_string(
                                comment=f"Energy: {conf_energy:.3f} kJ/mol"))

                self.ostream.print_info(
                    f"{len(conformers_dict['energies'])} conformers with " +
                    f"the lowest energies are saved in folder {str(save_path)}"
                )
            else:
                self.ostream.print_info(
                    f"{len(conformers_dict['energies'])} conformers remain " +
                    "after removal of duplicate conformers.")

            self.ostream.print_info(
                "Total time spent in generating conformers: " +
                f"{time.time() - conf_gen_t0:.2f} sec")
            self.ostream.flush()

            self.conformer_dict = conformers_dict

            return conformers_dict
        else:
            return None

    def show_global_minimum(self, atom_indices=False, atom_labels=False):

        if self._rank == mpi_master():
            print(
                f"Global minimum conformer with energy {self.global_minimum_energy:.3f} kJ/mol"
            )
            self.global_minimum_conformer.show(atom_indices=atom_indices,
                                               atom_labels=atom_labels)

    def show_conformers(self, number=1, atom_indices=False, atom_labels=False):

        if self._rank == mpi_master():
            if number > len(self.conformer_dict["energies"]):
                number = len(self.conformer_dict["energies"])
                print(f"Only {number} conformers available, showing all.")
            for i in range(number):
                print(
                    f"Conformer {i + 1} with energy {self.conformer_dict['energies'][i]:.3f} kJ/mol"
                )
                self.conformer_dict["molecules"][i].show(
                    atom_indices=atom_indices, atom_labels=atom_labels)

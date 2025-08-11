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
import numpy as np
import time
import json
import h5py
import sys
import copy

from .veloxchemlib import mpi_master
from .molecule import Molecule
from .outputstream import OutputStream
from .mmforcefieldgenerator import MMForceFieldGenerator
from .evbsystembuilder import EvbSystemBuilder
from .evbfepdriver import EvbFepDriver
from .evbffbuilder import EvbForceFieldBuilder
from .evbdataprocessing import EvbDataProcessing
from .evbsystembuilder import EvbForceGroup
from .solvationbuilder import SolvationBuilder
from .errorhandler import assert_msg_critical
from .sanitychecks import molecule_sanity_check

try:
    import openmm as mm
    import openmm.app as mmapp
    import openmm.unit as mmunit
except ImportError:
    pass


class EvbDriver():

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

        # mpi information
        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()

        self.temperature: float = 300
        self.Lambda: list[float] = None

        self.reactant: MMForceFieldGenerator = None
        self.product: MMForceFieldGenerator = None

        self.name: str = None
        self.results = None
        self.system_confs: list[dict] = []
        self.mute_scf = True

        self.t_label = int(time.time())
        self.water_model = 'spce'

    def build_and_run_default_water_EVB(
        self,
        reactant: Molecule | list[Molecule],
        product: Molecule | list[Molecule],
        barrier,
        free_energy,
    ):
        self.ostream.print_blank()
        self.ostream.print_header("Building forcefields")
        self.ostream.flush()

        self.build_ff_from_molecules(
            reactant,
            product,
            ordered_input=True,
            optimize_mol=True,
        )
        self.ostream.print_blank()
        self.ostream.print_header("Building systems")
        self.ostream.flush()
        self.build_systems(configurations=["vacuum", "water"])

        self.ostream.print_blank()
        self.ostream.print_header("Running FEP")
        self.ostream.flush()
        self.run_FEP()

        self.ostream.print_blank()
        self.ostream.print_header("Computing energy profiles")
        self.ostream.flush()
        self.compute_energy_profiles(barrier, free_energy)

        self.ostream.flush()

    def build_ff_from_molecules(
            self,
            reactant: Molecule | list[Molecule],
            product: Molecule | list[Molecule],
            reactant_partial_charges: list[float] | list[list[float]] = None,
            product_partial_charges: list[float] | list[list[float]] = None,
            reactant_total_multiplicity: int = -1,
            product_total_multiplicity: int = -1,
            reparameterize: bool = True,
            optimize_mol: bool = False,
            optimize_ff: bool = True,
            mm_opt_constrain_bonds: bool = True,
            breaking_bonds: set[tuple[int, int]] | tuple = set(),
            reactant_hessians: np.ndarray | list[np.ndarray|None] | None = None,
            product_hessians: list[np.ndarray|None] | None = None,
            mute_scf: bool = True,
    ):
        
        ffbuilder = EvbForceFieldBuilder(ostream=self.ostream)
        ffbuilder.reactant_partial_charges = reactant_partial_charges
        ffbuilder.product_partial_charges = product_partial_charges
        ffbuilder.reactant_total_multiplicity = reactant_total_multiplicity
        ffbuilder.product_total_multiplicity = product_total_multiplicity
        ffbuilder.reparameterize = reparameterize
        ffbuilder.optimize_ff = optimize_ff
        ffbuilder.optimize_mol = optimize_mol
        ffbuilder.mm_opt_constrain_bonds = mm_opt_constrain_bonds
        ffbuilder.breaking_bonds = breaking_bonds
        ffbuilder.reactant_hessians = reactant_hessians
        ffbuilder.product_hessians = product_hessians
        ffbuilder.mute_scf = mute_scf

        ffbuilder.water_model = self.water_model

        self.reactant, self.product, self.forming_bonds, self.breaking_bonds, self.reactants, self.products, self.product_mapping = ffbuilder.build_forcefields(
            reactant=reactant,
            product=product,
        )

    def build_systems(
        self,
        configurations: list[str] | list[dict],  # type: ignore
        Lambda: list[float] | np.ndarray = None,
        constraints: dict | list[dict] | None = None,
    ):
        """Build OpenMM systems for the given configurations with interpolated forcefields for each lambda value. Saves the systems as xml files, the topology as a pdb file and the options as a json file to the disk.

        Args:
            configurations (list[str] | list[dict]): The given configurations for which to perform an FEP. The first configuration will be regarded as the reference configuration. 
            Lambda (list[float] | np.ndarray): The Lambda vector to be used for the FEP. Should start with 0, end with 1 and be monotonically increasing. 
                Defaults to None, in which case default values will be assigned depending on if debugging is enabled or not.
                If a string is given, the return value of default_system_configurations() will be used. See this function for default configurations.
            constraints (dict | list[dict] | None, optional): Dictionary of harmonic bond, angle or (improper) torsion forces to apply over in every FEP frame. Defaults to None.
        """

        assert_msg_critical('openmm' in sys.modules,
                            'openmm is required for EvbDriver.')

        if all(isinstance(conf, str) for conf in configurations):
            configurations = [
                self.default_system_configurations(conf)
                for conf in configurations
            ]

        assert all(
            isinstance(conf, dict) for conf in configurations
        ), "Configurations must be a list of strings or a list of dictionaries"
        self.configurations: list[dict] = configurations  # type: ignore
        if constraints is None:
            constraints = []
        if isinstance(constraints, dict):
            constraints = [constraints]

        if Lambda is None:
            if configurations[0].get("debug", False):
                Lambda = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
            else:
                Lambda = np.linspace(0, 0.1, 6)
                Lambda = np.append(Lambda[:-1], np.linspace(0.1, 0.9, 21))
                Lambda = np.append(Lambda[:-1], np.linspace(0.9, 1, 6))
                Lambda = np.round(Lambda, 3)
                self.ostream.print_info(
                    f"Using default lambda vector: {list(Lambda)}")
        assert (Lambda[0] == 0 and Lambda[-1]
                == 1), f"Lambda must start at 0 and end at 1. Lambda = {Lambda}"
        assert np.all(
            np.diff(Lambda) >
            0), f"Lambda must be monotonically increasing. Lambda = {Lambda}"
        Lambda = [round(lam, 3) for lam in Lambda]
        self.Lambda = Lambda

        #Per configuration
        for conf in self.configurations:
            #create folders,
            data_folder = f"EVB_{self.name}_{conf['name']}_data_{self.t_label}"
            while Path(data_folder).exists():
                self.t_label += 1
                data_folder = f"EVB_{self.name}_{conf['name']}_data_{self.t_label}"

            run_folder = str(Path(data_folder) / "run")
            conf["data_folder"] = data_folder
            conf["run_folder"] = run_folder

            cwd = Path().cwd()
            data_folder_path = cwd / data_folder
            run_folder_path = cwd / run_folder

            data_folder_path.mkdir(parents=True)
            run_folder_path.mkdir(parents=True)

            #
            self.reactant.molecule.write_xyz_file(
                str(data_folder_path / "reactant_struct.xyz"))
            self.product.molecule.write_xyz_file(
                str(data_folder_path / "product_struct.xyz"))

            MMForceFieldGenerator.save_forcefield(
                self.reactant, str(data_folder_path / f"reactant_ff_data.json"))
            MMForceFieldGenerator.save_forcefield(
                self.product, str(data_folder_path / f"product_ff_data.json"))

            if conf.get('solvent', None) is None and conf.get('pressure',
                                                              -1) > 0:
                self.ostream.print_warning(
                    f"A pressure is defined for {conf['name']}, but no solvent is defined. Removing pressure definition."
                )
                conf.pop("pressure")
            # build the system
            system_builder = EvbSystemBuilder(ostream=self.ostream)
            system_builder.water_model = self.water_model
            self.ostream.print_blank()
            self.ostream.print_header(f"Building systems for {conf['name']}")
            self.ostream.flush()
            systems, topology, initial_positions = system_builder.build_systems(
                reactant=self.reactant,
                product=self.product,
                Lambda=self.Lambda,
                configuration=conf,
                constraints=constraints,
            )

            conf["systems"] = systems
            conf["topology"] = topology
            conf["initial_positions"] = initial_positions
            conf['forming_bonds'] = list(self.forming_bonds)
            conf['breaking_bonds'] = list(self.breaking_bonds)

            self.ostream.print_info(f"Saving files to {data_folder_path}")
            self.ostream.flush()
            self.save_systems_as_xml(systems, conf["run_folder"])

            top_path = cwd / data_folder / "topology.pdb"

            mmapp.PDBFile.writeFile(
                topology,
                initial_positions,  # positions are handled in nanometers, but pdb's should be in angstroms
                open(top_path, "w"),
            )

            dump_conf = copy.copy(conf)
            dump_conf.pop('systems')
            dump_conf.pop('topology')
            dump_conf.pop('initial_positions')
            self.update_options_json(dump_conf, conf)
            self.update_options_json(
                {
                    "Lambda": Lambda,
                },
                conf,
            )

        self.system_confs = configurations

        self.create_viamd_environment_files()
        self.ostream.flush()

    def load_initialisation(self,
                            data_folder: str,
                            name: str,
                            load_systems=False,
                            load_pdb=False):
        """Load a configuration from a data folder for which the systems have already been generated, such that an FEP can be performed. 
        The topology, initial positions, temperature and Lambda vector will be loaded from the data folder.

        Args:
            data_folder (str): The folder to load the data from
            name (str): The name of the configuration. Can be arbitrary, but should be unique.
            load_systems (bool, optional): If set to true, the systems will be loaded from the xml files. Used for debugging. Defaults to False.
            lead_pdb (bool, optional): If set to true, the topology will be loaded from the pdb file. Used for debugging. Defaults to False.
        """

        assert_msg_critical('openmm' in sys.modules,
                            'openmm is required for EvbDriver.')

        with open(str(Path(data_folder) / "options.json"), "r") as file:
            options = json.load(file)
            temperature = options["temperature"]
            Lambda = options["Lambda"]
        if self.Lambda != Lambda and self.Lambda is not None:
            self.ostream.print_warning(
                f"Lambda vector in {data_folder}/options.json does not match the current Lambda vector. Overwriting current Lambda vector with the one from the file."
            )

        self.Lambda = Lambda

        conf = {
            "name": name,
            "data_folder": data_folder,
            "run_folder": str(Path(data_folder) / "run"),
            "temperature": temperature,
            "Lambda": Lambda
        }
        if load_systems:
            systems = self.load_systems_from_xml(str(Path(data_folder) / "run"))
            conf["systems"] = systems
        else:
            systems = []

        if load_pdb:
            pdb = mmapp.PDBFile(str(Path(data_folder) / "topology.pdb"))
            conf["topology"] = pdb.getTopology()
            conf["initial_positions"] = pdb.getPositions(
                asNumpy=True).value_in_unit(mmunit.nanometers)

        self.system_confs.append(conf)
        self.ostream.print_info(
            f"Initialised configuration with {len(systems)} systems, topology, initial positions, temperatue {temperature} and Lambda vector {Lambda} from {data_folder}"
        )
        self.ostream.print_info(
            f"Current configurations: {[conf['name'] for conf in self.system_confs]}"
        )
        self.ostream.flush()

    def save_systems_as_xml(self, systems: dict, folder: str):
        """Save the systems as xml files to the given folder.

        Args:
            systems (dict): The systems to save
            folder (str): The folder relative to the current working directory to save the systems to.
        """

        assert_msg_critical('openmm' in sys.modules,
                            'openmm is required for EvbDriver.')

        path = Path().cwd() / folder
        self.ostream.print_info(f"Saving systems to {path}")
        self.ostream.flush()
        for name, system in systems.items():
            if isinstance(name, float) or isinstance(name, int):
                filename = f"{name:.3f}_sys.xml"
            else:
                filename = f"{name}_sys.xml"
            with open(path / filename, mode="w", encoding="utf-8") as output:
                output.write(mm.XmlSerializer.serialize(system))

    def load_systems_from_xml(self, folder: str):
        """Load the systems from xml files in the given folder.

        Args:
            folder (str): The folder relative to the current working directory to load the systems from.
        Returns:
            dict: The loaded systems
        """

        assert_msg_critical('openmm' in sys.modules,
                            'openmm is required for EvbDriver.')

        systems = {}
        path = Path().cwd() / folder
        for lam in self.Lambda:
            with open(path / f"{lam:.3f}_sys.xml", mode="r",
                      encoding="utf-8") as input:
                systems[lam] = mm.XmlSerializer.deserialize(input.read())
        return systems

    def run_FEP(
        self,
        platform=None,
        platform_properties=None,
    ):
        """Run the the FEP calculations for all configurations in self.system_confs.

        Args:
            equil_steps (int, optional): The amount of timesteps to equilibrate at the beginning af each Lambda frame. Equilibration is done with frozen H-bonds. Defaults to 5000.
            sample_steps (int, optional): The amount of steps to sample. Defaults to 100000.
            write_step (int, optional): Per how many steps to take a sample and save its data as well as the trajectory point. Defaults to 1000.
            initial_equil_steps (int, optional): The amount of timesteps to add to the equilibration at the first Lambda frame. Defaults to 5000.
            step_size (float, optional): The step size during the sampling in picoseconds. Defaults to 0.001.
            equil_step_size (float, optional): The step size during the equilibration in picoseconds. Is typically larger then step_size as equilibration is done with frozen H-bonds. Defaults to 0.002.
            initial_equil_step_size (float, optional): The step size during initial equilibration in picoseconds. Defaults to 0.002.
        """

        for conf in self.system_confs:
            self.ostream.print_blank()
            self.ostream.print_header(f"Running FEP for {conf['name']}")
            self.ostream.flush()
            FEP = EvbFepDriver(ostream=self.ostream)
            FEP.run_FEP(
                Lambda=self.Lambda,
                configuration=conf,
                platform=platform,
                platform_properties=platform_properties,
            )

    def update_options_json(self, dict, conf):

        cwd = Path().cwd()
        path = cwd / conf["data_folder"] / "options.json"
        if not path.exists():
            with open(path, "w") as file:
                json.dump(dict, file, indent=4)
        else:
            with open(path, "r") as file:
                options = json.load(file)
            options.update(dict)
            with open(path, "w") as file:
                json.dump(options, file, indent=4)

    def compute_energy_profiles(
        self,
        barrier,
        free_energy,
        lambda_sub_sample=1,
        lambda_sub_sample_ends=False,
        time_sub_sample=1,
        dE_range=None,
        alpha=None,
        H12=None,
        alpha_guess=None,
        H12_guess=None,
    ):
        """Compute the EVB energy profiles using the FEP results, print the results and save them to an h5 file

        Args:
            barrier (float): the reaction barrier in kJ/mol of the reference system
            free_energy (float): the reaction free energy in kJ/mol of the reference system
            lambda_sub_sample (int, optional): Factor with which the lambda vector will be subsampled. Setting this to two will discard every other lambda frame. Defaults to 1.
            lambda_sub_sample_ends (bool, optional): If set to False, the lambda frames up to 0.1 and from 0.9 will not be subsampled. Defaults to False.
            time_sub_sample (int, optional): Factor with which the time vector will be subsampled. Setting this to two will discard every other snapshot. Defaults to 1.
        """
        dp = EvbDataProcessing(ostream=self.ostream)
        results = self._load_output_from_folders(lambda_sub_sample,
                                                 lambda_sub_sample_ends,
                                                 time_sub_sample)
        self.ostream.flush()

        if alpha is not None: dp.alpha = alpha
        if H12 is not None: dp.H12 = H12
        if alpha_guess is not None: dp.alpha_guess = alpha_guess
        if H12_guess is not None: dp.H12_guess = H12_guess
        if dE_range is not None:
            dp.coordinate_bins = np.linspace(dE_range[0], dE_range[1], 200)

        self.dataprocessing = dp
        results = dp.compute(results, barrier, free_energy)
        self.results = results
        self.print_results()
        # self._save_dict_as_h5(results, f"results_{self.name}")
        self.ostream.flush()
        return self.results

    def print_results(self, results: dict = None, file_name: str = None):
        """Print EVB results. Uses the provided dictionary first, then tries to load it from the disk, and last it uses the results attribute of this object.

        Args:
            results (dict, optional): A dictionary with EVB results. Defaults to None.
            file_name (str, optional): Filename of an h5 file containing EVB results. Defaults to None.
        """
        if results is None:
            if file_name is None:
                assert self.results is not None, "No results known, and none provided"

            if self.results is None:
                self.results = self._load_dict_from_h5(file_name)
            else:
                results = self.results

        dp = EvbDataProcessing()
        dp.print_results(self.results, self.ostream)
        self.ostream.flush()

    def plot_results(self,
                     results: dict = None,
                     file_name: str = None,
                     **kwargs):
        """Plot EVB results. Uses the provided dictionary first, then tries to load it from the disk, and last it uses the results attribute of this object.

        Args:
            results (dict, optional): A dictionary with EVB results. Defaults to None.
            file_name (str, optional): Filename of an h5 file containing EVB results. Defaults to None.
        """
        if results is None:
            if file_name is None:
                assert self.results is not None, "No results known, and none provided"

            if self.results is None:
                self.results = self._load_dict_from_h5(file_name)
            else:
                results = self.results
        dp = EvbDataProcessing()
        dp.plot_results(self.results, **kwargs)
        self.ostream.flush()

    def _load_output_from_folders(
        self,
        lambda_sub_sample,
        lambda_sub_sample_ends,
        time_sub_sample,
    ) -> dict:
        reference_folder = self.system_confs[0]["data_folder"]
        target_folders = [conf["data_folder"] for conf in self.system_confs[1:]]

        reference_name = self.system_confs[0]["name"]
        target_names = [conf["name"] for conf in self.system_confs[1:]]

        folders = [reference_folder] + target_folders
        results = {}
        cwd = Path().cwd()

        common_results = []
        specific_results = {}
        for name, folder in zip([reference_name] + target_names, folders):
            E_file = str(cwd / folder / "Energies.csv")
            data_file = str(cwd / folder / "Data_combined.csv")
            options_file = str(cwd / folder / "options.json")
            fg_file = str(cwd / folder / "ForceGroups.csv")
            rea_fg_file = str(cwd / folder / "ForceGroups_rea.csv")
            pro_fg_file = str(cwd / folder / "ForceGroups_pro.csv")
            decomp_file = str(cwd / folder / "NB_decompositions.csv")
            specific, common = self._load_output_files(
                E_file,
                data_file,
                options_file,
                fg_file,
                rea_fg_file,
                pro_fg_file,
                decomp_file,
                lambda_sub_sample,
                lambda_sub_sample_ends,
                time_sub_sample,
            )
            specific_results.update({name: specific})
            common_results.append(common)

        results.update({"configuration_results": specific_results})
        for common in common_results[1:]:
            for key, val in common.items():
                if isinstance(common[key], list) or isinstance(
                        common[key], np.ndarray):
                    assert np.all(
                        common[key] == common_results[0][key]
                    ), f"Common results are not the same for all configurations. Key: {key}, value: {val}"
                else:
                    assert common[key] == common_results[0][
                        key], f"Common results are not the same for all configurations. Key: {key}, value: {val}"

        for key, val in common_results[0].items():
            results.update({key: val})
        return results

    def _load_output_files(
        self,
        E_file,
        data_file,
        options_file,
        fg_file=None,
        fg_rea_file=None,
        fg_pro_file=None,
        decomp_file=None,
        lambda_sub_sample=1,
        lambda_sub_sample_ends=False,
        time_sub_sample=1,
    ):
        with open(options_file, "r") as file:
            options = json.load(file)
        Lambda = options["Lambda"]
        Temp_set = options["temperature"]
        options.pop("Lambda")

        if lambda_sub_sample > 1:
            if lambda_sub_sample_ends:
                Lambda = Lambda[::lambda_sub_sample]
            else:
                arg01 = np.where(np.array(Lambda) <= 0.1)[0][-1] + 1
                arg09 = np.where(np.array(Lambda) >= 0.9)[0][0] - 1
                Lambda = Lambda[:arg01] + Lambda[
                    arg01:arg09:lambda_sub_sample] + Lambda[arg09:]
                # Lambda_middle =
            self.ostream.print_info(
                f"Subsampling Lambda vector with factor {lambda_sub_sample}. New Lambda vector: {Lambda}"
            )

        if Lambda[-1] != 1:
            self.ostream.print_info(
                "Lambda vector does not end at 1. Appending 1 to the Lambda vector"
            )
            Lambda = np.append(Lambda, 1)

        E_data = np.loadtxt(E_file, skiprows=1, delimiter=',').T
        fg_data = []
        rea_fg_data = []
        pro_fg_data = []

        l_sub_indices = np.where([lf in Lambda for lf in E_data[0]])[0]

        sub_indices = l_sub_indices[::time_sub_sample]

        Lambda_frame = E_data[0, sub_indices]
        E1_pes = E_data[1, sub_indices]
        E2_pes = E_data[2, sub_indices]
        E1_int = E_data[3, sub_indices]
        E2_int = E_data[4, sub_indices]
        E_m_pes = E_data[5, sub_indices]
        # E_m_int = E_data[6, sub_indices]

        step, Ep, Ek, Temp, Vol, Dens = np.loadtxt(
            data_file,
            skiprows=1,
            delimiter=',',
        ).T[:, sub_indices]

        specific_result = {
            "E1_pes": E1_pes,
            "E2_pes": E2_pes,
            "E1_int": E1_int,
            "E2_int": E2_int,
            "E_m_pes": E_m_pes,
            "Ep": Ep,
            "Ek": Ek,
            "Temp_step": Temp,
            "Vol": Vol,
            "Dens": Dens,
            "options": options,
            "Temp_set": Temp_set,
        }

        if fg_file is not None and Path(fg_file).is_file():
            fg_data = np.loadtxt(fg_file, skiprows=1, delimiter=',').T
            specific_result.update({"E_m_fg": fg_data})
        if fg_rea_file is not None and Path(fg_rea_file).is_file():
            rea_fg_data = np.loadtxt(fg_rea_file, skiprows=1, delimiter=',').T
            specific_result.update({"E1_fg": rea_fg_data})
        if fg_pro_file is not None and Path(fg_pro_file).is_file():
            pro_fg_data = np.loadtxt(fg_pro_file, skiprows=1, delimiter=',').T
            specific_result.update({"E2_fg": pro_fg_data})
        if decomp_file is not None and Path(decomp_file).exists():
            decomp_data = np.loadtxt(decomp_file, skiprows=1, delimiter=',').T
            decomp_rea = decomp_data[decomp_data.shape[0] // 2:, :]
            decomp_pro = decomp_data[:decomp_data.shape[0] // 2, :]
            with open(decomp_file, "r") as file:
                decomp_names = file.readline().strip().split(",")
            decomp_names = [name.replace("_rea", "") for name in decomp_names]
            decomp_names = decomp_names[:len(decomp_names) // 2]
            specific_result.update({
                "decompositions": {
                    "E1": decomp_rea,
                    "E2": decomp_pro,
                    "names": decomp_names,
                }
            })

        lambda_indices = [
            np.where(np.round(Lambda, 3) == L)[0][0] for L in Lambda_frame
        ]
        common_result = {
            "Lambda": Lambda,
            "Lambda_frame": Lambda_frame,
            "Lambda_indices": lambda_indices,
        }
        return specific_result, common_result

    @staticmethod
    def _load_dict_from_h5(file):
        """Load a dictionary from from an h5 file

        Args:
            file (path): The file to load the results from.

        Returns:
            dict: Dictionary with the results
        """
        with h5py.File(file, "r") as f:

            def load_group(group):
                data = {}
                for k, v in group.items():
                    if isinstance(v, h5py.Group):
                        data[k] = load_group(v)
                    elif isinstance(v, h5py.Dataset):
                        data[k] = v[()]
                    else:
                        data[k] = v
                return data

            data = load_group(f)
        return data

    def _save_dict_as_h5(self, data: dict, file_name: str, overwrite=True):
        """Save the provided dictionary to an h5 file

        Args:
            results (dict): Dictionary to be saved.
        """
        cwd = Path.cwd()

        file_path = str(cwd / f"{file_name}.h5")

        with h5py.File(file_path, "w") as file:
            self.ostream.print_info(f"Saving results to {file_path}")

            def save_group(data, group):
                for k, v in data.items():
                    if isinstance(v, dict):
                        subgroup = group.create_group(k)
                        save_group(v, subgroup)
                    elif isinstance(v, np.ndarray) or isinstance(v, list):
                        group.create_dataset(k, data=v)
                    elif isinstance(v, set):
                        group.create_dataset(k, data=list(v))
                    else:
                        group[k] = v

            save_group(data, file)

    def create_viamd_environment_files(self):
        for conf in self.system_confs:
            base = ("[Files]\n"
                    "MoleculeFile=./topology.pdb\n"
                    "TrajectoryFile=./trajectory.xtc\n"
                    "CoarseGrained=0\n"
                    "\n"
                    "[RenderSettings]\n"
                    "SsaoEnabled=0\n"
                    "DofEnabled=0\n"
                    "\n"
                    "[Representation]\n"
                    "Name=Reaction\n"
                    'Filter=resname("REA")\n'
                    "Enabled=1\n"
                    "Type=2\n"
                    "ColorMapping=1\n"
                    "Saturation=1.000000\n"
                    "Param=1.000000,1.000000,1.000000,1.000000\n"
                    "DynamicEval=0\n")

            script = ("[Script]\n"
                      'Text="""\n')
            rea_script = 'rea = resname("REA");'
            sol_script = ""
            if conf.get("solvent", None) is not None:

                sol_script = (
                    f'sol = resname("SOL");\n'
                    'close_sol = (within(5, rea) and resname("SOL"));\n')
            pdb_script = ""
            if conf.get('pdb', None) is not None:
                resids = [
                    res['residue'] for res in conf.get("pdb_active_res", [])
                ]

                if len(resids) > 0:
                    s = "".join([f" or resid({id})" for id in resids])
                    rea_script = rea_script[:-1] + s + ";"
                pdb_script = "pocket = residue(protein and within(3,rea)) and not element('H');\n"
            script += rea_script + "\n"
            script += sol_script + "\n"
            script += pdb_script + "\n"

            script += '"""'

            solvent_rep = ("[Representation]\n"
                           "Name=Solvent\n"
                           "Filter=close_sol\n"
                           "Enabled=1\n"
                           "Type=1\n"
                           "ColorMapping=1\n"
                           "Saturation=1.000000\n"
                           "Param=0.354000,1.000000,1.000000,1.000000\n"
                           "DynamicEval=1\n")

            protein_rep = ("[Representation]\n"
                           "Name=Protein\n"
                           "Filter=protein\n"
                           "Enabled=1\n"
                           "Type=4\n"
                           "ColorMapping=8\n"
                           "StaticColor=1.000000,1.000000,1.000000,1.000000\n"
                           "Saturation=1.000000\n"
                           "Param=1.000000,1.000000,1.000000,1.000000\n"
                           "DynamicEval=0\n"
                           "\n"
                           "[Representation]\n"
                           "Name=pocket\n"
                           "Filter=pocket\n"
                           "Enabled=1\n"
                           "Type=0\n"
                           "ColorMapping=1\n"
                           "StaticColor=1.000000,1.000000,1.000000,1.000000\n"
                           "Saturation=0.570000\n"
                           "Param=1.000000,1.000000,1.000000,1.000000\n"
                           "DynamicEval=0\n")

            carbon_rep = ("[Representation]\n"
                          "Name=Carbon\n"
                          'Filter=resname("CCC")\n'
                          "Enabled=1\n"
                          "Type=2\n"
                          "ColorMapping=1\n"
                          "Saturation=1.000000\n"
                          "Param=1.000000,1.000000,1.000000,1.000000\n"
                          "DynamicEval=0\n")

            string = base + "\n"
            if conf.get("solvent", None) is not None:
                string += solvent_rep + "\n"
            if conf.get('pdb', None) is not None:
                string += protein_rep + "\n"
            if conf.get('CNT', False) or conf.get('graphene', False):
                string += carbon_rep + "\n"

            string += script + "\n"

            with open(f"{conf['data_folder']}/workspace.via", "w") as file:
                file.write(string)

    def default_system_configurations(self, name: str) -> dict:
        """Return a dictionary with a default configuration. Options not given in the dictionary will be set to default values in the build_systems function.

        Args:
            name (string): The name of the configuration to be used. Options are "vacuum", "water", "CNT", "graphene", "E_field", "no_reactant"
        """
        if name == "vacuum" or name == "vacuum_NVT":
            conf = {
                "name": name,
                "temperature": self.temperature,
            }
        elif name == "vacuum_NVE":
            conf = {
                "name": "vacuum_NVE",
            }
        elif name == "debug":
            conf = {
                "name": "debug",
                "debug": True,
                "temperature": self.temperature,
                "pressure": 1,
                "equil_NVT_steps": 100,
                "equil_NPT_steps": 100,
                "sample_steps": 1000,
                "write_step": 1,
                "initial_equil_NVT_steps": 0,
                "initial_equil_NPT_steps": 0,
            }
        elif name == "water" or name == "water_NPT":
            conf = {
                "name": f"water_{self.water_model}_NPT",
                "solvent": self.water_model,
                "temperature": self.temperature,
                "pressure": 1,
                "padding": 1.5,
                "ion_count": 0,
                "neutralize": False
            }
        elif name == "water_NVT":
            conf = {
                "name": f"water_{self.water_model}_NVT",
                "solvent": self.water_model,
                "temperature": self.temperature,
                "padding": 1.5,
                "ion_count": 0,
                "neutralize": False
            }
        elif name == "E_field":
            conf = {
                "name": f"water_E_field_{self.water_model}",
                "solvent": self.water_model,
                "temperature": self.temperature,
                "pressure": 1,
                "padding": 1.5,
                "ion_count": 0,
                "E_field": [0, 0, 10],
            }
        elif name == "no_reactant":
            conf = {
                "name": "no_reactant",
                "solvent": self.water_model,
                "temperature": self.temperature,
                "pressure": 1,
                "padding": 1.5,
                "ion_count": 0,
                "no_reactant": True,
            }
        elif name == "ts_guesser":
            conf = {
                "name": "vacuum",
                "temperature": self.temperature,
                "bonded_integration": True,
                "soft_core_coulomb_pes": True,
                "soft_core_lj_pes": True,
                "soft_core_coulomb_int": False,
                "soft_core_lj_int": False,
            }
        else:
            try:
                solvent = SolvationBuilder()._solvent_properties(name)
                conf = {
                    "name": name,
                    "solvent": name,
                    "temperature": self.temperature,
                    "pressure": 1,
                    "padding": 1.5,
                    "ion_count": 0,
                }
            except:
                raise ValueError(f"Unknown system configuration {name}")

        return conf

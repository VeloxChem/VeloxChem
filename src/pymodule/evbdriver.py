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

        self.input_folder: str

        self.reactant: MMForceFieldGenerator = None
        self.product: MMForceFieldGenerator = None
        self.input_folder: str = "input_files"

        self.name: str = None
        self.results = None
        self.system_confs: list[dict] = []
        self.fast_run = False

        self.t_label = int(time.time())
        self.water_model = 'spce'

    def build_and_run_default_water_EVB(
        self,
        reactant: str | Molecule,
        product: str | list[str] | Molecule | list[Molecule],
        barrier,
        free_energy,
        ordered_input=False,
    ):
        """Automatically perform an EVB calculation using a vacuum system as reference and a system solvated in water as target system.

        Args:
            reactant (str | Molecule): The reactant. If a string is given, the corresponding xyz file must be present in the input_files folder.
            product (str | list[str] | Molecule | list[Molecule]): A list of products. If a (list of) string(s) is given, the corresponding xyz file(s) must be present in the input_files folder.
            barrier (float): the reaction barrier in kJ/mol of the vacuum system
            free_energy (float): the reaction free energy in kJ/mol of the vacuum system
            ordered_input (bool, optional): If set to true, assumes that the reactant and product have the same ordering of atoms, and thus will not attempt to generate a mapping. Defaults to False.
        """
        self.ostream.print_blank()
        self.ostream.print_header("Building forcefields")
        self.ostream.flush()
        if isinstance(reactant, str) and isinstance(product, str):
            self.build_ff_from_files(
                reactant,
                product,
                ordered_input=ordered_input,
                optimize=True,
            )
        elif isinstance(reactant, Molecule) and isinstance(product, Molecule):
            self.build_ff_from_molecules(
                reactant,
                product,
                ordered_input=True,
                optimize=True,
            )
        else:
            assert_msg_critical(
                False,
                "Either both reactant and product should be strings or both should be Molecule objects.",
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
        reparameterize: bool = True,
        optimize: bool = False,
        ordered_input: bool = False,
        breaking_bonds: list[tuple[int, int]] = [],
        name=None,
    ):
        if self.name is None:
            self.name = name
        cwd = Path().cwd()
        input_path = cwd / self.input_folder
        if not input_path.exists():
            input_path.mkdir(parents=True, exist_ok=True)

        rea_input, reactant_total_charge = self._process_molecule_input(
            reactant,
            reactant_partial_charges,
        )
        pro_input, product_total_charge = self._process_molecule_input(
            product,
            product_partial_charges,
        )
        assert reactant_total_charge == product_total_charge, f"Total charge of reactants {reactant_total_charge} and products {product_total_charge} must match"

        ffbuilder = EvbForceFieldBuilder(ostream=self.ostream)
        ffbuilder.water_model = self.water_model
        ffbuilder.reparameterize = reparameterize
        ffbuilder.optimize = optimize

        if isinstance(breaking_bonds, tuple):
            breaking_bonds = [breaking_bonds]

        self.reactant, self.product, self.formed_bonds, self.broken_bonds, self.reactants, self.products = ffbuilder.build_forcefields(
            rea_input,
            pro_input,
            ordered_input,
            breaking_bonds,
        )

    @staticmethod
    def _process_molecule_input(molecules, partial_charges):
        if isinstance(molecules, Molecule):
            molecules = [molecules]

        if partial_charges is None:
            partial_charges = [None] * len(molecules)
        elif isinstance(partial_charges[0], float) or isinstance(
                partial_charges[0], int):
            partial_charges = [partial_charges]

        assert len(molecules) == len(
            partial_charges
        ), "Amount of input molecules and lists of partial charges must match"

        for i, (molecule,
                partial_charge) in enumerate(zip(molecules, partial_charges)):
            charge = molecule.get_charge()
            molecule_sanity_check(molecule)
            if partial_charge is not None:
                assert abs(
                    sum(partial_charge) - charge
                ) < 0.001, f"Sum of partial charges of reactant {sum(partial_charge)} must match the total foral charge of the system {charge} for input {i+1}"

        input = [{
            "molecule": mol,
            "optimize": None,
            "forcefield": None,
            "hessian": None,
            "charges": charge
        } for mol, charge in zip(molecules, partial_charges)]
        total_charge = sum([mol.get_charge() for mol in molecules])
        return input, total_charge

    def build_ff_from_files(
        self,
        reactant: str | list[str],
        product: str | list[str],
        reactant_charge: int | list[int] = 0,
        product_charge: int | list[int] = 0,
        reactant_multiplicity: int | list[int] = 1,
        product_multiplicity: int | list[int] = 1,
        reparameterize: bool = True,
        optimize: bool = False,
        ordered_input: bool = False,
        breaking_bonds: list[tuple[int, int]] = [],
        save_output: bool = True,
        force_recalculation: bool = False,
    ):
        if isinstance(reactant, list):
            combined_reactant_name = '_'.join(reactant)
        else:
            combined_reactant_name = reactant
        if self.name is None:
            self.name = combined_reactant_name
        combined_rea_input = self._get_input_files(combined_reactant_name)

        if isinstance(product, list):
            combined_product_name = '_'.join(product)
        else:
            combined_product_name = product
        combined_pro_input = self._get_input_files(combined_product_name)
        # combined_pro_input = self._process_file_input(
        #     combined_product_name,
        #     product_charge,product_multiplicity,
        # )[0]

        cwd = Path().cwd()
        mapped_product_path = cwd / self.input_folder / f"{combined_product_name}_mapped.xyz"
        if (combined_rea_input['forcefield'] is not None
                and combined_rea_input['molecule'] is not None
                and combined_pro_input['forcefield'] is not None
                and mapped_product_path.exists() and not force_recalculation):
            mapped_product_molecule = Molecule.read_xyz_file(
                str(mapped_product_path))
            combined_pro_input['forcefield'].molecule = mapped_product_molecule
            combined_pro_input['molecule'] = mapped_product_molecule
            self.ostream.print_info(
                f"Found both reactant and product forcefield data. Not generating new forcefields"
            )
            self.reactant = combined_rea_input["forcefield"]
            self.product = combined_pro_input["forcefield"]
            self.ostream.flush()
        else:
            rea_input = self._process_file_input(
                reactant,
                reactant_charge,
                reactant_multiplicity,
            )

            pro_input = self._process_file_input(
                product,
                product_charge,
                product_multiplicity,
            )

            ffbuilder = EvbForceFieldBuilder(ostream=self.ostream)
            ffbuilder.water_model = self.water_model
            ffbuilder.reparameterize = reparameterize
            ffbuilder.optimize = optimize

            if isinstance(breaking_bonds, tuple):
                breaking_bonds = [breaking_bonds]

            if force_recalculation:
                self.ostream.print_warning(
                    f"Forcing recalculation of forcefields, even though they might all be present"
                )
                for rea in rea_input:
                    rea['forcefield'] = None
                for pro in pro_input:
                    pro['forcefield'] = None
            self.reactant, self.product, self.formed_bonds, self.broken_bonds, self.reactants, self.products = ffbuilder.build_forcefields(
                rea_input,
                pro_input,
                ordered_input,
                breaking_bonds,
            )
            if save_output:
                self.ostream.print_info(
                    f"Saving forcefield and structure data to {self.input_folder} folder"
                )
                cwd = Path().cwd()
                reactant_ff_path = cwd / self.input_folder / f"{combined_reactant_name}_ff_data.json"
                product_ff_path = cwd / self.input_folder / f"{combined_product_name}_ff_data.json"
                self.save_forcefield(self.reactant, str(reactant_ff_path))
                self.save_forcefield(self.product, str(product_ff_path))

                reactant_mol_path = cwd / self.input_folder / f"{combined_reactant_name}.xyz"
                self.reactant.molecule.write_xyz_file(str(reactant_mol_path))
                self.product.molecule.write_xyz_file(str(mapped_product_path))

    def _process_file_input(self, filenames, charge, multiplicity):
        if isinstance(filenames, str):
            filenames = [filenames]

        # by default, give all molecules 0 charge and 1 multiplicity
        if isinstance(charge, int):
            if charge == 0:
                charge = [charge] * len(filenames)
            else:
                charge = [charge]

        if isinstance(multiplicity, int):
            if multiplicity == 1:
                multiplicity = [multiplicity] * len(filenames)
            else:
                multiplicity = [multiplicity]

        assert len(filenames) == len(
            charge), "Number of reactants and charges must match"
        assert len(filenames) == len(
            multiplicity), "Number of reactants and multiplicities must match"

        input = []
        for rea, charge, mult in zip(filenames, charge, multiplicity):
            input.append(self._get_input_files(rea))
            if input[-1]['molecule'] is not None:
                input[-1]['molecule'].set_charge(charge)
                input[-1]['molecule'].set_multiplicity(mult)
                molecule_sanity_check(input[-1]['molecule'])

            partial_charge = input[-1]['charges']
            if partial_charge is not None:
                assert abs(
                    sum(partial_charge) - charge
                ) < 0.001, f"Sum of partial charges of reactant {sum(partial_charge)} must match the total foral charge of the system {charge} for input {i+1}"

        for inp, filename in zip(input, filenames):
            assert inp['molecule'] is not None or inp[
                'forcefield'] is not None, f"Could not load {filename} file. Check if a corresponding xyz or json file exists and if the input folder is set correct"

        return input

    def _get_input_files(self, filename: str):
        # Build a molecule from a (possibly optimized) geometry
        optimize = True
        cwd = Path().cwd()

        opt_path = cwd / self.input_folder / f"{filename}_xtb_opt.xyz"
        struct_path = cwd / self.input_folder / f"{filename}.xyz"
        if opt_path.exists():
            self.ostream.print_info(
                f"Loading optimized geometry from {opt_path}")
            molecule = Molecule.read_xyz_file(str(opt_path))
            optimize = False
        elif struct_path.exists():
            self.ostream.print_info(
                f"Loading (possibly unoptimized) geometry from {struct_path}")
            molecule = Molecule.read_xyz_file(str(struct_path))
        else:
            molecule = None
            self.ostream.print_info(
                f"Could not load {'/'.join(opt_path.parts[-2:])} and {'/'.join(struct_path.parts[-2:])}"
            )
        self.ostream.flush()

        charges = None
        charge_path = cwd / self.input_folder / f"{filename}_charges.txt"
        if charge_path.exists():
            with open(charge_path, "r", encoding="utf-8") as file:
                charges = []
                for line in file:
                    try:
                        charges.append(float(line))
                    except ValueError:
                        self.ostream.print_info(
                            f"Could not read line {line} from {charge_path}. Continuing"
                        )
            print_charge = sum([round(charge, 3) for charge in charges])
            self.ostream.print_info(
                f"Loading charges from {charge_path} file, total charge: {print_charge}"
            )
        else:
            self.ostream.print_info(
                f"Could not load {'/'.join(charge_path.parts[-2:])} file. Continuing without charges"
            )
        self.ostream.flush()

        forcefield = None
        json_path = cwd / self.input_folder / f"{filename}_ff_data.json"
        if json_path.exists():
            self.ostream.print_info(
                f"Loading force field data from {json_path}")
            forcefield = self.load_forcefield_from_json(str(json_path))
            if charges is not None:
                forcefield.partial_charges = charges
            if molecule is not None:
                forcefield.molecule = molecule
        else:
            self.ostream.print_info(
                f"Could not find force field data file {self.input_folder}/{filename}_ff_data.json."
            )
        self.ostream.flush()

        hessian = None
        hessian_path = cwd / self.input_folder / f"{filename}_hess.np"
        if hessian_path.exists():
            self.ostream.print_info(
                f"Found hessian file at {hessian_path}, using it to reparameterize."
            )
            hessian = np.loadtxt(hessian_path)
        else:
            self.ostream.print_info(
                f"Could not find hessian file at {hessian_path}")
        self.ostream.flush()

        return {
            "molecule": molecule,
            "optimize": optimize,
            "forcefield": forcefield,
            "hessian": hessian,
            "charges": charges
        }

    #todo, should be moved to forcefieldgenerator class
    @staticmethod
    def load_forcefield_from_json(path: str) -> MMForceFieldGenerator:
        """
        Load forcefield data from a JSON file.

        Args:
            Path (str): The path to the JSON file containing the forcefield data.

        Returns:
            MMForceFieldGenerator: The updated forcefield object with the loaded data.
        """
        with open(path, "r", encoding="utf-8") as file:
            forcefield = MMForceFieldGenerator()
            ff_data = json.load(file)

            forcefield.atoms = EvbDriver._str_to_tuple_key(ff_data["atoms"])
            forcefield.bonds = EvbDriver._str_to_tuple_key(ff_data["bonds"])
            forcefield.angles = EvbDriver._str_to_tuple_key(ff_data["angles"])
            forcefield.dihedrals = EvbDriver._str_to_tuple_key(
                ff_data["dihedrals"])
            forcefield.impropers = EvbDriver._str_to_tuple_key(
                ff_data["impropers"])
        return forcefield

    #todo, should be moved to forcefieldgenerator class
    @staticmethod
    def save_forcefield(forcefield: MMForceFieldGenerator, path: str):
        """
        Save the forcefield data of the forcefieldgenerator to a JSON file, converting all tuples to strings

        Args:
            forcefield (MMForceFieldGenerator): The forcefield object containing the data to be saved.
            filename (str): The name of the file to save the forcefield data to.

        Returns:
            None
        """
        ff_data = {
            "atoms": forcefield.atoms,
            "bonds": EvbDriver._tuple_to_str_key(forcefield.bonds),
            "angles": EvbDriver._tuple_to_str_key(forcefield.angles),
            "dihedrals": EvbDriver._tuple_to_str_key(forcefield.dihedrals),
            "impropers": EvbDriver._tuple_to_str_key(forcefield.impropers),
        }
        with open(path, "w", encoding="utf-8") as file:
            json.dump(ff_data, file, indent=4)

    @staticmethod
    def _str_to_tuple_key(dictionary: dict) -> dict:
        """
        Converts the keys of a dictionary from string to tuple.

        Args:
            dictionary (dict): The dictionary to convert.

        Returns:
            dict: The dictionary with keys converted to tuple.
        """
        str_keys = list(dictionary.keys())
        tup_keys = []
        for str_key in str_keys:
            tuple = ()
            for item in str_key.split(","):
                item = item.replace("(", "")
                item = item.replace(")", "")
                item = item.replace(" ", "")
                tuple += (int(item), )
            if len(tuple) == 1:
                tuple = tuple[0]
            tup_keys.append(tuple)
        return {key: value for key, value in zip(tup_keys, dictionary.values())}

    @staticmethod
    def _tuple_to_str_key(dictionary: dict) -> dict:
        """
        Converts the keys of a dictionary from tuples to strings.

        Args:
            dictionary (dict): The dictionary to be converted.

        Returns:
            dict: The dictionary with string keys.

        """
        return {str(key): value for key, value in dictionary.items()}

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

        if Lambda is None:
            if configurations[0] == "debug":
                Lambda = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
            else:
                Lambda = np.linspace(0, 0.1, 11)
                Lambda = np.append(Lambda[:-1], np.linspace(0.1, 0.9, 41))
                Lambda = np.append(Lambda[:-1], np.linspace(0.9, 1, 11))
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

            self.ostream.print_info(f"Saving files to {data_folder_path}")
            self.ostream.flush()
            self.save_systems_as_xml(systems, conf["run_folder"])

            top_path = cwd / data_folder / "topology.pdb"

            mmapp.PDBFile.writeFile(
                topology,
                initial_positions *
                10,  # positions are handled in nanometers, but pdb's should be in angstroms
                open(top_path, "w"),
            )

            dump_conf = copy.copy(conf)
            dump_conf.pop('systems')
            dump_conf.pop('topology')
            dump_conf.pop('initial_positions')
            self.update_options_json(dump_conf, conf)
            self.update_options_json(
                {
                    "Lambda":
                    Lambda,
                    "integration forcegroups":
                    list(EvbForceGroup.integration_force_groups()),
                    "pes forcegroups":
                    list(EvbForceGroup.pes_force_groups()),
                },
                conf,
            )

        self.system_confs = configurations
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
        for lam in self.Lambda:
            file_path = str(path / f"{lam:.3f}_sys.xml")
            with open(file_path, mode="w", encoding="utf-8") as output:
                output.write(mm.XmlSerializer.serialize(systems[lam]))

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
        saved_frames_on_crash=None,
        platform=None,
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
            if saved_frames_on_crash is not None:
                FEP.save_frames = saved_frames_on_crash
            FEP.run_FEP(
                Lambda=self.Lambda,
                configuration=conf,
                platform=platform,
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
        self._save_dict_as_h5(results, f"results_{self.name}")
        self.results = results
        self.print_results()
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
            specific, common = self._load_output_files(
                E_file,
                fg_file,
                rea_fg_file,
                pro_fg_file,
                data_file,
                options_file,
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
        fg_file,
        fg_rea_file,
        fg_pro_file,
        data_file,
        options_file,
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
        if Path(fg_file).exists():
            fg_data = np.loadtxt(fg_file, skiprows=1, delimiter=',').T
        if Path(fg_rea_file).exists():
            rea_fg_data = np.loadtxt(fg_rea_file, skiprows=1, delimiter=',').T
        if Path(fg_pro_file).exists():
            pro_fg_data = np.loadtxt(fg_pro_file, skiprows=1, delimiter=',').T
        l_sub_indices = np.where([lf in Lambda for lf in E_data[0]])[0]

        sub_indices = l_sub_indices[::time_sub_sample]

        Lambda_frame = E_data[0, sub_indices]
        E1_pes = E_data[1, sub_indices]
        E2_pes = E_data[2, sub_indices]
        E1_int = E_data[3, sub_indices]
        E2_int = E_data[4, sub_indices]
        E_m_pes = E_data[5, sub_indices]
        E_m_int = E_data[6, sub_indices]

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
            "E_m_int": E_m_int,
            "E1_fg": rea_fg_data,
            "E2_fg": pro_fg_data,
            "E_m_fg": fg_data,
            "Ep": Ep,
            "Ek": Ek,
            "Temp_step": Temp,
            "Vol": Vol,
            "Dens": Dens,
            "options": options,
            "Temp_set": Temp_set,
        }

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
                    else:
                        group[k] = v

            save_group(data, file)

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

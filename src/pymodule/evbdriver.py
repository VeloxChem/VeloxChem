#
#                              VELOXCHEM
#         ----------------------------------------------------
#                     An Electronic Structure Code
#
#  Copyright © 2018-2024 by VeloxChem developers. All rights reserved.
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

from mpi4py import MPI
from pathlib import Path
import numpy as np
import time
import json
import h5py
import sys

from .veloxchemlib import mpi_master
from .molecule import Molecule
from .outputstream import OutputStream
from .mmforcefieldgenerator import MMForceFieldGenerator
from .evbsystembuilder import EvbSystemBuilder
from .evbfepdriver import EvbFepDriver
from .evbffbuilder import EvbForceFieldBuilder
from .evbdataprocessing import EvbDataProcessing
from .solvationbuilder import SolvationBuilder
from .errorhandler import assert_msg_critical

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
        self.debug = False
        self.fast_run = False

        self.t_label = int(time.time())

    def build_and_run_default_water_EVB(self,
                                        reactant: str | Molecule,
                                        product: str | list[str] | Molecule
                                        | list[Molecule],
                                        barrier,
                                        free_energy,
                                        ordered_input=False):
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
        self.build_forcefields(reactant,
                               product,
                               ordered_input=ordered_input,
                               optimise=True)
        self.ostream.print_blank()
        self.ostream.print_header("Building systems")
        self.ostream.flush()
        self.build_systems(configurations=["vacuum", "water"])

        self.ostream.print_blank()
        self.ostream.print_header("Running FEP")
        self.ostream.flush()
        self.run_FEP()
        if not self.debug:
            self.ostream.print_blank()
            self.ostream.print_header("Computing energy profiles")
            self.ostream.flush()
            self.compute_energy_profiles(barrier, free_energy)
        else:
            self.ostream.print_info(
                "Debugging option enabled. Skipping energy profile calculation because recalculation is necessary."
            )

        self.ostream.flush()

    def build_forcefields(
        self,
        reactant: str | Molecule,
        product: str | list[str] | Molecule | list[Molecule],
        reactant_partial_charges: list[float] = None,
        product_partial_charges: list[float] | list[list[float]] = None,
        product_charge: int | list[int] = 0,  # type: ignore
        reactant_multiplicity: int = 1,
        product_multiplicity: int | list[int] = 1,
        reparameterize: bool = True,
        optimise: bool = False,
        ordered_input: bool = False,
        breaking_bonds: tuple[int, int] | list[tuple[int, int]] = None,
        save_output: bool = True,
    ):
        """Build forcefields for the reactant and products, and set self.reactant and self.product to the respective forcefields as well as saving them as json to the input files folder. 
        Will calculate RESP charges and use an xtb Hessian for any necessary reparameterisation. If these files are already present in the input_files folder, they will be loaded instead of recalculated.

        Args:
            reactant (str | Molecule): The reactant. If a string is given, the corresponding xyz file must be present in the input_files folder.
            product (str | list[str] | Molecule | list[Molecule]): A list of products. If a (list of) string(s) is given, the corresponding xyz file(s) must be present in the input_files folder.
            reactant_partial_charges (list[float], optional): The partial charges of the reactant. If not provided, the charges will be calculated using the RESP method. Defaults to None.
            product_partial_charges (list[float] | list[list[float]], optional): The partial charges of each provided product. If not provided, the charges will be calculated using the RESP method.
            product_charge (int | list[int], optional): The nominal charge of each provided product. List should have the same length as the amount of products provided. 
                The reactant will be assigned the sum of the product charges. Defaults to 0.
            reactant_multiplicity (int, optional): The multiplicity of the reactant. Defaults to 1.
            product_multiplicity (int | list[int], optional): The multiplicity of each provided product. List should have the same length as the amount of products provided. Defaults to 1.
            reparameterize (bool, optional): If unknown parameters in the forcefield should be reparameterized. Defaults to True.
            optimise (bool, optional): If the provided structure should be optimised before the forcefield is generated. Defaults to False.
            ordered_input (bool, optional): If set to true, assumes that the reactant and product have the same ordering of atoms, and thus will not attempt to generate a mapping. Defaults to False.
            breaking_bonds (list[tuple[int, int]], optional): A list of tuples of atom-indices of breaking bonds. 
                The atom indices are 0-indexed with respect to the reactant structure, and not all breaking bonds have to be provided. Defaults to None.

        Raises:
            ValueError: If the reactant and product are not given both as a molecule or both as a file.
        """

        if isinstance(reactant, Molecule):
            assert isinstance(product, Molecule) or all(
                isinstance(pro, Molecule) for pro in product
            ), "All products must be Molecule objects if the reactant is a Molecule object"
            if not isinstance(product, list):
                product = [product]

            reactant_name = "reactant"
            combined_product_name = "product"

            reactant_charge = reactant.get_charge()
            reactant_multiplicity = reactant.get_multiplicity()
            if isinstance(product, list):
                product_charge = [pro.get_charge() for pro in product]
                assert reactant_charge == sum(
                    product_charge
                ), "Total charge of reactant and products must match"

                product_multiplicity = [
                    pro.get_multiplicity() for pro in product
                ]
            else:
                product_charge = product.get_charge()
                assert reactant_charge == product_charge, "Total charge of reactant and products must match"
                product_multiplicity = [product.get_multiplicity()]

            rea_input = {
                "molecule": reactant,
                "optimise": None,
                "forcefield": None,
                "hessian": None,
                "charges": None
            }
            pro_input = [{
                "molecule": pro,
                "optimise": None,
                "forcefield": None,
                "hessian": None,
                "charges": None
            } for pro in product]

        elif isinstance(reactant, str):
            assert isinstance(product, str) or all(
                isinstance(pro, str) for pro in product
            ), "All products must be strings if the reactant is a string"
            if not isinstance(product, list):
                product = [product]

            reactant_name = reactant
            product_names = product
            combined_product_name = "_".join(product)
            rea_input = self._get_input_files(reactant_name)
            pro_input = [self._get_input_files(file) for file in product_names]

            if isinstance(product_charge, int):
                if len(product) != 1:
                    assert product_charge == 0, "A charge should be provided for every provided product"

                product_charge: list[int] = [product_charge] * len(product)
                pro_input[0]["molecule"].set_charge(product_charge[0])

                reactant_charge: int = product_charge[0]  # type: ignore
                rea_input["molecule"].set_charge(reactant_charge)

            else:
                reactant_charge = sum(product_charge)
                for pro, charge in zip(pro_input, product_charge):
                    pro["molecule"].set_charge(charge)
                rea_input["molecule"].set_charge(reactant_charge)

            if isinstance(product_multiplicity, int):
                product_multiplicity = [product_multiplicity] * len(product)

            assert len(product) == len(
                product_charge), "Number of products and charges must match"
            assert len(product) == len(
                product_multiplicity
            ), "Number of products and multiplicities must match"

        else:
            raise ValueError(
                "Reactant and product must be either a both string or a Molecule object"
            )

        if reactant_partial_charges is not None:
            rea_input["charges"] = reactant_partial_charges
        if product_partial_charges is not None:
            if isinstance(product_partial_charges[0], float):
                assert len(
                    product
                ) == 1, "Provide a list of seperate partial charges for every separate product"
                pro_input[0]["charges"] = product_partial_charges
            else:
                assert len(product) == len(
                    product_partial_charges
                ), "Amount of products and lists of partial charges must match"
                for pro, charges in zip(pro_input, product_partial_charges):
                    pro["charges"] = charges

        cwd = Path().cwd()
        reactant_path = cwd / self.input_folder / f"{reactant_name}_ff_data.json"
        combined_product_path = cwd / self.input_folder / f"{combined_product_name}_ff_data.json"

        rea_atoms = rea_input['molecule'].number_of_atoms()
        pro_atoms = [pro['molecule'].number_of_atoms() for pro in pro_input]
        assert rea_atoms == sum(
            pro_atoms
        ), f"Number of atoms in reactant ({rea_atoms}) and products ({pro_atoms}, sum={sum(pro_atoms)}) must match"
        if self.name is None:
            self.name = reactant_name

        if rea_input["forcefield"] is not None and combined_product_path.exists(
        ):
            self.ostream.print_info(
                f"Loading combined forcefield data from {combined_product_path}"
            )
            self.ostream.print_info(
                "Found both reactant and product forcefield data. Not generating new forcefields"
            )
            self.reactant = rea_input["forcefield"]
            self.product = self.load_forcefield_from_json(
                str(combined_product_path))
        else:
            ffbuilder = EvbForceFieldBuilder()
            ffbuilder.reparameterize = reparameterize
            ffbuilder.optimise = optimise

            if isinstance(breaking_bonds, tuple):
                breaking_bonds = [breaking_bonds]

            self.reactant, self.product, self.formed_bonds, self.broken_bonds = ffbuilder.build_forcefields(
                rea_input,
                pro_input,
                reactant_charge,
                product_charge,
                reactant_multiplicity,
                product_multiplicity,
                ordered_input,
                breaking_bonds,
            )
        if save_output:
            self.save_forcefield(self.reactant, str(reactant_path))
            self.save_forcefield(self.product, str(combined_product_path))
        self.ostream.flush()

    def _get_input_files(self, filename: str):
        # Build a molecule from a (possibly optimised) geometry
        optimise = True
        cwd = Path().cwd()

        opt_path = cwd / self.input_folder / f"{filename}_xtb_opt.xyz"
        if opt_path.exists():
            self.ostream.print_info(
                f"Loading optimised geometry from {opt_path}")
            molecule = Molecule.read_xyz_file(str(opt_path))
            optimise = False
        else:
            struct_path = cwd / self.input_folder / f"{filename}.xyz"
            self.ostream.print_info(
                f"Loading (possibly unoptimised) geometry from {struct_path}")
            molecule = Molecule.read_xyz_file(str(struct_path))

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

        forcefield = None
        json_path = cwd / self.input_folder / f"{filename}_ff_data.json"
        if json_path.exists():
            self.ostream.print_info(
                f"Loading force field data from {json_path}")
            forcefield = self.load_forcefield_from_json(str(json_path))
            forcefield.molecule = molecule
            if charges is not None:
                forcefield.partial_charges = charges
        else:
            self.ostream.print_info(
                f"Could not find force field data file {self.input_folder}/{filename}_ff_data.json."
            )

        hessian = None
        hessian_path = cwd / self.input_folder / f"{filename}_hess.np"
        if hessian_path.exists():
            self.ostream.print_info(
                f"Found hessian file at {hessian_path}, using it to reparameterize."
            )
            hessian = np.loadtxt(hessian_path)
        else:
            self.ostream.print_info(
                f"Could not find hessian file at {hessian_path}, calculating hessian with xtb and saving it"
            )

        return {
            "molecule": molecule,
            "optimise": optimise,
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
        save_output=True,
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
            if not self.debug:
                if self.fast_run:
                    Lambda = np.linspace(0, 0.1, 6)
                    Lambda = np.append(Lambda[:-1], np.linspace(0.1, 0.9, 21))
                    Lambda = np.append(Lambda[:-1], np.linspace(0.9, 1, 6))
                else:
                    Lambda = np.linspace(0, 0.1, 11)
                    Lambda = np.append(Lambda[:-1], np.linspace(0.1, 0.9, 41))
                    Lambda = np.append(Lambda[:-1], np.linspace(0.9, 1, 11))
                Lambda = np.round(Lambda, 3)
            else:
                Lambda = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
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
        configurations: list[dict] = configurations  # type: ignore
        if constraints is None:
            constraints = []
        if isinstance(constraints, dict):
            constraints = [constraints]

        #Per configuration
        for conf in configurations:
            #create folders,
            if save_output:
                data_folder = f"EVB_{self.name}_{conf['name']}_data_{self.t_label}"
                conf["data_folder"] = data_folder
                run_folder = str(Path(data_folder) / "run")
                conf["run_folder"] = run_folder
                cwd = Path().cwd()
                data_folder_path = cwd / data_folder
                run_folder_path = cwd / run_folder

                data_folder_path.mkdir(parents=True, exist_ok=True)
                run_folder_path.mkdir(parents=True, exist_ok=True)

            # build the system
            system_builder = EvbSystemBuilder()
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

            if save_output:
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

                options_path = cwd / data_folder / "options.json"
                with open(options_path, "w") as file:
                    json.dump(
                        {
                            "temperature": conf.get("temperature",
                                                    self.temperature),
                            "Lambda": Lambda,
                        },
                        file,
                    )

        self.system_confs = configurations
        self.ostream.flush()

    def default_system_configurations(self, name: str) -> dict:
        """Return a dictionary with a default configuration. Options not given in the dictionary will be set to default values in the build_systems function.

        Args:
            name (string): The name of the configuration to be used. Options are "vacuum", "water", "CNT", "graphene", "E_field", "no_reactant"
        """
        if name == "vacuum":
            conf = {
                "name": "vacuum",
                "temperature": self.temperature,
            }
        elif name == "water":
            conf = {
                "name": "water",
                "solvent": "spce",
                "temperature": self.temperature,
                "NPT": True,
                "pressure": 1,
                "padding": 1,
                "ion_count": 0,
                "neutralize": False
            }
        # elif name == "CNT":
        #     conf = {
        #         "name": "water_CNT",
        #         "solvent": "spce",
        #         "temperature": self.temperature,
        #         "NPT": True,
        #         "pressure": 1,
        #         "ion_count": 0,
        #         "CNT": True,
        #         "CNT_radius": 0.5,
        #     }
        # elif name == "graphene":
        #     conf = {
        #         "name": "water_graphene",
        #         "solvent": "spce",
        #         "temperature": self.temperature,
        #         "NPT": True,
        #         "pressure": 1,
        #         "ion_count": 0,
        #         "graphene": True,
        #         "graphene_size": 2,
        #     }
        elif name == "E_field":
            conf = {
                "name": "water_E_field",
                "solvent": "spce",
                "temperature": self.temperature,
                "NPT": True,
                "pressure": 1,
                "padding": 1,
                "ion_count": 0,
                "E_field": [0, 0, 10],
            }
        elif name == "no_reactant":
            conf = {
                "name": "no_reactant",
                "solvent": "spce",
                "temperature": self.temperature,
                "NPT": True,
                "pressure": 1,
                "padding": 1,
                "ion_count": 0,
                "no_reactant": True,
            }
        elif SolvationBuilder()._solvent_properties(name) is not None:
            conf = {
                "name": name,
                "solvent": name,
                "temperature": self.temperature,
                "NPT": True,
                "pressure": 1,
                "padding": 1,
                "ion_count": 0,
            }
        else:
            raise ValueError(f"Unknown system configuration {name}")

        return conf

    def load_initialisation(self,
                            data_folder: str,
                            name: str,
                            skip_systems=False,
                            skip_pdb=False):
        """Load a configuration from a data folder for which the systems have already been generated, such that an FEP can be performed. 
        The topology, initial positions, temperature and Lambda vector will be loaded from the data folder.

        Args:
            data_folder (str): The folder to load the data from
            name (str): The name of the configuration. Can be arbitrary, but should be unique.
            skip_systems (bool, optional): If set to true, the systems will not be loaded from the xml files. Used for debugging. Defaults to False.
            skip_topology (bool, optional): If set to true, the topology will not be loaded from the pdb file. Used for debugging. Defaults to False.
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
        if not skip_systems:
            systems = self.load_systems_from_xml(str(Path(data_folder) / "run"))
            conf["systems"] = systems
        else:
            systems = []

        if not skip_pdb:
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
        equil_steps=10000,
        sample_steps=100000,
        write_step=1000,
        initial_equil_steps=10000,
        step_size=0.001,
        equil_step_size=0.001,
        initial_equil_step_size=0.001,
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
        # if self.debug:
        #     self.ostream.print_warning(
        #         "Debugging enabled, using low number of steps. Do not use for production"
        #     )
        #     self.ostream.flush()
        #     equil_steps = 100
        #     sample_steps = 200
        #     write_step = 5
        #     initial_equil_steps = 100
        #     step_size = 0.001
        #     equil_step_size = 0.001
        #     initial_equil_step_size = 0.001

        if self.fast_run:
            self.ostream.print_warning(
                "Fast run enabled, using modest number of steps. Be careful with using results"
            )
            sample_steps = 25000

        for conf in self.system_confs:
            cwd = Path().cwd()
            options_path = cwd / conf["data_folder"] / "options.json"
            with open(options_path, "r") as file:
                options = json.load(file)

            options.update({
                "equil_steps": equil_steps,
                "sample_steps": sample_steps,
                "write_step": write_step,
                "initial_equil_steps": initial_equil_steps,
                "step_size": step_size,
                "equil_step_size": equil_step_size,
                "initial_equil_step_size": initial_equil_step_size,
            })

            with open(options_path, "w") as file:
                json.dump(options, file, indent=4)

            self.ostream.print_blank()
            self.ostream.print_header(f"Running FEP for {conf['name']}")
            self.ostream.flush()
            FEP = EvbFepDriver()
            FEP.debug = self.debug
            FEP.run_FEP(
                equilibration_steps=equil_steps,
                total_sample_steps=sample_steps,
                write_step=write_step,
                lambda_0_equilibration_steps=initial_equil_steps,
                step_size=step_size,
                equil_step_size=equil_step_size,
                initial_equil_step_size=initial_equil_step_size,
                Lambda=self.Lambda,
                configuration=conf,
            )

    def compute_energy_profiles(
        self,
        barrier,
        free_energy,
        lambda_sub_sample=1,
        lambda_sub_sample_ends=False,
        time_sub_sample=1,
        alpha=None,
        H12=None,
        alpha_guess=None,
        H12_guess=None,
        coordinate_bins=None,
    ):
        """Compute the EVB energy profiles using the FEP results, print the results and save them to an h5 file

        Args:
            barrier (float): the reaction barrier in kJ/mol of the reference system
            free_energy (float): the reaction free energy in kJ/mol of the reference system
            lambda_sub_sample (int, optional): Factor with which the lambda vector will be subsampled. Setting this to two will discard every other lambda frame. Defaults to 1.
            lambda_sub_sample_ends (bool, optional): If set to False, the lambda frames up to 0.1 and from 0.9 will not be subsampled. Defaults to False.
            time_sub_sample (int, optional): Factor with which the time vector will be subsampled. Setting this to two will discard every other snapshot. Defaults to 1.
        """
        dp = EvbDataProcessing()
        results = self._load_output_from_folders(lambda_sub_sample,
                                                 lambda_sub_sample_ends,
                                                 time_sub_sample)
        self.ostream.flush()

        if alpha is not None: dp.alpha = alpha
        if H12 is not None: dp.H12 = H12
        if alpha_guess is not None: dp.alpha_guess = alpha_guess
        if H12_guess is not None: dp.H12_guess = H12_guess
        if coordinate_bins is not None: dp.coordinate_bins = coordinate_bins

        results = dp.compute(results, barrier, free_energy)
        self._save_dict_as_h5(results, f"results_{self.name}_{self.t_label}")
        self.results = results
        self.dataprocessing = dp
        self.print_results()
        self.ostream.flush()

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
        dp.print_results(results, self.ostream)
        self.ostream.flush()

    def plot_results(self, results: dict = None, file_name: str = None):
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
        dp.plot_results(results)
        self.ostream.flush()

    def _load_output_from_folders(self, lambda_sub_sample,
                                  lambda_sub_sample_ends,
                                  time_sub_sample) -> dict:
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
            specific, common = self._load_output_files(E_file, data_file,
                                                       options_file,
                                                       lambda_sub_sample,
                                                       lambda_sub_sample_ends,
                                                       time_sub_sample)
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

    def _load_output_files(self,
                           E_file,
                           data_file,
                           options_file,
                           lambda_sub_sample=1,
                           lambda_sub_sample_ends=False,
                           time_sub_sample=1):
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
            "Ep": Ep,
            "Ek": Ek,
            "Temp_step": Temp,
            "Vol": Vol,
            "Dens": Dens,
            "options": options,
            "Temp_set": Temp_set,
        }

        if len(E_data) > 7:
            E_m_pes = E_data[6, sub_indices]
            E_m_int = E_data[7, sub_indices]
            specific_result.update({"E_m_pes": E_m_pes, "E_m_int": E_m_int})

        lambda_indices = [
            np.where(np.round(Lambda, 3) == L)[0][0] for L in Lambda_frame
        ]
        common_result = {
            "step": step,
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

    def _save_dict_as_h5(self, data: dict, file_name: str):
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

import sys
import time
import json
import h5py
from pathlib import Path

import numpy as np
import openmm as mm
import openmm.app as mmapp
import openmm.unit as mmunit

from mpi4py import MPI

from .outputstream import OutputStream
from .forcefieldgenerator import ForceFieldGenerator

from .evbsystembuilder import EvbSystemBuilder
from .evbfepdriver import FepDriver
from .evbffbuilder import EvbForceFieldBuilder
from .evbdataprocessing import EvbDataProcessing
from .veloxchemlib import mpi_master, Molecule


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

        self.reactant: ForceFieldGenerator = None
        self.product: ForceFieldGenerator = None
        self.input_folder: str = "input_files"

        self.name: str = None
        self.system_confs: list[dict] = []
        self.debug = False

    def build_and_run_default_water_EVB(self, reactant: str | Molecule, product: str | list[str] | Molecule | list[Molecule], barrier, free_energy, ordered_input=False):
        """Automatically perform an EVB calculation using a vacuum system as reference and a system solvated in water as target system.

        Args:
            reactant (str | Molecule): The reactant. If a string is given, the corresponding xyz file must be present in the input_files folder.
            product (str | list[str] | Molecule | list[Molecule]): A list of products. If a (list of) string(s) is given, the corresponding xyz file(s) must be present in the input_files folder.
            barrier (float): the reaction barrier in kcal/mol of the vacuum system
            free_energy (float): the reaction free energy in kcal/mol of the vacuum system
            ordered_input (bool, optional): If set to true, assumes that the reactant and product have the same ordering of atoms, and thus will not attempt to generate a mapping. Defaults to False.
        """
        if not self.debug:
            Lambda = np.linspace(0,0.1,11)
            Lambda = np.append(Lambda[:-1],np.linspace(0.1,0.9,41))
            Lambda = np.append(Lambda[:-1],np.linspace(0.9,1,11))
            Lambda = np.round(Lambda,3)
        else:
            Lambda = [0, 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]
        self.ostream.print_info("Building forcefields")

        self.build_forcefields(reactant, product, ordered_input=ordered_input,optimise=True)
        self.ostream.print_blank()
        self.ostream.print_info("Building systems")
        self.ostream.flush()
        self.build_systems(Lambda=Lambda, configurations=["vacuum", "water"])

        self.ostream.print_blank()
        self.ostream.print_info("Running FEP")
        self.ostream.flush()
        self.run_FEP()
        if not self.debug:
            self.ostream.print_blank()
            self.ostream.print_info("Computing energy profiles")
            self.ostream.flush()
            self.compute_energy_profiles(barrier, free_energy)
        else:
            self.ostream.print_info("Debugging option enabled. Skipping energy profile calculation because recalculation is necessary.")

        self.ostream.flush()

    def build_forcefields(
        self,
        reactant: str | Molecule,
        product: str | list[str] | Molecule | list[Molecule],
        product_charge: int | list[int] = 0,  # type: ignore
        reactant_multiplicity: int = 1,
        product_multiplicity: int | list[int] = 1,
        reparameterise: bool = True,
        optimise: bool = False,
        ordered_input: bool = False,
    ):
        """Build forcefields for the reactant and products, and set self.reactant and self.product to the respective forcefields as well as saving them as json to the input files folder. 
        Will calculate RESP charges and use an xtb Hessian for any necessary reparameterisation. If these files are already present in the input_files folder, they will be loaded instead of recalculated.

        Args:
            reactant (str | Molecule): The reactant. If a string is given, the corresponding xyz file must be present in the input_files folder.
            product (str | list[str] | Molecule | list[Molecule]): A list of products. If a (list of) string(s) is given, the corresponding xyz file(s) must be present in the input_files folder.
            product_charge (int | list[int], optional): The nominal charge of each provided product. List should have the same length as the amount of products provided. 
                The reactant will be assigned the sum of the product charges. Defaults to 0.
            reactant_multiplicity (int, optional): The multiplicity of the reactant. Defaults to 1.
            product_multiplicity (int | list[int], optional): The multiplicity of each provided product. List should have the same length as the amount of products provided. Defaults to 1.
            reparameterise (bool, optional): If unknown parameters in the forcefield should be reparameterised. Defaults to True.
            optimise (bool, optional): If the provided structure should be optimised before the forcefield is generated. Defaults to False.
            ordered_input (bool, optional): If set to true, assumes that the reactant and product have the same ordering of atoms, and thus will not attempt to generate a mapping. Defaults to False.

        Raises:
            ValueError: If the reactant and product are not given both as a molecule or both as a file.
        """
        loaded_forcefield = False
        if isinstance(reactant, Molecule):
            assert isinstance(product, Molecule) or all(isinstance(pro, Molecule) for pro in product), "All products must be Molecule objects if the reactant is a Molecule object"
            if not isinstance(product, list):
                product = [product]

            reactant_name = "reactant"
            combined_product_name = "psroduct"

            reactant_charge = reactant.get_charge()
            reactant_multiplicity = reactant.get_multiplicity()
            if isinstance(product, list):
                product_charge = [pro.get_charge() for pro in product]
                assert reactant_charge == sum(product_charge), "Total charge of reactant and products must match"

                product_multiplicity = [pro.get_multiplicity() for pro in product]
            else:
                product_charge = product.get_charge()
                assert reactant_charge == product_charge, "Total charge of reactant and products must match"
                product_multiplicity = product.get_multiplicity()

            rea_input = {"molecule": reactant, "optimise": None, "forcefield": None, "hessian": None, "charges": None}
            pro_input = [{"molecule": pro, "optimise": None, "forcefield": None, "hessian": None, "charges": None} for pro in product]

        elif isinstance(reactant, str):
            assert isinstance(product, str) or all(isinstance(pro, str) for pro in product), "All products must be strings if the reactant is a string"
            if not isinstance(product, list):
                product = [product]

            reactant_name = reactant
            product_names = product
            combined_product_name = "_".join(product)
            rea_input = self._get_input_files(reactant_name)
            pro_input = [self._get_input_files(file) for file in product_names]

            if isinstance(product_charge, int):
                product_charge: list[int] = [product_charge]
                reactant_charge: int = product_charge[0]  # type: ignore
            else:
                reactant_charge = sum(product_charge)

            if isinstance(product_multiplicity, int):
                product_multiplicity = [product_multiplicity] * len(product)
            
            assert len(product) == len(product_charge), "Number of products and charges must match"
            assert len(product) == len(product_multiplicity), "Number of products and multiplicities must match"

        else:
            raise ValueError("Reactant and product must be either a both string or a Molecule object")

        cwd = Path().cwd()
        reactant_path = cwd / self.input_folder / f"{reactant_name}_ff_data.json"
        combined_product_path = cwd / self.input_folder / f"{combined_product_name}_ff_data.json"

        rea_atoms = rea_input['molecule'].number_of_atoms()
        pro_atoms = [pro['molecule'].number_of_atoms() for pro in pro_input]
        assert rea_atoms == sum(pro_atoms), f"Number of atoms in reactant ({rea_atoms}) and products ({pro_atoms}, sum={sum(pro_atoms)}) must match"
        self.name = reactant_name

        if rea_input["forcefield"] is not None and combined_product_path.exists():
            self.ostream.print_info(f"Loading combined forcefield data from {combined_product_path}")
            self.ostream.print_info("Found both reactant and product forcefield data. Not generating new forcefields")
            self.reactant = rea_input["forcefield"]
            self.product = self.load_forcefield_from_json(str(combined_product_path))
        else:
            ffbuilder = EvbForceFieldBuilder()
            ffbuilder.reparameterise = reparameterise
            ffbuilder.optimise = optimise

            self.reactant, self.product = ffbuilder.build_forcefields(
                rea_input,
                pro_input,
                reactant_charge,
                product_charge,
                reactant_multiplicity,
                product_multiplicity,
                ordered_input,
            )
        self.save_forcefield(self.reactant, str(reactant_path))
        self.save_forcefield(self.product, str(combined_product_path))

    def _get_input_files(self, filename: str):
        # Build a molecule from a (possibly optimised) geometry
        optimise = True
        cwd = Path().cwd()

        opt_path = cwd / self.input_folder / f"{filename}_xtb_opt.xyz"
        if opt_path.exists():
            self.ostream.print_info(f"Loading optimised geometry from {opt_path}")
            molecule = Molecule.read_xyz_file(str(opt_path))
            optimise = False
        else:
            struct_path = cwd / self.input_folder / f"{filename}.xyz"
            self.ostream.print_info(f"Loading (possibly unoptimised) geometry from {struct_path}")
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
                        self.ostream.print_info(f"Could not read line {line} from {charge_path}. Continuing")
            print_charge = [round(charge, 3) for charge in charges]
            self.ostream.print_info(
                f"Loading charges from {charge_path} file, total charge = {print_charge}"
            )

        forcefield = None
        json_path = cwd / self.input_folder / f"{filename}_ff_data.json"
        if json_path.exists():
            self.ostream.print_info(f"Loading force field data from {json_path}")
            forcefield = self.load_forcefield_from_json(str(json_path))
            forcefield.molecule = molecule
            if charges is not None:
                forcefield.partial_charges = charges
        else:
            self.ostream.print_info(f"Could not find force field data file {self.input_folder}/{filename}_ff_data.json.")

        hessian = None
        hessian_path = cwd / self.input_folder / f"{filename}_hess.np"
        if hessian_path.exists():
            self.ostream.print_info(
                f"Found hessian file at {hessian_path}, using it to reparameterise.")
            hessian = np.loadtxt(hessian_path)
        else:
            self.ostream.print_info(
                f"Could not find hessian file at {hessian_path}, calculating hessian with xtb and saving it"
            )

        return {"molecule": molecule, "optimise": optimise, "forcefield": forcefield, "hessian": hessian, "charges": charges}

    #todo, should be moved to forcefieldgenerator class
    @staticmethod
    def load_forcefield_from_json(path: str) -> ForceFieldGenerator:
        """
        Load forcefield data from a JSON file.

        Args:
            Path (str): The path to the JSON file containing the forcefield data.

        Returns:
            ForceFieldGenerator: The updated forcefield object with the loaded data.
        """
        with open(path, "r", encoding="utf-8") as file:
            forcefield = ForceFieldGenerator()
            ff_data = json.load(file)

            forcefield.atoms = EvbDriver._str_to_tuple_key(ff_data["atoms"])
            forcefield.bonds = EvbDriver._str_to_tuple_key(ff_data["bonds"])
            forcefield.angles = EvbDriver._str_to_tuple_key(ff_data["angles"])
            forcefield.dihedrals = EvbDriver._str_to_tuple_key(ff_data["dihedrals"])
            forcefield.impropers = EvbDriver._str_to_tuple_key(ff_data["impropers"])
        return forcefield

    #todo, should be moved to forcefieldgenerator class
    @staticmethod
    def save_forcefield(forcefield: ForceFieldGenerator, path: str):
        """
        Save the forcefield data of the forcefieldgenerator to a JSON file, converting all tuples to strings

        Args:
            forcefield (ForceFieldGenerator): The forcefield object containing the data to be saved.
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
                tuple += (int(item),)
            if len(tuple) == 1:
                tuple = tuple[0]
            tup_keys.append(tuple)
        return {key: value for key, value in zip(tup_keys,dictionary.values())}

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
        Lambda: list[float]|np.ndarray,
        configurations: list[str] | list[dict],  # type: ignore
        constraints: dict | list[dict] | None = None,
    ):
        """Build OpenMM systems for the given configurations with interpolated forcefields for each lambda value. Saves the systems as xml files, the topology as a pdb file and the options as a json file to the disk.

        Args:
            Lambda (list[float] | np.ndarray): The Lambda vector to be used for the FEP. Should start with 0, end with 1 and be monotonically increasing.
            configurations (list[str] | list[dict]): The given configurations for which to perform an FEP. The first configuration will be regarded as the reference configuration. 
                If a string is given, the return value of default_system_configurations() will be used. See this function for default configurations.
            constraints (dict | list[dict] | None, optional): Dictionary of harmonic bond, angle or (improper) torsion forces to apply over in every FEP frame. Defaults to None.
        """
        assert (Lambda[0] == 0 and Lambda[-1] == 1), f"Lambda must start at 0 and end at 1. Lambda = {Lambda}"
        assert np.all(np.diff(Lambda) > 0), f"Lambda must be monotonically increasing. Lambda = {Lambda}"
        Lambda = [round(lam, 3) for lam in Lambda]
        self.Lambda = Lambda

        if all(isinstance(conf, str) for conf in configurations):
            configurations = [self.default_system_configurations(conf) for conf in configurations]

        assert all(isinstance(conf, dict)
                   for conf in configurations), "Configurations must be a list of strings or a list of dictionaries"
        configurations: list[dict] = configurations  # type: ignore
        if constraints is None:
            constraints = []
        if isinstance(constraints, dict):
            constraints = [constraints]

        #Per configuration
        for conf in configurations:
            t = int(time.time())
            #create folders,
            
            data_folder = f"EVB_{self.name}_{conf["name"]}_data_{t}"
            conf["data_folder"] = data_folder
            run_folder = f"{data_folder}/run"
            conf["run_folder"] = run_folder
            cwd = Path().cwd()
            data_folder_path = cwd / data_folder
            run_folder_path = cwd / run_folder
            self.ostream.print_info(f"Saving files to {data_folder_path} and {run_folder_path}")

            data_folder_path.mkdir(parents=True, exist_ok=True)
            run_folder_path.mkdir(parents=True, exist_ok=True)

            # build the system
            system_builder = EvbSystemBuilder()

            systems, topology, initial_positions = system_builder.build_systems(
                reactant=self.reactant,
                product=self.product,
                Lambda=self.Lambda,
                configuration=conf,
                constraints=constraints,
            )

            self.save_systems_as_xml(systems, conf["run_folder"])

            top_path = cwd / data_folder / "topology.pdb"
            mmapp.PDBFile.writeFile(
                topology,
                initial_positions * 10,  # positions are handled in nanometers, but pdb's should be in angstroms
                open(top_path, "w"),
            )

            options_path = cwd / data_folder / "options.json"
            with open(options_path, "w") as file:
                json.dump(
                    {
                        "temperature": conf.get("temperature", self.temperature),
                        "Lambda": Lambda,
                    },
                    file,
                )

            conf["systems"] = systems
            conf["topology"] = topology
            conf["initial_positions"] = initial_positions

        self.system_confs = configurations

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
                "padding": 1.2,
                "ion_count": 0,
            }
        elif name == "CNT":
            conf = {
                "name": "water_CNT",
                "solvent": "spce",
                "temperature": self.temperature,
                "NPT": True,
                "pressure": 1,
                "ion_count": 0,
                "CNT": True,
                "M": 5,
                "N": 9,
            }
        elif name == "graphene":
            conf = {
                "name": "water_graphene",
                "solvent": "spce",
                "temperature": self.temperature,
                "NPT": True,
                "pressure": 1,
                "ion_count": 0,
                "graphene": True,
                "M": 5,
                "N": 9,
            }
        elif name == "E_field":
            conf = {
                "name": "water_E_field",
                "solvent": "spce",
                "temperature": self.temperature,
                "NPT": True,
                "pressure": 1,
                "padding": 1.5,
                "ion_count": 0,
                "E_field": [0, 0, 10],
            }

        elif name == "no_reactant":
            conf = {
                "name": "no_reactant",
                "solvent": "spce",
                "temperature": temperature,
                "NPT": True,
                "pressure": 1,
                "padding": 1.5,
                "ion_count": 0,
                "no_reactant": True,
            }
        else:
            raise ValueError(f"Unknown system configuration {name}")

        return conf

    # def load_initialisation(self, data_folder: str):
    #     self.data_folder = data_folder
    #     self.run_folder = f"{data_folder}/run"

    #     pdb = mmapp.PDBFile(f"{data_folder}/topology.pdb")
    #     self.topology = pdb.getTopology()
    #     self.initial_positions = pdb.getPositions(asNumpy=True).value_in_unit(mmunit.nanometers)

    #     with open(f"{data_folder}/options.json", "r") as file:
    #         options = json.load(file)
    #         self.temperature = options["temperature"]
    #         self.Lambda = options["Lambda"]

    #     self.systems = self.load_systems_from_xml(self.run_folder)

    #     self.ostream.print_info(f"Loaded systems, topology, initial positions, temperatue and Lambda from {data_folder}")

    def save_systems_as_xml(self, systems: dict, folder: str):
        """Save the systems as xml files to the given folder.

        Args:
            systems (dict): The systems to save
            folder (str): The folder relative to the current working directory to save the systems to.
        """
        path = Path().cwd() / folder
        self.ostream.print_info(f"Saving systems to {path}")
        for lam in self.Lambda:
            file_path = str(path / f"{lam:.3f}_sys.xml")
            with open(file_path, mode="w", encoding="utf-8") as output:
                output.write(mm.XmlSerializer.serialize(systems[lam]))

        file_path = str(path / "reactant.xml")
        with open(file_path, mode="w", encoding="utf-8") as output:
            output.write(mm.XmlSerializer.serialize(systems["reactant"]))

        file_path = str(path / "product.xml")
        with open(file_path, mode="w", encoding="utf-8") as output:
            output.write(mm.XmlSerializer.serialize(systems["product"]))

    def load_systems_from_xml(self, folder: str):
        """Load the systems from xml files in the given folder.

        Args:
            folder (str): The folder relative to the current working directory to load the systems from.
        Returns:
            dict: The loaded systems
        """
        systems = {}
        path = Path().cwd() / folder
        for lam in self.Lambda:
            with open(path / f"{lam:.3f}_sys.xml", mode="r", encoding="utf-8") as input:
                systems[lam] = mm.XmlSerializer.deserialize(input.read())
        with open(path / "reactant.xml", mode="r", encoding="utf-8") as input:
            systems["reactant"] = mm.XmlSerializer.deserialize(input.read())
        with open(path / "product.xml", mode="r", encoding="utf-8") as input:
            systems["product"] = mm.XmlSerializer.deserialize(input.read())
        return systems

    def run_FEP(
        self,
        equil_steps=5000,
        sample_steps=100000,
        write_step=1000,
        initial_equil_steps=5000,
        step_size=0.001,
        equil_step_size=0.002,
        initial_equil_step_size=0.002,
    ):
        """Run the the FEP calculations for all configurations in self.system_confs.

        Args:
            equil_steps (int, optional): The amount of timesteps to equiliberate at the beginning af each Lambda frame. Equiliberation is done with frozen H-bonds. Defaults to 5000.
            sample_steps (int, optional): The amount of steps to sample. Defaults to 100000.
            write_step (int, optional): Per how many steps to take a sample and save its data as well as the trajectory point. Defaults to 1000.
            initial_equil_steps (int, optional): The amount of timesteps to add to the equiliberation at the first Lambda frame. Defaults to 5000.
            step_size (float, optional): The step size during the sampling in picoseconds. Defaults to 0.001.
            equil_step_size (float, optional): The step size during the equiliberation in picoseconds. Is typically larger then step_size as equilliberation is done with frozen H-bonds. Defaults to 0.002.
            initial_equil_step_size (float, optional): The step size during initial equiliberation in picoseconds. Defaults to 0.002.
        """
        if self.debug:
            self.ostream.print_info("Debugging enabled, using low number of steps. Do not use for production")
            equil_steps = 100
            sample_steps = 500
            write_step = 1
            initial_equil_steps = 100
            step_size = 0.001
            equil_step_size = 0.001
            initial_equil_step_size = 0.001

        for conf in self.system_confs:
            self.ostream.print_info(f"Running FEP for {conf['name']}")
            FEP = FepDriver()
            # FEP.constrain_H = False
            FEP.run_FEP(
                equilliberation_steps=equil_steps,
                total_sample_steps=sample_steps,
                write_step=write_step,
                lambda_0_equilliberation_steps=initial_equil_steps,
                step_size=step_size,
                equil_step_size=equil_step_size,
                initial_equil_step_size=initial_equil_step_size,
                Lambda=self.Lambda,
                configuration=conf
            )

            if self.debug:
                self.ostream.print_info("Debugging option enabled.Skipping recalculation.")
            else:
                FEP.recalculate(interpolated_potential=True, force_contributions=True)

    def compute_energy_profiles(self, barrier, free_energy):
        """Compute the EVB energy profiles using the FEP results, print the results and save them to an h5 file

        Args:
            barrier (float): the reaction barrier in kcal/mol of the reference system
            free_energy (float): the reaction free energy in kcal/mol of the reference system
        """
        reference_folder = self.system_confs[0]["data_folder"]
        target_folders = [conf["data_folder"] for conf in self.system_confs[1:]]
        dp = EvbDataProcessing()
        results = self.load_output_from_folders(reference_folder, target_folders)
        dp.compute(results, barrier, free_energy)
        self.save_results(dp.results)
        dp.print_results()

    def load_output_from_folders(self, reference_folder, target_folders) -> dict:
        """Load results from the output of the FEP calculations, and return a dictionary. Looks for Energies.dat, Data_combined.dat and options.json in each folder.

        Args:
            reference_folder (str): The folder containing the output files of the reference system.
            target_folders (str): A list of folders containing the output files of the target systems.

        Returns:
            dict: A dictionary containing the results of the FEP calculations for the reference and target systems.
        """
        reference = reference_folder.split("/")[-1]
        targets = []
        for target in target_folders:
            targets.append(target.split("/")[-1])

        folders = [reference_folder] + target_folders
        results = {}
        cwd = Path().cwd()
        for name, folder in zip([reference] + targets, folders):
            E_file = str(cwd / folder / "Energies.dat")
            data_file = str(cwd / folder / "Data_combined.dat")
            options_file = str(cwd / folder / "options.json")
            result = self.load_output_files(E_file, data_file, options_file)
            results.update({name: result})
        return results

    @staticmethod
    def load_output_files(E_file, data_file, options_file):
        """Load the output from one FEP calculation

        Args:
            E_file (path): The location of the energies file
            data_file (path): The location of the data file
            options_file (path): The location of the options json file.

        Returns:
            dict: A dictionary containing the results of the FEP calculation
        """
        E = np.loadtxt(E_file)
        joule_to_cal = 1 / 4.184
        E *= joule_to_cal
        E1_ref, E2_ref, E1_run, E2_run, E_m = E.T

        Data = np.loadtxt(data_file)
        step, Ep, Ek, Temp, Vol, Dens, Lambda_frame = Data.T

        with open(options_file, "r") as file:
            options = json.load(file)
        Lambda = options["Lambda"]
        Temp_set = options["temperature"]
        options.pop("Lambda")
        result = {
            "E1_ref": E1_ref,
            "E2_ref": E2_ref,
            "E1_run": E1_run,
            "E2_run": E2_run,
            "E_m": E_m,
            "step": step,
            "Ep": Ep,
            "Ek": Ek,
            "Temp_step": Temp,
            "Temp_set": Temp_set,
            "Vol": Vol,
            "Dens": Dens,
            "Lambda": Lambda,
            "Lambda_frame": Lambda_frame,
            "Lambda_indices": [np.where(np.round(Lambda, 3) == L)[0][0] for L in Lambda_frame],
            "options": options,
        }
        return result

    @staticmethod
    def load_result(file):
        """Load results from an h5 file

        Args:
            file (path): The file to load the results from.

        Returns:
            dict: Dictionary with the results
        """
        with h5py.File(file, "r") as f:
            result = {}
            for key, value in f.items():
                if isinstance(value, h5py.Dataset):
                    result[key] = value[()]
                elif isinstance(value, h5py.Group):
                    group = {}
                    for k, v in value.items():
                        group[k] = v[()]
                    result[key] = group
                else:
                    result[key] = value
        return result

    def save_results(self, results: dict[str, dict]):
        """Save the results to seperate h5 files

        Args:
            results (dict[str, dict]): Dictionary of dictionaries containing the results to save. The key string will be used as filename
        """
        cwd = Path.cwd()
        for name, result in results.items():
            file_path = str(cwd / f"{name}.h5")
            with h5py.File(file_path, "w") as file:
                self.ostream.print_info(f"Saving results to {file_path}")
                for key, value in result.items():
                    if isinstance(value, np.ndarray) or isinstance(value, list):
                        file.create_dataset(key, data=value)
                    elif isinstance(value, dict):
                        group = file.create_group(key)
                        for k, v in value.items():
                            group.create_dataset(k, data=v)
                    else:
                        file[key] = value
                pass
        pass

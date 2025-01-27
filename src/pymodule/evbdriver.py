import sys
import os
import json
import h5py
from pathlib import Path

import numpy as np
import openmm as mm
import openmm.app as mmapp
import openmm.unit as mmunit

from mpi4py import MPI

from .veloxchemlib import mpi_master
from .outputstream import OutputStream
from .forcefieldgenerator import ForceFieldGenerator

from .evbsystembuilder import EvbSystemBuilder
from .evbfepdriver import FepDriver
from .evbffbuilder import EvbForceFieldBuilder
from .evbdataprocessing import EvbDataProcessing
from .veloxchemlib import Molecule


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

        # MPI information
        self.comm = comm
        self.rank = comm.Get_rank()
        self.size = comm.Get_size()
        self.temperature: float = 300
        self.Lambda: list[float]

        self.input_folder: str

        self.reactant: ForceFieldGenerator
        self.product: ForceFieldGenerator
        self.products: list[ForceFieldGenerator]
        self.input_folder: str = "./input_files"

        self.name: str
        self.system_confs: list[dict] = []
        self.debug = False

    def build_and_run_default_water_EVB(self, reactant: str, product: str | list[str], barrier, free_energy, ordered_input=False):

        if not self.debug:
            Lambda = np.linspace(0,0.1,11)
            Lambda = np.append(Lambda[:-1],np.linspace(0.1,0.9,41))
            Lambda = np.append(Lambda[:-1],np.linspace(0.9,1,11))
            Lambda = np.round(Lambda,3)
        else:
            Lambda = [0, 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]

        self.build_forcefields(reactant, product, ordered_input=ordered_input,optimise=True)
        self.build_systems(Lambda=Lambda, configurations=["vacuum", "water"])

        self.run_FEP()
        if not self.debug:
            self.compute_energy_profiles(barrier, free_energy)
        else:
            print("Debugging option enabled. Skipping energy profile calculation because recalculation is necessary.")

    def build_forcefields(
        self,
        reactant_file: str,
        product_file: str | list[str],
        product_charge: int | list[int] | None = None,  # type: ignore
        reactant_multiplicity=1,
        product_multiplicity: int | list[int] | None = None,
        reparameterise: bool = True,
        optimise: bool = False,
        ordered_input: bool = False,
    ):

        if isinstance(product_file, str):
            product_file = [product_file]

        if product_charge is None:
            product_charge = [0] * len(product_file)
            reactant_charge = 0

        elif isinstance(product_charge, int):
            product_charge: list[int] = [product_charge]
            reactant_charge: int = product_charge[0]  # type: ignore
        else:
            reactant_charge = sum(product_charge)

        if isinstance(product_multiplicity, int):
            product_multiplicity = [product_multiplicity]
        if product_multiplicity is None:
            product_multiplicity = [1] * len(product_file)

        assert len(product_file) == len(product_charge), "Number of products and charges must match"
        assert len(product_file) == len(product_multiplicity), "Number of products and multiplicities must match"
        self.name = reactant_file
        
        ffbuilder = EvbForceFieldBuilder()

        ffbuilder.reparameterise = reparameterise
        ffbuilder.optimise = optimise

        rea_input = self._get_input_files(reactant_file)
        pro_input = [self._get_input_files(file) for file in product_file]

        reactant_path = f"{self.input_folder}/{reactant_file}_ff_data.json"

        combined_product_name = "_".join(product_file)
        combined_product_path = f"{self.input_folder}/{combined_product_name}_ff_data.json"
        combined_product_exists = os.path.exists(combined_product_path)
        if rea_input["forcefield"] is not None and combined_product_exists:
            print(f"Loading combined forcefield data from {combined_product_path}")
            print("Found both reactant and product forcefield data. Not generating new forcefields")
            self.reactant = rea_input["forcefield"]
            self.product = self.load_forcefield_from_json(combined_product_path)
        else:
            if not combined_product_exists:
                print(f"Could not find combined forcefield data file {combined_product_path}. Generating new forcefields")

            self.reactant, self.product = ffbuilder.build_forcefields(
                rea_input, 
                pro_input,
                ordered_input,
                reactant_charge=reactant_charge,
                product_charge=product_charge,
                reactant_multiplicity=reactant_multiplicity,
                product_multiplicity=product_multiplicity,
            )
            self.save_forcefield(self.reactant, reactant_path)
            self.save_forcefield(self.product, combined_product_path)
        

    def _get_input_files(self, filename):
        # Build a molecule from a (possibly optimised) geometry
        optimise = True
        if os.path.exists(f"{self.input_folder}/{filename}_xtb_opt.xyz"):
            print(f"Loading optimised geometry from {self.input_folder}/{filename}_xtb_opt.xyz")
            molecule = Molecule.read_xyz_file(f"{self.input_folder}/{filename}_xtb_opt.xyz")
            optimise = False
        else:
            print(f"Loading (possibly unoptimised) geometry from {self.input_folder}/{filename}.xyz")
            molecule = Molecule.read_xyz_file(f"{self.input_folder}/{filename}.xyz")

        json_path = f"{self.input_folder}/{filename}_ff_data.json"

        charges = None
        if os.path.exists(f"{self.input_folder}/{filename}_charges.txt"):
            with open(f"{self.input_folder}/{filename}_charges.txt", "r", encoding="utf-8") as file:
                charges = self._load_charges(f"{self.input_folder}/{filename}_charges.txt")
            print(
                f"Loading charges from {self.input_folder}/{filename}_charges.txt file, total charge = {charges}"
            )

        forcefield = None
        if os.path.exists(json_path):
            print(f"Loading force field data from {json_path}")
            forcefield = self.load_forcefield_from_json(json_path)
            forcefield.molecule = molecule
            if charges is not None:
                forcefield.partial_charges = charges
        else:
            print(f"Could not find force field data file {self.input_folder}/{filename}_ff_data.json.")

        hessian = None
        if os.path.exists(f"{self.input_folder}/{filename}_hess.np"):
            print(
                f"Found hessian file at {self.input_folder}/{filename}_hess.np, using it to reparameterise.")
            hessian = np.loadtxt(f"{self.input_folder}/{filename}_hess.np")
        else:
            print(
                f"Could not find hessian file at {self.input_folder}/{filename}_hess.np, calculating hessian with xtb and saving it"
            )

        return {"molecule": molecule, "optimise": optimise, "forcefield": forcefield, "hessian": hessian, "charges": charges}

    #todo, should be moved to forcefieldgenerator class
    @staticmethod
    def load_forcefield_from_json(path: str) -> ForceFieldGenerator:
        """
        Load forcefield data from a JSON file.

        Args:
            forcefield (ForceFieldGenerator): The forcefield object to load the data into.
            filename (str): The name of the JSON file, without _ff_data.json.

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
    def _load_charges(path: str) -> list:
        """
        Load a list of charges from a file.

        Args:
            filename (str): The name of the file to load charges from.

        Returns:
            list: A list of charges read from the file.
        """
        with open(path, "r", encoding="utf-8") as file:
            charges = []
            for line in file:
                try:
                    charges.append(float(line))
                except ValueError:
                    print(rf"Could not read line {line} from {filename}_charges.txt. Continuing")
        return charges

    # def _save_charges(self, charges: list, filename: str):
    #     """
    #     Save the given list of charges to a text file.

    #     Args:
    #         charges (list): A list of charges to be saved.
    #         filename (str): The name of the file to save the charges to.

    #     Returns:
    #         None
    #     """
    #     with open(f"{self.input_folder}/{filename}_charges.txt", "w", encoding="utf-8") as file:
    #         for charge in charges:
    #             file.write(f"{charge}\n")

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
            #create folders,
            data_folder = f"EVB_{self.name}_{conf["name"]}_data"
            conf["data_folder"] = data_folder
            run_folder = f"{data_folder}/run"
            conf["run_folder"] = run_folder
            print(f"Saving files to {data_folder}")

            if not os.path.exists(data_folder):
                os.makedirs(data_folder)
            else:
                import shutil

                folder = f"./{data_folder}"
                for filename in os.listdir(folder):
                    file_path = os.path.join(folder, filename)
                    try:
                        if os.path.isfile(file_path) or os.path.islink(file_path):
                            os.unlink(file_path)
                        elif os.path.isdir(file_path):
                            shutil.rmtree(file_path)
                    except Exception as e:
                        print("Failed to delete %s. Reason: %s" % (file_path, e))

            if not os.path.exists(run_folder):
                os.makedirs(run_folder)
            else:
                import shutil

                folder = f"./{run_folder}"
                for filename in os.listdir(folder):
                    file_path = os.path.join(folder, filename)
                    try:
                        if os.path.isfile(file_path) or os.path.islink(file_path):
                            os.unlink(file_path)
                        elif os.path.isdir(file_path):
                            shutil.rmtree(file_path)
                    except Exception as e:
                        print("Failed to delete %s. Reason: %s" % (file_path, e))

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

            mmapp.PDBFile.writeFile(
                topology,
                initial_positions * 10,  # positions are handled in nanometers, but pdb's should be in angstroms
                open(f"{data_folder}/topology.pdb", "w"),
            )

            with open(f"{data_folder}/options.json", "w") as file:
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

    def default_system_configurations(self, name):
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

    #     print(f"Loaded systems, topology, initial positions, temperatue and Lambda from {data_folder}")

    def save_systems_as_xml(self, systems: dict, folder: str):
        path = Path().cwd() / folder
        print(f"Saving systems to {path}")
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
        systems = {}
        for lam in self.Lambda:
            with open(f"{folder}/{lam:.3f}_sys.xml", mode="r", encoding="utf-8") as input:
                systems[lam] = mm.XmlSerializer.deserialize(input.read())
        with open(f"{folder}/reactant.xml", mode="r", encoding="utf-8") as input:
            systems["reactant"] = mm.XmlSerializer.deserialize(input.read())
        with open(f"{folder}/product.xml", mode="r", encoding="utf-8") as input:
            systems["product"] = mm.XmlSerializer.deserialize(input.read())
        return systems

    def run_FEP(
        self,
        equil_steps=5000,
        total_sample_steps=100000,
        write_step=1000,
        initial_equil_steps=5000,
        step_size=0.001,
        equil_step_size=0.002,
        initial_equil_step_size=0.002,
    ):

        if self.debug:
            print("Debugging enabled, using low number of steps. Do not use for production")
            equil_steps = 100
            total_sample_steps = 500
            write_step = 1
            initial_equil_steps = 100
            step_size = 0.001
            equil_step_size = 0.001
            initial_equil_step_size = 0.001

        for conf in self.system_confs:
            print(f"Running FEP for {conf['name']}")
            FEP = FepDriver()
            # FEP.constrain_H = False
            FEP.run_FEP(
                equilliberation_steps=equil_steps,
                total_sample_steps=total_sample_steps,
                write_step=write_step,
                lambda_0_equilliberation_steps=initial_equil_steps,
                step_size=step_size,
                equil_step_size=equil_step_size,
                initial_equil_step_size=initial_equil_step_size,
                Lambda=self.Lambda,
                configuration=conf
            )

            if self.debug:
                print("Debugging option enabled.Skipping recalculation.")
            else:
                FEP.recalculate(interpolated_potential=True, force_contributions=True)

    def compute_energy_profiles(self, barrier, free_energy):
        reference_folder = self.system_confs[0]["data_folder"]
        target_folders = [conf["data_folder"] for conf in self.system_confs[1:]]
        dp = EvbDataProcessing()
        results = self.load_output_from_folders(reference_folder, target_folders)
        dp.compute(results, barrier, free_energy)
        self.save_results(dp.results)
        dp.print_results()

    def load_output_from_folders(self, reference_folder, target_folders):
        reference = reference_folder.split("/")[-1]
        targets = []
        for target in target_folders:
            targets.append(target.split("/")[-1])

        folders = [reference_folder] + target_folders
        results = {}
        for name, folder in zip([reference] + targets, folders):
            E_file = f"{folder}/Energies.dat"
            data_file = f"{folder}/Data_combined.dat"
            options_file = f"{folder}/options.json"
            result = self.load_output_files(E_file, data_file, options_file)
            results.update({name: result})
        return results

    @staticmethod
    def load_output_files(E_file, data_file, options_file):
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

    @staticmethod
    def save_results(results):
        cwd = Path.cwd()
        for name, result in results.items():
            file_path = str(cwd / f"{name}.h5")
            with h5py.File(file_path, "w") as file:
                print(f"Saving results to {file_path}")
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

    @staticmethod
    def show_snapshots(folder):
        import py3Dmol
        # Read the pdb file and split it into models
        with open(f"{folder}/traj_combined.pdb", "r") as file:
            models = file.read().split("ENDMDL")

        # Extract the first and last model
        first_model = models[0] + "ENDMDL"
        last_model = models[-2] + "ENDMDL"  # -2 because the last element is an empty string

        # Display the first model
        view = py3Dmol.view(width=400, height=300)
        view.addModel(first_model, "pdb", {"keepH": True})
        view.setStyle({}, {"stick": {}, "sphere": {"scale": 0.25}})
        view.zoomTo()
        view.show()

        # Display the last model
        view = py3Dmol.view(width=400, height=300)
        view.addModel(last_model, "pdb", {"keepH": True})
        view.setStyle({}, {"stick": {}, "sphere": {"scale": 0.25}})
        view.zoomTo()
        view.show()

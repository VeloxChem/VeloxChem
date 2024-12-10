import sys
import os
import json

import numpy as np
import openmm as mm
import openmm.app as mmapp
import openmm.unit as mmunit

from mpi4py import MPI

from .veloxchemlib import mpi_master
from .outputstream import OutputStream
from .forcefieldgenerator import ForceFieldGenerator

from .evbsystembuilder import EvbSystemBuilder
from .fepdriver import FepDriver
from .evbffbuilder import EvbForceFieldBuilder
from .evbdataprocessing import EvbDataProcessing


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
        
        Lambda = list(np.linspace(0, 1, 51))
        if self.debug:
            Lambda = [0, 0.5, 1]
        self.build_forcefields(reactant, product, ordered_input=ordered_input)
        self.build_systems(Lambda=Lambda, configurations=["vacuum", "water"])
        self.run_FEP(test_run=self.debug)
        self.compute_energy_profiles(barrier, free_energy)

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
        self.name = reactant_file
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
        ffbuilder = EvbForceFieldBuilder()
        ffbuilder.reactant_charge = reactant_charge
        ffbuilder.reactant_multiplicity = reactant_multiplicity
        ffbuilder.product_charge = product_charge
        ffbuilder.product_multiplicity = product_multiplicity

        ffbuilder.reparameterise = reparameterise
        ffbuilder.optimise = optimise

        self.reactant, self.product = ffbuilder.build_forcefields(reactant_file, product_file, ordered_input)

    def build_systems(
        self,
        Lambda: list[float],
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
                temperature=conf.get("temperature", self.temperature),
                NPT=conf.get("NPT", False),
                pressure=conf.get("pressure", 1),
                solvent=conf.get("solvent", None),
                padding=conf.get("padding", 1.2),
                CNT=conf.get("CNT", False),
                Graphene=conf.get("graphene", False),
                M=conf.get("M", 5),
                N=conf.get("N", 9),
                ion_count=conf.get("ion_count", 0),
                no_reactant=conf.get("no_reactant", False),
                E_field=conf.get("E_field", [0, 0, 0]),
            )

            self.save_systems_as_xml(systems, conf["run_folder"])

            mmapp.PDBFile.writeFile(
                topology,
                initial_positions * 10,  # positions are handled in nanometers, but pdb's should be in angstroms
                open(f"{data_folder}/topology.pdb", "w"),
            )

            #todo save all options from configuration
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
                "padding": 1.5,
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
        for lam in self.Lambda:
            with open(f"{folder}/{lam:.3f}_sys.xml", mode="w", encoding="utf-8") as output:
                output.write(mm.XmlSerializer.serialize(systems[lam]))
        with open(f"{folder}/reactant.xml", mode="w", encoding="utf-8") as output:
            output.write(mm.XmlSerializer.serialize(systems["reactant"]))
        with open(f"{folder}/product.xml", mode="w", encoding="utf-8") as output:
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
        platform="CPU",
        test_run=False
    ):
                
        if test_run:
            print("Running a test run")
            equil_steps = 100
            total_sample_steps = 1000
            write_step = 100
            initial_equil_steps = 100
            step_size = 0.001
            equil_step_size = 0.002
            initial_equil_step_size = 0.002
        
        #todo print out details of run

        for conf in self.system_confs:
            print(f"Running FEP for {conf['name']}")
            FEP = FepDriver()
            FEP.run_FEP(
                equilliberation_steps=equil_steps,
                total_sample_steps=total_sample_steps,
                write_step=write_step,
                lambda_0_equilliberation_steps=initial_equil_steps,
                step_size=step_size,
                equil_step_size=equil_step_size,
                initial_equil_step_size=initial_equil_step_size,
                Lambda=self.Lambda,
                systems=conf["systems"],
                topology=conf["topology"],
                temperature=conf["temperature"],
                initial_positions=conf["initial_positions"],
                run_folder=conf["run_folder"],
                data_folder=conf["data_folder"],
            )

            FEP.recalculate(interpolated_potential=True, force_contributions=True)

    def compute_energy_profiles(self, barrier, free_energy):
        reference_folder = self.system_confs[0]["data_folder"]
        target_folders = [conf["data_folder"] for conf in self.system_confs[1:]]
        dp = EvbDataProcessing()
        dp.compute(reference_folder, target_folders, barrier, free_energy)
        dp.print_results()

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

import json
import os
import numpy as np
import openmm as mm
import openmm.app as mmapp
import openmm.unit as mmunit

from .evbsystembuilder import EvbSystemBuilder
from .fepdriver import FepDriver


class EvbRunDriver:

    def __init__(self, reactant, product, data_folder=None, run_folder=None):
        self.reactant = reactant
        self.product = product
        self.data_folder = data_folder
        self.run_folder = run_folder

        self.M = 0
        self.N = 0
        self.CNT = False
        self.graphene = False
        self.E_field: list[float] = [0., 0., 0.]
        self.ion_count = 10

        self.Lambda: list[float]
        self.temperature: float
        self.constraints: dict | list[dict] | None = None
        self.pressure = 1
        self.solvent: str | None = None
        self.NPT: bool = False
        self.padding = 1.0  # in nm

        self.no_reactant: bool = False  #for debugging purposes

    def set_graphene(self, M=5, N=9, CNT=False):
        self.M = M
        self.N = N
        if CNT:
            self.CNT = CNT
        else:
            self.graphene = True

    def build_systems(self, Lambda):

        if self.constraints is None:
            self.constraints = []
        if isinstance(self.constraints, dict):
            self.constraints = [self.constraints]

        assert (Lambda[0] == 0 and Lambda[-1] == 1), f"Lambda must start at 0 and end at 1. Lambda = {Lambda}"
        assert np.all(np.diff(Lambda) > 0), f"Lambda must be monotonically increasing. Lambda = {Lambda}"
        Lambda = [round(lam, 3) for lam in Lambda]
        self.Lambda = Lambda

        folder_postfix = ""

        if self.solvent is None and self.graphene is False and self.CNT is False:
            folder_postfix = "_vacuum"
        else:
            if self.solvent is not None:
                folder_postfix += "_" + self.solvent
            if self.graphene:
                folder_postfix += "_graphene"
            if self.CNT:
                folder_postfix += "_CNT"

        if self.data_folder is None:
            self.data_folder = f"EVB_data{folder_postfix}"

        if self.run_folder is None:
            self.run_folder = f"{self.data_folder}/run"

        print(f"Saving files to {self.data_folder}")

        if not os.path.exists(self.data_folder):
            os.makedirs(self.data_folder)

        if not os.path.exists(self.run_folder):
            os.makedirs(self.run_folder)
        else:
            import shutil

            folder = f"./{self.run_folder}"
            for filename in os.listdir(folder):
                file_path = os.path.join(folder, filename)
                try:
                    if os.path.isfile(file_path) or os.path.islink(file_path):
                        os.unlink(file_path)
                    elif os.path.isdir(file_path):
                        shutil.rmtree(file_path)
                except Exception as e:
                    print("Failed to delete %s. Reason: %s" % (file_path, e))

        system_builder = EvbSystemBuilder()
        system_builder.constraints = self.constraints
        # self.system_builder = EvbSystemBuilder(self.reactant, self.product, l,
        #                                        temperature, pressure, NPT)

        self.systems, self.topology, self.initial_positions = system_builder.build_systems(
            reactant=self.reactant,
            product=self.product,
            Lambda=self.Lambda,
            NPT=self.NPT,
            temperature=self.temperature,
            pressure=self.pressure,
            solvent=self.solvent,
            padding=self.padding,
            CNT=self.CNT,
            Graphene=self.graphene,
            M=self.M,
            N=self.N,
            ion_count=self.ion_count,
            no_reactant=self.no_reactant,
            E_field=self.E_field,
        )
        self.save_systems_as_xml(self.run_folder)

        mmapp.PDBFile.writeFile(
            self.topology,
            self.initial_positions * 10,  # positions are handled in nanometers, but pdb's should be in angstroms
            open(f"{self.data_folder}/topology.pdb", "w"),
        )

        #todo add other options as well
        with open(f"{self.data_folder}/options.json", "w") as file:
            json.dump({"temperature": self.temperature, "Lambda": Lambda}, file, indent=4)

        self.temperature = self.temperature

    def load_initialisation(self, data_folder: str):
        self.data_folder = data_folder
        self.run_folder = f"{data_folder}/run"

        pdb = mmapp.PDBFile(f"{data_folder}/topology.pdb")
        self.topology = pdb.getTopology()
        self.initial_positions = pdb.getPositions(asNumpy=True).value_in_unit(mmunit.nanometers)

        with open(f"{data_folder}/options.json", "r") as file:
            options = json.load(file)
            self.temperature = options["temperature"]
            self.Lambda = options["Lambda"]

        self.systems = self.load_systems_from_xml(self.run_folder)

        print(f"Loaded systems, topology, initial positions, temperatue and Lambda from {data_folder}")

    def save_systems_as_xml(self, folder: str):
        for lam in self.Lambda:
            with open(f"{folder}/{lam:.3f}_sys.xml", mode="w", encoding="utf-8") as output:
                output.write(mm.XmlSerializer.serialize(self.systems[lam]))
        with open(f"{folder}/reactant.xml", mode="w", encoding="utf-8") as output:
            output.write(mm.XmlSerializer.serialize(self.systems["reactant"]))
        with open(f"{folder}/product.xml", mode="w", encoding="utf-8") as output:
            output.write(mm.XmlSerializer.serialize(self.systems["product"]))

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
        equil_steps,
        total_sample_steps,
        write_step,
        initial_equil_steps=0,
        step_size=0.001,
        equil_step_size=0.0005,
        initial_equil_step_size=0.0001,
        platform="CPU",
    ):
        self.FEP = FepDriver()

        self.FEP.run_FEP(equilliberation_steps=equil_steps,
                         total_sample_steps=total_sample_steps,
                         write_step=write_step,
                         lambda_0_equilliberation_steps=initial_equil_steps,
                         step_size=step_size,
                         Lambda=self.Lambda,
                         systems=self.systems,
                         topology=self.topology,
                         temperature=self.temperature,
                         initial_positions=self.initial_positions,
                         run_folder=self.run_folder,
                         data_folder=self.data_folder,
                         equil_step_size=equil_step_size,
                         initial_equil_step_size=initial_equil_step_size)

        self.FEP.recalculate(interpolated_potential=True, force_contributions=True)

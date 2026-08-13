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
from .reactionsystembuilder import ReactionSystemBuilder
from .evbfepdriver import EvbFepDriver
from .reaffbuilder import ReactionForceFieldBuilder
from .evbdataprocessing import EvbDataProcessing
from .solvationbuilder import SolvationBuilder
from .errorhandler import assert_msg_critical

try:
    import openmm as mm
    import openmm.app as mmapp
    import openmm.unit as mmunit
except ImportError:
    pass


class EvbDriver:

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
        self.lambda_vec: list[float] = None

        # The states of the reaction path, in order. A reactant -> product
        # reaction is simply the shortest path, with two states and one step, so
        # there is no separate single-step representation anywhere below.
        # reactive_bonds and bond_registry describe the pairs that react
        # anywhere along the path, which is what keeps the potential energy
        # surfaces E_1 ... E_N consistent across all of the individual FEPs.
        self.states: list[MMForceFieldGenerator] = []
        self.reactive_bonds: set = set()
        self.bond_registry: dict = {}
        # Per step, one entry each; a one-step reaction has one entry.
        self.forming_bonds: list[set] = []
        self.breaking_bonds: list[set] = []

        self.results = None
        self.system_confs: list[dict] = []
        self.mute_scf = True

        self.t_label = int(time.time())
        self.water_model = 'cspce'
        self.data_folder_override = None

        self.ffbuilder = ReactionForceFieldBuilder(ostream=self.ostream)
        self.dataprocessing = EvbDataProcessing(ostream=self.ostream)

    @property
    def reactant(self):
        """The first state of the path. Read-only: assign self.states."""
        return self.states[0] if self.states else None

    @property
    def product(self):
        """The last state of the path. Read-only: assign self.states."""
        return self.states[-1] if self.states else None

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

    def build_ff_from_molecules(self, reactant: Molecule | list[Molecule],
                                product: Molecule | list[Molecule], **kwargs):
        """_summary_

        Args:
            reactant (Molecule | list[Molecule]): The reactant molecule or a list of reactant molecules.
            product (Molecule | list[Molecule]): The product molecule or a list of product molecules.
            reactant_partial_charges (list[float], list[list[float]]): Partial charges for the reactant. Will be calculated if not provided. Defaults to None.
            product_partial_charges (list[float], list[list[float]]): Partial charges for the product. Will be calculated if not provided. Defaults to None.
            reparameterize (bool): If True, reparameterizes unknown force constants with the Seminario method. Defaults to True
            reactant_hessians (np.ndarray, list[np.ndarray]): Hessians for the reactant for the Seminario method. Will be calculated if not provided. Defaults to None.
            product_hessians (np.ndarray, list[np.ndarray]): Hessians for the product for the Seminario method. Will be calculated if not provided. Defaults to None.
            mm_opt_constrain_bonds (list[tuple[int, int]]): Bonds to constrain during MM optimization.
            reactant_total_multiplicity (int): Total multiplicity for the reactant to override calculated value. Defaults to -1.
            product_total_multiplicity (int): Total multiplicity for the product to override calculated value. Defaults to -1.
            breaking_bonds (list[tuple[int, int]]): List of bond(s) that is forced to break and is not allowed to recombine over the reaction. Defaults to None.
            mute_ff_scf (bool): If True, mutes SCF output from RESP calculations. Has no effect if mute_ff_build is True. Defaults to True.
            optimize_mol (bool): If True, does an xtb optimization of every provided molecule object before reparameterisation. Defaults to False.
            optimize_ff (bool): If True, does an mm optimization of the combined reactant and product after reparameterisation. Defaults to True.
        """

        # A reactant -> product reaction is the shortest reaction path, so this
        # only translates the two-state keyword names into the per-state and
        # per-step lists the general builder takes.
        per_state = {
            'partial_charges':
                ('reactant_partial_charges', 'product_partial_charges'),
            'hessians': ('reactant_hessians', 'product_hessians'),
            'total_multiplicities':
                ('reactant_total_multiplicity', 'product_total_multiplicity'),
        }
        per_step = {
            'forced_breaking_bonds': 'forced_breaking_bonds',
            'forced_forming_bonds': 'forced_forming_bonds',
            'mappings': 'product_mapping',
        }

        arguments = {}
        for name, (rea_key, pro_key) in per_state.items():
            if rea_key in kwargs or pro_key in kwargs:
                default = -1 if name == 'total_multiplicities' else None
                arguments[name] = [
                    kwargs.pop(rea_key, default),
                    kwargs.pop(pro_key, default),
                ]
        for name, key in per_step.items():
            if key in kwargs:
                arguments[name] = [kwargs.pop(key)]

        results = self.build_ffs_from_molecules([reactant, product],
                                                **arguments,
                                                **kwargs)
        self.product_mapping = results['mappings'][0]
        return results

    def build_ffs_from_molecules(self, states: list, **kwargs):
        """Build matched force fields for a reaction path of arbitrary length.

        Every state is brought into one and the same atom ordering, and the
        bonds that react anywhere along the path are recorded so that each
        state's potential energy surface is described identically in every step
        it takes part in.

        Args:
            states (list): The states of the reaction path, in order. Each entry
                is a Molecule or a list of Molecule fragments.
            partial_charges (list): One entry per state; each None, a flat list
                of charges, or one list of charges per fragment of that state.
            hessians (list): One entry per state, for the Seminario method.
            total_multiplicities (list): One entry per state; -1 to calculate.
            forced_breaking_bonds (list): One entry per step (one fewer than the
                number of states), each a set of 1-indexed atom pairs.
            forced_forming_bonds (list): One entry per step.
            mappings (list): One entry per step; a 1-indexed mapping from that
                step's state to the next to skip reaction matching.
        """

        self.ffbuilder.water_model = self.water_model

        results = self.ffbuilder.build_many_force_fields(states=states,
                                                         **kwargs)

        self.states = results['states']
        self.forming_bonds = results['forming_bonds']
        self.breaking_bonds = results['breaking_bonds']
        self.reactive_bonds = results['reactive_bonds']
        self.bond_registry = results['registry']
        self.state_fragments = results['fragments']
        self.state_mappings = results['mappings']

        self.ostream.print_info(
            f"Built force fields for {len(self.states)} states "
            f"({len(self.states) - 1} steps) with "
            f"{len(self.reactive_bonds)} reacting bonds in total.")
        self.ostream.flush()
        return results

    def build_systems(
        self,
        configurations: list[str] | list[dict],  # type: ignore
        Lambda: list[float] | np.ndarray = None,
        constraints: dict | list[dict] | None = None,
    ):
        """Build OpenMM systems for the given configurations with interpolated forcefields for each lambda value. Saves the systems as xml files, the topology as a cif file and the options as a json file to the disk.

        Args:
            configurations (list[str] | list[dict]): The given configurations for which to perform an FEP. The first configuration will be regarded as the reference configuration.
            Lambda (list[float] | np.ndarray): The Lambda vector to be used for the FEP. Should start with 0, end with 1 and be monotonically increasing.
                Defaults to None, in which case default values will be assigned depending on if debugging is enabled or not.
                If a string is given, the return value of default_system_configurations() will be used. See this function for default configurations.
            constraints (dict | list[dict] | None, optional): Dictionary of harmonic bond, angle or (improper) torsion forces to apply over in every FEP frame. Defaults to None.
        """

        assert_msg_critical('openmm' in sys.modules,
                            'openmm is required for EvbDriver.')

        # System building and the associated folder / file creation happen only
        # on the master rank. In async-reporter mode (nodes > 1) the reporter
        # worker reconstructs the systems from the master's folder in run_FEP,
        # so only one output folder is ever created and all ranks agree on it.
        if self.rank != mpi_master():
            return

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
        self.lambda_vec = Lambda

        # Per configuration
        for conf in self.configurations:
            # create folders,
            if self.data_folder_override is not None:
                data_folder = self.data_folder_override
            else:
                data_folder = f"{conf['name']}_{self.t_label}"
                while Path(data_folder).exists():
                    self.t_label += 1
                    data_folder = f"{conf['name']}_{self.t_label}"

            run_folder = str(Path(data_folder) / "run")
            conf["data_folder"] = data_folder
            conf["run_folder"] = run_folder

            cwd = Path.cwd()
            data_folder_path = cwd / data_folder
            run_folder_path = cwd / run_folder

            data_folder_path.mkdir(parents=True)
            run_folder_path.mkdir(parents=True)

            # The ends of the path keep the reactant / product names, so a
            # one-step reaction writes exactly the files it always has.
            self.states[0].molecule.write_xyz_file(
                str(data_folder_path / "reactant_struct.xyz"))
            self.states[-1].molecule.write_xyz_file(
                str(data_folder_path / "product_struct.xyz"))

            MMForceFieldGenerator.save_forcefield_as_json(
                self.states[0],
                str(data_folder_path / "reactant_ff_data.json"))
            MMForceFieldGenerator.save_forcefield_as_json(
                self.states[-1],
                str(data_folder_path / "product_ff_data.json"))

            if len(self.states) > 2:
                for k, state in enumerate(self.states):
                    state.molecule.write_xyz_file(
                        str(data_folder_path / f"state_{k}_struct.xyz"))
                    MMForceFieldGenerator.save_forcefield_as_json(
                        state,
                        str(data_folder_path / f"state_{k}_ff_data.json"))

            if conf.get('solvent', None) is None and conf.get('pressure',
                                                              -1) > 0:
                self.ostream.print_warning(
                    f"A pressure is defined for {conf['name']}, but no solvent is defined. Removing pressure definition."
                )
                conf.pop("pressure")
            # build the system
            system_builder = ReactionSystemBuilder(ostream=self.ostream)
            system_builder.water_model = self.water_model
            self.ostream.print_blank()
            self.ostream.print_header(f"Building systems for {conf['name']}")
            self.ostream.flush()

            topology, initial_positions = self._build_step_confs(
                conf,
                system_builder,
                constraints,
                data_folder,
            )

            top_path = cwd / data_folder / "topology.cif"

            mmapp.PDBxFile.writeFile(
                topology,
                initial_positions,  # positions are handled in nanometers, but pdb's should be in angstroms
                open(top_path, "w"),
            )

            dump_conf = copy.copy(conf)
            dump_conf.pop('systems', None)
            dump_conf.pop('topology', None)
            dump_conf.pop('initial_positions', None)
            if 'steps' in dump_conf:
                # The per-step dicts carry the same non-serialisable entries.
                dump_conf['steps'] = [{
                    key: value
                    for key, value in step.items()
                    if key not in ('systems', 'topology', 'initial_positions')
                } for step in dump_conf['steps']]
            self.update_options_json(dump_conf, conf)
            self.update_options_json(
                {
                    "Lambda": Lambda,
                },
                conf,
            )

        self.system_confs = configurations

        self.ostream.flush()

    def _build_step_confs(self, conf, system_builder, constraints,
                          data_folder):
        """Build the systems of one configuration, step by step.

        The environment is built once, so all steps share one topology, one
        particle ordering and one run folder; each step then gets its own output
        folder and its own configuration dict, which is what lets the FEP driver
        run a step without knowing that it is part of a longer path.
        """

        cwd = Path.cwd()
        n_states = len(self.states)
        single_step = n_states == 2

        step_systems, pes_systems, topology, initial_positions = (
            system_builder.build_multi_step_systems(
                states=self.states,
                lambda_vec=self.lambda_vec,
                configuration=conf,
                constraints=constraints,
                reactive_bonds=self.reactive_bonds,
                bond_registry=self.bond_registry,
            ))

        self.ostream.print_info(f"Saving files to {cwd / data_folder}")
        self.ostream.flush()

        # A one-step reaction keeps the names it always had: no step prefix on
        # the lambda systems, and 'reactant' / 'product' for the two surfaces.
        # A longer path needs both, because the lambda values of different steps
        # would otherwise collide in the one shared run folder.
        def prefix_of(step):
            return "" if single_step else f"step{step}_"

        if single_step:
            pes_by_name = {'reactant': pes_systems[0], 'product': pes_systems[1]}
        else:
            pes_by_name = {
                f"state_{k}": system for k, system in pes_systems.items()
            }

        for step, systems in enumerate(step_systems):
            system_builder.save_systems_as_xml(systems,
                                               conf["run_folder"],
                                               prefix=prefix_of(step))
        system_builder.save_systems_as_xml(pes_by_name, conf["run_folder"])

        steps = []
        for step, systems in enumerate(step_systems):
            if single_step:
                step_folder = data_folder
            else:
                step_folder = str(Path(data_folder) / f"step_{step:02d}")
                (cwd / step_folder).mkdir(parents=True, exist_ok=True)

            step_conf = copy.copy(conf)
            step_conf['step'] = step
            step_conf['n_states'] = n_states
            step_conf['state_left'] = step
            step_conf['state_right'] = step + 1
            step_conf['data_folder'] = step_folder
            step_conf['run_folder'] = conf["run_folder"]
            step_conf['system_prefix'] = prefix_of(step)
            step_conf['forming_bonds'] = sorted(self.forming_bonds[step])
            step_conf['breaking_bonds'] = sorted(self.breaking_bonds[step])
            # Every state's PES system is visible to every step, so a frame
            # sampled anywhere along the path can be scored on all of E_1..E_N.
            step_conf['systems'] = {**systems, **pes_by_name}
            step_conf['topology'] = topology
            step_conf.pop('steps', None)
            if not single_step:
                step_conf['name'] = f"{conf['name']}_step_{step + 1}"
            steps.append(step_conf)

        # Only the first step starts from the built geometry; every later step
        # is seeded with the equilibrated lambda=1 state of the one before it.
        steps[0]['initial_positions'] = initial_positions

        conf['n_states'] = n_states
        conf['reactive_bonds'] = sorted(
            [list(pair) for pair in self.reactive_bonds])
        conf['steps'] = steps
        conf['topology'] = topology
        conf['initial_positions'] = initial_positions
        conf['systems'] = steps[0]['systems'] if single_step else None
        conf['forming_bonds'] = [sorted(bonds) for bonds in self.forming_bonds]
        conf['breaking_bonds'] = [sorted(bonds) for bonds in self.breaking_bonds]

        return topology, initial_positions

    def load_initialisation(self,
                            data_folder: str,
                            load_systems=False,
                            load_top=False,
                            restart_pdb: str = None,
                            rename_data=False):
        """Load a configuration from a data folder for which the systems have already been generated, such that an FEP can be performed.
        The topology, initial positions, temperature and Lambda vector will be loaded from the data folder.

        Args:
            data_folder (str): The folder to load the data from
            name (str): The name of the configuration. Can be arbitrary, but should be unique.
            load_systems (bool, optional): If set to true, the systems will be loaded from the xml files. Used for debugging. Defaults to False.
            load_top (bool, optional): If set to true, the topology will be loaded from the cif file. Used for debugging. Defaults to False.
            restart_pdb (str, optional): If given, the topology and initial positions will be loaded from the given pdb file. Used for restarting an FEP. Defaults to None.
            rename_data (bool, optional): If set to true, any csv or xtc files in the data folder will be renamed with the prefix OLD_ to avoid overwriting them during the new FEP run. Defaults to False.
        """

        assert_msg_critical('openmm' in sys.modules,
                            'openmm is required for EvbDriver.')

        # Loading / folder preparation happens only on the master rank; in
        # async-reporter mode the reporter worker reconstructs what it needs
        # from the shared folder in run_FEP.
        if self.rank != mpi_master():
            return

        options_path = Path(data_folder) / "options.json"
        with options_path.open("r") as file:
            conf = json.load(file)

        existing_names = {
            other_conf["name"]
            for other_conf in self.system_confs
        }

        if conf['name'] in existing_names:
            name = conf['name']
            i = 1
            unique_name = f"{name}_{i}"
            while unique_name in existing_names:
                i += 1
                unique_name = f"{name}_{i}"
            self.ostream.print_warning(
                f"Name {name} is already in use. Using {unique_name} instead.")
            self.ostream.print_info(
                f"The name can be edited in options.json before loading")
            name = unique_name
            conf["name"] = name

        if self.lambda_vec != conf["Lambda"] and self.lambda_vec is not None:
            self.ostream.print_warning(
                f"Lambda vector in {options_path} does not match the current Lambda vector. Overwriting current Lambda vector with the one from the file."
            )
        self.temperature = conf['temperature']
        self.lambda_vec = conf['Lambda']
        conf["data_folder"] = str(data_folder)
        conf["run_folder"] = str(Path(data_folder) / "run")

        if load_systems:
            systems = self._load_step_systems(conf)
            conf["systems"] = conf['steps'][0]['systems'] if conf.get(
                'steps') else systems
        else:
            systems = []

        if restart_pdb is not None:
            try:
                pdb = mmapp.PDBxFile(str(Path(data_folder) / restart_pdb))
            except Exception as e:
                raise ValueError(
                    f"Could not load restart pdb file {restart_pdb} from {data_folder}. Error: {e}"
                )
            self.ostream.print_info(
                f"Loading topology and initial positions from restart pdb file {restart_pdb} from {data_folder}."
            )
            self.ostream.print_info(
                "Turning on skip_initial_equil, so the FEP will start directly with sampling instead of equilibration."
            )
            self.ostream.flush()
            conf["topology"] = pdb.getTopology()
            conf["initial_positions"] = pdb.getPositions(
                asNumpy=True).value_in_unit(mmunit.nanometers)
            conf['skip_initial_equil'] = True
        elif load_top:
            pdb = mmapp.PDBxFile(str(Path(data_folder) / "topology.cif"))
            conf["topology"] = pdb.getTopology()
            conf["initial_positions"] = pdb.getPositions(
                asNumpy=True).value_in_unit(mmunit.nanometers)

        # The steps recorded in options.json only carry names and folders; give
        # them back the topology and the starting configuration that run_FEP
        # needs, and point each one at its own output folder.
        self._restore_step_confs(conf, data_folder, rename_data)

        self.system_confs.append(conf)
        self.ostream.print_info(
            f"Initialised configuration with {len(systems)} systems, temperatue {self.temperature} and Lambda vector {self.lambda_vec} from {data_folder}"
        )
        self.ostream.print_info(
            f"Current configurations: {[conf['name'] for conf in self.system_confs]}"
        )

        self.ostream.flush()

    def _load_step_systems(self, conf):
        """Load every step's lambda systems plus the shared PES systems.

        The steps share one run folder, so a step's lambda systems are selected
        by the prefix they were saved under. A one-step reaction was saved
        without a prefix and with the reactant / product names, exactly as it
        always has been.
        """

        sysbuilder = ReactionSystemBuilder(ostream=self.ostream)
        run_folder = conf['run_folder']
        steps = conf.get('steps') or []
        single_step = len(steps) <= 1

        if single_step:
            # One flat load: lambda systems and the two surfaces sit together
            # under their own names.
            systems = sysbuilder.load_systems_from_xml(run_folder)
        else:
            pes_systems = {
                f"state_{name}": system
                for name, system in sysbuilder.load_systems_from_xml(
                    run_folder, prefix="state_", parse_lambda=False).items()
            }
            systems = pes_systems

        for step_conf in steps:
            if single_step:
                step_conf['systems'] = systems
            else:
                step_conf['systems'] = {
                    **sysbuilder.load_systems_from_xml(
                        run_folder, prefix=step_conf['system_prefix']),
                    **systems,
                }
        return systems

    def _restore_step_confs(self, conf, data_folder, rename_data):
        """Give the steps read back from options.json what run_FEP needs."""

        steps = conf.get('steps')
        if not steps:
            return

        data_folder = Path(data_folder)
        for step, step_conf in enumerate(steps):
            step_conf['run_folder'] = conf['run_folder']
            if len(steps) == 1:
                step_conf['data_folder'] = str(data_folder)
            else:
                step_conf['data_folder'] = str(data_folder /
                                               f"step_{step:02d}")
            if 'topology' in conf:
                step_conf['topology'] = conf['topology']
            if step > 0:
                # Seeded from the preceding step's endpoint, written when that
                # step runs; only step 0 starts from a stored configuration.
                previous = (Path(step_conf['data_folder']).parent /
                            f"step_{step - 1:02d}" / "final_state.xml")
                if previous.is_file():
                    step_conf['initial_state'] = str(previous)
                    step_conf['skip_initial_equil'] = True
            elif 'initial_positions' in conf:
                step_conf['initial_positions'] = conf['initial_positions']
            if 'skip_initial_equil' in conf and step == 0:
                step_conf['skip_initial_equil'] = conf['skip_initial_equil']

            if rename_data:
                self._rename_existing_data(Path(step_conf['data_folder']))

    def _rename_existing_data(self, folder: Path):
        """Move any csv / xtc left over from an earlier run out of the way."""
        if not folder.is_dir():
            return
        if not (any(folder.glob("*.csv")) or any(folder.glob("*.xtc"))):
            return
        self.ostream.print_warning(
            f"Found csv or xtc files in {folder}. These might be from a "
            "previous FEP run using the same folder. The files will be renamed "
            "with the prefix OLD_ to avoid overwriting them during the new FEP "
            "run.")
        for file in folder.iterdir():
            if file.is_file() and (file.name.endswith(".csv")
                                   or file.name.endswith(".xtc")):
                new_name = file.parent / f"OLD_{file.name}"
                file.rename(new_name)
                self.ostream.print_info(
                    f"Renamed {file.name} to {new_name.name}")

    def load_results(
        self,
        data_folder: str,
    ):
        assert_msg_critical('openmm' in sys.modules,
                            'openmm is required for EvbDriver.')

        # Loading / folder preparation happens only on the master rank; in
        # async-reporter mode the reporter worker reconstructs what it needs
        # from the shared folder in run_FEP.
        if self.rank != mpi_master():
            return

        options_path = Path(data_folder) / "options.json"
        with options_path.open("r") as file:
            conf = json.load(file)

        if self.lambda_vec != conf["Lambda"] and self.lambda_vec is not None:
            self.ostream.print_warning(
                f"Lambda vector in {options_path} does not match the current Lambda vector. Overwriting current Lambda vector with the one from the file."
            )
        self.temperature = conf['temperature']
        self.lambda_vec = conf['Lambda']
        conf["data_folder"] = str(data_folder)
        conf["run_folder"] = str(Path(data_folder) / "run")
        self.system_confs.append(conf)
        self.ostream.print_info(
            f"Loaded configuration with temperatue {self.temperature} and Lambda vector {self.lambda_vec} from {data_folder}"
        )
        self.ostream.print_info(
            f"Current configurations: {[conf['name'] for conf in self.system_confs]}"
        )

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

        # In async-reporter mode only the master built the systems and owns the
        # output folder. Broadcast the lightweight per-configuration metadata
        # (everything except the non-picklable OpenMM systems / topology /
        # positions) plus the lambda vector so every rank iterates the same
        # configurations and agrees on the folder the reporter worker writes to.
        if self.nodes > 1:
            if self.rank == mpi_master():
                meta = [{
                    key: value
                    for key, value in conf.items()
                    if key not in ("systems", "topology", "initial_positions")
                } for conf in self.system_confs]
            else:
                meta = None
            meta = self.comm.bcast(meta, root=mpi_master())
            self.lambda_vec = self.comm.bcast(self.lambda_vec,
                                              root=mpi_master())
            if self.rank != mpi_master():
                self.system_confs = meta

        for conf in self.system_confs:
            self.ostream.print_blank()
            self.ostream.print_header(f"Running FEP for {conf['name']}")
            self.ostream.flush()

            # Always a list of steps; a reactant -> product reaction has one.
            step_confs = conf.get('steps', [conf])
            for index, step_conf in enumerate(step_confs):
                self.ostream.print_blank()
                if len(step_confs) > 1:
                    self.ostream.print_header(
                        f"Step {step_conf['step'] + 1} of "
                        f"{len(step_confs)}: state "
                        f"{step_conf['state_left'] + 1} -> state "
                        f"{step_conf['state_right'] + 1}")
                    self.ostream.flush()

                FEP = EvbFepDriver(ostream=self.ostream)
                FEP.run_replicas(
                    Lambda=self.lambda_vec,
                    configuration=step_conf,
                    platform=platform,
                    platform_properties=platform_properties,
                )

                # Each step is sampled on its own, seeded by the equilibrated
                # lambda=1 configuration the previous step ended its first
                # forward sweep on. That configuration is a valid starting point
                # because the two Hamiltonians are identical there by
                # construction: lambda=1 of one step and lambda=0 of the next
                # both describe the state they share.
                if index + 1 < len(step_confs):
                    final_state = (Path(step_conf['data_folder']) /
                                   "final_state.xml")
                    assert_msg_critical(
                        final_state.is_file(),
                        f"Step {step_conf['step'] + 1} did not write "
                        f"{final_state}, so the next step has nothing to start "
                        "from.")
                    next_conf = step_confs[index + 1]
                    next_conf['initial_state'] = str(final_state)
                    next_conf['skip_initial_equil'] = True
                    self.ostream.print_info(
                        f"Seeding step {next_conf['step'] + 1} with the "
                        f"lambda=1 configuration from {final_state}")
                    self.ostream.flush()

    def compute_force_groups(self,
                             platform=None,
                             platform_properties=None,
                             decompose_nb=None,
                             decompose_bonded=False):
        """Compute reactant/product OpenMM force-group energy decompositions
        for every configuration in self.system_confs, from the trajectory
        already sampled by run_FEP() (or reloaded via
        load_initialisation(..., load_systems=True, load_top=True)). This is
        a standalone post-processing step: it is never run automatically and
        does not affect sampling.

        Writes ForceGroups_rea.csv / ForceGroups_pro.csv into each
        configuration's data_folder.

        The two finer decompositions are exclusively options of this method (no
        other code path touches them). When requested, the required decomposition systems are
        regenerated on the fly, saved to the run folder as unique {name}_sys.xml
        files (so a later load_initialisation picks them up), and then evaluated
        against the existing trajectory - no re-sampling:

        - decompose_nb: list of reactant-atom-index groups; writes
          NB_decompositions.csv (solute-solvent Coulomb / LJ contributions).
          Requires a solvated configuration.
        - decompose_bonded: writes the per-bonded-term decomposition
          (bonded_E1_decomp.csv / bonded_E2_decomp.csv / bonded_params.csv),
          which uses the configuration's forming_bonds/breaking_bonds.

        Args:
            platform (str, optional): OpenMM platform name to use for the
                recalculation (e.g. "CUDA", "CPU"). Defaults to None, which
                lets OpenMM auto-select the fastest available platform.
            platform_properties (dict, optional): OpenMM platform properties.
            decompose_nb (list, optional): reactant-atom-index groups to
                decompose the solute-solvent nonbonded energy over. Defaults to
                None (no nonbonded decomposition).
            decompose_bonded (bool, optional): also compute the per-bonded-term
                energy decomposition. Defaults to False.
        """
        # Post-processing reads the trajectory/Energies.csv that only the
        # master rank produced; the reporter / helper ranks have nothing to
        # recompute.
        if self.rank != mpi_master():
            return

        assert_msg_critical(
            len(self.states) <= 2,
            "compute_force_groups is not available for a multi-step reaction "
            "path: it writes one reactant/product decomposition per "
            "configuration, which has no meaning once a configuration holds "
            "several steps and several states.")

        for conf in self.system_confs:
            self.ostream.print_blank()
            self.ostream.print_header(
                f"Computing force-group contributions for {conf['name']}")
            self.ostream.flush()

            FEP = EvbFepDriver(ostream=self.ostream)
            FEP.compute_force_group_contributions(
                configuration=conf,
                platform=platform,
                platform_properties=platform_properties,
                decompose_nb=decompose_nb,
                decompose_bonded=decompose_bonded,
            )

    def update_options_json(self, dict, conf):

        cwd = Path.cwd()
        path = cwd / conf["data_folder"] / "options.json"
        if not path.exists():
            with path.open("w") as file:
                json.dump(dict, file, indent=4)
        else:
            with path.open("r") as file:
                options = json.load(file)
            options.update(dict)
            with path.open("w") as file:
                json.dump(options, file, indent=4)

    def compute_energy_profiles(self,
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
                                results_filename="results"):
        """Compute the EVB energy profiles using the FEP results, print the results and save them to an h5 file

        Args:
            barrier (float): the reaction barrier in kJ/mol of the reference system
            free_energy (float): the reaction free energy in kJ/mol of the reference system
            lambda_sub_sample (int, optional): Factor with which the lambda vector will be subsampled. Setting this to two will discard every other lambda frame. Defaults to 1.
            lambda_sub_sample_ends (bool, optional): If set to False, the lambda frames up to 0.1 and from 0.9 will not be subsampled. Defaults to False.
            time_sub_sample (int, optional): Factor with which the time vector will be subsampled. Setting this to two will discard every other snapshot. Defaults to 1.
        """
        # Post-processing reads the FEP output files that only the master rank
        # produced; the reporter / helper ranks have no data to analyse.
        if self.rank != mpi_master():
            return

        results = self._load_output_from_folders(lambda_sub_sample,
                                                 lambda_sub_sample_ends,
                                                 time_sub_sample)
        self.ostream.flush()

        multi_step = any(
            "steps" in result
            for result in results["configuration_results"].values())
        assert_msg_critical(
            not multi_step,
            "compute_energy_profiles cannot yet fit a reaction path of more "
            "than two states: EvbDataProcessing solves one two-state EVB "
            "Hamiltonian per configuration, and a path has one per step. The "
            "sampled output itself loads fine - call "
            "_load_output_from_folders directly to get the per-step energies.")

        if alpha is not None:
            self.dataprocessing.alpha = alpha
        if H12 is not None:
            self.dataprocessing.H12 = H12
        if alpha_guess is not None:
            self.dataprocessing.alpha_guess = alpha_guess
        if H12_guess is not None:
            self.dataprocessing.H12_guess = H12_guess
        if dE_range is not None:
            self.dataprocessing.coordinate_bins = np.linspace(
                dE_range[0], dE_range[1], 200)

        results = self.dataprocessing.compute(results, barrier, free_energy)
        self.results = results
        self.print_results()
        self._save_dict_as_h5(results, results_filename)
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

    def print_uncertainties(self, results: dict = None, file_name: str = None):
        """Print total, per-replica/direction, and hysteresis uncertainty diagnostics. Uses the provided dictionary first, then tries to load it from the disk, and last it uses the results attribute of this object.

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
        dp.print_uncertainties(self.results, self.ostream)
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
            if file_name is not None:
                results = self._load_dict_from_h5(file_name)
            else:
                assert self.results is not None, "No results known, and none provided"
                results = self.results
        dp = EvbDataProcessing()
        dp.plot_results(results, **kwargs)
        self.ostream.flush()

    def _load_output_from_folders(
        self,
        lambda_sub_sample,
        lambda_sub_sample_ends,
        time_sub_sample,
    ) -> dict:
        """Read every configuration's sampled output back off disk.

        Each step of a reaction path writes into its own folder, so the files
        are read per step; options.json describes the whole configuration and is
        read from the configuration folder. A one-step reaction keeps the flat
        result shape it has always had, since that is what the analysis
        consumes; a longer path nests its steps under 'steps', because there is
        then no single unambiguous set of energies to put at the top level.
        """

        results = {}
        cwd = Path.cwd()

        common_results = []
        specific_results = {}
        for conf in self.system_confs:
            options_file = str(cwd / conf["data_folder"] / "options.json")
            # Always a list of steps; a reactant -> product reaction has one.
            step_confs = conf.get("steps") or [conf]

            step_results = []
            for step_conf in step_confs:
                folder = cwd / step_conf["data_folder"]
                specific, common = self._load_output_files(
                    str(folder / "Energies.csv"),
                    str(folder / "Data_combined.csv"),
                    options_file,
                    str(folder / "ForceGroups.csv"),
                    str(folder / "ForceGroups_rea.csv"),
                    str(folder / "ForceGroups_pro.csv"),
                    str(folder / "NB_decompositions.csv"),
                    lambda_sub_sample,
                    lambda_sub_sample_ends,
                    time_sub_sample,
                    state_left=step_conf.get("state_left"),
                    state_right=step_conf.get("state_right"),
                )
                step_results.append(specific)
                common_results.append(common)

            if len(step_results) == 1:
                specific_results[conf["name"]] = step_results[0]
            else:
                specific_results[conf["name"]] = {
                    "steps": step_results,
                    "n_states": conf.get("n_states", len(step_results) + 1),
                    "options": step_results[0]["options"],
                    "Temp_set": step_results[0]["Temp_set"],
                }

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
        state_left=None,
        state_right=None,
    ):
        """Read one step's sampled output.

        Energies.csv is read by column name rather than by position, so a file
        holding any number of potential energy surfaces is understood: a
        reactant -> product run names its two 'reactant' / 'product', a longer
        path names them 'state_0' ... 'state_N-1' and adds a 'step' column.
        Every surface is returned in E_pes / E_pes_names, and E1_pes / E2_pes
        are the two states flanking this step, which is what the two-state
        analysis works with.
        """
        with Path(options_file).open("r") as file:
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
        columns = self._parse_energies_header(E_file, E_data.shape[0])
        fg_data = []
        rea_fg_data = []
        pro_fg_data = []

        l_sub_indices = np.where([lf in Lambda for lf in E_data[0]])[0]

        sub_indices = l_sub_indices[::time_sub_sample]

        def column(name, default=None):
            index = columns.get(name)
            return default if index is None else E_data[index, sub_indices]

        Lambda_frame = E_data[columns["Lambda"], sub_indices]

        # Every potential energy surface in the file, in the order it was
        # written. On a reaction path these are the whole path's states, so a
        # frame sampled in one step carries the energies of all of them.
        pes_names = [
            name[:-len(" PES")] for name in columns if name.endswith(" PES")
        ]
        E_pes = np.array([column(f"{name} PES") for name in pes_names])

        # The two states flanking this step. Named per state on a reaction path,
        # 'reactant' / 'product' for a plain reactant -> product run; and with
        # pes_states_to_report='adjacent' only these two are in the file at all,
        # so they are picked by name rather than by position.
        left_name, right_name = pes_names[0], pes_names[-1]
        if state_left is not None and f"state_{state_left}" in pes_names:
            left_name = f"state_{state_left}"
        if state_right is not None and f"state_{state_right}" in pes_names:
            right_name = f"state_{state_right}"

        E1_pes = column(f"{left_name} PES")
        E2_pes = column(f"{right_name} PES")

        integration = [name for name in columns if name.endswith(" integration")]
        E1_int = column(integration[0])
        E2_int = column(integration[-1])

        # 'Em' in files written by EvbReporter.format_energies_row; 'E_m_pes' in
        # older ones, which also carry an E_m_int column that nothing reads.
        E_m_pes = column("Em") if "Em" in columns else column("E_m_pes")

        # replica / direction were added later; older Energies.csv files (e.g.
        # canned test fixtures) don't have them. direction: 0 = forward sweep
        # (l: 0 -> 1), 1 = backward (l: 1 -> 0), see EvbFepDriver.run_replicas.
        replica_frame = column("replica")
        direction_frame = column("direction")

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

        # Only added when it says something E1_pes / E2_pes do not: with two
        # surfaces in the file those are already both of them, and repeating
        # them here would double the arrays in every saved h5.
        if len(pes_names) > 2:
            specific_result.update({
                "E_pes": E_pes,
                # Comma-joined rather than a list, for the same reason as
                # E1_fg_names below: a numpy unicode array has no native h5
                # representation.
                "E_pes_names": ",".join(pes_names),
            })

        if "step" in columns:
            specific_result["Step_frame"] = column("step").astype(int)
            specific_result["step"] = int(specific_result["Step_frame"][0])

        if replica_frame is not None:
            specific_result.update({
                "Replica_frame": replica_frame.astype(int),
                "Direction_frame": direction_frame.astype(int),
            })

        if fg_file is not None and Path(fg_file).is_file():
            fg_data = np.loadtxt(fg_file, skiprows=1, delimiter=',').T
            specific_result.update({"E_m_fg": fg_data})
        if fg_rea_file is not None and Path(fg_rea_file).is_file():
            rea_fg_data = np.loadtxt(fg_rea_file, skiprows=1, delimiter=',').T
            specific_result.update({
                "E1_fg":
                rea_fg_data,
                # Comma-joined, not a list: the generic h5 dict saver
                # (_save_dict_as_h5) can't natively store a numpy unicode
                # array (h5py has no conversion path for '<U' dtype) and
                # would silently fall back to a repr() string on save; a
                # plain str round-trips through h5 either way, so this is the
                # one representation that works identically whether read
                # fresh this session or reloaded from an h5 file.
                "E1_fg_names":
                ",".join(self._parse_force_group_header(fg_rea_file)),
            })
        if fg_pro_file is not None and Path(fg_pro_file).is_file():
            pro_fg_data = np.loadtxt(fg_pro_file, skiprows=1, delimiter=',').T
            specific_result.update({
                "E2_fg":
                pro_fg_data,
                "E2_fg_names":
                ",".join(self._parse_force_group_header(fg_pro_file)),
            })
        if decomp_file is not None and Path(decomp_file).exists():
            decomp_data = np.loadtxt(decomp_file, skiprows=1, delimiter=',').T
            decomp_rea = decomp_data[decomp_data.shape[0] // 2:, :]
            decomp_pro = decomp_data[:decomp_data.shape[0] // 2, :]
            with Path(decomp_file).open("r") as file:
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
    def _parse_energies_header(path, n_columns):
        """Column name -> index for an Energies.csv.

        Reading by name is what lets one loader handle a file with any number
        of potential energy surfaces, and keeps files written before a column
        was added readable. A file whose header does not match its width is
        rejected rather than silently mis-indexed.
        """
        with Path(path).open("r") as file:
            header = [cell.strip() for cell in file.readline().split(",")]

        assert_msg_critical(
            len(header) == n_columns,
            f"{path} has {n_columns} columns but a header naming "
            f"{len(header)}: {header}")

        columns = {name: index for index, name in enumerate(header)}
        assert_msg_critical(
            "Lambda" in columns,
            f"{path} has no Lambda column; header: {header}")
        assert_msg_critical(
            any(name.endswith(" PES") for name in columns),
            f"{path} names no potential energy surface; header: {header}")
        return columns

    @staticmethod
    def _parse_force_group_header(path):
        """Column names, in order, from a ForceGroups_rea/pro.csv header
        line ("NAME(value), NAME(value), ..."). The numeric values are
        per-run (see ReactionSystemBuilder.save_systems_as_xml /
        force_group_name.json) and not meaningful across runs, so only the
        names - which are what E1_fg/E2_fg should always be indexed by - are
        kept here.
        """
        with Path(path).open("r") as file:
            header = file.readline().strip()
        return [cell.strip().split("(")[0] for cell in header.split(",")]

    @staticmethod
    def _load_dict_from_h5(file):
        """Load a dictionary from an h5 file

        Args:
            file (path): The file to load the results from.

        Returns:
            dict: Dictionary with the results
        """

        def decode_bytes(v):
            # h5py reads back a scalar string dataset/attr (saved from a
            # plain Python str, e.g. group[k] = v in _save_dict_as_h5) as
            # bytes, not str - silently changing the value's type across a
            # save/load round trip. Undo that so callers get back exactly
            # what they saved.
            if isinstance(v, bytes):
                return v.decode('utf-8')
            if isinstance(v, np.ndarray) and v.dtype == object:
                # String arrays are saved via h5py's variable-length string
                # dtype (see _save_dict_as_h5) and read back as an object
                # array of bytes - decode each element to undo that.
                return np.array([
                    x.decode('utf-8') if isinstance(x, bytes) else x for x in v
                ])
            return v

        with h5py.File(file, "r") as f:

            def load_group(group):
                data = {}
                for k, v in group.items():
                    if isinstance(v, h5py.Group):
                        data[k] = load_group(v)
                    elif isinstance(v, h5py.Dataset):
                        data[k] = decode_bytes(v[()])
                    else:
                        data[k] = v
                # Values that Cn't be stored natively (see
                # _save_dict_as_h5) are saved as group attrs instead of
                # datasets/subgroups; without this they'd silently disappear.
                for k, v in group.attrs.items():
                    data[k] = None if v == 'None' else decode_bytes(v)
                return data

            data = load_group(f)
        return data

    def _save_dict_as_h5(self, data: dict, file_name: str, overwrite=True):
        """Save the provided dictionary to an h5 file

        Any value that can't be represented natively in HDF5 (heterogeneous
        lists, unsupported object attributes, circular references, etc.) is
        stored as a ``repr()`` string attribute instead of aborting the save.

        Args:
            results (dict): Dictionary to be saved.
            file_name (str): Name (without extension) of the h5 file to write.
            overwrite (bool): If False and the target file already exists, raise
                instead of silently overwriting it.
        """
        cwd = Path.cwd()

        file_path = str(cwd / f"{file_name}.h5")

        if not overwrite and Path(file_path).exists():
            raise FileExistsError(
                f"{file_path} already exists and overwrite=False")

        # Tracks ids of custom objects currently being recursed into, to avoid
        # infinite recursion on circular references (e.g. an object holding a
        # reference back to something already being saved).
        seen_object_ids = set()

        def sanitize_key(k):
            # '/' is the HDF5 path separator: an unsanitized key would silently
            # create nested groups (or collide with an existing dataset).
            k = str(k)
            return k.replace('/', '_') if '/' in k else k

        def store_as_repr(group, key, value, reason):
            group.attrs[key] = repr(value)
            self.ostream.print_warning(
                f"Key '{key}': {reason}, stored {type(value).__name__} as string repr"
            )

        def save_group(data, group):
            for raw_k, v in data.items():
                k = sanitize_key(raw_k)
                if k == 'pdb_active_res':
                    continue
                try:
                    save_item(k, v, group)
                except Exception as e:
                    store_as_repr(group, k, v,
                                  f"failed with {type(e).__name__}: {e}")

        def save_item(k, v, group):
            if v is None:
                group.attrs[k] = 'None'
            elif isinstance(v, dict):
                # Recurse into nested dicts as HDF5 subgroups
                subgroup = group.create_group(k)
                save_group(v, subgroup)
            elif isinstance(v, (np.ndarray, list, set, tuple)):
                # sets/tuples are converted to list first; np.array handles both
                arr = np.array(list(v) if isinstance(v, (set, tuple)) else v)
                if arr.dtype == object:
                    # Heterogeneous/ragged data (e.g. a list of dicts or custom
                    # objects) has no native HDF5 representation.
                    store_as_repr(
                        group, k, v,
                        "no native HDF5 representation for object dtype")
                elif arr.dtype.kind == 'U':
                    # h5py has no conversion path for numpy's fixed-width
                    # unicode dtype ('<U..'); store as variable-length
                    # UTF-8 strings instead.
                    group.create_dataset(
                        k, data=arr.astype(h5py.string_dtype(encoding='utf-8')))
                else:
                    group.create_dataset(k, data=arr)
            elif isinstance(v, (bool, int, float, str, bytes, np.generic)):
                # np.generic covers numpy scalars (np.float64, np.int32, etc.)
                group[k] = v
            elif hasattr(v, '__dict__'):
                # Custom objects: recurse into their attributes as a subgroup
                if id(v) in seen_object_ids:
                    store_as_repr(group, k, v, "circular reference detected")
                    return
                seen_object_ids.add(id(v))
                try:
                    subgroup = group.create_group(k)
                    save_group(vars(v), subgroup)
                finally:
                    seen_object_ids.discard(id(v))
            else:
                # Last resort: let h5py try; fall back to repr string if it fails
                try:
                    group[k] = v
                except TypeError:
                    store_as_repr(group, k, v, "unsupported type for h5py")

        with h5py.File(file_path, "w") as file:
            self.ostream.print_info(f"Saving results to {file_path}")
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
        elif name == "water" or name == "water_NVT":
            conf = {
                "name": f"water_{self.water_model}_NVT",
                "solvent": self.water_model,
                "temperature": self.temperature,
                "pressure": 1,
                "isothermal": True,
                "padding": 1.5,
                "ion_count": 0,
                "neutralize": False
            }
        elif name == "water_NVE":
            conf = {
                "name": f"water_{self.water_model}_NVT",
                "solvent": self.water_model,
                "temperature": self.temperature,
                "pressure": 1,
                "isothermal": False,
                "padding": 1.5,
                "ion_count": 0,
                "neutralize": False
            }
        elif name == "water_NPT":
            conf = {
                "name": f"water_{self.water_model}_NPT",
                "solvent": self.water_model,
                "temperature": self.temperature,
                "pressure": 1,
                "isobaric": True,
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
                # "bonded_integration": True,
                "soft_core_coulomb_pes": True,
                "soft_core_lj_pes": True,
                "soft_core_coulomb_int": False,
                "soft_core_lj_int": False,
            }
        else:
            try:
                solvent_prop_not_used = SolvationBuilder(
                    ostream=self.ostream)._solvent_properties(name)
                conf = {
                    "name": name,
                    "solvent": name,
                    "temperature": self.temperature,
                    "pressure": 1,
                    "padding": 1.5,
                    "ion_count": 0,
                }
            except ValueError:
                # _solvent_properties raises ValueError for unrecognized names
                raise ValueError(f"Unknown system configuration {name}")

        return conf

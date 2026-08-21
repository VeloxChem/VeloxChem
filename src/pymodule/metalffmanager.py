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
import sys

from networkx.algorithms.isomorphism import GraphMatcher
import networkx as nx

from .veloxchemlib import mpi_master
from .molecule import Molecule
from .outputstream import OutputStream
from .metalsiteffbuilder import MetalSiteForceFieldBuilder
from .optimizationdriver import OptimizationDriver
from .superimpose import svd_superimpose
from .errorhandler import assert_msg_critical


class MetalForceFieldManager:
    """
    Keeps a set of metal force fields as templates and transfers them onto new
    enzymes.

    The parameters of a metal site are what a run of MetalSiteForceFieldBuilder
    pays QM for: the Seminario metal-ligand terms and the RESP charges. A site
    built from the same residues in the same coordination is described by the
    same numbers, so a folder left behind by an earlier run can stand in for the
    QM of the next one. A template is matched to a new structure by graph
    isomorphism of the active site and by the RMSD between the two geometries,
    and the mapping the isomorphism gives is what carries the parameters across.

    Only settings are kept on the object, apart from the templates themselves.
    Everything a match produces is returned by match_enzyme_to_template.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - templates: The loaded templates, keyed by name.
        - builder: The MetalSiteForceFieldBuilder whose steps are used to
          extract and parameterize an active site. Its settings apply here,
          so set them on it.
        - match_criterion: What a geometry has to satisfy to count as the same
          site as a template: 'xyz_rmsd' for the cartesian RMSD alone,
          'ic_rmsd' for the internal coordinate deviations alone, 'ic_or_xyz'
          for either and 'ic_and_xyz' for both.
        - rmsd_region: What every RMSD is measured over, cartesian and
          internal alike: 'active_site' for all of it, 'metal_shell' for the
          metals and everything within metal_shell_bonds bonds of them, or
          'metal_beta_carbons' for the metals and the beta carbons alone. The
          atom mapping is chosen on the whole site whichever is set, since a
          correspondence between two sites is not a question about a region.
        - metal_shell_bonds: How many bonds out from a metal the metal_shell
          region reaches.
        - geometry_rmsd_threshold: The RMSD in Angstrom below which an active
          site counts as the same site as a template.
        - ic_rmsd_thresholds: The internal coordinate deviations a match may
          have, as {'bonds': {'rms': , 'max': }, 'angles': ..., 'dihedrals':
          ...}, in Angstrom for the bonds and degrees for the rest. A type set
          to None, or a missing rms or max, is not checked.
        - metal_bond_tolerance: How far a metal-ligand bond may differ from the
          template, in Angstrom, before the match is rejected however good the
          RMSD is.
        - rmsd_heavy_atoms_only: The flag for leaving the hydrogens out of
          every RMSD, cartesian and internal alike. On by default: protonate
          places the hydrogens through Modeller.addHydrogens, which does not
          do it reproducibly, so two runs over one structure differ by around
          0.3 A over all atoms and by nothing at all over the heavy ones. Both
          cartesian numbers are printed either way.
        - max_mappings: The limit on the number of isomorphisms evaluated for
          one template.
        - do_prepare_protein: The flag for repairing the query structure before
          anything is extracted from it. Needed for an enzyme system, and it
          needs pdbfixer.
        - do_mm_fallback: The flag for relaxing the query active site on a
          crude force field and comparing again when nothing matched.
        - mm_fallback_literature_bonds: The flag for pulling the metal bonds of
          that relaxation toward LITERATURE_METAL_BONDS, which normalizes the
          query toward the geometry a template was optimized to. Ignored when
          the builder already carries metal bond equilibria of its own.
        - build_enzyme_system: The flag for going on to build a protein system
          from a matched force field. Off by default, since it needs pdbfixer
          and a structure a protein force field can template.
        - comm: The MPI communicator.
        - rank: The rank of the MPI process.
        - nodes: The number of MPI processes.
        - ostream: The output stream.
    """

    # The geometry a run leaves behind, in the order it is looked for. Which
    # of them a template is allowed to be built from is what the fallback
    # argument of load_ff_from_folder decides.
    GEOMETRY_KINDS = ('qm_opt', 'mm_opt')

    # What a geometry has to satisfy to count as a match. The cartesian RMSD
    # answers whether the whole site sits where the template has it; the
    # internal coordinates answer whether it is built the same way, which is
    # the question a force field is actually transferred on.
    MATCH_CRITERIA = ('xyz_rmsd', 'ic_rmsd', 'ic_or_xyz', 'ic_and_xyz')

    # Which atoms an RMSD is measured over. The whole active site answers
    # whether two sites are the same site; the metals with everything within
    # metal_shell_bonds bonds of them answers whether the coordination sphere
    # is the same, which is what the transferred parameters describe and all
    # they describe; the metals with the beta carbons answers whether the
    # residues are anchored in the same places, which is the frame of the site
    # with every sidechain conformation left out of it.
    RMSD_REGIONS = ('active_site', 'metal_shell', 'metal_beta_carbons')

    # The internal coordinate types get_ic_rmsd reports, with their units.
    IC_TYPES = {
        'bonds': 'A',
        'angles': 'deg',
        'dihedrals': 'deg',
    }

    def __init__(self, comm=None, ostream=None):
        """
        Initializes the metal force field manager.
        """

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

        self.templates = {}

        # the builder carries every setting of the structural pipeline, so it
        # is exposed rather than wrapped
        self.builder = MetalSiteForceFieldBuilder(comm, ostream)

        # matching
        self.match_criterion = 'xyz_rmsd'
        self.geometry_rmsd_threshold = 0.5
        self.rmsd_region = 'active_site'
        self.metal_shell_bonds = 2

        # Internal coordinate deviations, through
        # OptimizationDriver.get_ic_rmsd. The dihedrals are left unchecked:
        # the ones about a metal center are floppy, and on a binuclear zinc
        # site two geometries whose heavy atoms are identical still differ by
        # tens of degrees over them. Set a threshold here to check them
        # anyway. These numbers assume rmsd_heavy_atoms_only; with the
        # hydrogens in, the placement noise alone reaches 3.6 degrees rms and
        # 16.8 degrees max over the angles of the same structure compared to
        # itself, and the thresholds have to clear that.
        self.ic_rmsd_thresholds = {
            'bonds': {'rms': 0.05, 'max': 0.15},
            'angles': {'rms': 5.0, 'max': 20.0},
            'dihedrals': None,
        }
        self.metal_bond_tolerance = 0.25
        # the hydrogens carry the noise of Modeller.addHydrogens rather than
        # anything about the site; see the class docstring
        self.rmsd_heavy_atoms_only = True
        self.max_mappings = 10000

        # what to run on the query structure
        self.do_prepare_protein = True
        self.do_mm_fallback = True
        self.mm_fallback_literature_bonds = True
        self.build_enzyme_system = False

    def load_ff_from_folder(self, folder, name=None, fallback=None):
        """
        Loads one template from a folder left behind by a builder run.

        Two files are read: the force field, which carries the fitted metal
        terms, the charges and the bonding topology of the active site, and the
        geometry the fit was done on. Nothing else of the folder is needed, and
        the atom indices of the force field are not assumed to mean anything
        outside it: the mapping onto a new site is solved when a match is made.

        The geometry a run writes to opt_active_site.xyz is whatever the active
        site ended up as, which is only a QM optimized geometry if the run
        actually optimized. That is checked against mm_active_site.xyz where the
        crude MM pass left one: two identical files mean no QM optimization ever
        ran. Where the folder holds no MM geometry to compare against, the file
        name is taken at its word.

        :param folder:
            The folder of an earlier run.
        :param name:
            The name to store the template under. Defaults to the name of the
            folder.
        :param fallback:
            The lowest geometry kind that is acceptable. None accepts a QM
            optimized geometry only; 'mm_opt' also accepts one that was merely
            relaxed on the crude force field.

        :return:
            The template dictionary.
        """

        folder = Path(folder)
        assert_msg_critical(
            folder.is_dir(),
            f'MetalForceFieldManager: template folder {folder} not found')

        assert_msg_critical(
            fallback is None or fallback in self.GEOMETRY_KINDS,
            'MetalForceFieldManager: fallback must be None or one of '
            f'{list(self.GEOMETRY_KINDS)}, got {fallback}')

        if name is None:
            name = folder.resolve().name

        ff_path = folder / MetalSiteForceFieldBuilder.FORCEFIELD_FILE
        assert_msg_critical(
            ff_path.is_file(),
            f'MetalForceFieldManager: {ff_path} not found, so {folder} holds '
            'no force field to use as a template. It is written by a run with '
            'save_output on, and by save_results.')

        forcefield = self.builder.load_forcefield(ff_path)
        geometry, kind = self._load_geometry(folder, fallback)
        forcefield.molecule = geometry

        template = self._build_template(name, forcefield, geometry, kind,
                                        folder)

        if name in self.templates:
            self.ostream.print_warning(
                f'A template named {name} is already loaded; replacing it '
                f'with the one from {folder}')

        self.templates[name] = template

        self._print_template(template)
        self._print_templates()

        return template

    def load_folders(self, folders):
        """
        Loads several template folders at once.

        Not implemented yet. It is where a shared naming policy and a single
        summary over a whole series of runs will live; until then, call
        load_ff_from_folder once per folder.

        :param folders:
            The folders of earlier runs.
        """

        raise NotImplementedError(
            'MetalForceFieldManager.load_folders is not implemented yet; call '
            'load_ff_from_folder once per folder.')

    def _load_geometry(self, folder, fallback):
        """
        Reads the active site geometry of a template folder and works out what
        kind of geometry it is.

        :param folder:
            The folder of an earlier run.
        :param fallback:
            The lowest acceptable kind, as given to load_ff_from_folder.

        :return:
            The tuple of the molecule and its kind.
        """

        opt_path = folder / MetalSiteForceFieldBuilder.GEOMETRY_FILE
        mm_path = folder / MetalSiteForceFieldBuilder.MM_GEOMETRY_FILE

        allowed = self.GEOMETRY_KINDS[:1 if fallback is None else 2]

        if opt_path.is_file():
            geometry = Molecule.read_xyz_file(str(opt_path))

            if mm_path.is_file() and self._same_geometry(
                    geometry, Molecule.read_xyz_file(str(mm_path))):
                # save_results writes whatever the active site ended up as
                # under the name of the optimized geometry, so an untouched
                # copy of the MM geometry means the QM optimization never ran
                kind = 'mm_opt'
                assert_msg_critical(
                    kind in allowed,
                    f'MetalForceFieldManager: {opt_path} is identical to '
                    f'{mm_path.name}, so it is the crude MM geometry and no '
                    'QM optimization ran. Pass fallback="mm_opt" to use it '
                    'anyway.')
            else:
                kind = 'qm_opt'

            return geometry, kind

        assert_msg_critical(
            mm_path.is_file(),
            f'MetalForceFieldManager: neither {opt_path.name} nor '
            f'{mm_path.name} found in {folder}, so there is no geometry to '
            'build a template from.')

        assert_msg_critical(
            'mm_opt' in allowed,
            f'MetalForceFieldManager: {folder} holds only {mm_path.name}, a '
            'geometry relaxed on the crude force field. Pass '
            'fallback="mm_opt" to use it anyway.')

        return Molecule.read_xyz_file(str(mm_path)), 'mm_opt'

    @staticmethod
    def _same_geometry(first, second):
        """
        Whether two molecules hold the same coordinates.

        :param first:
            The first molecule.
        :param second:
            The second molecule.

        :return:
            True when they match to within the precision of an xyz file.
        """

        if first.number_of_atoms() != second.number_of_atoms():
            return False

        return np.allclose(first.get_coordinates_in_angstrom(),
                           second.get_coordinates_in_angstrom(),
                           atol=1.0e-6)

    def _build_template(self, name, forcefield, geometry, kind, folder):
        """
        Assembles a template from a loaded force field and geometry.

        The bonding topology is not a separate file: the bonds of the force
        field are the edges of the active site, metal-ligand bonds included,
        since build_forcefield builds them from the connectivity matrix. The
        capping hydrogens are the ones left at exactly zero charge by
        redistribute_cap_charges, which is how a cap is told apart from the
        hydrogens it is symmetric with.

        :param name:
            The name of the template.
        :param forcefield:
            The force field generator loaded from the folder.
        :param geometry:
            The active site geometry.
        :param kind:
            The geometry kind.
        :param folder:
            The folder the template came from.

        :return:
            The template dictionary.
        """

        labels = list(geometry.get_labels())
        elements = self.builder._forcefield_elements(forcefield)

        assert_msg_critical(
            len(labels) == len(elements),
            f'MetalForceFieldManager: the force field of {folder} has '
            f'{len(elements)} atoms but its geometry has {len(labels)}')

        assert_msg_critical(
            labels == elements,
            f'MetalForceFieldManager: the elements of the force field of '
            f'{folder} do not match its geometry, so the two files describe '
            'different structures')

        metal_indices = [
            index for index, label in enumerate(labels)
            if label in self.builder.metal_elements
        ]

        assert_msg_critical(
            len(metal_indices) > 0,
            f'MetalForceFieldManager: the template of {folder} holds no metal '
            f'center. Recognized elements: {self.builder.metal_elements}')

        charges = np.array(
            [atom['charge'] for atom in forcefield.atoms.values()])
        cap_indices = [
            index for index, label in enumerate(labels)
            if label == 'H' and charges[index] == 0.0
        ]

        # A force field built without charges has every atom at zero, which
        # makes every hydrogen look like a cap and would put nothing but
        # zeroes on whatever it is transferred to.
        if not np.any(charges):
            self.ostream.print_warning(
                f'The force field of {folder} carries no charges at all, so '
                'it was built without RESP. Its capping hydrogens cannot be '
                'told from the others, and a site matched to it would be '
                'given a charge of zero on every atom.')
            self.ostream.flush()

        graph = self._graph(labels, forcefield.bonds.keys(), cap_indices)

        return {
            'name': name,
            'forcefield': forcefield,
            'geometry': geometry,
            'geometry_kind': kind,
            'labels': labels,
            'metal_indices': metal_indices,
            'cap_indices': cap_indices,
            'charges': charges,
            'graph': graph,
            'folder': str(folder),
        }

    @staticmethod
    def _graph(labels, edges, cap_indices):
        """
        Builds the graph an isomorphism is solved on.

        :param labels:
            The element of every atom.
        :param edges:
            The bonds, as index pairs.
        :param cap_indices:
            The indices of the capping hydrogens.

        :return:
            The graph, with an elem and an is_cap attribute per node.
        """

        caps = set(cap_indices)
        graph = nx.Graph()

        for index, label in enumerate(labels):
            graph.add_node(index, elem=label, is_cap=index in caps)

        for i, j in edges:
            graph.add_edge(i, j)

        return graph

    def _active_site_graph(self, active_site, connectivity_matrix):
        """
        Builds the graph of an extracted active site.

        :param active_site:
            The active site.
        :param connectivity_matrix:
            Its connectivity.

        :return:
            The graph.
        """

        matrix = np.asarray(connectivity_matrix)
        edges = [(int(i), int(j))
                 for i, j in zip(*np.triu_indices_from(matrix, k=1))
                 if matrix[i, j]]

        return self._graph(list(active_site['labels']), edges,
                           active_site['cap_indices'])

    def match_enzyme_to_template(self, structure):
        """
        Matches the metal site of a structure against the loaded templates and
        builds its force field from the one that fits.

        The active site of the structure is extracted exactly as a builder run
        would extract it, and compared to every template in two steps. The
        bonding topology decides first: an isomorphism between the two active
        sites is what says the same residues coordinate the same metals in the
        same way, and it is solved rather than assumed, so the residues of the
        structure may come in any order. The geometry decides second, on the
        RMSD of the best of the isomorphisms.

        A structure that matches nothing on its input geometry is relaxed on a
        crude force field of its own and compared again, which is what takes
        the slack out of an unrelaxed design model.

        :param structure:
            The path to a PDB or mmCIF file.

        :return:
            The results dictionary, or None when no template matched.
        """

        assert_msg_critical(
            len(self.templates) > 0,
            'MetalForceFieldManager.match_enzyme_to_template: no templates '
            'loaded. Call load_ff_from_folder first.')

        assert_msg_critical(
            self.match_criterion in self.MATCH_CRITERIA,
            'MetalForceFieldManager: match_criterion must be one of '
            f'{list(self.MATCH_CRITERIA)}, got {self.match_criterion}')

        assert_msg_critical(
            self.rmsd_region in self.RMSD_REGIONS,
            'MetalForceFieldManager: rmsd_region must be one of '
            f'{list(self.RMSD_REGIONS)}, got {self.rmsd_region}')

        query = self._prepare_query(structure)
        active_site = query['active_site']
        self.last_mol = active_site['molecule']
        graph = self._active_site_graph(active_site,
                                        query['connectivity_matrix'])

        self._print_query(structure, active_site)

        mappings = {}
        report = []
        hit = None

        coordinates = active_site['molecule'].get_coordinates_in_angstrom()

        # the relaxed geometry is appended to this list while it is being
        # walked, so the fallback is a second pass of the same loop rather
        # than a copy of it
        stages = [('input', coordinates)]

        for stage, coords in stages:
            candidates = []

            for name, template in self.templates.items():
                if name not in mappings:
                    # the isomorphisms depend on the two graphs alone, so a
                    # second stage on a relaxed geometry reuses them
                    mappings[name] = self._find_mappings(
                        template, graph, active_site)

                if not mappings[name]:
                    report.append((name, stage, None, None, 'no isomorphism'))
                    continue

                mapping, rmsd, heavy_rmsd, deviation = self._best_mapping(
                    template, mappings[name], coords)

                scored = heavy_rmsd if self.rmsd_heavy_atoms_only else rmsd
                ic_rmsd = None

                if deviation > self.metal_bond_tolerance:
                    verdict = f'metal bond off by {deviation:.2f} A'
                else:
                    matched, verdict, ic_rmsd = self._apply_criterion(
                        template, mapping, coords, scored)
                    if matched:
                        candidates.append((scored, name, mapping, rmsd,
                                           heavy_rmsd, ic_rmsd))

                report.append((name, stage, rmsd, heavy_rmsd, verdict))

            if candidates:
                scored, name, mapping, rmsd, heavy_rmsd, ic_rmsd = min(
                    candidates, key=lambda row: row[0])

                if len(candidates) > 1:
                    self._print_multiple_hits(candidates, name, stage)

                hit = {
                    'name': name,
                    'mapping': mapping,
                    'rmsd': scored,
                    'rmsd_all_atoms': rmsd,
                    'rmsd_heavy_atoms': heavy_rmsd,
                    'ic_rmsd': ic_rmsd,
                    'stage': stage,
                }
                break

            if stage == 'input' and self.do_mm_fallback:
                self.ostream.print_info(
                    'No template matched the input geometry; relaxing the '
                    'active site on a crude force field and comparing again.')
                self.ostream.flush()
                relaxed = self._mm_relax(query)
                self.last_relaxed = relaxed
                stages.append(
                    ('mm_relaxed', relaxed.get_coordinates_in_angstrom()))

        self._print_report(report)

        if hit is None:
            self.ostream.print_info(
                'No template matched. Returning None; a metal site the '
                'templates do not cover has to be built with '
                'MetalSiteForceFieldBuilder.')
            self.ostream.flush()
            return None

        template = self.templates[hit['name']]
        self.ostream.print_info(
            f'Matched {hit["name"]} on the {hit["stage"]} geometry with an '
            f'RMSD of {hit["rmsd"]:.3f} A. Transferring its metal terms and '
            'charges.')
        self.ostream.flush()

        forcefield = self._transfer(template, hit['mapping'], active_site,
                                    query['connectivity_matrix'])

        enzyme_system = None
        if self.build_enzyme_system:
            enzyme_system, _ = self.builder.create_enzyme_system(
                query['topology'], active_site, forcefield,
                forcefield.partial_charges)

        self.builder._print_metal_parameters(active_site, forcefield)

        return {
            'template': hit['name'],
            'forcefield': forcefield,
            'mapping': hit['mapping'],
            'rmsd': hit['rmsd'],
            'rmsd_all_atoms': hit['rmsd_all_atoms'],
            'rmsd_heavy_atoms': hit['rmsd_heavy_atoms'],
            'ic_rmsd': hit['ic_rmsd'],
            'stage': hit['stage'],
            'geometry_kind': template['geometry_kind'],
            'topology': query['topology'],
            'positions': query['positions'],
            'binding_modes': query['binding_modes'],
            'active_site': active_site,
            'connectivity_matrix': query['connectivity_matrix'],
            'enzyme_system': enzyme_system,
        }

    def _prepare_query(self, structure):
        """
        Runs the structural half of the builder pipeline on a structure.

        :param structure:
            The path to a PDB or mmCIF file.

        :return:
            A dictionary with the topology, the positions, the binding modes,
            the active site and its connectivity.
        """

        assert_msg_critical('openmm' in sys.modules,
                            'MetalForceFieldManager: openmm is required')

        builder = self.builder

        topology, positions = builder.load_structure(structure)

        if self.do_prepare_protein:
            # every index downstream refers to this topology, so the repair
            # has to happen before anything is extracted from it
            topology, positions = builder.prepare_protein(topology, positions)

        binding_modes = builder.suggest_binding_modes(topology, positions)
        topology, positions, binding_modes = builder.protonate(
            topology, positions, binding_modes)
        active_site = builder.extract_active_site(topology, positions,
                                                  binding_modes)
        connectivity_matrix = builder.build_connectivity(
            topology, active_site, binding_modes)

        return {
            'topology': topology,
            'positions': positions,
            'binding_modes': binding_modes,
            'active_site': active_site,
            'connectivity_matrix': connectivity_matrix,
        }

    def _mm_relax(self, query):
        """
        Relaxes the query active site on a crude force field of its own.

        The metal bonds of that force field are seeded from the geometry it is
        handed, so on their own they would keep whatever the structure file
        happened to hold. Pulling them to the literature distances instead
        moves the query toward the geometry a template was optimized to, which
        is the whole point of comparing again after relaxing.

        :param query:
            The prepared query.

        :return:
            The relaxed molecule.
        """

        builder = self.builder
        equilibria = builder.metal_bond_equilibria
        override = (self.mm_fallback_literature_bonds and equilibria is None)

        if override:
            builder.metal_bond_equilibria = builder.LITERATURE_METAL_BONDS

        # this geometry is a way of comparing, not a result of a run, so it
        # does not get a folder of its own
        save_output = builder.save_output
        builder.save_output = False

        try:
            return builder.mm_optimize_active_site(query['active_site'],
                                                   query['connectivity_matrix'])
        finally:
            builder.save_output = save_output
            if override:
                builder.metal_bond_equilibria = equilibria

    def _find_mappings(self, template, graph, active_site):
        """
        Finds every atom mapping from a template onto an active site.

        A truncated site is full of symmetry - the hydrogens of a CB, the two
        oxygens of a chelating carboxylate, two histidines in the same role -
        so there is rarely one mapping. All of them are kept and the geometry
        picks between them.

        Capping hydrogens are pinned to capping hydrogens where the template
        says which of its atoms they are. Without that a cap is free to map
        onto a real hydrogen of the same CB, which would put the charge of an
        alpha carbon on a hydrogen that has to keep its own.

        :param template:
            The template.
        :param graph:
            The graph of the active site.
        :param active_site:
            The active site, for the number of capping hydrogens.

        :return:
            The list of mappings from template index to active site index.
        """

        if not self._invariants_match(template, graph):
            return []

        pin_caps = len(template['cap_indices']) == len(
            active_site['cap_indices'])

        if not pin_caps:
            self.ostream.print_warning(
                f'Template {template["name"]} has '
                f'{len(template["cap_indices"])} capping hydrogen(s) at zero '
                f'charge but the active site has '
                f'{len(active_site["cap_indices"])}; matching on the elements '
                'alone, so a cap may end up mapped onto an ordinary hydrogen')

        def node_match(a, b):
            if a['elem'] != b['elem']:
                return False
            if pin_caps:
                return a['is_cap'] == b['is_cap']
            return True

        matcher = GraphMatcher(template['graph'], graph, node_match=node_match)

        mappings = []
        for mapping in matcher.isomorphisms_iter():
            mappings.append(mapping)
            if len(mappings) >= self.max_mappings:
                self.ostream.print_warning(
                    f'Template {template["name"]} reached the limit of '
                    f'{self.max_mappings} mappings; the best of the ones '
                    'found is used, which need not be the best there is')
                break

        return mappings

    @staticmethod
    def _invariants_match(template, graph):
        """
        Rejects a template that cannot possibly be isomorphic, before the
        isomorphism itself is solved.

        :param template:
            The template.
        :param graph:
            The graph of the active site.

        :return:
            Whether the two graphs share their cheap invariants.
        """

        first, second = template['graph'], graph

        if first.number_of_nodes() != second.number_of_nodes():
            return False

        if first.number_of_edges() != second.number_of_edges():
            return False

        def signature(one):
            elements = sorted(data['elem'] for _, data in one.nodes(data=True))
            degrees = sorted(degree for _, degree in one.degree())
            return elements, degrees

        return signature(first) == signature(second)

    def _best_mapping(self, template, mappings, coordinates):
        """
        Picks the mapping that superimposes the two geometries best.

        :param template:
            The template.
        :param mappings:
            The mappings to choose between.
        :param coordinates:
            The coordinates of the active site, in Angstrom.

        The mapping is chosen on the whole active site whatever rmsd_region
        says: which atom of one site is which atom of the other is not a
        question about a region, and deciding it on a fragment would leave the
        atoms outside that fragment to be paired up arbitrarily. The RMSDs
        that come back are measured over the region.

        :return:
            The tuple of the best mapping, its RMSD over all atoms of the
            region, its RMSD over the heavy atoms of the region and the
            largest metal-ligand bond difference under it.
        """

        reference = template['geometry'].get_coordinates_in_angstrom()
        region = self._region_indices(template)
        heavy = [
            index for index in region if template['labels'][index] != 'H'
        ]

        best = None

        for mapping in mappings:
            order = [mapping[index] for index in range(len(reference))]
            moved = coordinates[order]

            # the whole site decides the mapping
            whole, _, _ = svd_superimpose(moved, reference)

            rmsd, _, _ = svd_superimpose(moved[region], reference[region])
            heavy_rmsd, _, _ = svd_superimpose(moved[heavy], reference[heavy])

            if best is None or whole < best[0]:
                best = (whole, mapping, rmsd, heavy_rmsd)

        _, mapping, rmsd, heavy_rmsd = best
        deviation = self._metal_bond_deviation(template, mapping, coordinates)

        return mapping, rmsd, heavy_rmsd, deviation

    def _rmsd_indices(self, template):
        """
        Returns the template indices every RMSD is measured over, which is the
        region less the hydrogens when they are being left out.

        :param template:
            The template.

        :return:
            The indices, in order.
        """

        region = self._region_indices(template)

        if not self.rmsd_heavy_atoms_only:
            return region

        return [
            index for index in region if template['labels'][index] != 'H'
        ]

    def _beta_carbons(self, template, graph):
        """
        Finds the beta carbons of a template.

        A template is a force field and a geometry, and neither of them says
        which atom is which residue's CB. The truncation does: it cuts at the
        CA-CB bond and puts a capping hydrogen where the alpha carbon was, so
        the one atom a cap is bonded to is the beta carbon it was cut from.

        :param template:
            The template.
        :param graph:
            Its graph.

        :return:
            The indices of the beta carbons.
        """

        assert_msg_critical(
            len(template['cap_indices']) > 0,
            f'MetalForceFieldManager: template {template["name"]} has no '
            'capping hydrogens, so its beta carbons cannot be found')

        beta_carbons = set()

        for cap in template['cap_indices']:
            neighbors = list(graph.neighbors(cap))

            assert_msg_critical(
                len(neighbors) == 1 and template['labels'][neighbors[0]] == 'C',
                f'MetalForceFieldManager: in template {template["name"]} the '
                f'capping hydrogen {cap} is not bonded to exactly one carbon, '
                'so it is not a cap. The capping hydrogens are read off the '
                'charges, which a force field built without RESP does not '
                'have.')

            beta_carbons.add(neighbors[0])

        return beta_carbons

    def _region_indices(self, template):
        """
        Returns the template indices of the region an RMSD is measured over.

        :param template:
            The template.

        :return:
            The indices, in order.
        """

        if self.rmsd_region == 'active_site':
            return list(range(len(template['labels'])))

        graph = template['graph']

        if self.rmsd_region == 'metal_beta_carbons':
            return sorted(
                set(template['metal_indices'])
                | self._beta_carbons(template, graph))
        region = set()

        for metal in template['metal_indices']:
            region.update(
                nx.single_source_shortest_path_length(
                    graph, metal, cutoff=self.metal_shell_bonds))

        assert_msg_critical(
            len(region) > len(template['metal_indices']),
            'MetalForceFieldManager: the metal centers of template '
            f'{template["name"]} have nothing bonded to them within '
            f'{self.metal_shell_bonds} bond(s)')

        return sorted(region)

    def _apply_criterion(self, template, mapping, coordinates, scored):
        """
        Decides whether a geometry counts as the template's site.

        Which of the two RMSDs has to pass is what match_criterion says. The
        cartesian one asks whether the site sits where the template has it;
        the internal one asks whether it is built the same way, bond by bond
        and angle by angle, which is closer to the question a transferred
        force field turns on.

        :param template:
            The template being compared to.
        :param mapping:
            The mapping from template index to active site index.
        :param coordinates:
            The coordinates of the active site, in Angstrom.
        :param scored:
            The cartesian RMSD already measured under that mapping.

        :return:
            The tuple of whether it matched, the verdict for the report and
            the internal coordinate deviations, or None when they were not
            needed.
        """

        criterion = self.match_criterion
        xyz_ok = scored <= self.geometry_rmsd_threshold

        if criterion == 'xyz_rmsd':
            return xyz_ok, 'match' if xyz_ok else 'geometry off', None

        ic_rmsd = self._ic_rmsd(template, mapping, coordinates)
        violation = self._ic_violation(ic_rmsd)
        ic_ok = violation is None

        if criterion == 'ic_rmsd':
            matched = ic_ok
            verdict = violation
        elif criterion == 'ic_and_xyz':
            matched = ic_ok and xyz_ok
            # name whichever failed first, so the report says what to loosen
            verdict = 'geometry off' if not xyz_ok else violation
        else:
            matched = ic_ok or xyz_ok
            verdict = 'geometry and IC off'

        return matched, 'match' if matched else verdict, ic_rmsd

    def _ic_rmsd(self, template, mapping, coordinates):
        """
        Measures the internal coordinate deviations from a template.

        The active site is reordered onto the atoms of the template first,
        since get_ic_rmsd pairs the two geometries by position and refuses two
        molecules whose elements do not line up. Only the atoms an RMSD is
        measured over are handed across, so that restricting one to the
        coordination sphere or to the heavy atoms restricts the internal
        coordinates that get built out of it in the same way.

        :param template:
            The template being compared to.
        :param mapping:
            The mapping from template index to active site index.
        :param coordinates:
            The coordinates of the active site, in Angstrom.

        :return:
            The deviations, as get_ic_rmsd reports them, or None when they
            could not be measured.
        """

        region = self._rmsd_indices(template)
        labels = [template['labels'][index] for index in region]

        reference = template['geometry'].get_coordinates_in_angstrom()
        order = [mapping[index] for index in region]

        mapped = Molecule(labels, coordinates[order], 'angstrom')
        wanted = Molecule(labels, reference[region], 'angstrom')

        ic_rmsd = OptimizationDriver.get_ic_rmsd(mapped, wanted)

        if isinstance(ic_rmsd, str):
            # get_ic_rmsd reports its refusals as a string rather than raising
            self.ostream.print_warning(
                f'Could not measure the internal coordinates against template '
                f'{template["name"]}: {ic_rmsd}')
            self.ostream.flush()
            return None

        return ic_rmsd

    def _ic_violation(self, ic_rmsd):
        """
        Checks internal coordinate deviations against ic_rmsd_thresholds.

        A type left at None, or a threshold left out of one, is not checked.
        The dihedrals are None by default: the ones about a metal center are
        floppy enough that they say nothing about whether two sites are the
        same.

        :param ic_rmsd:
            The deviations, as get_ic_rmsd reports them.

        :return:
            A description of the first threshold that was exceeded, or None
            when every one of them holds.
        """

        if ic_rmsd is None:
            return 'no internal coords'

        for name, unit in self.IC_TYPES.items():
            thresholds = self.ic_rmsd_thresholds.get(name)
            if not thresholds:
                continue

            found = ic_rmsd.get(name)
            if found is None:
                # the region holds no internal coordinate of this type
                continue

            for measure in ('rms', 'max'):
                limit = thresholds.get(measure)
                if limit is None:
                    continue
                if found[measure] > limit:
                    return (f'{name} {measure} {found[measure]:.2f} > '
                            f'{limit:.2f} {unit}')

        return None

    def _metal_bond_deviation(self, template, mapping, coordinates):
        """
        Returns the largest difference between the metal-ligand bonds of a
        template and those of an active site under a mapping.

        A whole site superimposes well while one metal contact is a long way
        off, and that contact is the parameter being transferred, so it is
        checked on its own.

        :param template:
            The template.
        :param mapping:
            The mapping from template index to active site index.
        :param coordinates:
            The coordinates of the active site, in Angstrom.

        :return:
            The largest difference in Angstrom, or 0.0 without metal bonds.
        """

        reference = template['geometry'].get_coordinates_in_angstrom()
        bonds, _ = self._metal_keys(template)

        deviation = 0.0

        for i, j in bonds:
            expected = np.linalg.norm(reference[i] - reference[j])
            found = np.linalg.norm(coordinates[mapping[i]] -
                                   coordinates[mapping[j]])
            deviation = max(deviation, abs(found - expected))

        return deviation

    def _metal_keys(self, template):
        """
        Returns the metal bond and angle keys of a template.

        :param template:
            The template.

        :return:
            The tuple of the bond key list and the angle key list.
        """

        # get_metal_keys reads nothing of the active site but the metal
        # indices, and a template knows those from its elements
        return self.builder.get_metal_keys(
            template['forcefield'],
            {'metal_indices': template['metal_indices']})

    def _transfer(self, template, mapping, active_site, connectivity_matrix):
        """
        Builds a force field for an active site out of a template.

        Only what a builder run pays QM for is taken from the template: the
        fitted metal bonds and angles, and the charges. Everything else is
        built for the site in front of us, so the atom types and the bonded
        terms of the residues come from the structure rather than from
        somewhere else.

        :param template:
            The template that matched.
        :param mapping:
            The mapping from template index to active site index.
        :param active_site:
            The active site of the query.
        :param connectivity_matrix:
            Its connectivity.

        :return:
            The force field generator.
        """

        template_ff = template['forcefield']
        charges = np.zeros(len(active_site['labels']))

        for template_index, site_index in mapping.items():
            charges[site_index] = template['charges'][template_index]

        total = float(np.sum(charges))
        if abs(total - active_site['charge']) > 1.0e-3:
            self.ostream.print_warning(
                f'The transferred charges sum to {total:+.3f}, but the active '
                f'site charge is {active_site["charge"]:+d}')

        # build_forcefield writes the charges onto the atoms and redistributes
        # the caps; the metal terms it seeds here are overwritten below
        forcefield = self.builder.build_forcefield(active_site,
                                                   connectivity_matrix, None,
                                                   charges)

        bonds, angles = self._metal_keys(template)

        for key in bonds:
            target = self._map_key(key, mapping, forcefield.bonds, 'bond')
            forcefield.bonds[target] = self._transferred(
                template_ff.bonds[key], template['name'])

        for key in angles:
            target = self._map_key(key, mapping, forcefield.angles, 'angle')
            forcefield.angles[target] = self._transferred(
                template_ff.angles[key], template['name'])

        self.ostream.print_info(
            f'Transferred {len(bonds)} metal bond(s), {len(angles)} metal '
            f'angle(s) and {len(charges)} charge(s) from {template["name"]}.')
        self.ostream.flush()

        return forcefield

    @staticmethod
    def _map_key(key, mapping, table, kind):
        """
        Maps a force field key onto the active site.

        A bond and an angle read the same forwards and backwards, and the
        generator stores only one of the two orders, so the reverse is tried
        before giving up. Silently dropping a term that came out the wrong way
        round would leave a metal site half parameterized.

        :param key:
            The key in the template.
        :param mapping:
            The mapping from template index to active site index.
        :param table:
            The bonds or angles of the force field being built.
        :param kind:
            The name of the term, for the error message.

        :return:
            The key in the force field being built.
        """

        mapped = tuple(mapping[index] for index in key)

        if mapped in table:
            return mapped

        if mapped[::-1] in table:
            return mapped[::-1]

        assert_msg_critical(
            False, f'MetalForceFieldManager: the template {kind} {key} maps '
            f'onto {mapped}, which the active site force field does not have')

    @staticmethod
    def _transferred(params, name):
        """
        Copies one set of parameters, recording where it came from.

        :param params:
            The parameters of the template.
        :param name:
            The name of the template.

        :return:
            The copied parameters.
        """

        params = dict(params)
        comment = params.get('comment', '')
        params['comment'] = f'{comment} (template {name})'.strip()

        return params

    @staticmethod
    def _param(label, value, label_width=26, value_width=20):
        """
        Formats one parameter line with fixed label and value widths.

        print_header centers what it is given, so every line has to be the
        same width to come out left aligned against the others.

        :param label:
            The label.
        :param value:
            The value.

        :return:
            The formatted line.
        """

        return f'{label:<{label_width}} : {str(value):>{value_width}}'

    def _print_template(self, template):
        """
        Prints what one template holds.

        :param template:
            The template.
        """

        bonds, angles = self._metal_keys(template)
        metals = ', '.join(template['labels'][index]
                           for index in template['metal_indices'])

        self.ostream.print_blank()
        self.ostream.print_header(f'Template {template["name"]}')
        self.ostream.print_header((9 + len(template['name'])) * '-')
        self.ostream.print_header(
            self._param('geometry', template['geometry_kind']))
        self.ostream.print_header(self._param('atoms', len(template['labels'])))
        self.ostream.print_header(self._param('metal centers', metals))
        self.ostream.print_header(
            self._param('capping hydrogens', len(template['cap_indices'])))
        self.ostream.print_header(self._param('metal bonds', len(bonds)))
        self.ostream.print_header(self._param('metal angles', len(angles)))
        self.ostream.print_header(
            self._param('total charge',
                        f'{float(np.sum(template["charges"])):+.3f}'))
        self.ostream.print_blank()
        self.ostream.print_info(f'Loaded from {template["folder"]}')
        self.ostream.flush()

    def _print_templates(self):
        """
        Prints every template that is loaded.
        """

        self.ostream.print_header(f'Loaded templates ({len(self.templates)})')
        self.ostream.print_header(60 * '-')
        valstr = '{:>24} | {:>7} | {:>7} | {:>13}'.format(
            'name', 'atoms', 'metals', 'geometry')
        self.ostream.print_header(valstr)
        self.ostream.print_header(60 * '-')

        for name, template in self.templates.items():
            valstr = '{:>24} | {:>7} | {:>7} | {:>13}'.format(
                name[:24], len(template['labels']),
                len(template['metal_indices']), template['geometry_kind'])
            self.ostream.print_header(valstr)

        self.ostream.print_blank()
        self.ostream.flush()

    def _print_query(self, structure, active_site):
        """
        Prints the active site a structure was reduced to.

        :param structure:
            The structure file.
        :param active_site:
            The extracted active site.
        """

        metals = ', '.join(active_site['labels'][index]
                           for index in active_site['metal_indices'])

        self.ostream.print_blank()
        self.ostream.print_header('Matching against the templates')
        self.ostream.print_header(30 * '-')
        self.ostream.print_header(self._param('structure',
                                              Path(structure).name))
        self.ostream.print_header(
            self._param('active site atoms', len(active_site['labels'])))
        self.ostream.print_header(self._param('metal centers', metals))
        self.ostream.print_header(
            self._param('match criterion', self.match_criterion))

        if self.rmsd_region == 'metal_shell':
            region = f'metals + {self.metal_shell_bonds} bond(s)'
        elif self.rmsd_region == 'metal_beta_carbons':
            region = 'metals + beta C'
        else:
            region = 'whole active site'
        self.ostream.print_header(self._param('RMSD over', region))
        self.ostream.print_header(
            self._param(
                'RMSD atoms',
                'heavy atoms' if self.rmsd_heavy_atoms_only else 'all atoms'))

        if self.match_criterion != 'ic_rmsd':
            self.ostream.print_header(
                self._param('RMSD threshold',
                            f'{self.geometry_rmsd_threshold:.2f} A'))

        if self.match_criterion != 'xyz_rmsd' and not self.rmsd_heavy_atoms_only:
            self.ostream.print_warning(
                'The internal coordinates are being measured with the '
                'hydrogens in. protonate does not place them reproducibly, so '
                'a site compared to a template built from the same structure '
                'reaches several degrees rms over its angles on that account '
                'alone; ic_rmsd_thresholds has to clear that.')

        if self.match_criterion != 'xyz_rmsd':
            for name, unit in self.IC_TYPES.items():
                # the unit goes in the label, so that the values stay inside
                # the width that keeps print_header looking left aligned
                label = f'IC {name} ({unit})'
                thresholds = self.ic_rmsd_thresholds.get(name)
                if not thresholds:
                    self.ostream.print_header(
                        self._param(label, 'not checked'))
                    continue
                limits = ', '.join(
                    f'{measure} {thresholds[measure]:.2f}'
                    for measure in ('rms', 'max') if thresholds.get(measure))
                self.ostream.print_header(self._param(label, limits))

        self.ostream.print_header(
            self._param('metal bond tolerance',
                        f'{self.metal_bond_tolerance:.2f} A'))
        self.ostream.print_blank()
        self.ostream.print_info(
            f'Residues: {", ".join(active_site["residues"])}')
        self.ostream.flush()

    def _print_multiple_hits(self, candidates, chosen, stage):
        """
        Reports that more than one template fits the same site.

        The closest of them is taken, which is the best that can be done from
        the geometry alone. Whether the others are copies of one run or
        genuinely different sites that the thresholds cannot tell apart is a
        question about the templates rather than about the structure, so the
        whole field is printed rather than only the winner.

        :param candidates:
            The templates that matched, as collected while matching.
        :param chosen:
            The name of the template that was taken.
        :param stage:
            The geometry they matched on.
        """

        self.ostream.print_info(
            f'{len(candidates)} templates match the {stage} geometry. Taking '
            f'{chosen}, the closest of them.')
        self.ostream.print_blank()

        valstr = '{:>24} | {:>10} | {:>10} | {:>8}'.format(
            'template', 'RMSD (A)', 'heavy (A)', 'taken')
        self.ostream.print_header(valstr)
        self.ostream.print_header(62 * '-')

        for row in sorted(candidates, key=lambda row: row[0]):
            name, rmsd, heavy_rmsd = row[1], row[3], row[4]
            valstr = '{:>24} | {:>10.3f} | {:>10.3f} | {:>8}'.format(
                name[:24], rmsd, heavy_rmsd, 'yes' if name == chosen else '')
            self.ostream.print_header(valstr)

        self.ostream.print_blank()
        self.ostream.flush()

    def _print_report(self, report):
        """
        Prints how every template fared, so that a near miss is visible.

        :param report:
            The rows collected while matching.
        """

        self.ostream.print_blank()
        # wide enough for the verdicts that carry a number and a threshold
        row = '{:>20} | {:>11} | {:>10} | {:>10} | {:>26}'
        valstr = row.format('template', 'geometry', 'RMSD (A)', 'heavy (A)',
                            'verdict')
        self.ostream.print_header(valstr)
        self.ostream.print_header(89 * '-')

        for name, stage, rmsd, heavy_rmsd, verdict in report:
            rmsd_str = '' if rmsd is None else f'{rmsd:.3f}'
            heavy_str = '' if heavy_rmsd is None else f'{heavy_rmsd:.3f}'
            valstr = row.format(name[:20], stage, rmsd_str, heavy_str,
                                verdict[:26])
            self.ostream.print_header(valstr)

        self.ostream.print_blank()
        self.ostream.flush()

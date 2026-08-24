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
from collections import Counter
from itertools import permutations, product
from pathlib import Path
import numpy as np
import math
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
        - selection_criteria: What build_ff_from_comparison holds a template to
          before its parameters are transferred: 'tight', 'loose' or a set of
          the same shape as those in SELECTION_CRITERIA. Unrelated to
          match_criterion, which is what the older matching path gates on.
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

    # What a template has to be within, region by region, before its
    # parameters are put on a new site. The bonds are what is checked: they
    # are the internal coordinates a transferred force field is written in,
    # and they carry none of the noise the floppy dihedrals about a metal do.
    # A region is read as {ic type: {rms: , max: }}, the shape
    # ic_rmsd_thresholds has, so that angles can be added to a set without
    # anything here changing. The tight set is what two runs of the same site
    # look like; the loose one is what an unrelaxed design model looks like.
    SELECTION_CRITERIA = {
        'tight': {
            'active_site': {
                'bonds': {
                    'rms': 0.5,
                    'max': 1.0
                }
            },
            'metal_shell': {
                'bonds': {
                    'rms': 0.25,
                    'max': 0.5
                }
            },
            'metal_beta_carbons': {
                'bonds': {
                    'rms': 0.5,
                    'max': 1.0
                }
            },
        },
        'loose': {
            'active_site': {
                'bonds': {
                    'rms': 1.0,
                    'max': 2.0
                }
            },
            'metal_shell': {
                'bonds': {
                    'rms': 0.5,
                    'max': 1.0
                }
            },
            'metal_beta_carbons': {
                'bonds': {
                    'rms': 1.0,
                    'max': 1.5
                }
            },
        },
    }

    # Which region decides between several templates that all pass. The
    # cartesian RMSD is not what is ranked on: the whole site is what a
    # transferred force field describes, and its bonds are what it is written
    # in, so the site that is built the most like the template wins.
    SELECTION_RANKED_ON = ('active_site', 'bonds', 'rms')

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
            'bonds': {
                'rms': 0.05,
                'max': 0.15
            },
            'angles': {
                'rms': 5.0,
                'max': 20.0
            },
            'dihedrals': None,
        }
        self.metal_bond_tolerance = 0.25
        # what a template has to be within before build_ff_from_comparison
        # puts it on a site; the call may override it
        self.selection_criteria = self.SELECTION_CRITERIA['tight']
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

    def _build_template(self, name, forcefield, molecule, kind, folder):
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
        :param molecule:
            The active site geometry.
        :param kind:
            The geometry kind.
        :param folder:
            The folder the template came from.

        :return:
            The template dictionary.
        """

        labels = list(molecule.get_labels())
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

        marker = MetalSiteForceFieldBuilder.BETA_CARBON_COMMENT
        beta_carbon_indices = [
            index for index, atom in forcefield.atoms.items()
            if marker in (atom.get('comment', '') or '')
        ]

        graph = self._graph(labels, forcefield.bonds.keys(), cap_indices)

        # A template is the description of a site with the two files it was
        # loaded from added to it, rather than a dictionary holding a
        # description: the labels, the graph and the metal and cap indices
        # would otherwise be in both of them, and a template could not be
        # handed to the matching steps that take a description.
        template = self._describe_site(labels, graph, metal_indices,
                                       cap_indices)

        template.update({
            'name': name,
            'forcefield': forcefield,
            'molecule': molecule,
            'geometry_kind': kind,
            'beta_carbon_indices': beta_carbon_indices,
            'charges': charges,
            'folder': str(folder),
        })

        return template

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

    def _describe_site(self, labels, graph, metal_indices, cap_indices):
        """
        Describes a site at the two levels a comparison is made on.

        The fine level is the atoms of each residue, where the bonding is
        ordinary covalent chemistry. The coarse level is which residue
        coordinates which metal, and nothing about how many atoms do the
        coordinating: a carboxylate that grips one metal with one oxygen and
        one that grips it with two are the same residue on the same metal, and
        only the coarse level is allowed to have an opinion about metal bonds.

        The residues are the connected components left once the metals are
        taken out of the graph. That decomposition is untouched by how the
        metal is gripped, which is the whole reason for splitting the
        comparison in two.

        :param labels:
            The element of every atom.
        :param graph:
            The graph of the site.
        :param metal_indices:
            The indices of the metal centers.
        :param cap_indices:
            The indices of the capping hydrogens.

        :return:
            The description: the residues, their keys, the coarse graph and
            the printable site spec.
        """

        sidechains = graph.copy()
        sidechains.remove_nodes_from(metal_indices)

        residues = []
        for component in nx.connected_components(sidechains):
            nodes = sorted(component)
            residues.append({
                'formula':
                self._formula(labels, nodes),
                'key':
                self._fragment_key(graph, nodes),
                'metals':
                sorted(metal for metal in metal_indices if any(
                    graph.has_edge(metal, node) for node in nodes)),
                'heavy':
                self._heavy_subgraph(labels, graph, nodes),
            })

        coarse = nx.Graph()

        for metal in metal_indices:
            coarse.add_node(('metal', metal), kind='metal', key=labels[metal])

        for index, residue in enumerate(residues):
            node = ('residue', index)
            coarse.add_node(node, kind='residue', key=residue['key'])
            for metal in residue['metals']:
                coarse.add_edge(node, ('metal', metal))

        return {
            'labels': labels,
            'graph': graph,
            'metal_indices': list(metal_indices),
            'cap_indices': list(cap_indices),
            'residues': residues,
            'coarse': coarse,
            'spec': self._site_spec(labels, metal_indices, residues),
        }

    @staticmethod
    def _formula(labels, nodes):
        """
        Returns the formula of one residue, for reading rather than for
        comparing: two residues can share a formula and still be different
        residues, which is what the key is for.

        :param labels:
            The element of every atom.
        :param nodes:
            The atoms of the fragment.

        :return:
            The formula.
        """

        counts = Counter(labels[node] for node in nodes)
        order = [element for element in ('C', 'H') if element in counts]
        order += sorted(element for element in counts
                        if element not in ('C', 'H'))

        return ''.join(
            f'{element}{counts[element]}' if counts[element] > 1 else element
            for element in order)

    @staticmethod
    def _fragment_key(graph, nodes):
        """
        Returns a canonical key for one residue.

        The hash is taken over the whole fragment, hydrogens included, so that
        a protonated carboxylate and a deprotonated one are not the same
        residue.

        :param graph:
            The graph of the site.
        :param nodes:
            The atoms of the fragment.

        :return:
            The key.
        """

        return nx.weisfeiler_lehman_graph_hash(graph.subgraph(nodes),
                                               node_attr='elem')

    @staticmethod
    def _heavy_subgraph(labels, graph, nodes):
        """
        Returns the heavy atom graph of one residue.

        The hydrogens are the bulk of the atoms and almost all of the symmetry,
        so the atom mapping is enumerated over the heavy atoms alone. How many
        hydrogens each of them carries goes on as an attribute, since without
        it a CH2 would map onto a CH3 and a protonated oxygen onto a bare one.

        :param labels:
            The element of every atom.
        :param graph:
            The graph of the site.
        :param nodes:
            The atoms of the fragment.

        :return:
            The subgraph, with an elem and an h_count attribute per node.
        """

        heavy = nx.Graph()

        for node in nodes:
            if labels[node] == 'H':
                continue
            hydrogens = [
                other for other in graph.neighbors(node) if labels[other] == 'H'
            ]
            heavy.add_node(node, elem=labels[node], h_count=len(hydrogens))

        for node in heavy:
            for other in graph.neighbors(node):
                if other in heavy:
                    heavy.add_edge(node, other)

        return heavy

    @staticmethod
    def _site_spec(labels, metal_indices, residues):
        """
        Returns what a site is made of, in a form two sites can be compared
        and reported by.

        :param labels:
            The element of every atom.
        :param metal_indices:
            The indices of the metal centers.
        :param residues:
            The residues of the site.

        :return:
            One entry per metal center and one for the residues themselves.
        """

        spec = {}

        for metal in metal_indices:
            bound = sorted(residue['key'] for residue in residues
                           if metal in residue['metals'])
            spec[f'{labels[metal]}{metal}'] = bound

        spec['residues'] = sorted(residue['key'] for residue in residues)
        spec['bridging'] = sorted(residue['key'] for residue in residues
                                  if len(residue['metals']) > 1)

        return spec

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

        query = self._extract_active_site(structure)
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
                        candidates.append(
                            (scored, name, mapping, rmsd, heavy_rmsd, ic_rmsd))

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

    def compare_to_templates(self,
                             structure,
                             mm_opt=True,
                             include_hydrogens=False):
        """
        Measures a structure against every template and reports, deciding
        nothing.

        Where match_enzyme_to_template picks one template and builds a force
        field from it, this walks the whole set and says how each of them
        stands: whether it holds the same atoms at all, whether it is wired
        the same way, and how far away it is in every region under both
        measures. It is what to reach for when a match came back None and the
        question is which template was closest, and by what.

        Every number it returns is also printed as a table, and the dictionary
        holds one entry per template, so a decision the thresholds cannot make
        can be made on the numbers afterwards.

        The hydrogens are left out of every measurement by default, since
        protonate does not place them reproducibly; include_hydrogens puts
        them back for this call whatever rmsd_heavy_atoms_only says. One
        geometry is measured per call, so call it twice to see both.

        :param structure:
            The path to a PDB or mmCIF file.
        :param mm_opt:
            Whether to relax the active site on a crude force field first,
            which is what takes the slack out of an unrelaxed structure.
        :param include_hydrogens:
            Whether the hydrogens count toward the RMSDs and the internal
            coordinates.

        :return:
            The findings, with one entry per template under 'templates'.
        """

        assert_msg_critical(
            len(self.templates) > 0,
            'MetalForceFieldManager.compare_to_templates: no templates '
            'loaded. Call load_ff_from_folder first.')

        query = self._extract_active_site(structure)
        active_site = query['active_site']
        graph = self._active_site_graph(active_site,
                                        query['connectivity_matrix'])
        described = self._describe_site(list(active_site['labels']), graph,
                                        active_site['metal_indices'],
                                        active_site['cap_indices'])

        if mm_opt:
            molecule = self._mm_relax(query)
            geometry = 'mm_relaxed'
        else:
            molecule = active_site['molecule']
            geometry = 'input'

        coordinates = molecule.get_coordinates_in_angstrom()
        heavy_only = not include_hydrogens
        composition = sorted(active_site['labels'])

        findings = {}

        for name, template in self.templates.items():
            entry = {
                'status': 'measured',
                'mapping': None,
                'n_coarse_mappings': 0,
                'n_mappings': 0,
                'metal_bonds': None,
                'regions': {},
            }

            if sorted(template['labels']) != composition:
                # not the same atoms, so there is nothing to map onto
                entry['status'] = 'composition'
                findings[name] = entry
                continue

            # which residue coordinates which metal, with nothing said about
            # how many atoms of it do the coordinating
            coarse = self._coarse_mappings(template, described)

            if not coarse:
                # what it holds instead is printed from the template
                entry['status'] = 'spec'
                findings[name] = entry
                continue

            maps = []
            for coarse_mapping in coarse:
                maps.extend(
                    self._heavy_atom_maps(template, described, coarse_mapping))

            if not maps:
                # what it holds instead is printed from the template
                entry['status'] = 'spec'
                findings[name] = entry
                continue

            entry['n_coarse_mappings'] = len(coarse)
            entry['n_mappings'] = len(maps)

            heavy_map, rot, trans = self._best_heavy_map(
                template, maps, coordinates)
            mapping = self._complete_hydrogens(template, described, heavy_map,
                                               coordinates, rot, trans)

            entry['mapping'] = mapping
            entry['metal_bonds'] = self._metal_bond_summary(
                template, described, mapping, coordinates)

            for region in self.RMSD_REGIONS:
                entry['regions'][region] = self._measure_region(
                    template, mapping, coordinates, region, heavy_only)

            findings[name] = entry

        results = {
            'structure': str(structure),
            'geometry': geometry,
            'include_hydrogens': include_hydrogens,
            'active_site': active_site,
            'connectivity_matrix': query['connectivity_matrix'],
            'binding_modes': query['binding_modes'],
            'topology': query['topology'],
            'positions': query['positions'],
            # the geometry every number above was measured on, which is the
            # relaxed one rather than the site's own when mm_opt is set
            'molecule': molecule,
            'description': described,
            'templates': findings,
        }

        self._print_comparison(results)

        return results

    def select_template(self, comparison, template=None, criteria=None):
        """
        Picks the template a force field should be built from, and says why.

        A template is usable only if it maps onto every atom of the site

        If no template is specified, every template that maps completely is held to
        the the specified criteria. The passing templates are ranked and the best one is picked

        If, the criteria are measured and reported but do
        not decide, since naming a template is a statement that it is the
        right one.

        :param comparison:
            The findings of compare_to_templates, or a structure to run it on.
        :param template:
            The name of the template to use, or None to choose one.
        :param criteria:
            The criteria to hold a template to: 'tight', 'loose', a set of the
            same shape, or None for selection_criteria.

        :return:
            The decision, with the chosen name under 'name', None when nothing
            passed, and what was made of every template under 'verdicts'.
        """

        comparison = self._comparison(comparison)
        criteria_name, criteria = self._selection_criteria(criteria)

        decision = {
            'name': None,
            'entry': None,
            'forced': template is not None,
            'criteria': criteria,
            'criteria_name': criteria_name,
            'score': None,
            'verdicts': {},
            'scores': {},
            'candidates': [],
        }

        for name, entry in comparison['templates'].items():
            verdict = self._selection_verdict(comparison, entry, criteria)
            decision['verdicts'][name] = verdict
            decision['scores'][name] = self._selection_score(entry)

        if template is not None:
            assert_msg_critical(
                template in comparison['templates'],
                'MetalForceFieldManager.select_template: no template named '
                f'{template} was compared. Loaded: '
                f'{sorted(comparison["templates"])}')

            entry = comparison['templates'][template]
            verdict = decision['verdicts'][template]

            assert_msg_critical(
                self._maps_every_atom(comparison, entry),
                f'MetalForceFieldManager.select_template: template {template} '
                f'does not map onto every atom of the site ({verdict}), so '
                'its parameters cannot be transferred. Build this site with '
                'MetalSiteForceFieldBuilder.')

            if verdict is not None:
                self.ostream.print_warning(
                    f'Template {template} was named rather than chosen, and '
                    f'it is outside the {criteria_name} criteria: {verdict}. '
                    'Transferring it anyway.')
                self.ostream.flush()

            decision['name'] = template
            decision['entry'] = entry
            decision['score'] = decision['scores'][template]

        else:
            passed = [
                name for name, verdict in decision['verdicts'].items()
                if verdict is None
            ]
            decision['candidates'] = sorted(
                passed, key=lambda name: decision['scores'][name])

            if decision['candidates']:
                name = decision['candidates'][0]
                decision['name'] = name
                decision['entry'] = comparison['templates'][name]
                decision['score'] = decision['scores'][name]

        self._print_selection(comparison, decision)

        return decision

    def build_ff_from_comparison(self,
                                 comparison,
                                 template=None,
                                 criteria=None):
        """
        Builds a force field for a structure out of the template that fits it.

        This is the half of a builder run that QM was paid for: the fitted
        metal bonds and angles and the RESP charges are taken from a template
        that describes the same site, and everything else is built for the
        structure in front of us. Which template that is comes from
        select_template, on the numbers compare_to_templates measured.

        :param comparison:
            The findings of compare_to_templates, or a structure to run it on.
        :param template:
            The name of the template to use, or None to choose one.
        :param criteria:
            The criteria to hold a template to: 'tight', 'loose', a set of the
            same shape, or None for selection_criteria.

        :return:
            The results dictionary, or None when no template was good enough.
        """

        # todo accept pdb names?
        comparison = self._comparison(comparison)
        decision = self.select_template(comparison, template, criteria)

        # todo should this logic be moved into select template
        if decision['name'] is None:
            closest = min(decision['scores'],
                          key=lambda name: decision['scores'][name],
                          default=None)

            if closest is not None and math.isfinite(
                    decision['scores'][closest]):
                self.ostream.print_info(
                    f'No template is within the {decision["criteria_name"]} '
                    f'criteria. The closest is {closest}: '
                    f'{decision["verdicts"][closest]}.')
            else:
                self.ostream.print_info(
                    'No template describes this site: none of them maps onto '
                    'all of its atoms.')

            self.ostream.print_info(
                "Returning None. Try criteria='loose', or build this site "
                'with MetalSiteForceFieldBuilder.')
            self.ostream.flush()
            return None

        name = decision['name']
        template = self.templates[name]
        entry = decision['entry']
        active_site = comparison['active_site']

        self.ostream.print_info(
            f'Building the force field from {name}, measured on the '
            f'{comparison["geometry"]} geometry. Transferring its metal terms '
            'and charges.')
        self.ostream.flush()

        forcefield = self._transfer(template, entry['mapping'], active_site,
                                    comparison['connectivity_matrix'])

        enzyme_system = None
        if self.build_enzyme_system:
            enzyme_system, _ = self.builder.create_enzyme_system(
                comparison['topology'], active_site, forcefield,
                forcefield.partial_charges)

        self.builder._print_metal_parameters(active_site, forcefield)

        return {
            'template': name,
            'forcefield': forcefield,
            'mapping': entry['mapping'],
            'forced': decision['forced'],
            'criteria': decision['criteria_name'],
            'score': decision['score'],
            'verdicts': decision['verdicts'],
            'candidates': decision['candidates'],
            'regions': entry['regions'],
            'metal_bonds': entry['metal_bonds'],
            'geometry': comparison['geometry'],
            'geometry_kind': template['geometry_kind'],
            'topology': comparison['topology'],
            'positions': comparison['positions'],
            'binding_modes': comparison['binding_modes'],
            'active_site': active_site,
            'connectivity_matrix': comparison['connectivity_matrix'],
            'enzyme_system': enzyme_system,
        }

    def _comparison(self, comparison):
        """
        Takes either the findings of compare_to_templates or a structure to
        measure, and returns the findings.

        :param comparison:
            The findings, or the path to a PDB or mmCIF file.

        :return:
            The findings.
        """

        if isinstance(comparison, (str, Path)):
            return self.compare_to_templates(comparison)

        assert_msg_critical(
            isinstance(comparison, dict) and 'templates' in comparison,
            'MetalForceFieldManager: expected the findings of '
            'compare_to_templates or a structure to measure, got '
            f'{type(comparison).__name__}')

        return comparison

    def _selection_criteria(self, criteria):
        """
        Resolves the criteria a template is held to.

        :param criteria:
            A name in SELECTION_CRITERIA, a set of the same shape, or None for
            selection_criteria.

        :return:
            The tuple of the name to report it under and the set itself.
        """

        if criteria is None:
            criteria = self.selection_criteria

        if isinstance(criteria, str):
            assert_msg_critical(
                criteria in self.SELECTION_CRITERIA,
                'MetalForceFieldManager: the criteria must be one of '
                f'{sorted(self.SELECTION_CRITERIA)} or a set of thresholds, '
                f'got {criteria}')
            return criteria, self.SELECTION_CRITERIA[criteria]

        assert_msg_critical(
            isinstance(criteria, dict),
            'MetalForceFieldManager: the criteria must be one of '
            f'{sorted(self.SELECTION_CRITERIA)} or a set of thresholds, got '
            f'{type(criteria).__name__}')

        for region, thresholds in criteria.items():
            assert_msg_critical(
                region in self.RMSD_REGIONS,
                f'MetalForceFieldManager: {region} is not a region the '
                f'criteria can name; expected one of {list(self.RMSD_REGIONS)}')
            for name, limits in (thresholds or {}).items():
                assert_msg_critical(
                    name in self.IC_TYPES,
                    f'MetalForceFieldManager: {name} is not an internal '
                    'coordinate type the criteria can name; expected one of '
                    f'{list(self.IC_TYPES)}')
                for measure in (limits or {}):
                    assert_msg_critical(
                        measure in ('rms', 'max'),
                        f'MetalForceFieldManager: {measure} is not a measure '
                        "the criteria can name; expected 'rms' or 'max'")

        return 'custom', criteria

    @staticmethod
    def _maps_every_atom(comparison, entry):
        """
        Says whether a template covers every atom of the site.

        A template that holds other atoms, or coordinates them differently,
        never gets as far as a mapping; this also catches a mapping that came
        back incomplete, which would leave part of a site unparameterized.

        :param comparison:
            The findings of compare_to_templates.
        :param entry:
            What it measured for the template.

        :return:
            True when every atom of the site is mapped onto exactly once.
        """

        if entry['status'] != 'measured' or entry['mapping'] is None:
            return False

        atoms = len(comparison['active_site']['labels'])
        mapping = entry['mapping']

        return (len(mapping) == atoms
                and sorted(mapping.values()) == list(range(atoms)))

    def _selection_verdict(self, comparison, entry, criteria):
        """
        Holds one template to the criteria, region by region.

        :param comparison:
            The findings of compare_to_templates.
        :param entry:
            What it measured for the template.
        :param criteria:
            The criteria, as _selection_criteria resolves them.

        :return:
            A description of what stands in the way, or None when nothing
            does.
        """

        if entry['status'] == 'composition':
            return 'different atoms'

        if entry['status'] == 'spec':
            return 'different coordination'

        if not self._maps_every_atom(comparison, entry):
            return 'incomplete mapping'

        for region in self.RMSD_REGIONS:
            thresholds = criteria.get(region)
            if not thresholds:
                continue

            found = entry['regions'].get(region)
            if found is None:
                return f'{region} not measured'

            # a criterion that could not be evaluated is not one that was
            # passed, so the region is held to strictly here
            violation = self._ic_violation(found['ic_rmsd'],
                                           thresholds,
                                           strict=True)
            if violation is not None:
                return f'{region} {violation}'

        return None

    def _selection_score(self, entry):
        """
        Returns what several templates that all pass are ranked on.

        :param entry:
            What compare_to_templates measured for the template.

        :return:
            The measure named by SELECTION_RANKED_ON, or infinity where it was
            not measured.
        """

        region, ic_type, measure = self.SELECTION_RANKED_ON

        found = entry['regions'].get(region)
        if found is None or found['ic_rmsd'] is None:
            return math.inf

        found = found['ic_rmsd'].get(ic_type)
        if found is None:
            return math.inf

        return found[measure]

    def _extract_active_site(self, structure):
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

    def _coarse_mappings(self, template, query):
        """
        Matches the metals and the residues of two sites, ignoring how the
        residues grip the metals.

        The graph this runs on has one node per metal and one per residue and
        carries no denticity at all, so a monodentate carboxylate and a
        bidentate one are the same node with the same edge. It is also tiny -
        eight nodes for a binuclear site - which is what makes this cheap
        where an isomorphism of the whole site is not.

        :param template:
            The description of the template site.
        :param query:
            The description of the queried site.

        :return:
            One mapping of coarse nodes per way the two sites line up.
        """

        def node_match(a, b):
            return a['kind'] == b['kind'] and a['key'] == b['key']

        matcher = GraphMatcher(template['coarse'],
                               query['coarse'],
                               node_match=node_match)

        return list(matcher.isomorphisms_iter())

    def _heavy_atom_maps(self, template, query, coarse_mapping):
        """
        Builds the heavy atom mappings that one coarse mapping allows.

        Each residue is mapped onto the residue the coarse level paired it
        with, atom by atom and on its own, so the symmetry that survives is
        the symmetry of a sidechain: the two oxygens of a carboxylate, the two
        hydrogens of a CB. The hydrogens are left out here and put back once
        the geometry has chosen between what remains.

        :param template:
            The description of the template site.
        :param query:
            The description of the queried site.
        :param coarse_mapping:
            One mapping of coarse nodes.

        :return:
            The heavy atom mappings, from template index to query index.
        """

        def node_match(a, b):
            return a['elem'] == b['elem'] and a['h_count'] == b['h_count']

        metals = {}
        per_residue = []

        for node, image in coarse_mapping.items():
            if node[0] == 'metal':
                metals[node[1]] = image[1]
                continue

            first = template['residues'][node[1]]['heavy']
            second = query['residues'][image[1]]['heavy']

            matcher = GraphMatcher(first, second, node_match=node_match)
            found = list(matcher.isomorphisms_iter())

            if not found:
                # the coarse keys agreed, so this should not happen; a site
                # that manages it is not one to guess about
                return []

            per_residue.append(found)

        maps = []

        for combination in product(*per_residue):
            atom_map = dict(metals)
            for residue_map in combination:
                atom_map.update(residue_map)
            maps.append(atom_map)

            if len(maps) >= self.max_mappings:
                self.ostream.print_warning(
                    f'Reached the limit of {self.max_mappings} atom mappings; '
                    'the best of the ones built is used, which need not be '
                    'the best there is')
                self.ostream.flush()
                break

        return maps

    def _best_heavy_map(self, template, maps, coordinates):
        """
        Picks the heavy atom mapping that superimposes the two sites best.

        :param template:
            The template being compared to.
        :param maps:
            The heavy atom mappings to choose between.
        :param coordinates:
            The coordinates of the queried site, in Angstrom.

        :return:
            The tuple of the best mapping and the superposition it gives.
        """

        reference = template['molecule'].get_coordinates_in_angstrom()
        best = None

        for atom_map in maps:
            order = sorted(atom_map)
            moved = coordinates[[atom_map[index] for index in order]]

            rmsd, rot, trans = svd_superimpose(moved, reference[order])

            if best is None or rmsd < best[0]:
                best = (rmsd, atom_map, rot, trans)

        return best[1], best[2], best[3]

    def _complete_hydrogens(self, template, query, heavy_map, coordinates, rot,
                            trans):
        """
        Extends a heavy atom mapping over the hydrogens.

        The hydrogens of one heavy atom are equivalent to each other in
        everything the comparison measures, so which of them goes where is
        settled by proximity under the superposition the heavy atoms already
        found rather than by another round of enumeration. Capping hydrogens
        are kept apart from ordinary ones, since a cap stands for an alpha
        carbon and is not equivalent to anything else.

        :param template:
            The description of the template site.
        :param query:
            The description of the queried site.
        :param heavy_map:
            The mapping of the heavy atoms.
        :param coordinates:
            The coordinates of the queried site, in Angstrom.
        :param rot:
            The rotation of the superposition.
        :param trans:
            The translation of the superposition.

        :return:
            The mapping over every atom.
        """

        reference = template['molecule'].get_coordinates_in_angstrom()
        aligned = np.matmul(coordinates, rot) + trans

        template_caps = set(template['cap_indices'])
        query_caps = set(query['cap_indices'])

        atom_map = dict(heavy_map)

        def hydrogens(labels, graph, node):
            return [
                other for other in graph.neighbors(node) if labels[other] == 'H'
            ]

        for node, image in heavy_map.items():
            first = hydrogens(template['labels'], template['graph'], node)
            second = hydrogens(query['labels'], query['graph'], image)

            groups = [(first, second)]

            first_caps = [index for index in first if index in template_caps]
            second_caps = [index for index in second if index in query_caps]

            if first_caps and len(first_caps) == len(second_caps):
                # a cap stands for an alpha carbon, so it is not free to swap
                # with the ordinary hydrogens of the same carbon
                groups = [
                    (first_caps, second_caps),
                    ([index for index in first if index not in template_caps],
                     [index for index in second if index not in query_caps]),
                ]

            for ours, theirs in groups:
                if len(ours) != len(theirs):
                    # nothing sensible to pair up; leave them out rather than
                    # invent a correspondence
                    continue

                best = None
                for order in permutations(theirs):
                    distance = sum(
                        float(np.linalg.norm(reference[index] - aligned[other]))
                        for index, other in zip(ours, order))
                    if best is None or distance < best[0]:
                        best = (distance, order)

                atom_map.update(dict(zip(ours, best[1])))

        return atom_map

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

        reference = template['molecule'].get_coordinates_in_angstrom()
        region = self._region_indices(template)
        heavy = [index for index in region if template['labels'][index] != 'H']

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

    def _rmsd_indices(self, template, region=None, heavy_only=None):
        """
        Returns the template indices every RMSD is measured over, which is the
        region less the hydrogens when they are being left out.

        :param template:
            The template.
        :param region:
            The region to take, or None for rmsd_region.
        :param heavy_only:
            Whether to leave the hydrogens out, or None for
            rmsd_heavy_atoms_only.

        :return:
            The indices, in order.
        """

        if heavy_only is None:
            heavy_only = self.rmsd_heavy_atoms_only

        indices = self._region_indices(template, region)

        if not heavy_only:
            return indices

        return [index for index in indices if template['labels'][index] != 'H']

    def _region_indices(self, template, region=None):
        """
        Returns the template indices of the region an RMSD is measured over.

        :param template:
            The template.
        :param region:
            The region to take, or None for rmsd_region.

        :return:
            The indices, in order.
        """

        if region is None:
            region = self.rmsd_region

        if region == 'active_site':
            return list(range(len(template['labels'])))

        graph = template['graph']

        if region == 'metal_beta_carbons':
            assert_msg_critical(
                len(template['beta_carbon_indices']) > 0,
                f'MetalForceFieldManager: template {template["name"]} does '
                'not say which of its atoms are beta carbons. Its force field '
                'was written before annotate_atoms recorded them, so rebuild '
                'it or measure over another region.')

            return sorted(
                set(template['metal_indices'])
                | set(template['beta_carbon_indices']))

        shell = set()

        for metal in template['metal_indices']:
            shell.update(
                nx.single_source_shortest_path_length(
                    graph, metal, cutoff=self.metal_shell_bonds))

        assert_msg_critical(
            len(shell) > len(template['metal_indices']),
            'MetalForceFieldManager: the metal centers of template '
            f'{template["name"]} have nothing bonded to them within '
            f'{self.metal_shell_bonds} bond(s)')

        return sorted(shell)

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
        matched = self._criteria_flags(xyz_ok, violation is None)[criterion]

        if criterion == 'ic_rmsd':
            verdict = violation
        elif criterion == 'ic_and_xyz':
            # name whichever failed first, so the report says what to loosen
            verdict = 'geometry off' if not xyz_ok else violation
        else:
            verdict = 'geometry and IC off'

        return matched, 'match' if matched else verdict, ic_rmsd

    @staticmethod
    def _criteria_flags(xyz_ok, ic_ok):
        """
        Combines the two measurements into the four criteria.

        :param xyz_ok:
            Whether the cartesian RMSD is within its threshold.
        :param ic_ok:
            Whether the internal coordinates are within theirs.

        :return:
            What each criterion in MATCH_CRITERIA makes of the pair.
        """

        return {
            'xyz_rmsd': xyz_ok,
            'ic_rmsd': ic_ok,
            'ic_or_xyz': xyz_ok or ic_ok,
            'ic_and_xyz': xyz_ok and ic_ok,
        }

    def _measure_region(self,
                        template,
                        mapping,
                        coordinates,
                        region,
                        heavy_only=None):
        """
        Measures one region of an active site against a template.

        :param template:
            The template being compared to.
        :param mapping:
            The mapping from template index to active site index.
        :param coordinates:
            The coordinates of the active site, in Angstrom.
        :param region:
            The region to measure over.
        :param heavy_only:
            Whether to leave the hydrogens out, or None for
            rmsd_heavy_atoms_only.

        :return:
            The atom count, the cartesian RMSDs over the whole region and
            over its heavy atoms, the internal coordinate deviations and what
            each criterion makes of them.
        """

        if heavy_only is None:
            heavy_only = self.rmsd_heavy_atoms_only

        reference = template['molecule'].get_coordinates_in_angstrom()
        indices = self._region_indices(template, region)
        heavy = [index for index in indices if template['labels'][index] != 'H']

        order = [mapping[index] for index in range(len(reference))]
        moved = coordinates[order]

        rmsd, _, _ = svd_superimpose(moved[indices], reference[indices])
        heavy_rmsd, _, _ = svd_superimpose(moved[heavy], reference[heavy])

        ic_rmsd = self._ic_rmsd(template, mapping, coordinates, region,
                                heavy_only)

        scored = heavy_rmsd if heavy_only else rmsd
        flags = self._criteria_flags(scored <= self.geometry_rmsd_threshold,
                                     self._ic_violation(ic_rmsd) is None)

        return {
            'atoms': len(self._rmsd_indices(template, region, heavy_only)),
            'rmsd': rmsd,
            'rmsd_heavy': heavy_rmsd,
            'ic_rmsd': ic_rmsd,
            'criteria': flags,
        }

    def _ic_rmsd(self,
                 template,
                 mapping,
                 coordinates,
                 region=None,
                 heavy_only=None):
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
        :param region:
            The region to measure over, or None for rmsd_region.
        :param heavy_only:
            Whether to leave the hydrogens out, or None for
            rmsd_heavy_atoms_only.

        :return:
            The deviations, as get_ic_rmsd reports them, or None when they
            could not be measured.
        """

        indices = self._rmsd_indices(template, region, heavy_only)
        labels = [template['labels'][index] for index in indices]

        reference = template['molecule'].get_coordinates_in_angstrom()
        order = [mapping[index] for index in indices]

        mapped = Molecule(labels, coordinates[order], 'angstrom')
        wanted = Molecule(labels, reference[indices], 'angstrom')

        ic_rmsd = OptimizationDriver.get_ic_rmsd(mapped, wanted)

        if isinstance(ic_rmsd, str):
            # get_ic_rmsd reports its refusals as a string rather than raising
            self.ostream.print_warning(
                f'Could not measure the internal coordinates against template '
                f'{template["name"]}: {ic_rmsd}')
            self.ostream.flush()
            return None

        return ic_rmsd

    def _ic_violation(self, ic_rmsd, thresholds=None, strict=False):
        """
        Checks internal coordinate deviations against a set of thresholds.

        A type left at None, or a threshold left out of one, is not checked.
        The dihedrals are None by default: the ones about a metal center are
        floppy enough that they say nothing about whether two sites are the
        same.

        :param ic_rmsd:
            The deviations, as get_ic_rmsd reports them.
        :param thresholds:
            The thresholds to check against, or None for ic_rmsd_thresholds.
        :param strict:
            Whether a type that is asked for but was not measured counts as a
            violation. Off for the matching path, which reads a missing type
            as a region that holds none of them; on for the selection path,
            where a criterion that could not be evaluated is not one that was
            passed.

        :return:
            A description of the first threshold that was exceeded, or None
            when every one of them holds.
        """

        if thresholds is None:
            thresholds_for = self.ic_rmsd_thresholds.get
        else:
            thresholds_for = thresholds.get

        if ic_rmsd is None:
            return 'no internal coords'

        for name, unit in self.IC_TYPES.items():
            thresholds = thresholds_for(name)
            if not thresholds:
                continue

            found = ic_rmsd.get(name)
            if found is None:
                # the region holds no internal coordinate of this type
                if strict:
                    return f'no {name} measured'
                continue

            for measure in ('rms', 'max'):
                limit = thresholds.get(measure)
                if limit is None:
                    continue
                if found[measure] > limit:
                    return (f'{name} {measure} {found[measure]:.2f} > '
                            f'{limit:.2f} {unit}')

        return None

    def _metal_bond_summary(self, template, query, mapping, coordinates):
        """
        Compares the metal-ligand bonds of a template with those of a site.

        Only the bonds the two of them agree on are measured. Which atoms of a
        residue reach a metal is exactly what a comparison is meant not to
        turn on - it is a distance cutoff on an unrelaxed structure, not
        chemistry - so a contact one side makes and the other does not is
        counted and reported rather than allowed to fail the comparison.

        :param template:
            The template being compared to.
        :param query:
            The description of the queried site.
        :param mapping:
            The mapping from template index to active site index.
        :param coordinates:
            The coordinates of the queried site, in Angstrom.

        :return:
            The largest difference over the shared bonds, how many were
            shared, and how many either side makes alone.
        """

        reference = template['molecule'].get_coordinates_in_angstrom()
        bonds, _ = self._metal_keys(template)

        deviation = 0.0
        shared = 0
        template_only = 0
        mapped = set()

        for i, j in bonds:
            first, second = mapping[i], mapping[j]
            mapped.add(frozenset((first, second)))

            if not query['graph'].has_edge(first, second):
                template_only += 1
                continue

            shared += 1
            expected = np.linalg.norm(reference[i] - reference[j])
            found = np.linalg.norm(coordinates[first] - coordinates[second])
            deviation = max(deviation, abs(found - expected))

        metals = set(query['metal_indices'])
        query_only = sum(
            1 for first, second in query['graph'].edges()
            if ({first, second} & metals) and frozenset((first,
                                                         second)) not in mapped)

        return {
            'deviation': deviation,
            'shared': shared,
            'template_only': template_only,
            'query_only': query_only,
        }

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

        reference = template['molecule'].get_coordinates_in_angstrom()
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

    # todo come up with a more descriptive name for this
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
        expected = int(active_site['molecule'].get_charge())
        if abs(total - expected) > 1.0e-3:
            self.ostream.print_warning(
                f'The transferred charges sum to {total:+.3f}, but the active '
                f'site charge is {expected:+d}')

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
                    self.ostream.print_header(self._param(label, 'not checked'))
                    continue
                limits = ', '.join(f'{measure} {thresholds[measure]:.2f}'
                                   for measure in ('rms', 'max')
                                   if thresholds.get(measure))
                self.ostream.print_header(self._param(label, limits))

        self.ostream.print_header(
            self._param('metal bond tolerance',
                        f'{self.metal_bond_tolerance:.2f} A'))
        self.ostream.print_blank()
        self.ostream.print_info(
            f'Residues: {", ".join(active_site["residues"])}')
        self.ostream.flush()

    def _print_comparison(self, results):
        """
        Prints everything compare_to_templates measured.

        One table of numbers and one of verdicts per template that could be
        measured, and a closing summary ranking the templates by the region
        that is configured, so the closest one is visible without reading
        every table.

        :param results:
            The findings of compare_to_templates.
        """

        active_site = results['active_site']
        metals = ', '.join(active_site['labels'][index]
                           for index in active_site['metal_indices'])

        self.ostream.print_blank()
        self.ostream.print_header('Comparison against every template')
        self.ostream.print_header(33 * '-')
        self.ostream.print_header(
            self._param('structure',
                        Path(results['structure']).name))
        self.ostream.print_header(
            self._param('active site atoms', len(active_site['labels'])))
        self.ostream.print_header(self._param('metal centers', metals))
        self.ostream.print_header(self._param('geometry', results['geometry']))
        self.ostream.print_header(
            self._param(
                'measured over',
                'all atoms' if results['include_hydrogens'] else 'heavy atoms'))
        self.ostream.print_header(
            self._param('templates', len(results['templates'])))
        self.ostream.print_blank()
        self.ostream.print_info(
            f'Residues: {", ".join(active_site["residues"])}')
        self.ostream.print_info(
            f'Thresholds: RMSD {self.geometry_rmsd_threshold:.2f} A, '
            f'metal bond {self.metal_bond_tolerance:.2f} A, '
            f'IC {self._ic_threshold_summary()}')

        self._print_spec('the structure holds', results['description'])

        for name, entry in results['templates'].items():
            self._print_template_comparison(name, entry)

        self._print_comparison_summary(results)

    def _ic_threshold_summary(self):
        """
        Returns the internal coordinate thresholds on one line.

        :return:
            The summary.
        """

        parts = []

        for name, unit in self.IC_TYPES.items():
            thresholds = self.ic_rmsd_thresholds.get(name)
            if not thresholds:
                continue
            limits = '/'.join(f'{thresholds[measure]:g}'
                              for measure in ('rms', 'max')
                              if thresholds.get(measure))
            parts.append(f'{name} {limits} {unit}')

        return ', '.join(parts) if parts else 'none'

    def _print_spec(self, title, description):
        """
        Prints what a site is made of: which residues it holds and which of
        them coordinate which metal.

        The residues are named by their formula, which is for reading; two
        sites are compared on the keys behind them.

        :param title:
            What the block is describing.
        :param description:
            The description of the site.
        """

        labels = description['labels']

        def named(residues):
            return ', '.join(
                sorted(f'{residue["formula"]}/{residue["key"][:6]}'
                       for residue in residues))

        self.ostream.print_blank()
        self.ostream.print_info(f'Site spec, {title}:')

        for metal in description['metal_indices']:
            bound = [
                residue for residue in description['residues']
                if metal in residue['metals']
            ]
            self.ostream.print_info(f'  {labels[metal]}{metal}: {named(bound)}')

        bridging = [
            residue for residue in description['residues']
            if len(residue['metals']) > 1
        ]
        if bridging:
            self.ostream.print_info(f'  bridging: {named(bridging)}')

        self.ostream.flush()

    def _print_template_comparison(self, name, entry):
        """
        Prints the numbers and the verdicts of one template.

        :param name:
            The name of the template.
        :param entry:
            What compare_to_templates measured for it.
        """

        self.ostream.print_blank()

        if entry['status'] == 'composition':
            self.ostream.print_info(
                f'{name}: holds different atoms, so nothing was measured.')
            self.ostream.flush()
            return

        if entry['status'] == 'spec':
            self.ostream.print_info(
                f'{name}: coordinates a different set of residues, so '
                'nothing was measured.')
            self._print_spec(f'{name} holds', self.templates[name])
            self.ostream.flush()
            return

        bonds = entry['metal_bonds']
        summary = (f'{bonds["shared"]} metal bond(s) shared, within '
                   f'{bonds["deviation"]:.3f} A')
        if bonds['template_only']:
            summary += f', {bonds["template_only"]} only in the template'
        if bonds['query_only']:
            summary += f', {bonds["query_only"]} only in the structure'

        self.ostream.print_info(
            f'{name}: {entry["n_mappings"]} atom mapping(s) from '
            f'{entry["n_coarse_mappings"]} coarse mapping(s), {summary}')
        self.ostream.print_blank()

        row = '{:>19} | {:>5} | {:>7} | {:>7} | {:>14} | {:>14} | {:>14}'
        self.ostream.print_header(
            row.format('region', 'atoms', 'RMSD', 'heavy', 'bonds rms/max',
                       'angles rms/max', 'dihed rms/max'))
        self.ostream.print_header(98 * '-')

        for region, found in entry['regions'].items():
            self.ostream.print_header(
                row.format(region, found['atoms'], f'{found["rmsd"]:.3f}',
                           f'{found["rmsd_heavy"]:.3f}',
                           self._ic_cell(found['ic_rmsd'], 'bonds'),
                           self._ic_cell(found['ic_rmsd'], 'angles'),
                           self._ic_cell(found['ic_rmsd'], 'dihedrals')))

        self.ostream.print_blank()

        row = '{:>19} | {:>9} | {:>9} | {:>10} | {:>11}'
        self.ostream.print_header(row.format('region', *self.MATCH_CRITERIA))
        self.ostream.print_header(68 * '-')

        for region, found in entry['regions'].items():
            flags = found['criteria']
            self.ostream.print_header(
                row.format(
                    region, *[
                        'yes' if flags[criterion] else '-'
                        for criterion in self.MATCH_CRITERIA
                    ]))

        self.ostream.print_blank()
        self.ostream.flush()

    @staticmethod
    def _ic_cell(ic_rmsd, name):
        """
        Formats one internal coordinate type for a table cell.

        :param ic_rmsd:
            The deviations, as get_ic_rmsd reports them.
        :param name:
            The type to format.

        :return:
            The cell.
        """

        if ic_rmsd is None:
            return ''

        found = ic_rmsd.get(name)

        if found is None:
            return ''

        return f'{found["rms"]:.2f} / {found["max"]:.2f}'

    def _print_comparison_summary(self, results):
        """
        Ranks the templates by the region that is configured.

        :param results:
            The findings of compare_to_templates.
        """

        region = self.rmsd_region
        heavy = not results['include_hydrogens']

        def scored(entry):
            found = entry['regions'].get(region)
            if found is None:
                return math.inf
            return found['rmsd_heavy'] if heavy else found['rmsd']

        order = sorted(results['templates'].items(),
                       key=lambda item: scored(item[1]))

        self.ostream.print_blank()
        title = f'Ranked on the {region} region'
        self.ostream.print_header(title)
        self.ostream.print_header(len(title) * '-')

        row = '{:>24} | {:>14} | {:>9} | {:>10} | {:>10}'
        self.ostream.print_header(
            row.format('template', 'status', 'RMSD', 'metal bond', 'criteria'))
        self.ostream.print_header(78 * '-')

        for name, entry in order:
            found = entry['regions'].get(region)
            if found is None:
                self.ostream.print_header(
                    row.format(name[:24], entry['status'], '', '', ''))
                continue
            met = sum(1 for passed in found['criteria'].values() if passed)
            self.ostream.print_header(
                row.format(name[:24], entry['status'], f'{scored(entry):.3f}',
                           f'{entry["metal_bonds"]["deviation"]:.3f}',
                           f'{met} of {len(found["criteria"])}'))

        self.ostream.print_blank()
        self.ostream.flush()

    def _print_selection(self, comparison, decision):
        """
        Prints how every template stands against the criteria, and which one
        was taken.

        The whole field is printed rather than the winner alone: whether the
        others are near misses or a long way off is what says how much the
        chosen one is worth.

        :param comparison:
            The findings of compare_to_templates.
        :param decision:
            The decision, as select_template makes it.
        """

        regions = [
            region for region in self.RMSD_REGIONS
            if decision['criteria'].get(region)
        ]

        self.ostream.print_blank()
        title = f'Choosing a template on the {decision["criteria_name"]} criteria'
        self.ostream.print_header(title)
        self.ostream.print_header(len(title) * '-')
        self.ostream.print_blank()

        for region in regions:
            thresholds = decision['criteria'][region]
            # only the measures the set actually holds, since either of
            # them may be left out of one
            measures = {
                name:
                ' / '.join(f'{measure} {limit:.2f}'
                           for measure, limit in given.items()
                           if limit is not None)
                for name, given in thresholds.items() if given
            }
            limits = '; '.join(f'{name} {shown} {self.IC_TYPES[name]}'
                               for name, shown in measures.items())
            self.ostream.print_header(
                self._param(region, limits, value_width=44))

        self.ostream.print_blank()

        # one column per region the criteria name, so a custom set of them
        # prints as readably as the two that come with the class
        row = ' | '.join(['{:>22}'] + ['{:>13}'] * len(regions) +
                         ['{:>26}', '{:>5}'])
        header = row.format('template', *[region[:13] for region in regions],
                            'verdict', 'taken')
        self.ostream.print_header(header)
        self.ostream.print_header(len(header) * '-')

        order = sorted(comparison['templates'],
                       key=lambda name: (decision['verdicts'][name] is not None,
                                         decision['scores'][name], name))

        for name in order:
            entry = comparison['templates'][name]
            cells = []
            for region in regions:
                found = entry['regions'].get(region)
                cells.append('' if found is
                             None else self._ic_cell(found['ic_rmsd'], 'bonds'))

            verdict = decision['verdicts'][name] or 'within the criteria'
            self.ostream.print_header(
                row.format(name[:22], *cells, verdict[:26],
                           'yes' if name == decision['name'] else ''))

        self.ostream.print_blank()

        if decision['name'] is None:
            self.ostream.print_info('No template was taken.')
        elif decision['forced']:
            self.ostream.print_info(
                f'{decision["name"]} was named rather than chosen, so the '
                'criteria were measured but did not decide.')
        else:
            ranked = ' '.join(self.SELECTION_RANKED_ON)
            self.ostream.print_info(
                f'{len(decision["candidates"])} of '
                f'{len(comparison["templates"])} template(s) are within the '
                f'criteria. Taking {decision["name"]}, whose {ranked} of '
                f'{decision["score"]:.3f} is the lowest of them.')

        self.ostream.print_blank()
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

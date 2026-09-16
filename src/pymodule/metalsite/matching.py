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

"""
Describing a metal site, matching two of them, and measuring one against
the other.

The manager's algorithmic content, and pure computation over dictionaries
and graphs: nothing here holds state, touches the disk or knows what a
builder is.

Two rules that keep being worth restating, because both look like
something to simplify away:

  - Matching is two level -- residues onto residues first, atoms onto
    atoms second -- and denticity is deliberately not part of it. The
    whole-site isomorphism this replaced found 16 maps in 0.002 s against
    512 in 0.395 s, for the same best RMSD. Do not bring it back.

  - The atom mapping is chosen on the whole site whatever region is being
    measured. A region restricts what is measured, never what is matched;
    choosing a different mapping per region would compare each region
    against a differently oriented template.

The three settings the manager carries are keyword arguments here, with
its names and its defaults, so a run and a comparison cannot come to mean
different things by them.
"""

from collections import Counter
from pathlib import Path
from itertools import permutations, product
import numpy as np
import networkx as nx
from networkx.algorithms.isomorphism import GraphMatcher

from ..molecule import Molecule
from ..optimizationdriver import OptimizationDriver
from ..superimpose import svd_superimpose
from ..errorhandler import assert_msg_critical
from . import util
from .util import Shell, on_master, param, print_section, ic_cell
import math

# The defaults of the manager's matching settings, which reads them from
# here; the methods of this module take every setting they read as an
# argument and carry no default of their own.
DEFAULT_MAX_MAPPINGS = 10000
DEFAULT_METAL_SHELL_BONDS = 2
DEFAULT_RMSD_HEAVY_ATOMS_ONLY = True
DEFAULT_IC_TYPES = {
    'bonds': 'A',
    'angles': 'deg',
    'dihedrals': 'deg',
}


# ----------------------------------------------------------------------
# describing a site
# ----------------------------------------------------------------------


# ----------------------------------------------------------------------
# matching two sites
# ----------------------------------------------------------------------


# ----------------------------------------------------------------------
# measuring one against the other
# ----------------------------------------------------------------------


class SiteMatcher(Shell):
    """
    Measures a structure's active site against templates and decides
    between them: describes both as a coarse graph of metals and residues,
    maps the one onto the other atom by atom, reads the geometric
    agreement off the mapping (compare), holds every template to a set of
    criteria and ranks the passing ones (select_template), and prints the
    comparison and the decision. The criteria, the regions and what the
    ranking reads are the caller's: MetalForceFieldManager holds them as
    settings and class constants.

    Every method runs on the master rank and its result is broadcast. A
    node_match closure handed to networkx's GraphMatcher must not capture
    self: the matcher keeps the callable alive inside reference cycles of
    its own, and an instance caught in one is collected whenever the cycle
    collector gets to it rather than when the call returns.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.
    """

    @on_master
    def describe(self, active_site, edges):
        """
        Adds the two topologies a comparison is made on to an active site.

        The fine topology is the atoms and their bonds, metal-ligand bonds
        included, where the chemistry is ordinary covalent bonding. The
        coarse topology is one node per metal center and one per residue,
        with an edge wherever a residue coordinates a metal and nothing said
        about how many of its atoms do the coordinating: a carboxylate that
        grips a metal with one oxygen and one that grips it with two are the
        same residue on the same metal, and only the coarse level is allowed
        to have an opinion about metal bonds.

        The residues are the connected components left once the metals are
        taken out of the fine topology. That decomposition is untouched by
        how a metal is gripped, which is the whole reason for splitting the
        comparison in two. Each of them is a coarse node carrying its
        Weisfeiler-Lehman key, its formula, its atoms, its heavy atom subgraph
        and the family key of that subgraph, so the two levels are one object
        rather than a graph and a list that have to agree. The family key is hashed
        here, once per residue, rather than by each of the readers that
        compares one -- the isomorphism search asks for it per candidate
        node pair, which made it the most repeated hash in the module.

        The active site is not modified; a new dictionary is returned.

        :param active_site:
            The active site, as the builder extracts it.
        :param edges:
            Its bonds, as index pairs.

        :return:
            The active site with fine_topology, coarse_topology and
            composition added.
        """

        labels = active_site['molecule'].get_labels()
        metal_indices = list(active_site['metal_indices'])

        fine = self._graph(labels, edges, active_site['cap_indices'])

        sidechains = fine.copy()
        sidechains.remove_nodes_from(metal_indices)

        coarse = nx.Graph()

        for metal in metal_indices:
            coarse.add_node(('metal', metal), kind='metal', key=labels[metal])

        for index, component in enumerate(nx.connected_components(sidechains)):
            nodes = sorted(component)
            node = ('residue', index)
            heavy = self._heavy_subgraph(labels, fine, nodes)
            coarse.add_node(node,
                            kind='residue',
                            key=self._fragment_key(fine, nodes),
                            formula=self._formula(labels, nodes),
                            atoms=nodes,
                            heavy=heavy,
                            family=self.family_key(heavy))

            for metal in metal_indices:
                if any(fine.has_edge(metal, atom) for atom in nodes):
                    coarse.add_edge(node, ('metal', metal))

        return {
            **active_site,
            'fine_topology': fine,
            'coarse_topology': coarse,
            'composition': sorted(labels),
        }

    @on_master
    def _graph(self, labels, edges, cap_indices):
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

    @on_master
    def residue_nodes(self, coarse):
        """
        Returns the residue nodes of a coarse topology.

        :param coarse:
            The coarse topology.

        :return:
            The nodes, in the order they were found.
        """

        return [
            node for node, kind in coarse.nodes(data='kind')
            if kind == 'residue'
        ]

    @on_master
    def _formula(self, labels, nodes):
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

    @on_master
    def _fragment_key(self, graph, nodes):
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

    @on_master
    def _heavy_subgraph(self, labels, graph, nodes):
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

    @on_master
    def site_spec(self, described):
        """
        Returns what a site is made of, in a form two sites can be compared
        and reported by.

        The coarse topology answers this on its own: the residues on each
        metal are its neighbours, and a residue is named by its key.

        :param described:
            A described active site or a template.

        :return:
            One entry per metal center and one for the residues themselves.
        """

        coarse = described['coarse_topology']
        labels = described['molecule'].get_labels()
        spec = {}

        for node, kind in coarse.nodes(data='kind'):
            if kind != 'metal':
                continue
            metal = node[1]
            spec[f'{labels[metal]}{metal}'] = sorted(
                coarse.nodes[image]['key'] for image in coarse.neighbors(node))

        spec['residues'] = sorted(coarse.nodes[node]['key']
                                  for node in self.residue_nodes(coarse))
        spec['bridging'] = sorted(coarse.nodes[node]['key']
                                  for node in self.bridging_nodes(coarse))

        return spec

    @on_master
    def bridging_nodes(self, coarse):
        """
        Returns the residue nodes that coordinate more than one metal center.

        One definition, because the spec a site is compared by and the spec it
        is reported by have to agree on which residues bridge.

        :param coarse:
            The coarse topology.

        :return:
            The nodes, in the order they were found.
        """

        return [node for node in self.residue_nodes(coarse) if coarse.degree(node) > 1]

    @on_master
    def family_key(self, heavy):
        """
        Which amino acid a residue is, whatever it is protonated as.

        The key of a coarse node is hashed over the fragment with its
        hydrogens, so an ASP and an ASH are different residues -- which is
        what a comparison wants and what a shoehorning must see past. This
        hashes the heavy atoms alone, so the two come out the same and the
        question it answers is whether the site holds the residue at all.

        :param heavy:
            The heavy atom subgraph of one residue.

        :return:
            The key.
        """

        return nx.weisfeiler_lehman_graph_hash(heavy, node_attr='elem')

    @on_master
    def sidechain_heavy_atoms(self, residue):
        """
        The atoms of a residue the truncation keeps as heavy atoms.

        The cut is at the CA-CB bond and CA becomes a capping hydrogen, so
        what a residue contributes to an active site is its heavy atoms from
        CB outward. Derived on an untruncated topology, so that a residue can
        be identified before anything has been extracted.

        :param residue:
            The residue.

        :return:
            The atoms.
        """

        return [
            atom for atom in residue.atoms()
            if atom.name not in util.BACKBONE_ATOM_NAMES and atom.name != 'CA'
            and atom.element is not None and atom.element.symbol != 'H'
        ]

    @on_master
    def residue_family_key(self, topology, residue):
        """
        The family key of a residue that has not been extracted.

        Built to be the key _family_key gives the same residue once it is a
        fragment of an active site: the same atoms, the same bonds between
        them, and the hydrogens left out of both.

        :param topology:
            The topology the residue belongs to.
        :param residue:
            The residue.

        :return:
            The key.
        """

        atoms = {
            atom.index: atom for atom in self.sidechain_heavy_atoms(residue)
        }

        graph = nx.Graph()

        for index, atom in atoms.items():
            graph.add_node(index, elem=atom.element.symbol)

        for first, second in topology.bonds():
            if first.index in atoms and second.index in atoms:
                graph.add_edge(first.index, second.index)

        return self.family_key(graph)

    @on_master
    def coarse_mappings(self, template, query, match_protonation=True):
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
        :param match_protonation:
            Whether a residue has to be protonated like the one it is matched
            with. Always True for a comparison, since the key of a residue is
            hashed over its hydrogens and an ASP is not an ASH. The
            shoehorning turns it off for the stages that run before the
            protonation has been put right, where the question is only which
            amino acid coordinates which metal.

        :return:
            One mapping of coarse nodes per way the two sites line up.
        """

        # node_match must not capture self: GraphMatcher keeps the callable
        # alive inside reference cycles of its own, and an instance caught
        # in one is collected whenever the cycle collector gets to it rather
        # than when it is finished with. For a manager that means its output
        # stream -- and with it whatever sys.stdout was when it was made --
        # is closed in the middle of somebody else's work. Reading the
        # family key _describe already stored keeps it out of the closure
        # altogether, and off the hot path of the search.
        def node_match(a, b):
            if a['kind'] != b['kind']:
                return False
            if match_protonation or a['kind'] == 'metal':
                return a['key'] == b['key']
            return a['family'] == b['family']

        matcher = GraphMatcher(template['coarse_topology'],
                               query['coarse_topology'],
                               node_match=node_match)

        return list(matcher.isomorphisms_iter())

    @on_master
    def heavy_atom_maps(self, template, query, coarse_mapping, max_mappings,
                        match_h_count=True):
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
        :param match_h_count:
            Whether an atom has to carry the same number of hydrogens as the
            one it is mapped onto. Always True for a comparison: without it a
            CH2 maps onto a CH3 and a protonated oxygen onto a bare one. The
            shoehorning turns it off for the one pass that has to run before
            the protonation has been put right, since a residue that carries
            the wrong hydrogens is exactly what that pass is there to find.

        :return:
            The heavy atom mappings, from template index to query index.
        """

        def node_match(a, b):
            if match_h_count and a['h_count'] != b['h_count']:
                return False
            return a['elem'] == b['elem']

        metals = {}
        per_residue = []

        for node, image in coarse_mapping.items():
            if node[0] == 'metal':
                metals[node[1]] = image[1]
                continue

            first = template['coarse_topology'].nodes[node]['heavy']
            second = query['coarse_topology'].nodes[image]['heavy']

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

            if len(maps) >= max_mappings:
                self.ostream.print_warning(
                    f'Reached the limit of {max_mappings} atom mappings; '
                    'the best of the ones built is used, which need not be '
                    'the best there is')
                self.ostream.flush()
                break

        return maps

    @on_master
    def best_heavy_map(self, template, maps, coordinates):
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

        # every mapping in maps is keyed on the same template indices, so
        # the order and the slice it takes are loop invariants
        order = sorted(maps[0]) if maps else []
        reference_block = reference[order]

        for atom_map in maps:
            moved = coordinates[[atom_map[index] for index in order]]

            rmsd, rot, trans = svd_superimpose(moved, reference_block)

            if best is None or rmsd < best[0]:
                best = (rmsd, atom_map, rot, trans)

        return best[1], best[2], best[3]

    @on_master
    def complete_hydrogens(self, template, query, heavy_map, coordinates, rot, trans):
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

        template_labels = template['molecule'].get_labels()
        query_labels = query['molecule'].get_labels()

        def hydrogens(labels, graph, node):
            return [
                other for other in graph.neighbors(node) if labels[other] == 'H'
            ]

        for node, image in heavy_map.items():
            first = hydrogens(template_labels, template['fine_topology'], node)
            second = hydrogens(query_labels, query['fine_topology'], image)

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

    @on_master
    def rmsd_indices(self, template, region, heavy_only, metal_shell_bonds):
        """
        Returns the template indices every RMSD is measured over, which is the
        region less the hydrogens when they are being left out.

        :param template:
            The template.
        :param region:
            The region to take.
        :param heavy_only:
            Whether to leave the hydrogens out.
        :param metal_shell_bonds:
            How many bonds out from a metal the metal_shell region reaches.

        :return:
            The indices, in order.
        """

        indices = self.region_indices(template, region, metal_shell_bonds)

        if not heavy_only:
            return indices

        labels = template['molecule'].get_labels()

        return [index for index in indices if labels[index] != 'H']

    @on_master
    def region_indices(self, template, region, metal_shell_bonds):
        """
        Returns the template indices of the region an RMSD is measured over.

        :param template:
            The template.
        :param region:
            The region to take.
        :param metal_shell_bonds:
            How many bonds out from a metal the metal_shell region reaches.

        :return:
            The indices, in order.
        """

        if region == 'active_site':
            return list(range(template['molecule'].number_of_atoms()))

        graph = template['fine_topology']

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
                    graph, metal, cutoff=metal_shell_bonds))

        assert_msg_critical(
            len(shell) > len(template['metal_indices']),
            'MetalForceFieldManager: the metal centers of template '
            f'{template["name"]} have nothing bonded to them within '
            f'{metal_shell_bonds} bond(s)')

        return sorted(shell)

    @on_master
    def measure_region(self, template, mapping, coordinates, region,
                       heavy_only, metal_shell_bonds):
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
            Whether to leave the hydrogens out.
        :param metal_shell_bonds:
            How many bonds out from a metal the metal_shell region reaches.

        :return:
            The atom count, the cartesian RMSDs over the whole region and
            over its heavy atoms, and the internal coordinate deviations.
        """

        reference = template['molecule'].get_coordinates_in_angstrom()
        indices = self.region_indices(template, region, metal_shell_bonds)
        labels = template['molecule'].get_labels()
        heavy = [index for index in indices if labels[index] != 'H']

        order = [mapping[index] for index in range(len(reference))]
        moved = coordinates[order]

        rmsd, _, _ = svd_superimpose(moved[indices], reference[indices])
        heavy_rmsd, _, _ = svd_superimpose(moved[heavy], reference[heavy])

        ic_rmsd = self.measure_ic_rmsd(template, mapping, coordinates, region,
                                       heavy_only, metal_shell_bonds)

        return {
            # _rmsd_indices is the region less the hydrogens when they are
            # being left out, which is exactly the two lists already in hand
            'atoms': len(heavy) if heavy_only else len(indices),
            'rmsd': rmsd,
            'rmsd_heavy': heavy_rmsd,
            'ic_rmsd': ic_rmsd,
        }

    @on_master
    def measure_ic_rmsd(self, template, mapping, coordinates, region,
                        heavy_only, metal_shell_bonds):
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
            The region to measure over.
        :param heavy_only:
            Whether to leave the hydrogens out.
        :param metal_shell_bonds:
            How many bonds out from a metal the metal_shell region reaches.

        :return:
            The deviations, as get_ic_rmsd reports them, or None when they
            could not be measured.
        """

        indices = self.rmsd_indices(template, region, heavy_only,
                                    metal_shell_bonds)
        elements = template['molecule'].get_labels()
        labels = [elements[index] for index in indices]

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

    @on_master
    def ic_violation(self, ic_rmsd, thresholds, ic_types):
        """
        Checks internal coordinate deviations against a set of thresholds.

        A type left at None, or a threshold left out of one, is not checked.
        A type that is asked for but was not measured is a violation rather
        than a pass: a criterion that could not be evaluated is not one the
        template met.

        :param ic_rmsd:
            The deviations, as get_ic_rmsd reports them.
        :param thresholds:
            The thresholds to check against, as {type: {'rms': , 'max': }}.
        :param ic_types:
            The internal coordinate types to check, with their units.

        :return:
            A description of the first threshold that was exceeded, or None
            when every one of them holds.
        """

        if ic_rmsd is None:
            return 'no internal coords'

        for name, unit in ic_types.items():
            limits = thresholds.get(name)
            if not limits:
                continue

            found = ic_rmsd.get(name)
            if found is None:
                return f'no {name} measured'

            for measure in ('rms', 'max'):
                limit = limits.get(measure)
                if limit is None:
                    continue
                if found[measure] > limit:
                    return (f'{name} {measure} {found[measure]:.2f} > '
                            f'{limit:.2f} {unit}')

        return None

    @on_master
    def metal_bond_summary(self, template, query, mapping, coordinates):
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
        bonds, _ = self.metal_keys(template)

        deviation = 0.0
        shared = 0
        template_only = 0
        mapped = set()

        for i, j in bonds:
            first, second = mapping[i], mapping[j]
            mapped.add(frozenset((first, second)))

            if not query['fine_topology'].has_edge(first, second):
                template_only += 1
                continue

            shared += 1
            expected = np.linalg.norm(reference[i] - reference[j])
            found = np.linalg.norm(coordinates[first] - coordinates[second])
            deviation = max(deviation, abs(found - expected))

        metals = set(query['metal_indices'])
        query_only = sum(
            1 for first, second in query['fine_topology'].edges()
            if ({first, second} & metals) and frozenset((first,
                                                         second)) not in mapped)

        return {
            'deviation': deviation,
            'shared': shared,
            'template_only': template_only,
            'query_only': query_only,
        }

    @on_master
    def metal_keys(self, template):
        """
        Returns the metal bond and angle keys of a template.

        :param template:
            The template.

        :return:
            The tuple of the bond key list and the angle key list.
        """

        # get_metal_keys reads nothing of the active site but the metal
        # indices, and a template knows those from its elements
        return util.get_metal_keys(template['forcefield'],
                                   {'metal_indices': template['metal_indices']})

    @on_master
    def print_spec(self, title, described, bridging):
        """
        Prints what a site is made of: which residues it holds and which of them
        coordinate which metal.

        The residues are named by their formula, which is for reading; two sites
        are compared on the keys behind them.

        :param title:
            What the block is describing.
        :param described:
            A described active site or a template.
        :param bridging:
            Its residue nodes that coordinate more than one metal center, from
            matching.bridging_nodes -- the same answer the spec two sites are
            compared by is built from.
        """

        labels = described['molecule'].get_labels()
        coarse = described['coarse_topology']

        def named(nodes):
            return ', '.join(
                sorted(f'{coarse.nodes[node]["formula"]}/'
                       f'{coarse.nodes[node]["key"][:6]}' for node in nodes))

        self.ostream.print_blank()
        self.ostream.print_info(f'Site spec, {title}:')

        for metal in described['metal_indices']:
            node = ('metal', metal)
            self.ostream.print_info(
                f'  {labels[metal]}{metal}: {named(coarse.neighbors(node))}')

        if bridging:
            self.ostream.print_info(f'  bridging: {named(bridging)}')

        self.ostream.flush()

    @on_master
    def print_template_comparison(self, name, entry, spec=None):
        """
        Prints the numbers and the verdicts of one template.

        :param name:
            The name of the template.
        :param entry:
            What compare_active_site measured for it.
        :param spec:
            The template and its residue nodes, for the case where the site
            coordinates a different set of residues and the two specs are worth
            reading side by side. None when nothing was measured for another
            reason.
        """

        self.ostream.print_blank()

        if entry['status'] == 'composition':
            self.ostream.print_info(
                f'{name}: holds different atoms, so nothing was measured.')
            self.ostream.flush()
            return

        if entry['status'] == 'spec':
            self.ostream.print_info(f'{name}: coordinates a different set of residues, '
                                    'so nothing was measured.')
            if spec is not None:
                self.print_spec(f'{name} holds', *spec)
            self.ostream.flush()
            return

        bonds = entry['metal_bonds']
        summary = (f'{bonds["shared"]} metal bond(s) shared, within '
                   f'{bonds["deviation"]:.3f} A')
        if bonds['template_only']:
            summary += f', {bonds["template_only"]} only in the template'
        if bonds['query_only']:
            summary += f', {bonds["query_only"]} only in the structure'

        self.ostream.print_info(f'{name}: {entry["n_mappings"]} atom mapping(s) from '
                                f'{entry["n_coarse_mappings"]} coarse mapping(s), '
                                f'{summary}')
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
                           ic_cell(found['ic_rmsd'], 'bonds'),
                           ic_cell(found['ic_rmsd'], 'angles'),
                           ic_cell(found['ic_rmsd'], 'dihedrals')))

        self.ostream.print_blank()
        self.ostream.flush()

    @on_master
    def print_comparison_summary(self, results, ranked_on, scores):
        """
        Ranks the templates on what a selection is decided by.

        :param results:
            The last comparison, from compare_active_site.
        :param ranked_on:
            The (region, ic type, measure) a selection is ranked on.
        :param scores:
            That measure per template, computed by the caller.
        """

        region, ic_type, measure = ranked_on
        heavy = not results['include_hydrogens']

        def rmsd(entry):
            found = entry['regions'].get(region)
            if found is None:
                return None
            return found['rmsd_heavy'] if heavy else found['rmsd']

        order = sorted(results['templates'].items(),
                       key=lambda item: scores[item[0]])

        self.ostream.print_blank()
        print_section(f'Ranked on the {region} {ic_type} {measure}', self.ostream)

        row = '{:>24} | {:>14} | {:>11} | {:>9} | {:>10}'
        self.ostream.print_header(
            row.format('template', 'status', f'{ic_type} {measure}', 'RMSD',
                       'metal bond'))
        self.ostream.print_header(78 * '-')

        for name, entry in order:
            found = rmsd(entry)
            if found is None:
                self.ostream.print_header(
                    row.format(name[:24], entry['status'], '', '', ''))
                continue
            self.ostream.print_header(
                row.format(name[:24], entry['status'], f'{scores[name]:.3f}',
                           f'{found:.3f}',
                           f'{entry["metal_bonds"]["deviation"]:.3f}'))

        self.ostream.print_blank()
        self.ostream.flush()

    @on_master
    def print_comparison(self, results, specs, ranked_on, scores):
        """
        Prints everything compare_active_site measured.

        One table of numbers and one of verdicts per template that could be
        measured, and a closing summary ranking the templates by the region that
        is configured, so the closest one is visible without reading every table.

        :param results:
            The last comparison, from compare_active_site.
        :param specs:
            The (described, residue nodes) pair per site whose spec is printed:
            the structure under the key None, and a template under its name.
        :param ranked_on:
            The (region, ic type, measure) a selection is ranked on.
        :param scores:
            That measure per template, computed by the caller.
        """

        active_site = results['active_site']
        labels = active_site['molecule'].get_labels()
        metals = ', '.join(labels[index] for index in active_site['metal_indices'])

        self.ostream.print_blank()
        print_section('Comparison against every template', self.ostream)
        self.ostream.print_header(param('source', Path(results['source']).name))
        self.ostream.print_header(
            param('active site atoms', active_site['molecule'].number_of_atoms()))
        self.ostream.print_header(param('metal centers', metals))
        self.ostream.print_header(param('geometry', results['geometry']))
        self.ostream.print_header(
            param('measured over',
                  'all atoms' if results['include_hydrogens'] else 'heavy atoms'))
        self.ostream.print_header(param('templates', len(results['templates'])))
        self.ostream.print_blank()
        self.ostream.print_info(f'Residues: {", ".join(active_site["residues"])}')

        self.print_spec('the structure holds', *specs[None])

        for name, entry in results['templates'].items():
            self.print_template_comparison(name,
                                           entry,
                                           spec=specs.get(name))

        self.print_comparison_summary(results, ranked_on, scores)

    # ------------------------------------------------------------------
    # the comparison and the decision
    # ------------------------------------------------------------------

    @on_master
    def compare(self, templates, described, molecule, regions,
                include_hydrogens, max_mappings, metal_shell_bonds):
        """
        Measures a described site against every template, without an opinion.

        Per template: the same atoms or not ('composition'), the same coarse
        coordination or not ('spec'), and when both hold the best atom mapping
        and the geometric agreement of every region under it.

        :param templates:
            The templates, by name.
        :param described:
            The site, as describe returns it.
        :param molecule:
            The geometry every number is measured on.
        :param regions:
            The regions to measure.
        :param include_hydrogens:
            Whether the hydrogens take part in the measurements.
        :param max_mappings:
            The limit on how many atom mappings are built per template.
        :param metal_shell_bonds:
            How many bonds out from a metal the metal_shell region reaches.

        :return:
            One entry per template: 'status', 'mapping', 'n_coarse_mappings',
            'n_mappings', 'metal_bonds' and 'regions'.
        """

        coordinates = molecule.get_coordinates_in_angstrom()
        heavy_only = not include_hydrogens

        findings = {}

        for name, template in templates.items():
            entry = {
                'status': 'measured',
                'mapping': None,
                'n_coarse_mappings': 0,
                'n_mappings': 0,
                'metal_bonds': None,
                'regions': {},
            }

            if template['composition'] != described['composition']:
                # not the same atoms, so there is nothing to map onto
                entry['status'] = 'composition'
                findings[name] = entry
                continue

            # which residue coordinates which metal, with nothing said about
            # how many atoms of it do the coordinating
            coarse = self.coarse_mappings(template, described)
            if not coarse:
                entry['status'] = 'spec'
                findings[name] = entry
                continue

            maps = []
            for coarse_mapping in coarse:
                maps.extend(
                    self.heavy_atom_maps(template, described, coarse_mapping,
                                         max_mappings))
            if not maps:
                entry['status'] = 'spec'
                findings[name] = entry
                continue

            entry['n_coarse_mappings'] = len(coarse)
            entry['n_mappings'] = len(maps)

            heavy_map, rot, trans = self.best_heavy_map(template, maps,
                                                        coordinates)
            mapping = self.complete_hydrogens(template, described, heavy_map,
                                              coordinates, rot, trans)
            entry['mapping'] = mapping
            entry['metal_bonds'] = self.metal_bond_summary(
                template, described, mapping, coordinates)

            for region in regions:
                entry['regions'][region] = self.measure_region(
                    template, mapping, coordinates, region, heavy_only,
                    metal_shell_bonds)

            findings[name] = entry

        return findings

    @on_master
    def select_template(self, comparison,
                        criteria,
                        criteria_name,
                        regions,
                        ic_types,
                        ranked_on,
                        template):
        """
        Picks the template a force field should be built from, and says why.

        Reads the last comparison made by compare_active_site. A template is
        usable only if it maps onto every atom of the site. If none is
        named, every template that maps completely is held to
        selection_criteria and the passing ones are ranked, the best one
        winning. If one is named, its verdict is reported but does not
        decide whether it is picked here -- naming a template is a statement
        that it is the right one, and build_ff_from_template is what turns an
        outside-the-criteria verdict into a refusal.

        :param template:
            The name of the template to use, or None to choose one.

        :return:
            The decision, with the chosen name under 'name', None when
            nothing passed, and what was made of every template under
            'verdicts'.
        """

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
            verdict = self.selection_verdict(comparison, entry, criteria,
                                             regions, ic_types)
            decision['verdicts'][name] = verdict
            decision['scores'][name] = self.selection_score(entry, ranked_on)

        if template is not None:
            assert_msg_critical(
                template in comparison['templates'],
                'SiteMatcher.select_template: no template named '
                f'{template} was compared. Loaded: '
                f'{sorted(comparison["templates"])}')

            entry = comparison['templates'][template]
            verdict = decision['verdicts'][template]

            assert_msg_critical(
                self.maps_every_atom(comparison, entry),
                f'SiteMatcher.select_template: template {template} '
                f'does not map onto every atom of the site ({verdict}), so '
                'its parameters cannot be transferred. Build this site with '
                'MetalSiteForceFieldBuilder.')

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

        return decision

    @on_master
    def prefer_template(self, comparison, decision, name):
        """
        Takes the template the site was shoehorned into, over whatever the
        criteria came to.

        A shoehorning is the same statement naming a template is -- this
        site is to be built the way that one is -- so an unnamed call after
        one is not really unnamed, and what the criteria make of the field
        does not get to overrule it. Both ways they can differ are wrong on
        their own terms.

        With nothing within them, the criteria are the least able to judge:
        they measure a geometry whose coordination sphere is still open, and
        what would close it is the very force field being asked for, so the
        site cannot look like the template until after the transfer it is
        being refused.

        With something within them, the something is rarely alone. A
        shoehorned site passes against the template it was walked onto and
        against every sibling of that template's family too, and the ranking
        then separates them on an active-site bond rms that differs in
        thousandths of an Angstrom -- so the site gets walked onto one
        template and built from another, which is a geometry from one and
        parameters from the other and a description of neither. That is not
        a tie to be broken better: the shoehorning already said which one it
        is.

        What is not waived is the mapping. Every atom of the site has to
        land somewhere in the template or there is nothing to transfer onto
        it, and _build_ff_from_template fails outright on the first key it
        cannot map, so a template that maps incompletely is left refused and
        the criteria keep the decision.

        :param decision:
            The decision _select_template came to.

        :return:
            A decision naming the shoehorned template, or the one given when
            there is no shoehorning to take, when it is already what the
            criteria took, or when it does not map the site completely.
        """

        if name is None or name not in comparison['templates']:
            return decision

        if name == decision['name']:
            return decision

        entry = comparison['templates'][name]

        if not self.maps_every_atom(comparison, entry):
            self.ostream.print_warning(
                f'The site was shoehorned into {name}, but it does not map '
                'onto every atom of the site, so nothing can be transferred '
                'from it.')
            self.ostream.flush()
            return decision

        verdict = decision['verdicts'][name] or 'within the criteria'
        ranked = decision['name']

        decision = dict(decision)
        decision['name'] = name
        decision['entry'] = entry
        decision['forced'] = True
        decision['score'] = decision['scores'][name]

        if ranked is None:
            self.ostream.print_warning(
                f'No template is within the {decision["criteria_name"]} '
                f'criteria, but the site was shoehorned into {name}, which '
                'is the same statement as naming it. Building from it '
                f'anyway: {verdict}.')
            self.ostream.print_info(
                'The criteria measure a geometry whose coordination sphere '
                'the transferred parameters have not closed yet. Relax the '
                'site on the force field this returns '
                '(mm_optimize_active_site) and compare again to see whether '
                'it then matches on its own.')
        else:
            self.ostream.print_warning(
                f'The criteria ranked {ranked} first, but the site was '
                f'shoehorned into {name}, which is the same statement as '
                f'naming it. Building from {name} instead: {verdict}.')
            self.ostream.print_info(
                f'The site was walked onto {name}, so its parameters are the '
                f'ones that describe it. Name {ranked} in the call to build '
                'from that one instead.')

        self.ostream.flush()

        return decision

    @on_master
    def maps_every_atom(self, comparison, entry):
        """
        Says whether a template covers every atom of the site.

        A template that holds other atoms, or coordinates them differently,
        never gets as far as a mapping; this also catches a mapping that came
        back incomplete, which would leave part of a site unparameterized.

        :param comparison:
            The last comparison, from compare_active_site.
        :param entry:
            What it measured for the template.

        :return:
            True when every atom of the site is mapped onto exactly once.
        """

        if entry['status'] != 'measured' or entry['mapping'] is None:
            return False

        atoms = comparison['active_site']['molecule'].number_of_atoms()
        mapping = entry['mapping']

        return (len(mapping) == atoms
                and sorted(mapping.values()) == list(range(atoms)))

    @on_master
    def selection_verdict(self, comparison, entry, criteria, regions,
                          ic_types):
        """
        Holds one template to the criteria, region by region.

        :param comparison:
            The last comparison, from compare_active_site.
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

        if not self.maps_every_atom(comparison, entry):
            return 'incomplete mapping'

        for region in regions:
            thresholds = criteria.get(region)
            if not thresholds:
                continue

            found = entry['regions'].get(region)
            if found is None:
                return f'{region} not measured'

            # a criterion that could not be evaluated is not one that was
            # passed, so the region is held to strictly here
            violation = self.ic_violation(found['ic_rmsd'], thresholds,
                                          ic_types)
            if violation is not None:
                return f'{region} {violation}'

        return None

    @on_master
    def selection_score(self, entry, ranked_on):
        """
        Returns what several templates that all pass are ranked on.

        :param entry:
            What compare_active_site measured for the template.

        :return:
            The measure named by SELECTION_RANKED_ON, or infinity where it was
            not measured.
        """

        region, ic_type, measure = ranked_on

        found = entry['regions'].get(region)
        if found is None or found['ic_rmsd'] is None:
            return math.inf

        found = found['ic_rmsd'].get(ic_type)
        if found is None:
            return math.inf

        return found[measure]

    # ------------------------------------------------------------------
    # the decision, printed
    # ------------------------------------------------------------------

    @on_master
    def print_selection(self, comparison, decision, regions, ic_types,
                        ranked_on):
        """
        Prints how every template stands against the criteria, and which one was
        taken.

        The whole field is printed rather than the winner alone: whether the
        others are near misses or a long way off is what says how much the chosen
        one is worth.

        :param decision:
            The decision, as _select_template makes it.
        """

        regions = [
            region for region in regions
            if decision['criteria'].get(region)
        ]

        self.ostream.print_blank()
        print_section(
            f'Choosing a template on the {decision["criteria_name"]} criteria',
            self.ostream)
        self.ostream.print_blank()

        for region in regions:
            thresholds = decision['criteria'][region]
            # only the measures the set actually holds, since either of them may
            # be left out of one
            measures = {
                name:
                ' / '.join(f'{measure} {limit:.2f}'
                           for measure, limit in given.items() if limit is not None)
                for name, given in thresholds.items() if given
            }
            limits = '; '.join(f'{name} {shown} {ic_types[name]}'
                               for name, shown in measures.items())
            self.ostream.print_header(param(region, limits, value_width=44))

        self.ostream.print_blank()

        # one column per region the criteria name, so a custom set of them prints
        # as readably as the two that come with the class
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
                             None else ic_cell(found['ic_rmsd'], 'bonds'))

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
            ranked = ' '.join(ranked_on)
            self.ostream.print_info(
                f'{len(decision["candidates"])} of '
                f'{len(comparison["templates"])} template(s) are within the '
                f'criteria. Taking {decision["name"]}, whose {ranked} of '
                f'{decision["score"]:.3f} is the lowest of them.')

        self.ostream.print_blank()
        self.ostream.flush()

    @on_master
    def print_no_selection(self, decision):
        """
        Says which template came closest when none of them was good enough.

        :param decision:
            The decision, as _select_template makes it.
        """
        closest = min(decision['scores'],
                      key=lambda name: decision['scores'][name],
                      default=None)

        if closest is not None and math.isfinite(decision['scores'][closest]):
            self.ostream.print_info(
                f'No template is within the {decision["criteria_name"]} '
                f'criteria. The closest is {closest}: '
                f'{decision["verdicts"][closest]}.')
        else:
            self.ostream.print_info(
                'No template describes this site: none of them maps onto all '
                'of its atoms.')

        self.ostream.print_info(
            "Set selection_criteria to 'loose' to widen what counts as a "
            'match, or build this site with MetalSiteForceFieldBuilder.')
        self.ostream.flush()

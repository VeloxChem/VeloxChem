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
from itertools import permutations, product
import numpy as np
import networkx as nx
from networkx.algorithms.isomorphism import GraphMatcher

from ..molecule import Molecule
from ..optimizationdriver import OptimizationDriver
from ..superimpose import svd_superimpose
from ..errorhandler import assert_msg_critical
from . import core
from .printing import stream

# The manager's own defaults, repeated here so that a function of this
# module stands on its own; the manager always passes its settings in.
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


def describe(active_site, edges):
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

    fine = _graph(labels, edges, active_site['cap_indices'])

    sidechains = fine.copy()
    sidechains.remove_nodes_from(metal_indices)

    coarse = nx.Graph()

    for metal in metal_indices:
        coarse.add_node(('metal', metal), kind='metal', key=labels[metal])

    for index, component in enumerate(nx.connected_components(sidechains)):
        nodes = sorted(component)
        node = ('residue', index)
        heavy = _heavy_subgraph(labels, fine, nodes)
        coarse.add_node(node,
                        kind='residue',
                        key=_fragment_key(fine, nodes),
                        formula=_formula(labels, nodes),
                        atoms=nodes,
                        heavy=heavy,
                        family=family_key(heavy))

        for metal in metal_indices:
            if any(fine.has_edge(metal, atom) for atom in nodes):
                coarse.add_edge(node, ('metal', metal))

    return {
        **active_site,
        'fine_topology': fine,
        'coarse_topology': coarse,
        'composition': sorted(labels),
    }


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


def residue_nodes(coarse):
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


def site_spec(described):
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

    residues = residue_nodes(coarse)

    spec['residues'] = sorted(coarse.nodes[node]['key']
                              for node in residues)
    spec['bridging'] = sorted(coarse.nodes[node]['key'] for node in residues
                              if coarse.degree(node) > 1)

    return spec


def family_key(heavy):
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


def sidechain_heavy_atoms(residue):
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
        if atom.name not in core.BACKBONE_ATOM_NAMES and atom.name != 'CA'
        and atom.element is not None and atom.element.symbol != 'H'
    ]


def residue_family_key(topology, residue):
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
        atom.index: atom for atom in sidechain_heavy_atoms(residue)
    }

    graph = nx.Graph()

    for index, atom in atoms.items():
        graph.add_node(index, elem=atom.element.symbol)

    for first, second in topology.bonds():
        if first.index in atoms and second.index in atoms:
            graph.add_edge(first.index, second.index)

    return family_key(graph)

# ----------------------------------------------------------------------
# matching two sites
# ----------------------------------------------------------------------


def coarse_mappings(template, query, match_protonation=True):
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


def heavy_atom_maps(template,
                    query,
                    coarse_mapping,
                    match_h_count=True,
                    max_mappings=DEFAULT_MAX_MAPPINGS,
                    ostream=None):
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

    ostream = stream(ostream)

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
            ostream.print_warning(
                f'Reached the limit of {max_mappings} atom mappings; '
                'the best of the ones built is used, which need not be '
                'the best there is')
            ostream.flush()
            break

    return maps


def best_heavy_map(template, maps, coordinates):
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


def complete_hydrogens(template, query, heavy_map, coordinates, rot, trans):
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

# ----------------------------------------------------------------------
# measuring one against the other
# ----------------------------------------------------------------------


def rmsd_indices(template,
                 region,
                 heavy_only=None,
                 rmsd_heavy_atoms_only=DEFAULT_RMSD_HEAVY_ATOMS_ONLY,
                 metal_shell_bonds=DEFAULT_METAL_SHELL_BONDS):
    """
    Returns the template indices every RMSD is measured over, which is the
    region less the hydrogens when they are being left out.

    :param template:
        The template.
    :param region:
        The region to take.
    :param heavy_only:
        Whether to leave the hydrogens out, or None for
        rmsd_heavy_atoms_only.

    :return:
        The indices, in order.
    """

    if heavy_only is None:
        heavy_only = rmsd_heavy_atoms_only

    indices = region_indices(template,
                             region,
                             metal_shell_bonds=metal_shell_bonds)

    if not heavy_only:
        return indices

    labels = template['molecule'].get_labels()

    return [index for index in indices if labels[index] != 'H']


def region_indices(template, region, metal_shell_bonds=DEFAULT_METAL_SHELL_BONDS):
    """
    Returns the template indices of the region an RMSD is measured over.

    :param template:
        The template.
    :param region:
        The region to take.

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


def measure_region(template,
                   mapping,
                   coordinates,
                   region,
                   heavy_only=None,
                   rmsd_heavy_atoms_only=DEFAULT_RMSD_HEAVY_ATOMS_ONLY,
                   metal_shell_bonds=DEFAULT_METAL_SHELL_BONDS,
                   ostream=None):
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
        over its heavy atoms, and the internal coordinate deviations.
    """

    if heavy_only is None:
        heavy_only = rmsd_heavy_atoms_only

    reference = template['molecule'].get_coordinates_in_angstrom()
    indices = region_indices(template,
                             region,
                             metal_shell_bonds=metal_shell_bonds)
    labels = template['molecule'].get_labels()
    heavy = [index for index in indices if labels[index] != 'H']

    order = [mapping[index] for index in range(len(reference))]
    moved = coordinates[order]

    rmsd, _, _ = svd_superimpose(moved[indices], reference[indices])
    heavy_rmsd, _, _ = svd_superimpose(moved[heavy], reference[heavy])

    ic_rmsd = measure_ic_rmsd(template,
                              mapping,
                              coordinates,
                              region,
                              heavy_only,
                              rmsd_heavy_atoms_only=rmsd_heavy_atoms_only,
                              metal_shell_bonds=metal_shell_bonds,
                              ostream=ostream)

    return {
        # _rmsd_indices is the region less the hydrogens when they are
        # being left out, which is exactly the two lists already in hand
        'atoms': len(heavy) if heavy_only else len(indices),
        'rmsd': rmsd,
        'rmsd_heavy': heavy_rmsd,
        'ic_rmsd': ic_rmsd,
    }


def measure_ic_rmsd(template,
                    mapping,
                    coordinates,
                    region,
                    heavy_only=None,
                    rmsd_heavy_atoms_only=DEFAULT_RMSD_HEAVY_ATOMS_ONLY,
                    metal_shell_bonds=DEFAULT_METAL_SHELL_BONDS,
                    ostream=None):
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
        Whether to leave the hydrogens out, or None for
        rmsd_heavy_atoms_only.

    :return:
        The deviations, as get_ic_rmsd reports them, or None when they
        could not be measured.
    """

    ostream = stream(ostream)

    indices = rmsd_indices(template,
                           region,
                           heavy_only,
                           rmsd_heavy_atoms_only=rmsd_heavy_atoms_only,
                           metal_shell_bonds=metal_shell_bonds)
    elements = template['molecule'].get_labels()
    labels = [elements[index] for index in indices]

    reference = template['molecule'].get_coordinates_in_angstrom()
    order = [mapping[index] for index in indices]

    mapped = Molecule(labels, coordinates[order], 'angstrom')
    wanted = Molecule(labels, reference[indices], 'angstrom')

    ic_rmsd = OptimizationDriver.get_ic_rmsd(mapped, wanted)

    if isinstance(ic_rmsd, str):
        # get_ic_rmsd reports its refusals as a string rather than raising
        ostream.print_warning(
            f'Could not measure the internal coordinates against template '
            f'{template["name"]}: {ic_rmsd}')
        ostream.flush()
        return None

    return ic_rmsd


def ic_violation(ic_rmsd, thresholds, ic_types=DEFAULT_IC_TYPES):
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


def metal_bond_summary(template, query, mapping, coordinates):
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
    bonds, _ = metal_keys(template)

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


def metal_keys(template):
    """
    Returns the metal bond and angle keys of a template.

    :param template:
        The template.

    :return:
        The tuple of the bond key list and the angle key list.
    """

    # get_metal_keys reads nothing of the active site but the metal
    # indices, and a template knows those from its elements
    return core.get_metal_keys(template['forcefield'],
                               {'metal_indices': template['metal_indices']})

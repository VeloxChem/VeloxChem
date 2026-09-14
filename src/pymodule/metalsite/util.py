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
What every phase of the metal site force field code shares.

The constants, the Shell base class the stateless phase classes derive
from, the on_master / collective decorators that make their methods
MPI-safe, the formatting helpers, and the plain helpers more than one
phase reads. A plain helper here is deterministic and costs nothing, so it
is not decorated: broadcasting a string is not worth a collective.
Anything that reads a file, draws a random number or takes more than a few
milliseconds belongs on a shell, decorated.
"""

from mpi4py import MPI
from pathlib import Path
import functools
import numpy as np
import sys

from ..veloxchemlib import mpi_master
from ..outputstream import OutputStream
from ..mmforcefieldgenerator import MMForceFieldGenerator
from ..errorhandler import assert_msg_critical

try:
    import openmm.app as mmapp
except ImportError:
    pass

# Elements recognized when scanning a structure for metal centers. A
# center that is found but not supported is reported rather than ignored.
METAL_ELEMENTS = ('Zn', 'Fe', 'Cu', 'Mg', 'Mn', 'Co', 'Ni', 'Ca', 'Cd')

# Elements the builder is validated for. The literature distances, the
# formal charges and the coordination rules have only been checked against
# zinc sites, so anything else is rejected instead of silently treated as
# if it behaved the same way.
SUPPORTED_METAL_ELEMENTS = ('Zn', )

# Formal charges assumed for bare metal ions. Only used for the active site
# charge bookkeeping, which is checked rather than trusted.
METAL_FORMAL_CHARGES = {
    'Zn': 2,
    'Mg': 2,
    'Ca': 2,
    'Mn': 2,
    'Ni': 2,
    'Co': 2,
    'Cd': 2,
    'Fe': 2,
    'Cu': 2,
}

# Elements that can donate a lone pair to a metal center.
DONOR_ELEMENTS = ('N', 'O', 'S', 'Se')

# How much closer one carboxylate oxygen has to be than the other, in
# Angstrom, before the pair stops counting as a chelating bidentate and
# becomes a monodentate contact of the near oxygen alone.
BIDENTATE_ASYMMETRY = 0.75

# Distance, in Angstrom, within which a donor atom is taken to be bonded to
# a metal center. It is generous on purpose: the scan reads an unrelaxed
# structure, where a stretched bridging contact is still a bond.
METAL_BOND_CUTOFF = 3.0

# How much further than the bonding cutoff, in Angstrom, a contact is
# reported without being made a bond, so that a near miss is visible in the
# coordination table. It is a margin rather than a distance of its own
# because the two are not independent: a scan that stops before the bonding
# cutoff would drop contacts that are bonds, so raising metal_bond_cutoff
# past a fixed report_cutoff used to have no effect at all.
REPORT_CUTOFF_MARGIN = 0.5

# Distance, in Angstrom, out to which a contact is reported without being
# made a bond. Only the default pairing of the two; what a scan uses is
# resolved from the bonding cutoff it is given, in _collect_ligands.
REPORT_CUTOFF = METAL_BOND_CUTOFF + REPORT_CUTOFF_MARGIN

# Distance, in Angstrom, within which a donor atom is given Hessian blocks
# with a metal center whether or not it is bonded to one. The perception is
# a distance cutoff read on a single geometry, and the contact that falls
# just outside it is exactly the one a user adds by hand afterwards. Filling
# a block that was never computed costs the whole Hessian again, so what the
# partial Hessian covers is deliberately more forgiving than the bonding it
# is fitted to. It widens what is computed, never what is fitted.
PARTIAL_HESSIAN_CUTOFF = 3.5

# Length, in Angstrom, of the C-H bond of the hydrogen that caps a
# sidechain where the CA-CB bond was cut.
CAP_BOND_LENGTH = 1.09

# How much further out than its residue's shortest metal bond, in Angstrom,
# the weak arm of a bridging residue has to sit before the fit drops it.
WEAK_BRIDGE_TOLERANCE = 0.25

# Where _seed_metal_terms took a term's equilibrium from, written into the
# comment of every term it touches so that the force field says it without
# being asked again. printing.print_mm_optimization reads them back for the
# crude-pass table, so the two sides cannot drift on the wording.
SEEDED_FROM_REQUEST = 'asked for with the bond'
SEEDED_FROM_TABLE = 'given equilibrium'
SEEDED_FROM_GEOMETRY = 'measured on the input geometry'

# What the crude-pass table calls each of them.
SEEDED_EQUILIBRIUM_LABELS = {
    SEEDED_FROM_REQUEST: 'requested',
    SEEDED_FROM_TABLE: 'given',
    SEEDED_FROM_GEOMETRY: 'measured',
}

# Residues whose sidechain can bridge two metal centers. A carboxylate has
# two donor oxygens and a thiolate sulfur has several lone pairs, so both
# can reach two metals at once. An imidazole nitrogen has a single lone
# pair in the ring plane, so a histidine bound to two metals is a geometric
# artifact rather than a bridge, however short the second contact looks.
BRIDGING_RESIDUES = ('ASP', 'ASH', 'GLU', 'GLH', 'CYS', 'CYM', 'CYX')

# Residues whose sidechain ends in a carboxylate, whose two oxygens are
# interchangeable for as long as neither of them is bound.
CARBOXYLATE_RESIDUES = ('ASP', 'ASH', 'GLU', 'GLH')

# Atoms dropped when truncating a sidechain at the CA-CB bond. CA itself
# is replaced by a capping hydrogen rather than dropped.
BACKBONE_ATOM_NAMES = ('N', 'C', 'O', 'OXT', 'H', 'H2', 'H3', 'HA', 'HA2',
                       'HA3', 'HXT')

# Residues the CA-CB truncation cannot cut, and why. A residue that only
# coordinates a metal is always one of ASP/GLU/CYS/HIS and never lands
# here, but include_residue takes any residue at all, so the rule the
# truncation relies on is written down where it is enforced.
UNTRUNCATABLE_RESIDUES = {
    'GLY':
    'has no CB to cut at',
    'PRO':
    'has a sidechain that closes back onto the backbone nitrogen, '
    'so cutting at CA-CB leaves CD with a dangling valence and no cap',
}

# Net charge of each protonation variant, for the active site charge
# bookkeeping. Note that CYX here is Modeller's variant of that name - a
# cysteine with no HG, which for a metal-bound sidechain is a thiolate -
# and not Amber's disulfide-bridged CYX, which is neutral.
# Every residue the truncation can cut needs an entry, since
# include_residue can put any of them in the cluster.
VARIANT_CHARGES = {
    'ASP': -1,
    'ASH': 0,
    'GLU': -1,
    'GLH': 0,
    'CYS': 0,
    'CYX': -1,
    'HID': 0,
    'HIE': 0,
    'HIP': 1,
    'HIN': -1,
    'HIS': 0,
    'LYS': 1,
    'LYN': 0,
    'ARG': 1,
    'TYR': 0,
    'SER': 0,
    'THR': 0,
    'ASN': 0,
    'GLN': 0,
    'MET': 0,
    'TRP': 0,
    'ALA': 0,
    'GLY': 0,
    'ILE': 0,
    'LEU': 0,
    'PHE': 0,
    'PRO': 0,
    'VAL': 0,
}

# Names of the intermediates kept in the working folder. Each expensive
# step writes its result under one of these as soon as it has it, and picks
# it up again on a later run, so a failure part way through does not cost
# the steps that already succeeded.
GEOMETRY_FILE = 'opt_active_site.xyz'
# written for inspection only: the crude pass is cheap enough to repeat,
# so unlike the others this one is never read back
MM_GEOMETRY_FILE = 'mm_active_site.xyz'
HESSIAN_FILE = 'hessian.txt'
CHARGES_FILE = 'partial_charges.txt'
ENZYME_SYSTEM_FILE = 'enzyme_system.xml'
# written by compute() once the metal terms are fitted, for whatever
# comes after the run to read. Never read back by a run of its own: the
# fit is cheap next to everything that feeds it, and it belongs to the
# geometry it was made from. The crude pre-QM pass builds a force field
# of its own and is deliberately not written here, since it carries
# seeded force constants rather than fitted ones.
FORCEFIELD_FILE = 'forcefield.json'

# Written into the comment of an atom by annotate_atoms, so that a force
# field on its own still says which atoms these are. The truncation is
# what knows it, and nothing about a bare force field and a geometry
# recovers it afterwards without guessing.
BETA_CARBON_COMMENT = 'beta carbon'
CAP_COMMENT = 'capping hydrogen'

# What the generator writes on a term it found no parameters for, and on
# an atom GAFF could not type at all. A term marked this way carries a
# flat fallback constant rather than a tabulated or fitted one, which is
# what _check_atom_types looks for. Match these exactly: 'Guessed from
# Hessian' is a Seminario fit and 'X-c2-os-X(Guessed for ...)' is a
# wildcard match, and neither is a missing parameter.
UNPARAMETERIZED_COMMENT = 'Guessed'
UFF_TYPE_COMMENT = 'UFF'

# Literature equilibrium metal-ligand distances in nm. The crude pre-QM
# pass measures its equilibrium values on the input geometry instead;
# assign this to metal_bond_equilibria to impose these in their place.
LITERATURE_METAL_BONDS = {
    ('Zn', 'N'): 0.205,
    ('Zn', 'O'): 0.200,
    ('Zn', 'S'): 0.230,
}

# Force constants the crude pre-QM pass gives every metal term, in
# kJ/mol/nm^2 and kJ/mol/rad^2. Nothing at that stage says anything about the
# stiffness of a metal term, so they are flat: only the equilibrium values
# carry information until the Hessian is fitted.
DEFAULT_METAL_BOND_FORCE_CONSTANT = 100000.0
DEFAULT_METAL_ANGLE_FORCE_CONSTANT = 200.0

# Barrier, in kJ/mol, for the improper that nudges a metal into the plane
# of a coordinating histidine ring or a bidentate carboxylate
# (_add_metal_planarity_impropers). Deliberately weak -- of the same order
# as the generic sp2-planarity improper GAFF itself falls back on
# (1.1 kcal/mol, see MMForceFieldGenerator.populate_impropers) -- since it
# is a soft nudge towards the coordination geometry the lone pair actually
# favours, not a restraint the fit is meant to enforce.
DEFAULT_METAL_PLANARITY_FORCE_CONSTANT = 4.184

# ----------------------------------------------------------------------
# the shell and its decorators
# ----------------------------------------------------------------------


class Shell:
    """
    The base of the stateless phase classes: a communicator and an output
    stream, and nothing else. Settings arrive as keyword arguments and
    intermediates are returned, so a method takes what it uses and returns
    what it produces.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.
    """

    def __init__(self, comm=None, ostream=None):

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


# Per process, not per object: "am I inside a master-only section" is a
# fact about the process. A master-only method on one shell calling a
# master-only method on another must not broadcast from inside the
# master's body, and a per-instance counter could not see across objects.
_master_depth = 0


def on_master(method):
    """
    Runs the body on the master rank only and hands what it returned -- or
    the exception it raised -- to every rank. Nested calls run inline.

    The master keeps what it computed: bcast hands the root an unpickled
    copy of its own value as well, and a copy is not the same thing -- a
    force field generator comes back with a silent stream, and a large
    topology is pickled twice for nothing.

    :param method:
        The method to decorate.

    :return:
        The decorated method.
    """

    @functools.wraps(method)
    def wrapper(self, *args, **kwargs):
        global _master_depth

        if _master_depth > 0:
            return method(self, *args, **kwargs)

        _master_depth += 1
        try:
            outcome = None
            if self.rank == mpi_master():
                try:
                    outcome = ('value', method(self, *args, **kwargs))
                except Exception as error:
                    outcome = ('error', error)
            if self.nodes > 1:
                received = self.comm.bcast(outcome, root=mpi_master())
                if self.rank != mpi_master():
                    outcome = received
        finally:
            _master_depth -= 1

        kind, payload = outcome

        if kind == 'error':
            raise payload

        return payload

    return wrapper


def collective(method):
    """
    Marks a method every rank must enter, because it calls a VeloxChem
    driver on self.comm. Refuses to be called from inside an on_master
    body, which is the one way to deadlock under mpiexec.

    :param method:
        The method to decorate.

    :return:
        The decorated method.
    """

    @functools.wraps(method)
    def wrapper(self, *args, **kwargs):
        assert_msg_critical(
            _master_depth == 0, f'{method.__name__} is a collective step '
            'and was called from inside a master-only one')

        return method(self, *args, **kwargs)

    return wrapper


# ----------------------------------------------------------------------
# the shared line shapes
# ----------------------------------------------------------------------


def param(label, value, label_width=26, value_width=20):
    """
    Formats one parameter line with fixed label and value widths.

    print_header centers text, so all lines need the same total length to
    appear left-aligned relative to each other.

    :param label:
        The parameter name.
    :param value:
        The parameter value.
    :param label_width:
        The width of the label field.
    :param value_width:
        The width of the value field.

    :return:
        The formatted line.
    """

    return f'{label:<{label_width}} : {str(value):>{value_width}}'


def print_param_list(label, items, ostream, value_width=20):
    """
    Prints a list of values as parameter lines of uniform width.

    print_header centers each line, so a value that overflows the field
    would make its line start further left than the others. Long lists are
    therefore wrapped over several lines, with the label only on the
    first.

    :param label:
        The parameter name.
    :param items:
        The values to list.
    :param ostream:
        The output stream.
    :param value_width:
        The width of the value field.
    """

    chunks = []
    current = ''

    for item in items:
        candidate = item if not current else f'{current}, {item}'
        if len(candidate) > value_width and current:
            chunks.append(current + ',')
            current = item
        else:
            current = candidate

    if current:
        chunks.append(current)

    for i, chunk in enumerate(chunks):
        ostream.print_header(param(label if i == 0 else '', chunk))


def print_section(title, ostream):
    """
    Prints a title underlined to its own length.

    Hand-counting that length is how an underline ends up one character
    short of the title it sits under.

    :param title:
        The title.
    :param ostream:
        The output stream.
    """

    ostream.print_header(title)
    ostream.print_header(len(title) * '-')


# ----------------------------------------------------------------------
# plain helpers shared by more than one phase
# ----------------------------------------------------------------------


def _folder_file(name, folder=None):
    """
    Returns the path of an intermediate in the working folder, or None
    when it is not there.

    :param name:
        The file name, one of the class-level file name attributes.

    :return:
        The path, or None.
    """

    if folder is None:
        return None

    path = Path(folder) / name

    return path if path.is_file() else None


def residue_label(residue):
    """
    Returns the ASP130-style label a residue is named by.

    This label is load-bearing: it is what a user
    request is matched against, what update_protonation_state writes as a
    protonation_overrides key, and what active_site['residues'] holds. It
    is built in one place so that all of those agree on it.

    :param residue:
        The topology residue.

    :return:
        The label.
    """

    return f'{residue.name}{residue.id}'


def _site_index_map(active_site):
    """
    Maps a topology atom index back to the active site index that holds it.

    A capping hydrogen is mapped to the CA it replaces rather than to an
    atom of its own, so its entry is left out: the site holds no alpha
    carbon, and reading a cap's position back as that CA's is the one way
    this inversion goes wrong.

    :param active_site:
        The active site.

    :return:
        The dictionary from topology atom index to active site index.
    """

    caps = set(active_site['cap_indices'])

    return {
        top_index: site_index
        for site_index, top_index in active_site['atom_map'].items()
        if site_index not in caps
    }


def empty_request():
    """
    A new, empty record of what the user has decided about a metal site.

    - manual_bonds: the metal bonds added or removed by hand, as records
      of (residue index, atom name, metal residue index, action).
    - variants: the protonation variant of each residue, by residue index,
      as the protonation actually built it.
    - coordinating_residues: residues that are ligands whatever their
      distance.
    - extra_residues / excluded_residues: residues put into or taken out
      of the truncated cluster by hand.

    :return:
        The empty request.
    """

    return {
        'manual_bonds': [],
        'variants': {},
        'coordinating_residues': [],
        'extra_residues': [],
        'excluded_residues': [],
    }


def connectivity_bonds(connectivity_matrix):
    """
    Reads a connectivity matrix as the list of bonds Molecule.show draws.

    The pairs are zero-indexed and plain ints, since they are handed to
    RDKit, which does not take numpy integers.

    :param connectivity_matrix:
        The connectivity of an active site.

    :return:
        The bonds, as index pairs.
    """

    matrix = np.asarray(connectivity_matrix)

    return [(int(i), int(j)) for i, j in zip(*np.triu_indices_from(matrix, k=1))
            if matrix[i, j]]


def extract_pairs(connectivity_matrix, source_atoms, bond_count=2):
    """
    Finds the atom pairs needed for a pair-restricted Hessian.

    Walks the connectivity outwards from the source atoms and returns
    every bonded pair encountered within bond_count bonds. That is exactly
    what the Seminario method reads: a bond (i, j) uses the (i, j) block
    and an angle (i, j, k) uses the (i, j) and (j, k) blocks, never
    (i, k). A bond_count of one therefore covers the metal bonds and a
    bond_count of two additionally covers every angle involving a metal.

    :param connectivity_matrix:
        The connectivity matrix.
    :param source_atoms:
        The indices to walk out from, typically the metal centers.
    :param bond_count:
        The number of bonds to walk.

    :return:
        The tuple of the sorted pair list and the sorted atom list.
    """

    connectivity_matrix = np.asarray(connectivity_matrix)
    source_atoms = list(source_atoms)

    def neighbors(index):
        return set(np.where(connectivity_matrix[index])[0].tolist())

    visited = set(source_atoms)
    pairs = set()
    frontier = list(source_atoms)

    for _ in range(bond_count):
        next_frontier = []
        for index in frontier:
            for neighbor in neighbors(index):
                pairs.add((min(index, neighbor), max(index, neighbor)))
                if neighbor not in visited:
                    visited.add(neighbor)
                    next_frontier.append(neighbor)
        frontier = next_frontier

    pairs = sorted(pairs)
    atoms = sorted({index for pair in pairs for index in pair})

    return pairs, atoms


def get_metal_keys(forcefield, active_site):
    """
    Returns the bond and angle keys that involve a metal center.

    :param forcefield:
        The force field generator.
    :param active_site:
        The active site, for the indices of the metal centers.

    :return:
        The tuple of the bond key list and the angle key list.
    """

    metals = set(active_site['metal_indices'])

    bonds = [key for key in forcefield.bonds if metals & set(key)]
    angles = [key for key in forcefield.angles if metals & set(key)]

    return bonds, angles


def get_metal_impropers(forcefield, active_site):
    """
    Returns the improper keys that involve a metal center.

    Kept apart from get_metal_keys, whose two return values five callers
    unpack, and filtered to the metal the same way they are: every other
    improper of the site belongs to a residue the protein force field
    parameterizes itself, so transferring those would put a second copy of
    a term beside the one already there.

    :param forcefield:
        The force field generator.
    :param active_site:
        The active site, for the indices of the metal centers.

    :return:
        The improper key list.
    """

    metals = set(active_site['metal_indices'])

    return [key for key in forcefield.impropers if metals & set(key)]


def active_site_residues(binding_modes):
    """
    The residues the truncated active site holds.

    A residue is in the cluster when it coordinates a metal, or when
    include_residue asked for it, and out of it when remove_residue said
    so. Only those two requests are stored; the membership itself is
    worked out here every time it is needed, so that a re-detection which
    gains or loses a contact cannot leave a stored set behind.

    :param binding_modes:
        The binding modes.

    :return:
        The residue indices, sorted.
    """

    residues = {ligand['res_index'] for ligand in binding_modes['ligands']}
    residues |= set(binding_modes.get('extra_residues', []))
    residues -= set(binding_modes.get('excluded_residues', []))

    return sorted(residues)


def known_variants(res_name):
    """
    The protonation variants OpenMM will accept for a residue.

    Modeller keeps them in hydrogens.xml and refuses anything else from
    deep inside addHydrogens, with a bare ValueError. Asking it what it
    knows before it is called turns that into an error naming the legal
    set, and keeps the check and the thing it guards from drifting apart.

    :param res_name:
        The residue name, as the topology has it.

    :return:
        The variant names, as a tuple. Empty for a residue whose
        protonation Modeller does not offer a choice about.
    """

    assert_msg_critical('openmm' in sys.modules,
                        'known_variants: openmm is required')

    mmapp.Modeller._loadStandardHydrogenDefinitions()
    spec = mmapp.Modeller._residueHydrogens.get(res_name)

    if spec is None:
        return ()

    return tuple(spec.variants)


def check_variant(residue, variant):
    """
    Checks that a protonation variant can be asked for and paid for.

    Two things have to hold, and they come from different places: OpenMM
    has to know how to build the variant, and the charge bookkeeping of
    the active site has to know what it is worth. A variant that passes
    one and not the other fails much later and much less clearly - either
    inside addHydrogens or on the cluster charge - so both are checked
    here.

    :param residue:
        The residue the variant is for.
    :param variant:
        The variant name.
    """

    label = residue_label(residue)
    legal = known_variants(residue.name)

    assert_msg_critical(
        len(legal) > 0, 'check_variant: '
        f'{label} has no protonation variants to choose between; OpenMM '
        'builds it one way only')

    assert_msg_critical(
        variant in legal, 'check_variant: '
        f'{variant} is not a protonation variant of {label}. The variants '
        f'OpenMM knows for {residue.name} are {list(legal)}')

    assert_msg_critical(
        variant in VARIANT_CHARGES, 'check_variant: no charge known for '
        f'variant {variant}, so the charge of the active site could not be '
        'counted. Add it to VARIANT_CHARGES')


def _forcefield_elements(forcefield):
    """
    Returns the element of every atom of a force field generator.

    The generator names its atoms after the element followed by a counter,
    and reads the element back out of that name when it writes an OpenMM
    XML, so the name is where the element of a loaded force field lives.

    :param forcefield:
        The force field generator.

    :return:
        The element symbols in atom order.
    """

    elements = []

    for atom in forcefield.atoms.values():
        element = ''
        for character in atom['name']:
            if character.isdigit():
                break
            element += character
        elements.append(element)

    return elements


def _check_forcefield(forcefield, active_site, source=None):
    """
    Checks that a force field describes the extracted active site.

    Its bonds and angles are keyed by plain atom indices, which fit any
    cluster of the same size, so a force field belonging to another site
    would otherwise be applied to this one without a word.

    :param forcefield:
        The force field generator.
    :param active_site:
        The active site, to validate against.
    :param source:
        Where the pair came from, for the message. A caller loading one of
        many folders needs to be told which of them is the bad one; a caller
        with only the site in hand leaves it out.
    """

    labels = list(active_site['molecule'].get_labels())
    elements = _forcefield_elements(forcefield)
    named = f' of {source}' if source is not None else ''
    site = 'geometry' if source is not None else 'extracted active site'

    assert_msg_critical(
        len(elements) == len(labels),
        f'_check_forcefield: the force field{named} has {len(elements)} '
        f'atoms but the {site} has {len(labels)}')

    assert_msg_critical(
        elements == labels,
        f'_check_forcefield: the elements of the force field{named} do '
        f'not match the {site}, so it describes a different structure')


def save_forcefield(filename, forcefield):
    """
    Writes a force field to a JSON file.

    Only the parameters are written: the atoms, the bonds, the angles, the
    dihedrals and the impropers, with the redistributed charges already
    sitting on the atoms. The geometry is not part of it, which is why
    load_forcefield reads the molecule back off the active site.

    :param filename:
        The name of the JSON file.
    :param forcefield:
        The force field generator to write.
    """

    MMForceFieldGenerator.save_forcefield_as_json(forcefield, str(filename))


def load_forcefield(filename):
    """
    Reads a force field back from a JSON file.

    The file carries the parameters alone: no molecule, and nothing that
    says which site they belong to. A caller that wants to write OpenMM
    files or minimize on it attaches the molecule itself, and one loading
    from a folder it did not write checks the pair with _check_forcefield
    first, as templates.build does. partial_charges is restored from the
    charges of the atoms, which are the redistributed ones that
    build_forcefield wrote.

    :param filename:
        The name of the JSON file.

    :return:
        The force field generator.
    """

    assert_msg_critical(
        Path(filename).is_file(),
        f'load_forcefield: forcefield file {filename} not found')

    forcefield = MMForceFieldGenerator.load_forcefield_from_json_file(
        str(filename))
    forcefield.partial_charges = np.array(
        [atom['charge'] for atom in forcefield.atoms.values()])

    return forcefield


def _bond_separation(bonded, first, second, limit=3):
    """
    How many bonds apart two atoms are, up to a limit.

    :return:
        The number of bonds on the shortest path, or limit + 1 when there
        is none that short.
    """

    reached = {first}
    frontier = [first]
    for distance in range(1, limit + 1):
        frontier = [
            neighbour for atom in frontier
            for neighbour in bonded.get(atom, ()) if neighbour not in reached
        ]
        if not frontier:
            break
        if second in frontier:
            return distance
        reached.update(frontier)

    return limit + 1

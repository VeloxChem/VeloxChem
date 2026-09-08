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
from copy import deepcopy
from itertools import permutations, product
from pathlib import Path
import numpy as np
import math
import sys

from networkx.algorithms.isomorphism import GraphMatcher
import networkx as nx

try:
    from scipy.optimize import linear_sum_assignment
except ImportError:
    pass

from ..veloxchemlib import mpi_master
from ..molecule import Molecule
from ..outputstream import OutputStream
from .metalsiteffbuilder import MetalSiteForceFieldBuilder
from . import core
from ..optimizationdriver import OptimizationDriver
from ..superimpose import svd_superimpose
from ..errorhandler import assert_msg_critical


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
    How many atoms of a residue reach a metal is left out of the matching, and
    the site is wired the way the template is wired before its force field is
    built: a template that is monodentate builds a monodentate site out of a
    structure the cutoffs called bidentate, and the other way round.

    One workflow, three calls:

        manager.load_templates_from_folders([...])
        manager.compare_active_site(builder)   # builder has built_active_site
        manager.build_ff_from_template()

    At most one active site is tracked at a time -- the MetalSiteForceFieldBuilder
    handed to compare_active_site, read back through the active_site property.
    It can be edited in between (add_metal_bond, update_protonation_state, ...)
    and compared again with compare_active_site(), with no argument, to pick
    up the edit; build_ff_from_template hands the transferred force field back
    to that same builder so create_enzyme_system can be called on it directly.
    There is no other way to reach a result: a bare structure path is not
    accepted here, since building the active site is exactly what the builder
    already does.

    Where a comparison comes back a mismatch, shoehorn(name) makes those
    edits itself: it walks the site onto a named template's residues,
    coordination and protonation, so that comparing it again succeeds. It
    edits nothing else and builds nothing. A later unnamed
    build_ff_from_template then builds from that template rather than from
    whatever the criteria rank first, and built_from names the one it used.

    Settings are kept on the object, apart from the templates and the active
    site itself.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - templates: The loaded templates, keyed by name.
        - builder: A MetalSiteForceFieldBuilder held for its settings alone --
          never given an active site of its own. The steps themselves are the
          functions of the core module, which are called with the settings this
          carries, so set them on it. Distinct from active_site, which is the
          builder that holds the real, loaded site.
        - metal_shell_bonds: How many bonds out from a metal the metal_shell
          region reaches.
        - selection_criteria: What a template has to be within before
          build_ff_from_template puts its parameters on a site: 'tight',
          'loose', or a set of thresholds of the same shape as those in
          SELECTION_CRITERIA. This is the only thing a match is decided on.
        - rmsd_heavy_atoms_only: The flag for leaving the hydrogens out of
          every RMSD, cartesian and internal alike. On by default: protonate
          places the hydrogens through Modeller.addHydrogens, which does not
          do it reproducibly, so two runs over one structure differ by around
          0.3 A over all atoms and by nothing at all over the heavy ones. Both
          cartesian numbers are printed either way.
        - max_mappings: The limit on the number of isomorphisms evaluated for
          one template.
        - mm_fallback_literature_bonds: The flag for pulling the metal bonds
          of the crude relaxation toward LITERATURE_METAL_BONDS, which
          normalizes the query toward the geometry a template was optimized
          to. Ignored when the builder already carries metal bond equilibria
          of its own.
        - comm: The MPI communicator.
        - rank: The rank of the MPI process.
        - nodes: The number of MPI processes.
        - ostream: The output stream.
    """

    # ------------------------------------------------------------------
    # constants
    # ------------------------------------------------------------------

    # The geometry a run leaves behind, in the order it is looked for. Which
    # of them a template is allowed to be built from is what the fallback
    # argument of load_template_from_folder decides.
    GEOMETRY_KINDS = ('qm_opt', 'mm_opt')

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
    # A region is read as {ic type: {rms: , max: }}, so that angles can be
    # added to a set without anything here changing. The tight set is what
    # two runs of the same site look like; the loose one is what an
    # unrelaxed design model looks like.
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
        # is exposed rather than wrapped. It never gets an active site of its
        # own -- that is _active_site_builder, below.
        self.builder = MetalSiteForceFieldBuilder(comm, ostream)

        # the one active site the manager tracks, and the last comparison
        # made against it; see the active_site property and
        # compare_active_site
        self._active_site_builder = None
        self._comparison = None
        # the template the tracked site was last walked onto, which an
        # unnamed build_ff_from_template builds from whatever the criteria
        # come to. Cleared whenever a different site is handed in, since a
        # shoehorning is a statement about one site and not about the
        # manager.
        self._shoehorned = None
        # the template the last force field was built from, which is not
        # always the one the criteria rank first -- see built_from
        self._built_from = None

        # matching
        self.metal_shell_bonds = 2

        # What a template has to be within before its parameters are put on a
        # site. A name in SELECTION_CRITERIA or a set of thresholds of the
        # same shape; the dictionary is handed out rather than the name so
        # that one threshold can be changed without writing the rest.
        self.selection_criteria = self.SELECTION_CRITERIA['tight']
        # the hydrogens carry the noise of Modeller.addHydrogens rather than
        # anything about the site; see the class docstring
        self.rmsd_heavy_atoms_only = True
        self.max_mappings = 10000

        # what to run on the query structure
        self.mm_fallback_literature_bonds = True

    @property
    def active_site(self):
        """
        The MetalSiteForceFieldBuilder tracked as the current active site, or
        None before compare_active_site has been given one. Editing happens
        on this object directly (add_metal_bond, update_protonation_state,
        show_active_site, print_active_site, ...); create_enzyme_system is
        called on it once build_ff_from_template has adopted a force field
        onto it.
        """

        return self._active_site_builder

    @property
    def built_from(self):
        """
        The name of the template the last force field was built from, or None
        before build_ff_from_template has been called on the tracked site.

        Worth asking rather than assuming: an unnamed call is not always
        answered by the template the criteria rank first, since a shoehorning
        overrides them -- see _prefer_the_shoehorned_template.
        """

        return self._built_from

    # ------------------------------------------------------------------
    # templates
    # ------------------------------------------------------------------

    def load_template_from_folder(self, folder, name=None, fallback=None):
        """
        Loads one template from a folder left behind by a builder run.

        Two files are read: the force field, which carries the fitted metal
        terms, the charges and the bonding topology of the active site, and the
        geometry the fit was done on. Nothing else of the folder is needed, and
        the atom indices of the force field are not assumed to mean anything
        outside it: the mapping onto a new site is solved when a match is made.

        A residue of the folder that coordinates none of the metals is
        discarded here: it is a neighbour the truncation happened to keep,
        nothing is transferred from it, and every structure matched against
        the template would have to hold it too. The charge it takes with it
        is put back onto the atoms that stay.

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

        ff_path = folder / core.FORCEFIELD_FILE
        assert_msg_critical(
            ff_path.is_file(),
            f'MetalForceFieldManager: {ff_path} not found, so {folder} holds '
            'no force field to use as a template. It is written by '
            'MetalSiteForceFieldBuilder.build_forcefield.')

        forcefield = core.load_forcefield(ff_path)
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

    def load_templates_from_folders(self, folders):
        """
        Loads several template folders at once, each under the name
        load_template_from_folder would give it on its own.

        :param folders:
            The folders of earlier runs.
        """

        for folder in folders:
            self.load_template_from_folder(folder)

    def show_template(self, name, **kwargs):
        """
        Draws one loaded template, with its own explicit metal-ligand bonds.

        core.show_active_site is not used here: it reads
        active_site['labels'] and active_site['connectivity_matrix'], and a
        template dictionary (built by _build_template) has neither -- it
        carries forcefield.bonds instead, which is the same source every
        other template computation (_describe, _metal_keys, ...) reads its
        bonds from.

        :param name:
            The name of the template to show.
        :param kwargs:
            Further keyword arguments for Molecule.show. Passing bonds
            explicitly overrides what is worked out here.

        :return:
            Whatever Molecule.show returns.
        """

        known = name in self.templates
        assert_msg_critical(
            known, f'MetalForceFieldManager.show_template: no template '
            f'named {name}. Loaded: {sorted(self.templates)}')

        template = self.templates[name]
        kwargs.setdefault('bonds', list(template['forcefield'].bonds))

        return template['molecule'].show(**kwargs)

    def show_all_templates(self, **kwargs):
        """
        Draws every loaded template at once, in a grid, each with its own
        explicit metal-ligand bonds.

        :param kwargs:
            Further keyword arguments for Molecule.show_grid (e.g. grid,
            linked, width, height). Passing bonds explicitly overrides what
            is worked out here for every pane alike.

        :return:
            Whatever Molecule.show_grid returns.
        """

        names = list(self.templates)
        assert_msg_critical(
            len(names) > 0,
            'MetalForceFieldManager.show_all_templates: no templates '
            'loaded. Call load_template_from_folder first.')

        molecules = [self.templates[name]['molecule'] for name in names]
        bonds = [
            list(self.templates[name]['forcefield'].bonds) for name in names
        ]
        kwargs.setdefault('bonds', bonds)

        # print a grid of titles for reference

        return Molecule.show_grid(molecules, **kwargs)

    def _load_geometry(self, folder, fallback):
        """
        Reads the active site geometry of a template folder and works out what
        kind of geometry it is.

        :param folder:
            The folder of an earlier run.
        :param fallback:
            The lowest acceptable kind, as given to load_template_from_folder.

        :return:
            The tuple of the molecule and its kind.
        """

        opt_path = folder / core.GEOMETRY_FILE
        mm_path = folder / core.MM_GEOMETRY_FILE

        allowed = self.GEOMETRY_KINDS[:1 if fallback is None else 2]

        if opt_path.is_file():
            geometry = Molecule.read_xyz_file(str(opt_path))

            if mm_path.is_file() and self._same_geometry(
                    geometry, Molecule.read_xyz_file(str(mm_path))):
                # build_forcefield writes whatever the active site ended
                # up as under the name of the optimized geometry, so an
                # untouched copy of the MM geometry means the QM optimization
                # never ran
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
        capping hydrogens and the beta carbons come from the comments
        annotate_atoms wrote, which is how a cap is told apart from the
        hydrogens it is symmetric with.

        What the folder holds is not taken as it stands: the residues that
        coordinate no metal are discarded by _prune_unconnected_residues
        before the template is handed back.

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
        elements = core._forcefield_elements(forcefield)

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

        # The truncation records what it made in the comments, and those are
        # what a template reads. The capping hydrogens do sit at exactly zero
        # charge once redistribute_cap_charges has run, but that is the
        # correction doing its job rather than a marker, and it says nothing
        # in a force field whose charges are all zero.
        marked = {
            role: [
                index for index, atom in forcefield.atoms.items()
                if role in (atom.get('comment', '') or '')
            ]
            for role in (core.BETA_CARBON_COMMENT, core.CAP_COMMENT)
        }

        beta_carbon_indices = marked[core.BETA_CARBON_COMMENT]
        cap_indices = marked[core.CAP_COMMENT]

        assert_msg_critical(
            len(cap_indices) > 0,
            f'MetalForceFieldManager: the force field of {folder} does not '
            'say which of its atoms are capping hydrogens. It was written '
            'before annotate_atoms recorded them, so rebuild it.')

        # A template is an active site with the two files it was loaded from
        # added to it, rather than a dictionary holding a description of one:
        # every matching step then takes two dictionaries of one shape, and a
        # template is simply an active site that has no topology behind it.
        template = {
            'molecule': molecule,
            'metal_indices': metal_indices,
            'cap_indices': cap_indices,
            'beta_carbon_indices': beta_carbon_indices,
            'name': name,
            'forcefield': forcefield,
            'geometry_kind': kind,
            'charges': charges,
            'folder': str(folder),
        }

        described = self._describe(template, forcefield.bonds.keys())

        return self._prune_unconnected_residues(described)

    def _prune_unconnected_residues(self, template):
        """
        Discards the residues of a template that coordinate no metal.

        A template exists for the coordination it was fitted for, and a
        residue reaching none of the metals says nothing about that one: it
        is a neighbour the truncation happened to keep. Nothing is ever
        transferred from it -- the metal bonds and angles that are, are keyed
        on the residues that stay -- while every structure compared against
        the template has to hold a residue of its kind, in the right place,
        or map onto nothing at all. Dropping it therefore widens what the
        template matches and changes nothing it is used for.

        The atoms of a discarded residue take their charge with them, which
        leaves what is left off a whole number of electrons.
        _compensate_charges puts that back.

        :param template:
            The described template.

        :return:
            The template with those residues dropped and described again, or
            the object it was given when every residue coordinates a metal.
        """

        coarse = template['coarse_topology']
        unconnected = [
            node for node in self._residue_nodes(coarse)
            if coarse.degree(node) == 0
        ]

        if not unconnected:
            return template

        discard = sorted(atom for node in unconnected
                         for atom in coarse.nodes[node]['atoms'])

        pruned, shift = self._drop_atoms(template, discard)
        self._print_discarded_residues(template, unconnected, discard, shift)

        return self._describe(pruned, pruned['forcefield'].bonds.keys())

    def _drop_atoms(self, template, discard):
        """
        Removes atoms from a template and renumbers everything left.

        A template is indexed by atom throughout -- the geometry, the
        charges, the marked capping hydrogens and beta carbons, and every key
        of the force field -- so all of them are rewritten together against
        one map from the old index to the new. The force field generator is
        rewritten in place: it is loaded for this template alone and is
        shared with nothing.

        :param template:
            The described template.
        :param discard:
            The indices of the atoms to remove.

        :return:
            The tuple of a new template dictionary, without the topologies
            _describe adds to one, and the charge every remaining atom was
            shifted by.
        """

        molecule = template['molecule']
        labels = molecule.get_labels()
        dropped = set(discard)
        keep = [index for index in range(len(labels)) if index not in dropped]
        index_of = {old: new for new, old in enumerate(keep)}

        def renumbered(indices):
            return [index_of[index] for index in indices if index in index_of]

        forcefield = template['forcefield']
        forcefield.atoms = {
            index_of[index]: atom
            for index, atom in forcefield.atoms.items() if index in index_of
        }

        # a term of a discarded residue is internal to it, since a residue
        # that coordinates no metal shares no bond with the rest of the site;
        # the keys are filtered rather than assumed so anything else that
        # crosses the cut goes with it
        for name in ('bonds', 'angles', 'dihedrals', 'impropers'):
            table = getattr(forcefield, name)
            setattr(
                forcefield, name, {
                    tuple(index_of[index] for index in key): parameters
                    for key, parameters in table.items()
                    if all(index in index_of for index in key)
                })

        geometry = Molecule([labels[index] for index in keep],
                            molecule.get_coordinates_in_angstrom()[keep],
                            'angstrom')
        geometry.set_charge(molecule.get_charge())
        geometry.set_multiplicity(molecule.get_multiplicity())
        forcefield.molecule = geometry

        cap_indices = renumbered(template['cap_indices'])
        charges, shift = self._compensate_charges(template['charges'][keep],
                                                  cap_indices)

        # one array means the charges: what the template reports is what the
        # force field carries, here as everywhere else
        for index, atom in forcefield.atoms.items():
            atom['charge'] = float(charges[index])
        forcefield.partial_charges = charges

        pruned = {
            key: value
            for key, value in template.items()
            if key not in ('fine_topology', 'coarse_topology', 'composition')
        }
        pruned.update({
            'molecule': geometry,
            'metal_indices': renumbered(template['metal_indices']),
            'cap_indices': cap_indices,
            'beta_carbon_indices':
            renumbered(template['beta_carbon_indices']),
            'charges': charges,
        })

        return pruned, shift

    @staticmethod
    def _compensate_charges(charges, cap_indices):
        """
        Puts the charge the discarded residues took with them back onto the
        atoms that stay.

        What is left of a template is a whole number of electrons, and the
        atoms that went carried a share of the charge that is not one, so the
        difference from the nearest integer is added in equal parts to every
        remaining atom. The capping hydrogens are left out of it, for the
        reason redistribute_cap_charges empties them in the first place: a
        cap stands for a bond to a protein the site does not hold, and charge
        put on one is taken off again by the next redistribution.

        :param charges:
            The charges of the atoms that stay.
        :param cap_indices:
            Which of those are capping hydrogens.

        :return:
            The tuple of the compensated charges and the shift each atom was
            given.
        """

        charges = np.array(charges, dtype=float)
        caps = set(cap_indices)
        rest = [index for index in range(charges.size) if index not in caps]

        carriers = len(rest) > 0
        assert_msg_critical(
            carriers, 'MetalForceFieldManager: discarding the residues that '
            'coordinate no metal leaves nothing but capping hydrogens to '
            'carry the charge')

        total = float(np.sum(charges))
        shift = (round(total) - total) / len(rest)
        charges[rest] += shift

        return charges, shift

    def _print_discarded_residues(self, template, nodes, discard, shift):
        """
        Reports what loading a template threw away.

        :param template:
            The template as it stood before the atoms were dropped.
        :param nodes:
            The coarse nodes of the discarded residues.
        :param discard:
            The indices of the atoms that were dropped.
        :param shift:
            The charge every remaining atom was given to make up for them.
        """

        coarse = template['coarse_topology']
        charges = template['charges']

        named = ', '.join(
            sorted(f'{coarse.nodes[node]["formula"]}/'
                   f'{coarse.nodes[node]["key"][:6]}' for node in nodes))
        carried = float(sum(charges[index] for index in discard))

        self.ostream.print_info(
            f'Template {template["name"]}: discarding {len(nodes)} residue(s) '
            f'coordinating no metal ({named}), {len(discard)} atom(s) '
            f'carrying {carried:+.3f} e.')
        self.ostream.print_info(
            f'Shifting every remaining atom but the capping hydrogens by '
            f'{shift:+.4f} e to make the charge of the template whole again.')

        stray = carried - round(carried)
        if abs(stray) > 0.25:
            self.ostream.print_warning(
                f'Those residues carry {carried:+.3f} e, which is '
                f'{stray:+.3f} e away from a whole number of electrons, so '
                'what is left of the template was rounded onto the nearest '
                'integer and that may not be the charge of the site it '
                'describes.')

        self.ostream.flush()

    # ------------------------------------------------------------------
    # description
    #
    # A described site is the query and the template in one vocabulary:
    # a coarse residue graph, a heavy atom map and a formula per node.
    # Both the comparison and the shoehorning read it.
    # ------------------------------------------------------------------

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

    def _describe(self, active_site, edges):
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
                            family=self._family_key(heavy))

            for metal in metal_indices:
                if any(fine.has_edge(metal, atom) for atom in nodes):
                    coarse.add_edge(node, ('metal', metal))

        return {
            **active_site,
            'fine_topology': fine,
            'coarse_topology': coarse,
            'composition': sorted(labels),
        }

    @staticmethod
    def _residue_nodes(coarse):
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
    def _site_spec(described):
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

        residues = MetalForceFieldManager._residue_nodes(coarse)

        spec['residues'] = sorted(coarse.nodes[node]['key']
                                  for node in residues)
        spec['bridging'] = sorted(coarse.nodes[node]['key'] for node in residues
                                  if coarse.degree(node) > 1)

        return spec

    # ------------------------------------------------------------------
    # comparison
    # ------------------------------------------------------------------

    def compare_active_site(self,
                            active_site=None,
                            mm_opt=True,
                            include_hydrogens=False):
        """
        Compares the current active site against every loaded template,
        prints the full comparison, and reports whether an unforced match
        exists.

        An active site is a MetalSiteForceFieldBuilder on which
        build_active_site has already been called. Handing one in here
        replaces whatever active site was tracked before -- there is at most
        one at a time -- and it is then read back through the active_site
        property, so it can be edited (add_metal_bond,
        update_protonation_state, ...) and compared again by calling this
        with no argument, which describes it again from scratch and so picks
        up any edit made since the last call.

        The hydrogens are left out of every measurement by default, since
        protonate does not place them reproducibly; include_hydrogens puts
        them back for this call whatever rmsd_heavy_atoms_only says. One
        geometry is measured per call, so call it twice to see both.

        :param active_site:
            A MetalSiteForceFieldBuilder with a built active site, to become
            the one the manager tracks. None reuses the one already tracked.
        :param mm_opt:
            Whether to relax the active site on a crude force field first,
            which is what takes the slack out of an unrelaxed structure.
        :param include_hydrogens:
            Whether the hydrogens count toward the RMSDs and the internal
            coordinates.

        :return:
            True when at least one template is within selection_criteria,
            i.e. an unnamed build_ff_from_template() call would succeed.
        """

        assert_msg_critical(
            len(self.templates) > 0,
            'MetalForceFieldManager.compare_active_site: no templates '
            'loaded. Call load_template_from_folder first.')

        if active_site is not None:
            is_builder = isinstance(active_site, MetalSiteForceFieldBuilder)
            assert_msg_critical(
                is_builder,
                'MetalForceFieldManager.compare_active_site: active_site '
                'must be a MetalSiteForceFieldBuilder, got '
                f'{type(active_site).__name__}')
            has_site = active_site.active_site_molecule is not None
            assert_msg_critical(
                has_site,
                'MetalForceFieldManager.compare_active_site: this builder '
                'has no active site yet. Call build_active_site on it '
                'first.')
            self._active_site_builder = active_site
            # both belonged to the site being replaced
            self._shoehorned = None
            self._built_from = None

        have_builder = self._active_site_builder is not None
        assert_msg_critical(
            have_builder,
            'MetalForceFieldManager.compare_active_site: no active site '
            'loaded. Call compare_active_site with a '
            'MetalSiteForceFieldBuilder first.')

        builder = self._active_site_builder
        described = self._described_site(builder)

        if mm_opt:
            molecule = self._mm_relax({'active_site': described})
            geometry = 'mm_relaxed'
        else:
            molecule = described['molecule']
            geometry = 'input'

        coordinates = molecule.get_coordinates_in_angstrom()
        heavy_only = not include_hydrogens

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

            if template['composition'] != described['composition']:
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
            'source': str(builder.output_folder),
            'geometry': geometry,
            'include_hydrogens': include_hydrogens,
            'active_site': described,
            # the geometry every number above was measured on, which is the
            # relaxed one rather than the site's own when mm_opt is set
            'molecule': molecule,
            'templates': findings,
        }

        self._comparison = results
        self._print_comparison(results)

        return self._select_template(None)['name'] is not None

    # ------------------------------------------------------------------
    # selection
    # ------------------------------------------------------------------

    def _select_template(self, template=None):
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

        assert_msg_critical(
            self._comparison is not None,
            'MetalForceFieldManager._select_template: no comparison yet. '
            'Call compare_active_site first.')

        comparison = self._comparison
        criteria_name, criteria = self._selection_criteria()

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
                'MetalForceFieldManager._select_template: no template named '
                f'{template} was compared. Loaded: '
                f'{sorted(comparison["templates"])}')

            entry = comparison['templates'][template]
            verdict = decision['verdicts'][template]

            assert_msg_critical(
                self._maps_every_atom(comparison, entry),
                f'MetalForceFieldManager._select_template: template {template} '
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

        self._print_selection(comparison, decision)

        if decision['name'] is None:
            self._print_no_selection(decision)

        return decision

    def _print_no_selection(self, decision):
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

    def _prefer_the_shoehorned_template(self, decision):
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

        name = self._shoehorned

        if name is None or name not in self._comparison['templates']:
            return decision

        if name == decision['name']:
            return decision

        entry = self._comparison['templates'][name]

        if not self._maps_every_atom(self._comparison, entry):
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

    def build_ff_from_template(self, template=None):
        """
        Builds a force field for the current active site out of the template
        that fits it, and adopts it onto the active site so
        create_enzyme_system can be called on it directly.

        This is the half of a builder run that QM was paid for: the fitted
        metal bonds and angles and the RESP charges are taken from a template
        that describes the same site, and everything else is built for the
        site in front of us. Which template that is: with no name, the best
        match under selection_criteria, on the numbers compare_active_site
        measured; named, that template exactly. Either way a template that
        does not pass -- outside the criteria, or no match at all -- is a
        hard error naming why, since transferring parameters that were
        measured not to fit would silently mis-parameterize the site.

        The one exception is a site shoehorn walked onto a template: with no
        name given, that template is taken whatever the criteria came to --
        nothing within them, or another template within them and ranked
        ahead of it -- with a warning saying so and naming the verdict.
        Walking a site onto a template is already the statement that naming
        one makes, and the criteria cannot answer this case on their own:
        what they measure is a geometry whose coordination sphere is still
        open, since the parameters that would close it are the ones being
        asked for, and once it is open enough to pass at all it passes
        against the whole family the template belongs to. The mapping is
        still required to be complete: with an atom landing nowhere there is
        nothing to transfer, and that stays a hard error. Relaxing on the
        transferred force field and comparing once more is what turns this
        into an ordinary match.

        The coordination comes from the template too. Denticity is a distance
        cutoff on an unrelaxed structure rather than chemistry, so a site that
        matched a template on every other count is wired the way the template
        is wired before the force field is built, and what that changed is
        reported by _print_forced_bonds.

        :param template:
            The name of the template to use, or None to choose the best
            match.

        :return:
            The force field, carrying the template's transferred metal terms
            and charges.
        """

        assert_msg_critical(
            self._comparison is not None,
            'MetalForceFieldManager.build_ff_from_template: no comparison '
            'yet. Call compare_active_site first.')

        decision = self._select_template(template)

        if template is None:
            decision = self._prefer_the_shoehorned_template(decision)

        matched = decision['name'] is not None
        assert_msg_critical(
            matched,
            'MetalForceFieldManager.build_ff_from_template: no template is '
            f'within the {decision["criteria_name"]} criteria; see the '
            'comparison above, or loosen selection_criteria.')

        if template is not None:
            within_criteria = decision['verdicts'][template] is None
            assert_msg_critical(
                within_criteria,
                f'MetalForceFieldManager.build_ff_from_template: {template} '
                'does not match the active site: '
                f'{decision["verdicts"][template]}')

        name = decision['name']
        entry = decision['entry']
        self._built_from = name
        template_obj = self.templates[name]
        active_site = self._comparison['active_site']

        self.ostream.print_info(
            f'Building the force field from {name}, measured on the '
            f'{self._comparison["geometry"]} geometry. Transferring its '
            'metal terms and charges.')
        self.ostream.flush()

        forcefield, _, active_site = self._build_ff_from_template(
            template_obj, entry['mapping'], active_site)

        # strip the manager-only description keys before handing the site
        # back to the builder, which never produces them
        builder_active_site = {
            key: value
            for key, value in active_site.items()
            if key not in ('fine_topology', 'coarse_topology')
        }

        forcefield = self._active_site_builder.adopt_forcefield(
            forcefield, active_site=builder_active_site)

        core._print_metal_parameters(active_site,
                                     forcefield,
                                     ostream=self.ostream)

        return forcefield

    # ------------------------------------------------------------------
    # shoehorning
    #
    # Walks a site onto a template by editing it through the builder's
    # own public edit methods, and restores the request on any failure.
    # ------------------------------------------------------------------

    def shoehorn(self, template, max_include_radius=7.0):
        """
        Edits the current active site until it is built the way a template is
        built.

        compare_active_site refuses a template unless the site is made of the
        same residues, coordinated the same way and protonated the same way,
        and each of those three is something the cutoffs can get wrong on an
        unrelaxed structure rather than something about the chemistry. This
        walks the site onto the template's terms instead of refusing it:

        1. The residues. One the template holds and the site does not is
           looked for outward from the metal centers, out to
           max_include_radius, and the closest of the right kind is included.
           Which amino acid a residue is is decided on its heavy atoms alone,
           so a protonation that differs cannot hide it.
        2. The coarse topology -- which residue coordinates which metal. The
           metal centers are paired on what coordinates them, each pair is
           given the residues the template puts on it, and what is left over
           is dropped. Nothing is dropped until every metal has been dealt
           with: a residue on its way out of one center can be the one still
           holding another together.
        3. The protonation, and only then the denticity. In that order,
           because the atom mapping that says which of the site's atoms is
           which of the template's needs the hydrogens to agree before it can
           be solved, and that mapping is what the metal-ligand bonds are
           then forced through -- so a carboxylate the structure holds
           bidentate is made monodentate where the template is monodentate,
           and the other way round.

        Every edit goes through the builder's own edit methods, so what comes
        out is a site the builder could have been walked to by hand. Nothing
        is built here, but the site is measured again before returning: the
        edits are the whole reason the last comparison no longer describes
        it, and leaving that to the caller only invited
        build_ff_from_template to be answered out of stale numbers. The next
        call is build_ff_from_template.

        Which template was walked onto is remembered, and an unnamed
        build_ff_from_template falls back to it when the criteria pick
        nothing -- see there for why that is not the same as ignoring them.

        A failure leaves the active site as it was found -- the record of the
        edits is put back and the site rebuilt from it -- and says what stood
        in the way.

        :param template:
            The name of the template to shoehorn the active site into.
        :param max_include_radius:
            How far out from a metal center, in Angstrom, a residue may be
            picked up from.

        :return:
            True when the active site was walked onto the template, False
            when something the structure does not hold stood in the way.
        """

        have_site = self._active_site_builder is not None
        assert_msg_critical(
            have_site,
            'MetalForceFieldManager.shoehorn: no active site loaded. Call '
            'compare_active_site with a MetalSiteForceFieldBuilder first.')

        known = template in self.templates
        assert_msg_critical(
            known, f'MetalForceFieldManager.shoehorn: no template named '
            f'{template}. Loaded: {sorted(self.templates)}')

        builder = self._active_site_builder
        target = self.templates[template]

        # a run that fails leaves no record of one that succeeded
        self._shoehorned = None

        # everything the builder is is a function of the record of the edits
        # made on it, which is what its own edit methods write and what
        # build_active_site() rebuilds from. Keeping a copy of it is
        # therefore the whole of undoing a run that fails part way through.
        snapshot = deepcopy(builder._request)

        self._print_shoehorn_header(builder, target, max_include_radius)

        # Every edit rebuilds the site, and each rebuild reports the whole
        # cluster and relaxes it again. Neither says anything about what to
        # edit next -- the decisions are read off the protein, which the
        # crude relaxation does not touch -- so the builder is quietened and
        # its relaxation held back until there is a site worth relaxing. The
        # stream is swapped rather than muted because mute is reference
        # counted, and the manager's own stream may be this very object.
        stream = builder.ostream
        relaxing = builder._mm_opt

        builder.ostream = OutputStream(None)
        builder._mm_opt = False

        try:
            reason = self._shoehorn(builder, target, float(max_include_radius))
        except Exception:
            builder._mm_opt = relaxing
            self._restore(builder, snapshot)
            builder.ostream = stream
            raise

        builder._mm_opt = relaxing

        if reason is not None:
            self._restore(builder, snapshot)
            builder.ostream = stream
            self.ostream.print_warning(
                f'Could not shoehorn the site into {template}: {reason}. The '
                'active site is as it was.')
            self.ostream.flush()
            return False

        if relaxing:
            # what comes out is what a build_active_site would have left
            builder.build_active_site()

        builder.ostream = stream

        self._print_shoehorn_summary(builder, target)

        self._shoehorned = template

        # the edits are exactly what makes the last comparison stale, so the
        # site is measured again here rather than by a caller who has to
        # remember to
        self.compare_active_site()

        return True

    def _restore(self, builder, request):
        """
        Puts an active site back the way it was before a shoehorning.

        :param builder:
            The builder holding the site.
        :param request:
            The record of its edits, as it stood beforehand.
        """

        builder._request = deepcopy(request)
        builder.build_active_site()

    def _shoehorn(self, builder, template, max_include_radius):
        """
        The three stages of a shoehorning, in order.

        :param builder:
            The builder holding the site to edit.
        :param template:
            The template to edit it onto.
        :param max_include_radius:
            How far out from a metal center a residue may be picked up from.

        :return:
            What stood in the way, or None when nothing did.
        """

        reason = self._shoehorn_composition(builder, template,
                                            max_include_radius)
        if reason is not None:
            return reason

        # the protonation first: the mapping the denticity is forced through
        # cannot be solved while the hydrogens still differ
        for stage in (self._shoehorn_protonation, self._shoehorn_denticity):
            reason = stage(builder, template)
            if reason is not None:
                return reason

        return None

    def _described_site(self, builder):
        """
        Describes the active site of a builder as it now stands.

        Derived on every call rather than kept: every edit rebuilds the site
        and renumbers its atoms, so a description from before one says
        nothing about the site after it.

        :param builder:
            The builder holding the site.

        :return:
            The active site with its two topologies added.
        """

        active_site = builder.active_site

        return self._describe(
            active_site,
            core.connectivity_bonds(active_site['connectivity_matrix']))

    @staticmethod
    def _family_key(heavy):
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

    @staticmethod
    def _sidechain_heavy_atoms(residue):
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

    def _residue_family_key(self, topology, residue):
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
            atom.index: atom for atom in self._sidechain_heavy_atoms(residue)
        }

        graph = nx.Graph()

        for index, atom in atoms.items():
            graph.add_node(index, elem=atom.element.symbol)

        for first, second in topology.bonds():
            if first.index in atoms and second.index in atoms:
                graph.add_edge(first.index, second.index)

        return self._family_key(graph)

    def _candidate_residues(self, builder, max_radius):
        """
        Every residue a shoehorning is allowed to build the site out of.

        The site's own residues and everything else within max_radius of a
        metal center, each measured twice: how close its sidechain reaches a
        center at all, and how close its nearest donor atom comes, which is
        what a bond would be made through.

        :param builder:
            The builder holding the site.
        :param max_radius:
            How far out from a metal center, in Angstrom, to look.

        :return:
            The candidates, each naming its residue, its family and its
            distances to every metal center by that center's residue index.
        """

        topology = builder.enzyme_topology
        positions = np.asarray(builder.enzyme_positions)
        modes = builder.binding_modes

        members = set(core.active_site_residues(modes))
        metals = {
            entry['res_index']: positions[entry['index']]
            for entry in modes['metals']
        }

        candidates = []

        for residue in topology.residues():
            if residue.name in core.UNTRUNCATABLE_RESIDUES:
                continue

            atoms = self._sidechain_heavy_atoms(residue)
            if not atoms:
                continue

            reach = {}
            donor = {}

            for res_index, metal in metals.items():
                measured = [
                    (float(np.linalg.norm(positions[atom.index] - metal)),
                     atom) for atom in atoms
                ]
                reach[res_index] = min(distance for distance, _ in measured)
                donors = [(distance, atom) for distance, atom in measured
                          if atom.element.symbol in core.DONOR_ELEMENTS]
                if donors:
                    distance, atom = min(donors, key=lambda found: found[0])
                    donor[res_index] = (distance, atom.name)

            if residue.index not in members and min(
                    reach.values()) > max_radius:
                continue

            candidates.append({
                'res_index': residue.index,
                'resid': str(residue.id),
                'chain': str(residue.chain.id),
                'name': residue.name,
                'label': core.residue_label(residue),
                'key': self._residue_family_key(topology, residue),
                'reach': reach,
                'donor': donor,
                'member': residue.index in members,
            })

        return candidates

    def _template_slots(self, template):
        """
        What a template asks the site to be made of: one slot per residue of
        it, naming the kind of residue and the metal centers it is on.

        :param template:
            The template.

        :return:
            The slots, in the order the coarse topology holds them.
        """

        coarse = template['coarse_topology']
        slots = []

        for node in self._residue_nodes(coarse):
            slots.append({
                'key': coarse.nodes[node]['family'],
                'formula': coarse.nodes[node]['formula'],
                'metals': sorted(image[1] for image in coarse.neighbors(node)),
            })

        return slots

    def _metal_pairings(self, builder, template, described):
        """
        Every way the metal centers of a site can stand for those of a
        template.

        There are only ever a handful, so they are all tried and the one
        whose residues fit best is taken, rather than being guessed at from
        what coordinates them. The element has to agree outright: a zinc site
        is not a template for an iron one.

        The site's centers are named by the index of their residue, which no
        rebuild disturbs, unlike the atom indices every rebuild renumbers.

        :param builder:
            The builder holding the site.
        :param template:
            The template to match.
        :param described:
            The described active site.

        :return:
            One mapping of template metal index to site metal residue index
            per pairing.
        """

        template_labels = template['molecule'].get_labels()
        query_labels = described['molecule'].get_labels()

        template_metals = list(template['metal_indices'])
        query_metals = list(described['metal_indices'])

        if len(template_metals) != len(query_metals):
            return []

        pairings = []

        for order in permutations(query_metals):
            if any(template_labels[first] != query_labels[second]
                   for first, second in zip(template_metals, order)):
                continue
            pairings.append({
                first: self._metal_res_index(builder, described, second)
                for first, second in zip(template_metals, order)
            })

        return pairings

    @staticmethod
    def _slot_cost(slot, candidate, pairing, max_radius):
        """
        What it costs to build one of a template's residues out of one of the
        site's, as the total distance the bonds it asks for would span.

        :param slot:
            The template residue.
        :param candidate:
            The residue of the structure.
        :param pairing:
            Which site metal center stands for which of the template's.
        :param max_radius:
            The furthest a bond may be stretched to.

        :return:
            The cost, or infinity when this residue cannot fill this slot.
        """

        if candidate['key'] != slot['key']:
            return math.inf

        metals = [pairing[metal] for metal in slot['metals']]

        if len(metals) > 1 and candidate['name'] not in core.BRIDGING_RESIDUES:
            # an imidazole nitrogen has one lone pair in the ring plane, so a
            # second metal on it is an artifact however the template reads
            return math.inf

        if not metals:
            # a residue the template holds without coordinating anything;
            # what matters is only that it is in the same neighbourhood
            distance = min(candidate['reach'].values())
            return distance if distance <= max_radius else math.inf

        total = 0.0

        for metal in metals:
            found = candidate['donor'].get(metal)
            if found is None or found[0] > max_radius:
                return math.inf
            total += found[0]

        return total

    def _assign_residues(self, slots, candidates, pairing, max_radius):
        """
        Solves which residue of the structure fills which residue of the
        template, for one pairing of the metal centers.

        The whole site is assigned at once rather than one metal center at a
        time, which is what lets a residue move from one center to the other:
        a histidine the cutoffs put on the near metal is exactly the one the
        template wants on the far one, and no amount of filling one center's
        deficit from the outside will produce that.

        :param slots:
            The template's residues.
        :param candidates:
            The residues of the structure.
        :param pairing:
            Which site metal center stands for which of the template's.
        :param max_radius:
            The furthest a bond may be stretched to.

        :return:
            The tuple of the total cost and the slot to candidate mapping, or
            None when some residue of the template cannot be filled at all.
        """

        assert_msg_critical(
            'scipy' in sys.modules,
            'MetalForceFieldManager.shoehorn: scipy is required to solve '
            'which residue of the structure fills which of the template')

        if len(candidates) < len(slots):
            return None

        costs = np.array(
            [[
                self._slot_cost(slot, candidate, pairing, max_radius)
                for candidate in candidates
            ] for slot in slots])

        finite = np.isfinite(costs)

        if not finite.any(axis=1).all():
            # a slot nothing at all can fill; no assignment exists
            return None

        # linear_sum_assignment cannot be given infinities, so they are made
        # more expensive than any whole assignment of finite ones
        blocked = costs[finite].sum() + 1.0
        rows, columns = linear_sum_assignment(np.where(finite, costs, blocked))

        assignment = {}
        total = 0.0

        for row, column in zip(rows, columns):
            if not finite[row, column]:
                return None
            assignment[int(row)] = candidates[int(column)]
            total += float(costs[row, column])

        if len(assignment) != len(slots):
            return None

        return total, assignment

    def _assignment_failure(self, template, slots, candidates, pairings,
                            max_radius):
        """
        Says which of a template's residues the structure cannot supply.

        Reported off the pairing that got furthest, since that is the one the
        site came closest to being built like.

        :param template:
            The template that could not be matched.
        :param slots:
            Its residues.
        :param candidates:
            The residues of the structure.
        :param pairings:
            The pairings of the metal centers that were tried.
        :param max_radius:
            The furthest a bond may be stretched to.

        :return:
            What stood in the way.
        """

        best = None

        for pairing in pairings:
            missing = []
            for slot in slots:
                if any(
                        math.isfinite(
                            self._slot_cost(slot, candidate, pairing,
                                            max_radius))
                        for candidate in candidates):
                    continue
                missing.append(slot)
            if best is None or len(missing) < len(best):
                best = missing

        if best:
            formulas = ', '.join(sorted(slot['formula'] for slot in best))
            return (f'{template["name"]} is made of residues the structure '
                    f'does not have within {max_radius:.1f} A of the right '
                    f'metal center: {formulas}')

        # every residue of the template can be filled by something, but not
        # by enough different somethings at once: two of them are competing
        # for the one residue the structure has in reach
        counts = Counter(slot['key'] for slot in slots)
        available = Counter(candidate['key'] for candidate in candidates)
        short = [
            slot['formula'] for slot in slots
            if available[slot['key']] < counts[slot['key']]
        ]

        if short:
            return (f'{template["name"]} is made of more '
                    f'{", ".join(sorted(set(short)))} residues than the '
                    f'structure has within {max_radius:.1f} A of its metal '
                    'centers')

        return (f'the residues of the structure cannot be shared out over '
                f'those of {template["name"]}: two of them are wanted in '
                'places only one residue can reach')

    def _shoehorn_composition(self, builder, template, max_radius):
        """
        Makes the site hold the template's residues, coordinating the metal
        centers the template coordinates.

        Stages one and two of a shoehorning at once, because they are one
        question: which residue of the structure is to be which residue of
        the template. Answered as a single assignment over the whole site,
        for every way the metal centers can be paired up, and the cheapest
        answer is the one built -- so a residue is moved from one center to
        another where that is what the template says, rather than a center's
        deficit being filled from the outside while the residue that belongs
        in it sits on its neighbour.

        :param builder:
            The builder holding the site.
        :param template:
            The template to match.
        :param max_radius:
            How far out from a metal center a residue may be picked up from.

        :return:
            What stood in the way, or None when nothing did.
        """

        described = self._described_site(builder)

        # on the amino acids alone: the protonation is put right afterwards,
        # and until it is, an ASH is not the ASP a template holds
        if self._coarse_mappings(template, described, match_protonation=False):
            return None

        slots = self._template_slots(template)
        candidates = self._candidate_residues(builder, max_radius)
        pairings = self._metal_pairings(builder, template, described)

        if not pairings:
            return (f'the metal centers of {template["name"]} are not the '
                    'metal centers of the site')

        best = None

        for pairing in pairings:
            found = self._assign_residues(slots, candidates, pairing,
                                          max_radius)
            if found is None:
                continue
            if best is None or found[0] < best[0]:
                best = (found[0], pairing, found[1])

        if best is None:
            return self._assignment_failure(template, slots, candidates,
                                            pairings, max_radius)

        _, pairing, assignment = best

        self._apply_assignment(builder, template, slots, pairing, assignment)

        if not self._coarse_mappings(template,
                                     self._described_site(builder),
                                     match_protonation=False):
            return ('which residue coordinates which metal still differs from '
                    f'{template["name"]}')

        return None

    def _apply_assignment(self, builder, template, slots, pairing, assignment):
        """
        Edits the site into the assignment that was solved.

        The bonds the template does not make go first and the ones it makes
        second, so that a residue moving from one metal center to another is
        never bonded to both at once -- which for a histidine add_metal_bond
        refuses outright. What is left holding nothing is dropped last, once
        every center has what it is owed.

        :param builder:
            The builder holding the site.
        :param template:
            The template being matched.
        :param slots:
            Its residues.
        :param pairing:
            Which site metal center stands for which of the template's.
        :param assignment:
            Which residue of the structure fills which slot.
        """

        wanted = {}
        donors = {}

        for index, candidate in assignment.items():
            res_index = candidate['res_index']
            wanted.setdefault(res_index, set()).update(
                pairing[metal] for metal in slots[index]['metals'])
            donors[res_index] = candidate

        self._include_assigned(builder, template, donors)
        self._remove_unwanted_bonds(builder, wanted)
        self._add_wanted_bonds(builder, wanted, donors)
        self._drop_unassigned(builder, template, donors)

    def _include_assigned(self, builder, template, donors):
        """
        Puts every residue the assignment uses into the site before any bond
        is made.

        A residue remove_residue took out is remembered as excluded, and a
        bond to a residue the site excludes is refused by the extraction
        rather than quietly putting it back, so asking for it by name is
        what clears the way. It is also how a residue the template holds
        without coordinating anything gets in at all.

        :param builder:
            The builder holding the site.
        :param template:
            The template being matched.
        :param donors:
            The candidate of each assigned residue.
        """

        members = set(core.active_site_residues(builder.binding_modes))
        included = []

        for res_index, candidate in sorted(donors.items()):
            if res_index in members:
                continue
            builder.include_residue(candidate['resid'],
                                    chain=candidate['chain'])
            included.append(candidate['label'])

        if included:
            self.ostream.print_info(
                'Included ' + ', '.join(included) +
                f', which {template["name"]} is made of.')
            self.ostream.flush()

    def _remove_unwanted_bonds(self, builder, wanted):
        """
        Takes out every metal bond the assignment does not ask for.

        :param builder:
            The builder holding the site.
        :param wanted:
            The metal centers each residue is to coordinate, by residue
            index.
        """

        unwanted = []

        modes = builder.binding_modes
        metal_res = {
            entry['index']: entry['res_index'] for entry in modes['metals']
        }

        for ligand in modes['ligands']:
            for metal in ligand['metals']:
                if metal_res[metal] in wanted.get(ligand['res_index'], set()):
                    continue
                unwanted.append({
                    'resid': ligand['residue'][3:],
                    'chain': str(ligand['chain']),
                    'atom': ligand['atom'],
                    'metal': metal_res[metal],
                    'label': f'{ligand["residue"]} {ligand["atom"]}',
                })

        for entry in unwanted:
            metal = self._metal_entry(builder, entry['metal'])
            builder.remove_metal_bond(entry['resid'],
                                      metal=metal['index'],
                                      atom=entry['atom'],
                                      chain=entry['chain'])

        if unwanted:
            self.ostream.print_info(
                'Unbound ' + ', '.join(entry['label'] for entry in unwanted) +
                ', which the template does not coordinate.')
            self.ostream.flush()

    def _add_wanted_bonds(self, builder, wanted, donors):
        """
        Makes every metal bond the assignment asks for and the site does not
        already have.

        :param builder:
            The builder holding the site.
        :param wanted:
            The metal centers each residue is to coordinate, by residue
            index.
        :param donors:
            The candidate of each of those residues, which names the atom the
            bond is made through.
        """

        added = []

        for res_index in sorted(wanted):
            for metal_res_index in sorted(wanted[res_index]):
                if self._bonded(builder, res_index, metal_res_index):
                    continue

                candidate = donors[res_index]
                found = candidate['donor'].get(metal_res_index)

                metal = self._metal_entry(builder, metal_res_index)
                builder.add_metal_bond(candidate['resid'],
                                       metal['index'],
                                       atom=None if found is None else found[1],
                                       chain=candidate['chain'])
                added.append(f'{candidate["label"]} to '
                             f'{metal["element"]} {metal["index"]}')

        if added:
            self.ostream.print_info('Bonded ' + ', '.join(added) +
                                    ', which the template coordinates.')
            self.ostream.flush()

    @staticmethod
    def _bonded(builder, res_index, metal_res_index):
        """
        Whether one residue coordinates one metal center as things stand.

        :param builder:
            The builder holding the site.
        :param res_index:
            The residue.
        :param metal_res_index:
            The residue of the metal center.

        :return:
            True when the two are bonded.
        """

        modes = builder.binding_modes
        metal = [
            entry['index'] for entry in modes['metals']
            if entry['res_index'] == metal_res_index
        ]

        if not metal:
            return False

        return any(ligand['res_index'] == res_index and metal[0]
                   in ligand['metals'] for ligand in modes['ligands'])

    def _drop_unassigned(self, builder, template, donors):
        """
        Drops the residues the assignment had no use for.

        Last of the edits, so that a residue on its way out is never the one
        a metal center is still waiting for, and remove_residue's refusal to
        orphan a center cannot be tripped by the order things are done in.

        :param builder:
            The builder holding the site.
        :param template:
            The template being matched.
        :param donors:
            The candidate of each assigned residue.
        """

        residues = list(builder.enzyme_topology.residues())
        dropped = []

        for res_index in sorted(
                set(core.active_site_residues(builder.binding_modes)) -
                set(donors)):
            residue = residues[res_index]
            builder.remove_residue(str(residue.id),
                                   chain=str(residue.chain.id))
            dropped.append(core.residue_label(residue))

        if dropped:
            self.ostream.print_info(
                'Dropped ' + ', '.join(dropped) +
                f', which {template["name"]} is not made of.')
            self.ostream.flush()

    @staticmethod
    def _metal_res_index(builder, described, metal):
        """
        The residue index of one of the site's metal centers.

        :param builder:
            The builder holding the site.
        :param described:
            The described active site.
        :param metal:
            The index of the metal in that site.

        :return:
            The index of its residue in the topology.
        """

        index = described['atom_map'][metal]

        for entry in builder.binding_modes['metals']:
            if entry['index'] == index:
                return entry['res_index']

        assert_msg_critical(
            False, 'MetalForceFieldManager: the active site holds a metal '
            f'center at atom {index} that its binding modes do not')

    @staticmethod
    def _metal_entry(builder, res_index):
        """
        The metal center of one residue, as the coordination now has it.

        Looked up again before every edit rather than kept: an edit
        reprotonates the structure, which renumbers its atoms, so the atom
        index an edit takes is only good until the next one.

        :param builder:
            The builder holding the site.
        :param res_index:
            The index of the metal's residue.

        :return:
            The metal entry of the binding modes.
        """

        for entry in builder.binding_modes['metals']:
            if entry['res_index'] == res_index:
                return entry

        assert_msg_critical(
            False, 'MetalForceFieldManager: the structure no longer holds a '
            f'metal center in residue {res_index}')

    def _best_heavy_mapping(self, template, described, match_h_count=True):
        """
        Solves which of the site's atoms is which of the template's.

        Every way the two line up at the coarse level is taken down to the
        atoms and superimposed, and the one that fits best is the answer --
        which is the same thing compare_active_site measures a template by,
        run here to decide what to edit rather than what to report.

        :param template:
            The template to match.
        :param described:
            The described active site.
        :param match_h_count:
            Whether an atom has to carry as many hydrogens as the one it maps
            onto, at both levels. Off for the pass that runs before the
            protonation is put right, where no residue that differs in it
            would map at all.

        :return:
            The mapping over the heavy atoms, or None when there is none.
        """

        maps = []

        for coarse_mapping in self._coarse_mappings(
                template, described, match_protonation=match_h_count):
            maps.extend(
                self._heavy_atom_maps(template,
                                      described,
                                      coarse_mapping,
                                      match_h_count=match_h_count))

        if not maps:
            return None

        coordinates = described['molecule'].get_coordinates_in_angstrom()
        heavy_map, _, _ = self._best_heavy_map(template, maps, coordinates)

        return heavy_map

    @staticmethod
    def _hydrogen_count(described, index):
        """
        How many hydrogens one atom of a site carries.

        The capping hydrogen is one of them, on both sides of a comparison
        alike, so it cancels rather than having to be told apart here.

        :param described:
            A described active site or a template.
        :param index:
            The atom.

        :return:
            The count.
        """

        labels = described['molecule'].get_labels()

        return sum(1 for other in described['fine_topology'].neighbors(index)
                   if labels[other] == 'H')

    def _shoehorn_protonation(self, builder, template):
        """
        Protonates every residue of the site the way the template has it.

        Run before the denticity: an atom mapping that has to agree on the
        hydrogens cannot be solved while they differ, so this one is solved
        without them and is the last thing that has to be.

        :param builder:
            The builder holding the site.
        :param template:
            The template to match.

        :return:
            What stood in the way, or None when nothing did.
        """

        described = self._described_site(builder)
        heavy_map = self._best_heavy_mapping(template,
                                             described,
                                             match_h_count=False)

        if heavy_map is None:
            return ('the atoms of the site do not map onto those of '
                    f'{template["name"]}')

        changes = self._protonation_changes(builder, template, described,
                                            heavy_map)

        if isinstance(changes, str):
            return changes

        for change in changes:
            builder.update_protonation_state(change['resid'],
                                             change['variant'],
                                             chain=change['chain'])

        if changes:
            self.ostream.print_info(
                'Set ' + ', '.join(f'{change["label"]} to {change["variant"]}'
                                   for change in changes) +
                f', which is how {template["name"]} is protonated.')
            self.ostream.flush()

        return None

    def _protonation_changes(self, builder, template, described, heavy_map):
        """
        Works out which residues are protonated unlike the template.

        :param builder:
            The builder holding the site.
        :param template:
            The template to match.
        :param described:
            The described active site.
        :param heavy_map:
            The mapping from template index to site index.

        :return:
            What to hand to update_protonation_state, or what stood in the
            way as a string.
        """

        topology = builder.enzyme_topology
        atoms = list(topology.atoms())
        residues = list(topology.residues())
        variants = builder.binding_modes['variants']
        metals = set(described['metal_indices'])

        per_residue = {}

        for first, second in heavy_map.items():
            if second in metals:
                continue
            atom = atoms[described['atom_map'][second]]
            per_residue.setdefault(atom.residue.index, []).append(
                (first, atom.name))

        changes = []

        for res_index, pairs in sorted(per_residue.items()):
            residue = residues[res_index]
            current = variants.get(res_index, residue.name)

            wanted = {
                name: self._hydrogen_count(template, first)
                for first, name in pairs
            }
            delta = sum(wanted.values()) - sum(
                self._hydrogen_count(described, heavy_map[first])
                for first, _ in pairs)

            if delta == 0 and not self._tautomer_differs(
                    template, described, heavy_map, pairs):
                continue

            variant = self._target_variant(residue, current, delta, wanted)

            if variant is None:
                return (f'{core.residue_label(residue)} would have to be '
                        'protonated the way '
                        f'{template["name"]} has it, and there is no variant '
                        'of it that is')

            if variant == current:
                continue

            changes.append({
                'resid': str(residue.id),
                'chain': str(residue.chain.id),
                'label': core.residue_label(residue),
                'variant': variant,
            })

        return changes

    def _tautomer_differs(self, template, described, heavy_map, pairs):
        """
        Whether a residue carries its hydrogens on other atoms than the
        template does, with the same number of them.

        The histidines are what this is about: HID and HIE hold one hydrogen
        each and differ only in which ring nitrogen holds it.

        :param template:
            The template to match.
        :param described:
            The described active site.
        :param heavy_map:
            The mapping from template index to site index.
        :param pairs:
            The atoms of the residue, as template index and site atom name.

        :return:
            True when some atom of it disagrees.
        """

        return any(
            self._hydrogen_count(template, first) != self._hydrogen_count(
                described, heavy_map[first]) for first, _ in pairs)

    @staticmethod
    def _target_variant(residue, current, delta, wanted):
        """
        The protonation variant that gives a residue the template's
        hydrogens.

        A proton is a charge, so what the variant has to be is looked up by
        what the charge has to become: VARIANT_CHARGES is the table the
        active site charge is counted with and known_variants is what OpenMM
        will build. The histidines are the one family where two variants
        share a charge, and there the ring nitrogen that carries the hydrogen
        is what tells HID from HIE -- the same convention the automatic
        choice in _histidine_variant is made on.

        :param residue:
            The residue.
        :param current:
            The variant it is protonated as now.
        :param delta:
            How many hydrogens the template has on it more than the site
            does.
        :param wanted:
            How many hydrogens the template carries on each of its atoms, by
            atom name.

        :return:
            The variant, or None when no variant of the residue is it.
        """

        legal = core.known_variants(residue.name)
        charge = core.VARIANT_CHARGES.get(current)

        if not legal:
            # OpenMM builds it one way only, so as long as it carries as many
            # hydrogens as the template it is already protonated the only way
            # it can be. The relaxed mapping is free to swap two symmetric
            # heavy atoms, which is what puts one here with delta zero.
            return current if delta == 0 else None

        if charge is None:
            return None

        found = [
            variant for variant in legal
            if core.VARIANT_CHARGES.get(variant) == charge + delta
        ]

        if len(found) == 1:
            return found[0]

        if not found:
            return None

        # HID carries its hydrogen on ND1 and HIE on NE2
        tautomers = {'ND1': 'HID', 'NE2': 'HIE'}
        named = {
            tautomers[name]
            for name, count in wanted.items() if count and name in tautomers
        }
        named &= set(found)

        return named.pop() if len(named) == 1 else None

    def _shoehorn_denticity(self, builder, template):
        """
        Bonds every metal center to exactly the atoms the template bonds it
        to.

        How many atoms of a residue reach a metal is a distance cutoff on an
        unrelaxed structure rather than chemistry, which is why the matching
        refuses to have an opinion about it. What matched on those terms is
        built on the template's terms, and this is where the site is put onto
        them -- through the builder's own edits, so that the site itself says
        what it is bonded like rather than a force field built from it.

        :param builder:
            The builder holding the site.
        :param template:
            The template to match.

        :return:
            What stood in the way, or None when nothing did.
        """

        described = self._described_site(builder)
        heavy_map = self._best_heavy_mapping(template, described)

        if heavy_map is None:
            return ('the atoms of the site do not map onto those of '
                    f'{template["name"]} even once it is protonated like it')

        changes = self._denticity_changes(builder, template, described,
                                          heavy_map)

        for change in changes['added']:
            entry = self._metal_entry(builder, change['metal'])
            builder.add_metal_bond(change['resid'],
                                   entry['index'],
                                   atom=change['atom'],
                                   chain=change['chain'])

        for change in changes['removed']:
            entry = self._metal_entry(builder, change['metal'])
            builder.remove_metal_bond(change['resid'],
                                      metal=entry['index'],
                                      atom=change['atom'],
                                      chain=change['chain'])

        if changes['added'] or changes['removed']:
            self.ostream.print_info(
                f'Wired the metal centers as {template["name"]} wires its '
                f'own: {len(changes["added"])} metal bond(s) added, '
                f'{len(changes["removed"])} removed.')
            self.ostream.flush()

        return None

    def _denticity_changes(self, builder, template, described, heavy_map):
        """
        Works out which metal-ligand bonds the site makes and the template
        does not, and the other way round.

        Everything is named by residue id, atom name and the residue of the
        metal, none of which a rebuild disturbs, since applying the first of
        these changes renumbers the atoms the rest were found on.

        :param builder:
            The builder holding the site.
        :param template:
            The template to match.
        :param described:
            The described active site.
        :param heavy_map:
            The mapping from template index to site index.

        :return:
            What to add and what to remove.
        """

        metals = set(described['metal_indices'])
        bonds, _ = self._metal_keys(template)

        wanted = {
            frozenset((heavy_map[first], heavy_map[second]))
            for first, second in bonds
        }
        current = {
            frozenset((first, second))
            for first, second in core.connectivity_bonds(
                described['connectivity_matrix']) if metals & {first, second}
        }

        changes = {'added': [], 'removed': []}

        for kind, pairs in (('added', wanted - current), ('removed',
                                                          current - wanted)):
            for pair in sorted(pairs, key=sorted):
                changes[kind].append(
                    self._bond_record(builder, described, pair))

        return changes

    def _bond_record(self, builder, described, pair):
        """
        Names one metal-ligand bond the way an edit method takes it.

        :param builder:
            The builder holding the site.
        :param described:
            The described active site.
        :param pair:
            The two atoms of the bond, as active site indices.

        :return:
            The residue, its chain, the donor atom and the metal's residue.
        """

        metals = set(described['metal_indices'])
        metal, donor = sorted(pair, key=lambda index: index not in metals)

        atoms = list(builder.enzyme_topology.atoms())
        atom = atoms[described['atom_map'][donor]]

        return {
            'resid': str(atom.residue.id),
            'chain': str(atom.residue.chain.id),
            'atom': atom.name,
            'metal': self._metal_res_index(builder, described, metal),
            'label': f'{core.residue_label(atom.residue)} {atom.name}',
        }

    # ------------------------------------------------------------------
    # verdicts
    #
    # What the comparison is judged by, and the crude relaxation it is
    # judged after.
    # ------------------------------------------------------------------

    def _selection_criteria(self):
        """
        Resolves the criteria a template is held to, as selection_criteria
        holds them: a name in SELECTION_CRITERIA or a set of thresholds.

        A set that is one of the named ones is reported under its name, so
        that handing out SELECTION_CRITERIA['tight'] and naming 'tight' read
        the same way.

        :return:
            The tuple of the name to report it under and the set itself.
        """

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

        for name, known in self.SELECTION_CRITERIA.items():
            if criteria == known:
                return name, criteria

        return 'custom', criteria

    @staticmethod
    def _maps_every_atom(comparison, entry):
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

    def _selection_verdict(self, comparison, entry, criteria):
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
            violation = self._ic_violation(found['ic_rmsd'], thresholds)
            if violation is not None:
                return f'{region} {violation}'

        return None

    def _selection_score(self, entry):
        """
        Returns what several templates that all pass are ranked on.

        :param entry:
            What compare_active_site measured for the template.

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
        active_site = query['active_site']

        fit_kwargs = builder.fit_settings()
        if self.mm_fallback_literature_bonds and (
                fit_kwargs['metal_bond_equilibria'] is None):
            fit_kwargs['metal_bond_equilibria'] = core.LITERATURE_METAL_BONDS

        forcefield = core.build_forcefield(active_site,
                                           comm=MPI.COMM_SELF,
                                           ostream=self.ostream,
                                           **fit_kwargs)

        # this geometry is a way of comparing, not a result of a run, so
        # nothing about it is written to a folder
        return core.mm_optimize_active_site(
            active_site,
            forcefield,
            constrain_metals=builder.mm_constrain_metals,
            constrain_capping_hydrogens=builder.constrain_capping_hydrogens,
            max_iterations=builder.mm_max_iterations,
            bond_change_warning=builder.mm_bond_change_warning,
            ostream=self.ostream)

    # ------------------------------------------------------------------
    # matching
    # ------------------------------------------------------------------

    def _coarse_mappings(self, template, query, match_protonation=True):
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

    def _heavy_atom_maps(self,
                         template,
                         query,
                         coarse_mapping,
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

    # ------------------------------------------------------------------
    # measurement
    # ------------------------------------------------------------------

    def _rmsd_indices(self, template, region, heavy_only=None):
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
            heavy_only = self.rmsd_heavy_atoms_only

        indices = self._region_indices(template, region)

        if not heavy_only:
            return indices

        labels = template['molecule'].get_labels()

        return [index for index in indices if labels[index] != 'H']

    def _region_indices(self, template, region):
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
                    graph, metal, cutoff=self.metal_shell_bonds))

        assert_msg_critical(
            len(shell) > len(template['metal_indices']),
            'MetalForceFieldManager: the metal centers of template '
            f'{template["name"]} have nothing bonded to them within '
            f'{self.metal_shell_bonds} bond(s)')

        return sorted(shell)

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
            over its heavy atoms, and the internal coordinate deviations.
        """

        if heavy_only is None:
            heavy_only = self.rmsd_heavy_atoms_only

        reference = template['molecule'].get_coordinates_in_angstrom()
        indices = self._region_indices(template, region)
        labels = template['molecule'].get_labels()
        heavy = [index for index in indices if labels[index] != 'H']

        order = [mapping[index] for index in range(len(reference))]
        moved = coordinates[order]

        rmsd, _, _ = svd_superimpose(moved[indices], reference[indices])
        heavy_rmsd, _, _ = svd_superimpose(moved[heavy], reference[heavy])

        ic_rmsd = self._ic_rmsd(template, mapping, coordinates, region,
                                heavy_only)

        return {
            # _rmsd_indices is the region less the hydrogens when they are
            # being left out, which is exactly the two lists already in hand
            'atoms': len(heavy) if heavy_only else len(indices),
            'rmsd': rmsd,
            'rmsd_heavy': heavy_rmsd,
            'ic_rmsd': ic_rmsd,
        }

    def _ic_rmsd(self, template, mapping, coordinates, region, heavy_only=None):
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

        indices = self._rmsd_indices(template, region, heavy_only)
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

    def _ic_violation(self, ic_rmsd, thresholds):
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

        for name, unit in self.IC_TYPES.items():
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
        return core.get_metal_keys(template['forcefield'],
                                   {'metal_indices': template['metal_indices']})

    # ------------------------------------------------------------------
    # transfer
    # ------------------------------------------------------------------

    def _template_connectivity(self, template, mapping, active_site):
        """
        Wires a site's metal center exactly as the template wires its own.

        How many atoms of a residue reach a metal is a distance cutoff on an
        unrelaxed structure rather than chemistry, which is why the matching
        refuses to have an opinion about it: a carboxylate gripping a metal
        with one oxygen and one gripping it with two are the same residue on
        the same metal. What matches on those terms, though, has to be built
        on the template's terms as well. A template bond the site does not
        make has no atoms to land on, and a bond the site makes alone has no
        fitted parameters to land on it - it would keep the seeded guess and
        say nothing about it.

        So the template decides the coordination: its metal bonds are added
        where the site lacks them and the site's own are dropped where the
        template does not make them. Only bonds touching a metal are touched;
        the residues are wired as the structure has them.

        The active site is not modified; a new dictionary is returned, with
        the topologies rebuilt so they agree with the connectivity. When the
        two already agree the site is handed back as it came.

        :param template:
            The template that matched.
        :param mapping:
            The mapping from template index to active site index.
        :param active_site:
            The active site of the query.

        :return:
            The active site to build on, and what was added and removed.
        """

        metals = set(active_site['metal_indices'])

        assert_msg_critical(
            {mapping[index]
             for index in template['metal_indices']} == metals,
            'MetalForceFieldManager: the template maps its metal centers onto '
            'atoms of the site that are not metal centers')

        bonds, _ = self._metal_keys(template)
        wanted = {
            frozenset((mapping[first], mapping[second]))
            for first, second in bonds
        }

        matrix = np.array(active_site['connectivity_matrix'])
        changes = {'added': [], 'removed': []}

        for pair in wanted:
            first, second = sorted(pair)
            if not matrix[first, second]:
                matrix[first, second] = 1
                matrix[second, first] = 1
                changes['added'].append((first, second))

        for first, second in core.connectivity_bonds(matrix):
            if not ({first, second} & metals):
                continue
            if frozenset((first, second)) in wanted:
                continue
            matrix[first, second] = 0
            matrix[second, first] = 0
            changes['removed'].append((first, second))

        if not changes['added'] and not changes['removed']:
            return active_site, changes

        active_site = {
            **active_site,
            'connectivity_matrix': matrix,
        }

        return self._describe(active_site, core.connectivity_bonds(matrix)), changes

    def _print_forced_bonds(self, template, active_site, changes):
        """
        Reports the metal bonds the template decided against the site.

        A bond added over a long contact carries the template's equilibrium
        and will pull the two atoms together at the first minimization, and a
        short contact dropped is a pair left to the nonbonded terms alone.
        Both are worth seeing here rather than in a geometry afterwards, so
        anything on the wrong side of the coordination cutoff is a warning.

        :param template:
            The template that decided.
        :param active_site:
            The active site the distances are read off.
        :param changes:
            What _template_connectivity added and removed.
        """

        if not changes['added'] and not changes['removed']:
            return

        coordinates = active_site['molecule'].get_coordinates_in_angstrom()
        labels = active_site['molecule'].get_labels()
        cutoff = self.builder.metal_bond_cutoff

        self.ostream.print_info(
            f'The coordination of the site differs from {template["name"]}; '
            f'forcing it onto the template: {len(changes["added"])} metal '
            f'bond(s) added, {len(changes["removed"])} removed.')

        for kind, pairs in (('adding', changes['added']), ('removing',
                                                           changes['removed'])):
            for first, second in pairs:
                distance = np.linalg.norm(coordinates[first] -
                                          coordinates[second])
                line = (f'  {kind} {labels[first]}{first}-'
                        f'{labels[second]}{second}, {distance:.2f} A apart '
                        'in the site')

                far = (kind == 'adding' and distance > cutoff)
                near = (kind == 'removing' and distance <= cutoff)

                if far or near:
                    self.ostream.print_warning(line.strip())
                else:
                    self.ostream.print_info(line)

        self.ostream.flush()

    def _build_ff_from_template(self, template, mapping, active_site):
        """
        Builds a force field for an active site out of a template.

        Only what a builder run pays QM for is taken from the template: the
        fitted metal bonds and angles, and the charges. Everything else is
        built for the site in front of us, so the atom types and the bonded
        terms of the residues come from the structure rather than from
        somewhere else.

        The coordination itself is the template's, not the site's: the metal
        bonds are forced onto the template's by _template_connectivity before
        anything is built, so a residue the structure holds bidentate is built
        monodentate where the template is monodentate and the other way
        round. Everything the template was fitted for then has atoms to land
        on, and nothing is left carrying a seeded guess.

        :param template:
            The template that matched.
        :param mapping:
            The mapping from template index to active site index.
        :param active_site:
            The active site of the query, whose connectivity the force field
            is built on once the template has decided the metal bonds.

        :return:
            The force field generator, and the active site it was built on.
        """

        active_site, changes = self._template_connectivity(
            template, mapping, active_site)
        self._print_forced_bonds(template, active_site, changes)

        template_ff = template['forcefield']
        charges = np.zeros(active_site['molecule'].number_of_atoms())

        for template_index, site_index in mapping.items():
            charges[site_index] = template['charges'][template_index]

        total = float(np.sum(charges))
        expected = int(active_site['molecule'].get_charge())
        if abs(total - expected) > 1.0e-3:
            self.ostream.print_warning(
                f'The transferred charges sum to {total:+.3f}, but the active '
                f'site charge is {expected:+d}')

        # build_forcefield writes the charges onto the atoms and redistributes
        # the caps; the metal terms it seeds here are overwritten below.
        # Every rank builds its own: the work is cheap, and a communicator of
        # one keeps it clear of the collectives a shared one would invite.
        forcefield = core.build_forcefield(active_site,
                                           partial_charges=charges,
                                           comm=MPI.COMM_SELF,
                                           ostream=self.ostream,
                                           **self.builder.fit_settings())

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

        return forcefield, changes, active_site

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

    # ------------------------------------------------------------------
    # printing
    # ------------------------------------------------------------------

    def _print_shoehorn_header(self, builder, template, max_include_radius):
        """
        Says what a shoehorning is starting from.

        :param builder:
            The builder holding the site.
        :param template:
            The template it is being edited onto.
        :param max_include_radius:
            How far out from a metal center a residue may be picked up from.
        """

        name = template['name']

        self.ostream.print_blank()
        self.ostream.print_header(f'Shoehorning the site into {name}')
        self.ostream.print_header((26 + len(name)) * '-')
        self.ostream.print_blank()
        self.ostream.print_header(
            core._param('template residues',
                        len(self._residue_nodes(template['coarse_topology']))))
        self.ostream.print_header(
            core._param('site residues',
                        len(builder.active_site['residues'])))
        self.ostream.print_header(
            core._param('search radius', f'{max_include_radius:.1f} A'))
        self.ostream.print_blank()
        self.ostream.print_info(
            f'Site holds: {", ".join(builder.active_site["residues"])}')
        self.ostream.print_blank()
        self.ostream.flush()

    def _print_shoehorn_summary(self, builder, template):
        """
        Says what the site was made into.

        The coordination is what a shoehorning is for, so it is printed the
        way the builder prints it rather than left to be asked for: which
        residue ended up on which metal center, how far out, and what it is
        protonated as.

        :param builder:
            The builder holding the site.
        :param template:
            The template it was edited onto.
        """

        modes = builder.binding_modes
        residues = list(builder.enzyme_topology.residues())

        self.ostream.print_blank()
        self.ostream.print_info(
            f'The site is now built the way {template["name"]} is built.')
        self.ostream.print_blank()

        for metal in modes['metals']:
            bonds = []
            for ligand in modes['ligands']:
                for index, distance in zip(ligand['metals'],
                                           ligand['distances']):
                    if index != metal['index']:
                        continue
                    bonds.append(f'{ligand["residue"]} {ligand["atom"]} '
                                 f'{distance:.2f} A')
            self.ostream.print_info(
                f'  {metal["element"]} {metal["index"]}: ' +
                (', '.join(sorted(bonds)) if bonds else 'nothing'))

        protonation = ', '.join(
            f'{residues[res_index].name}{residues[res_index].id} {variant}'
            for res_index, variant in sorted(modes['variants'].items())
            if res_index in set(core.active_site_residues(modes)))

        self.ostream.print_info(f'  protonation: {protonation}')
        self.ostream.print_blank()
        self.ostream.print_info(
            'Call compare_active_site() to measure it, then '
            'build_ff_from_template.')
        self.ostream.print_blank()
        self.ostream.flush()

    def _print_template(self, template):
        """
        Prints what one template holds.

        :param template:
            The template.
        """

        bonds, angles = self._metal_keys(template)
        labels = template['molecule'].get_labels()
        metals = ', '.join(labels[index] for index in template['metal_indices'])

        self.ostream.print_blank()
        self.ostream.print_header(f'Template {template["name"]}')
        self.ostream.print_header((9 + len(template['name'])) * '-')
        self.ostream.print_header(
            core._param('geometry', template['geometry_kind']))
        self.ostream.print_header(
            core._param('atoms', template['molecule'].number_of_atoms()))
        self.ostream.print_header(core._param('metal centers', metals))
        self.ostream.print_header(
            core._param('capping hydrogens', len(template['cap_indices'])))
        self.ostream.print_header(core._param('metal bonds', len(bonds)))
        self.ostream.print_header(core._param('metal angles', len(angles)))
        self.ostream.print_header(
            core._param('total charge',
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
                name[:24], template['molecule'].number_of_atoms(),
                len(template['metal_indices']), template['geometry_kind'])
            self.ostream.print_header(valstr)

        self.ostream.print_blank()
        self.ostream.flush()

    def _print_comparison(self, results):
        """
        Prints everything compare_active_site measured.

        One table of numbers and one of verdicts per template that could be
        measured, and a closing summary ranking the templates by the region
        that is configured, so the closest one is visible without reading
        every table.

        :param results:
            The last comparison, from compare_active_site.
        """

        active_site = results['active_site']
        labels = active_site['molecule'].get_labels()
        metals = ', '.join(labels[index]
                           for index in active_site['metal_indices'])

        self.ostream.print_blank()
        self.ostream.print_header('Comparison against every template')
        self.ostream.print_header(33 * '-')
        self.ostream.print_header(
            core._param('source',
                        Path(results['source']).name))
        self.ostream.print_header(
            core._param('active site atoms',
                        active_site['molecule'].number_of_atoms()))
        self.ostream.print_header(core._param('metal centers', metals))
        self.ostream.print_header(core._param('geometry', results['geometry']))
        self.ostream.print_header(
            core._param(
                'measured over',
                'all atoms' if results['include_hydrogens'] else 'heavy atoms'))
        self.ostream.print_header(
            core._param('templates', len(results['templates'])))
        self.ostream.print_blank()
        self.ostream.print_info(
            f'Residues: {", ".join(active_site["residues"])}')

        self._print_spec('the structure holds', results['active_site'])

        for name, entry in results['templates'].items():
            self._print_template_comparison(name, entry)

        self._print_comparison_summary(results)

    def _print_spec(self, title, described):
        """
        Prints what a site is made of: which residues it holds and which of
        them coordinate which metal.

        The residues are named by their formula, which is for reading; two
        sites are compared on the keys behind them.

        :param title:
            What the block is describing.
        :param described:
            A described active site or a template.
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

        bridging = [
            node for node in self._residue_nodes(coarse)
            if coarse.degree(node) > 1
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
            What compare_active_site measured for it.
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
        Ranks the templates on what a selection is decided by.

        :param results:
            The last comparison, from compare_active_site.
        """

        region, ic_type, measure = self.SELECTION_RANKED_ON
        heavy = not results['include_hydrogens']

        def rmsd(entry):
            found = entry['regions'].get(region)
            if found is None:
                return None
            return found['rmsd_heavy'] if heavy else found['rmsd']

        order = sorted(results['templates'].items(),
                       key=lambda item: self._selection_score(item[1]))

        self.ostream.print_blank()
        title = f'Ranked on the {region} {ic_type} {measure}'
        self.ostream.print_header(title)
        self.ostream.print_header(len(title) * '-')

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
            score = self._selection_score(entry)
            self.ostream.print_header(
                row.format(name[:24], entry['status'], f'{score:.3f}',
                           f'{found:.3f}',
                           f'{entry["metal_bonds"]["deviation"]:.3f}'))

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
            The last comparison, from compare_active_site.
        :param decision:
            The decision, as _select_template makes it.
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
                core._param(region, limits, value_width=44))

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

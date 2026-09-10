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
import math
import sys

from ..veloxchemlib import mpi_master
from ..molecule import Molecule
from ..outputstream import OutputStream
from .metalsiteffbuilder import MetalSiteForceFieldBuilder
from . import core
from . import printing
from . import matching
from . import templates
from . import shoehorn
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
    # argument of load_template_from_folder decides. Defined in templates.py,
    # beside the loading that reads it.
    GEOMETRY_KINDS = templates.GEOMETRY_KINDS

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
                    'rms': 0.5,
                    'max': 1.0
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
                    'rms': 1.0,
                    'max': 1.5
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
        geometry, kind = templates.load_geometry(folder, fallback)
        forcefield.molecule = geometry

        template = templates.build(
            name,
            forcefield,
            geometry,
            kind,
            folder,
            metal_elements=self.builder.metal_elements,
            ostream=self.ostream)

        if name in self.templates:
            self.ostream.print_warning(
                f'A template named {name} is already loaded; replacing it '
                f'with the one from {folder}')

        self.templates[name] = template

        bonds, angles = matching.metal_keys(template)
        printing.print_template(template, bonds, angles, ostream=self.ostream)
        printing.print_templates(self.templates, ostream=self.ostream)

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
        template dictionary (built by templates.build) has neither -- it
        carries forcefield.bonds instead, which is the same source every
        other template computation (matching.describe, matching.metal_keys)
        reads its bonds from.

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
        described = shoehorn.described_site(builder)

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
            coarse = matching.coarse_mappings(template, described)

            if not coarse:
                # what it holds instead is printed from the template
                entry['status'] = 'spec'
                findings[name] = entry
                continue

            maps = []
            for coarse_mapping in coarse:
                maps.extend(
                    matching.heavy_atom_maps(template,
                                             described,
                                             coarse_mapping,
                                             max_mappings=self.max_mappings,
                                             ostream=self.ostream))

            if not maps:
                # what it holds instead is printed from the template
                entry['status'] = 'spec'
                findings[name] = entry
                continue

            entry['n_coarse_mappings'] = len(coarse)
            entry['n_mappings'] = len(maps)

            heavy_map, rot, trans = matching.best_heavy_map(
                template, maps, coordinates)
            mapping = matching.complete_hydrogens(template, described,
                                                  heavy_map, coordinates,
                                                  rot, trans)

            entry['mapping'] = mapping
            entry['metal_bonds'] = matching.metal_bond_summary(
                template, described, mapping, coordinates)

            for region in self.RMSD_REGIONS:
                entry['regions'][region] = matching.measure_region(
                    template,
                    mapping,
                    coordinates,
                    region,
                    heavy_only,
                    rmsd_heavy_atoms_only=self.rmsd_heavy_atoms_only,
                    metal_shell_bonds=self.metal_shell_bonds,
                    ostream=self.ostream)

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

        # the table this decision would print is printed by the one caller
        # that acts on it, build_ff_from_template
        return self._select_template(None)['name'] is not None

    def _print_comparison(self, results):
        """
        Gathers what the comparison tables are drawn from and prints them.

        The scores and the residue nodes are worked out here rather than in
        the printer, which computes nothing.

        :param results:
            The last comparison, from compare_active_site.
        """

        # a spec is printed for the structure, and for any template that
        # coordinates a different set of residues than it does
        specs = {None: self._spec_of(results['active_site'])}
        for name, entry in results['templates'].items():
            if entry['status'] == 'spec':
                specs[name] = self._spec_of(self.templates[name])

        scores = {
            name: self._selection_score(entry)
            for name, entry in results['templates'].items()
        }

        printing.print_comparison(results,
                                  specs,
                                  self.SELECTION_RANKED_ON,
                                  scores,
                                  ostream=self.ostream)

    def _spec_of(self, described):
        """
        Returns a described site paired with the residues that bridge it,
        which is what printing.print_spec is drawn from.

        :param described:
            A described active site or a template.

        :return:
            The pair.
        """

        return (described,
                matching.bridging_nodes(described['coarse_topology']))

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

        return decision

    def _print_decision(self, decision):
        """
        Prints how every template stood against the criteria, and which one
        was taken.

        Separate from _select_template because a decision is made twice --
        once to answer compare_active_site's boolean and once to act on in
        build_ff_from_template -- and printing the whole field twice for one
        comparison only invited the reader to look for what had changed
        between them. Nothing had.

        :param decision:
            The decision, as _select_template makes it.
        """

        printing.print_selection(self._comparison,
                                 decision,
                                 self.RMSD_REGIONS,
                                 self.IC_TYPES,
                                 self.SELECTION_RANKED_ON,
                                 ostream=self.ostream)

        if decision['name'] is None:
            self._print_no_selection(decision)

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

        self._print_decision(decision)

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

        # a run that fails leaves no record of one that succeeded
        self._shoehorned = None

        walked = shoehorn.run(self._active_site_builder,
                              self.templates[template],
                              max_include_radius,
                              max_mappings=self.max_mappings,
                              ostream=self.ostream)

        if not walked:
            return False

        self._shoehorned = template

        # the edits are exactly what makes the last comparison stale, so the
        # site is measured again here rather than by a caller who has to
        # remember to
        self.compare_active_site()

        return True

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
            violation = matching.ic_violation(found['ic_rmsd'], thresholds,
                                              ic_types=self.IC_TYPES)
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
            ostream=self.ostream,
            **builder.relax_settings())

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

        bonds, _ = matching.metal_keys(template)
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

        return matching.describe(active_site, core.connectivity_bonds(matrix)), changes

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

        bonds, angles = matching.metal_keys(template)

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

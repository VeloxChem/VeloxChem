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
Turning what a run left in a folder into a template.

A template is an active site with the two files it was loaded from added to
it, rather than a dictionary describing one: every matching step then takes
two dictionaries of the same shape, and a template is simply an active site
with no topology behind it.

The rules that live here: which geometry a folder is allowed to be built
from and how the two are told apart, that the bonding comes from the force
field's own bonds rather than a third file, that the caps and beta carbons
are read off the comments annotate_atoms wrote, and that a residue
coordinating no metal is discarded with its charge put back on what stays.
"""

import numpy as np

from ..molecule import Molecule
from ..errorhandler import assert_msg_critical
from . import util
from .util import Shell, on_master, param, print_section
from .matching import SiteMatcher
from .builder import ActiveSiteBuilder

# The geometry a run leaves behind, in the order it is looked for. Which of
# them a template is allowed to be built from is what the fallback argument
# of load_template_from_folder decides.
# The geometry a run leaves behind, in the order it is looked for. Which of
# them a template is allowed to be built from is what the fallback argument
# of load_template_from_folder decides.
GEOMETRY_KINDS = ('qm_opt', 'mm_opt')


class TemplateLoader(Shell):
    """
    Makes a template out of what a run left in its folder -- the force
    field and the geometry, described the way a site is described, with
    the residues that coordinate no metal discarded -- and puts a
    template's parameters onto a site that matched it: the metal bonds
    wired the way the template wires them, then the fitted metal terms and
    the charges transferred over the mapping.

    Every method runs on the master rank and its result is broadcast.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.
    """

    def _matcher(self):
        """
        The SiteMatcher a template is described with.

        :return:
            A SiteMatcher on this loader's communicator and stream.
        """

        return SiteMatcher(self.comm, self.ostream)

    def _sites(self):
        """
        The ActiveSiteBuilder a transferred force field is built with.

        :return:
            An ActiveSiteBuilder on this loader's communicator and stream.
        """

        return ActiveSiteBuilder(self.comm, self.ostream)

    @on_master
    def load_geometry(self, folder, fallback):
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

        opt_path = folder / util.GEOMETRY_FILE
        mm_path = folder / util.MM_GEOMETRY_FILE

        allowed = GEOMETRY_KINDS[:1 if fallback is None else 2]

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

    @on_master
    def _same_geometry(self, first, second):
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

    @on_master
    def build(self, name, forcefield, molecule, kind, folder, metal_elements):
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
        :param metal_elements:
            The elements treated as metal centers.

        :return:
            The template dictionary.
        """

        # the same check a run makes of its own force field, which is why it
        # is util's and not a second copy of it here. A template is an
        # active site with no topology behind it, so {'molecule': ...} is all
        # the adapter it needs.
        util._check_forcefield(forcefield, {'molecule': molecule}, folder)

        labels = list(molecule.get_labels())

        metal_indices = [
            index for index, label in enumerate(labels)
            if label in metal_elements
        ]

        assert_msg_critical(
            len(metal_indices) > 0,
            f'MetalForceFieldManager: the template of {folder} holds no metal '
            f'center. Recognized elements: {metal_elements}')

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
            for role in (util.BETA_CARBON_COMMENT, util.CAP_COMMENT)
        }

        beta_carbon_indices = marked[util.BETA_CARBON_COMMENT]
        cap_indices = marked[util.CAP_COMMENT]

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

        described = self._matcher().describe(template, forcefield.bonds.keys())

        return self._prune_unconnected_residues(described)

    @on_master
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
            node for node in self._matcher().residue_nodes(coarse)
            if coarse.degree(node) == 0
        ]

        if not unconnected:
            return template

        discard = sorted(atom for node in unconnected
                         for atom in coarse.nodes[node]['atoms'])

        pruned, shift = self._drop_atoms(template, discard)
        self._print_discarded_residues(template,
                                       unconnected,
                                       discard,
                                       shift)

        return self._matcher().describe(pruned, pruned['forcefield'].bonds.keys())

    @on_master
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
        charges, shift = self.compensate_charges(template['charges'][keep],
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

    @on_master
    def compensate_charges(self, charges, cap_indices):
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

    @on_master
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

    @on_master
    def print_template(self, template, bonds, angles):
        """
        Prints what one template holds.

        :param template:
            The template.
        :param bonds:
            Its metal bond keys.
        :param angles:
            Its metal angle keys.
        """

        labels = template['molecule'].get_labels()
        metals = ', '.join(labels[index] for index in template['metal_indices'])

        self.ostream.print_blank()
        print_section(f'Template {template["name"]}', self.ostream)
        self.ostream.print_header(param('geometry', template['geometry_kind']))
        self.ostream.print_header(param('atoms', template['molecule'].number_of_atoms()))
        self.ostream.print_header(param('metal centers', metals))
        self.ostream.print_header(
            param('capping hydrogens', len(template['cap_indices'])))
        self.ostream.print_header(param('metal bonds', len(bonds)))
        self.ostream.print_header(param('metal angles', len(angles)))
        self.ostream.print_header(
            param('total charge', f'{float(np.sum(template["charges"])):+.3f}'))
        self.ostream.print_blank()
        self.ostream.print_info(f'Loaded from {template["folder"]}')
        self.ostream.flush()

    @on_master
    def print_templates(self, templates):
        """
        Prints every template that is loaded.

        :param templates:
            The templates, by name.
        """

        self.ostream.print_header(f'Loaded templates ({len(templates)})')
        self.ostream.print_header(60 * '-')
        valstr = '{:>24} | {:>7} | {:>7} | {:>13}'.format('name', 'atoms', 'metals',
                                                          'geometry')
        self.ostream.print_header(valstr)
        self.ostream.print_header(60 * '-')

        for name, template in templates.items():
            valstr = '{:>24} | {:>7} | {:>7} | {:>13}'.format(
                name[:24], template['molecule'].number_of_atoms(),
                len(template['metal_indices']), template['geometry_kind'])
            self.ostream.print_header(valstr)

        self.ostream.print_blank()
        self.ostream.flush()

    # ------------------------------------------------------------------
    # the transfer
    # ------------------------------------------------------------------

    @on_master
    def template_connectivity(self, template, mapping, active_site):
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
            'TemplateLoader: the template maps its metal centers onto '
            'atoms of the site that are not metal centers')

        bonds, _ = self._matcher().metal_keys(template)
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

        for first, second in util.connectivity_bonds(matrix):
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

        described = self._matcher().describe(active_site,
                                             util.connectivity_bonds(matrix))

        return described, changes

    @on_master
    def print_forced_bonds(self, template, active_site, changes,
                           metal_bond_cutoff):
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
        cutoff = metal_bond_cutoff

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

    @on_master
    def build_forcefield_from_template(self, template, mapping, active_site,
                                       metal_bond_cutoff, **seed_settings):
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

        active_site, changes = self.template_connectivity(
            template, mapping, active_site)
        self.print_forced_bonds(template, active_site, changes,
                                metal_bond_cutoff)

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
        # the caps; the metal terms it seeds here are overwritten below
        forcefield = self._sites().build_forcefield(
            active_site, charges, **seed_settings)

        bonds, angles = self._matcher().metal_keys(template)

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

    @on_master
    def _map_key(self, key, mapping, table, kind):
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
            False, f'TemplateLoader: the template {kind} {key} maps '
            f'onto {mapped}, which the active site force field does not have')

    @on_master
    def _transferred(self, params, name):
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

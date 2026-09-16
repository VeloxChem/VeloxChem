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
The active site phase of the metal site force field builder: everything
between a structure file and a cluster with a seeded force field, which is
what the QM is then paid for.
"""

from mpi4py import MPI
from pathlib import Path
from copy import deepcopy
import numpy as np
import tempfile
import sys

from ..molecule import Molecule
from ..outputstream import OutputStream
from ..mmforcefieldgenerator import MMForceFieldGenerator
from ..errorhandler import assert_msg_critical
from .util import (
    Shell, on_master, param, print_param_list, SUPPORTED_METAL_ELEMENTS,
    DONOR_ELEMENTS, BIDENTATE_ASYMMETRY, SEEDED_FROM_REQUEST,
    SEEDED_FROM_TABLE, SEEDED_FROM_GEOMETRY, SEEDED_EQUILIBRIUM_LABELS,
    BRIDGING_RESIDUES, CARBOXYLATE_RESIDUES, BACKBONE_ATOM_NAMES,
    UNTRUNCATABLE_RESIDUES, VARIANT_CHARGES, BETA_CARBON_COMMENT, CAP_COMMENT,
    UNPARAMETERIZED_COMMENT, UFF_TYPE_COMMENT, residue_label, _site_index_map,
    connectivity_bonds, get_metal_keys, active_site_residues, check_variant,
    redistribute_cap_charges, constrained_indices)

try:
    import openmm as mm
    import openmm.app as mmapp
    import openmm.unit as mmunit
except ImportError:
    pass

try:
    from pdbfixer import PDBFixer
except ImportError:
    pass


class ActiveSiteBuilder(Shell):
    """
    Builds the active site of a metal-containing protein: finds the metal
    centers and what coordinates them, protonates the structure to match,
    truncates the cluster, relaxes it on a seeded force field and edits
    its coordination on request.

    Every method runs on the master rank and its result is broadcast: the
    protonation draws from Python's global random stream and the rest is
    cheap, so one rank does the work and every rank holds the same answer.
    A method takes what it uses and returns what it produces -- no
    intermediate is kept here -- and none of them modifies its arguments.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.
    """

    # ------------------------------------------------------------------
    # loading
    # ------------------------------------------------------------------

    @on_master
    def load_and_prepare_protein(self, structure, prepare):
        """
        Reads a structure and prepares it for a protein force field.

        Preparation adds missing heavy atoms so that a protein force field
        can match templates. Missing residues are deliberately not built.
        Necessary for building full enzymatic systems; can be skipped if the
        provided topology file is already correct.

        :param structure:
            The path to a .pdb, .cif or .pdbx file.
        :param prepare:
            Whether to run the preparation.

        :return:
            The tuple of the OpenMM topology and the positions as an (N, 3)
            numpy array in Angstrom.
        """

        assert_msg_critical('openmm' in sys.modules,
                            'load_and_prepare_protein: openmm is '
                            'required')

        if prepare:
            assert_msg_critical(
                'pdbfixer' in sys.modules or 'PDBFixer' in globals(),
                'prepare_protein: pdbfixer is require when preparing a protein')

        path = Path(structure)

        assert_msg_critical(path.is_file(), f'load_and_prepare_protein: {path} not '
                            'found')

        if path.suffix.lower() in ('.cif', '.pdbx', '.mmcif'):
            pdb = mmapp.PDBxFile(str(path))
        else:
            pdb = mmapp.PDBFile(str(path))

        positions = np.array(pdb.positions.value_in_unit(mmunit.angstrom))

        if not prepare:
            return pdb.topology, positions

        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / 'input.pdb'
            with path.open('w') as fh:
                mmapp.PDBFile.writeFile(pdb.topology,
                                        np.asarray(positions) * mmunit.angstrom,
                                        fh,
                                        keepIds=True)

            fixer = PDBFixer(filename=str(path))
            fixer.findMissingResidues()
            # do not build missing loops into a designed structure
            fixer.missingResidues = {}
            fixer.findMissingAtoms()
            fixer.addMissingAtoms()

        positions = np.array(fixer.positions.value_in_unit(mmunit.angstrom))

        return fixer.topology, positions

    # ------------------------------------------------------------------
    # detection
    # ------------------------------------------------------------------

    @on_master
    def derive_binding_modes(self,
                             topology,
                             positions,
                             request,
                             metal_elements,
                             metal_formal_charges,
                             metal_bond_cutoff,
                             report_cutoff,
                             coordinating_residues=None):
        """
        Derives the coordination topology of the metal centers from geometry.


        :param topology:
            The OpenMM topology.
        :param positions:
            The positions as an (N, 3) numpy array in Angstrom.
        :param request:
            The stored request -- manual metal bonds, protonation variants and
            the residue membership overrides. Replayed onto the detection.
        :param metal_elements:
            The elements treated as metal centers. Note that detecting one the
            builder is not validated for is a hard failure, so a zinc structure
            that also holds a calcium or magnesium ion has to narrow this to
            ('Zn',) by hand.
        :param metal_formal_charges:
            The formal charges assumed for the metal ions, by element.
        :param metal_bond_cutoff:
            The distance in Angstrom within which a donor atom is bonded to a
            metal.
        :param report_cutoff:
            The distance in Angstrom out to which a contact is reported
            without being made a bond.
        :param coordinating_residues:
            Residues that must be ligands whatever their distance, each given
            as a residue id ('130' or 130) or as a residue label ('ASP130').
            Merged into the request's own list.

        :return:
            The binding modes dictionary, freshly derived.
        """

        positions = np.asarray(positions)

        metals = []
        for atom in topology.atoms():
            if atom.element is None:
                continue
            if atom.element.symbol not in metal_elements:
                continue
            symbol = atom.element.symbol
            metals.append({
                'index': atom.index,
                'element': symbol,
                'chain': atom.residue.chain.id,
                'resid': atom.residue.id,
                'res_index': atom.residue.index,
                'formal_charge': metal_formal_charges.get(symbol, 2),
            })

        assert_msg_critical(
            len(metals) > 0, 'derive_binding_modes: no metal atom '
            f'found. Recognized elements: {metal_elements}')

        self._check_supported_metals(metals, 'derive_binding_modes')

        # Resolve which residues are forced to be ligands
        forced = self._resolve_residues(topology, coordinating_residues)
        forced |= set(request.get('coordinating_residues', []))

        atoms = list(topology.atoms())

        def position_of(index):
            return positions[index]

        notes = []
        # Collect all close-lying ligands
        ligands = self._collect_ligands(atoms, position_of, metals, notes,
                                        forced, metal_bond_cutoff,
                                        report_cutoff)

        # A bond decided by hand is not something the distances can be asked
        # about again, so the records are replayed onto every derivation. This
        # is what makes deriving twice give the same answer.
        records = request.get('manual_bonds', [])
        if records:
            self._apply_manual_bonds(ligands, records, metals, atoms, position_of, notes)
            self._assign_binding_modes(ligands, notes,
                                       self._manual_protected(records))

        binding_modes = {
            'metals': metals,
            'ligands': ligands,
            'variants': dict(request.get('variants', {})),
            'coordinating_residues': sorted(forced),
            'extra_residues': sorted(request.get('extra_residues', [])),
            'excluded_residues': sorted(request.get('excluded_residues', [])),
            'manual_bonds': deepcopy(records),
            'notes': notes,
        }

        return binding_modes

    @on_master
    def _collect_ligands(self, atoms, position_of, metals, notes, forced,
                         metal_bond_cutoff, report_cutoff):
        """
        Builds the classified ligand contact list of a set of candidate atoms.

        Shared by derive_binding_modes, which offers it every atom of the
        topology, and by update_binding_modes, which offers it the atoms of the
        truncated active site only, so that the two apply the same cutoffs and
        the same classification rules.

        :param atoms:
            The candidate atoms, as OpenMM atoms. Metals, non-donors and atoms
            beyond the secondary cutoff are skipped.
        :param position_of:
            The callable returning the position in Angstrom of an atom index of
            the topology. Called for the candidates and for the metals.
        :param metals:
            The metal entries of the binding modes.
        :param notes:
            The list of review notes. Appended to in place.
        :param forced:
            The residue indices that are ligands whatever their distance.
        :param metal_bond_cutoff:
            The distance in Angstrom within which a donor atom is bonded to a
            metal.
        :param report_cutoff:
            The distance in Angstrom out to which a contact is reported
            without being made a bond. A scan that stopped before the bonding
            cutoff would drop contacts that are bonds, so it is held to at
            least metal_bond_cutoff.

        :return:
            The list of ligand contacts, each carrying its mode.
        """

        report_cutoff = max(float(report_cutoff), float(metal_bond_cutoff))

        metal_indices = [metal['index'] for metal in metals]
        atoms = list(atoms)
        contacts = []

        # Find the contacting donor atom
        for atom in atoms:
            if atom.element is None or atom.index in metal_indices:
                continue
            if atom.element.symbol not in DONOR_ELEMENTS:
                continue

            position = position_of(atom.index)
            distances = {
                index: float(np.linalg.norm(position - position_of(index)))
                for index in metal_indices
            }
            closest = min(distances.values())
            if closest > report_cutoff:
                continue

            label = residue_label(atom.residue)

            if atom.name in BACKBONE_ATOM_NAMES:
                notes.append(
                    f'backbone atom {label} {atom.name} is {closest:.2f} A '
                    'from a metal; the truncation scheme currently is sidechain-only, '
                    'so it is not treated as a ligand')
                continue

            bonded_to = sorted([
                index for index, dist in distances.items()
                if dist <= metal_bond_cutoff
            ])

            contacts.append({
                'residue':
                label,
                'res_name':
                atom.residue.name,
                'res_index':
                atom.residue.index,
                'chain':
                atom.residue.chain.id,
                'atom':
                atom.name,
                'index':
                atom.index,
                'metals':
                bonded_to,
                'distances': [round(distances[i], 3) for i in bonded_to],
            })

            if not bonded_to:
                notes.append(
                    f'{label} {atom.name} at {closest:.2f} A is between the '
                    f'primary ({metal_bond_cutoff}) and secondary '
                    f'({report_cutoff}) cutoffs; review whether it '
                    'should be a ligand')

        self._force_ligands(atoms, position_of, metal_indices, contacts,
                            notes, forced, metal_bond_cutoff, report_cutoff)

        ligands = [contact for contact in contacts if contact['metals']]

        self._assign_binding_modes(ligands, notes)

        for metal in metals:
            n_ligands = sum(1 for ligand in ligands
                            if metal['index'] in ligand['metals'])
            if n_ligands < 3:
                notes.append(
                    f'metal {metal["element"]} (index {metal["index"]}) has '
                    f'only {n_ligands} ligand(s) within the primary cutoff; '
                    'check for a missing bridging ligand or water')

        return ligands

    @on_master
    def _resolve_residues(self, topology, coordinating_residues):
        """
        Turns the residues asked for into residue indices of the topology.

        A request is matched against the residue id ('130' or 130) and against
        the label the binding modes report ('ASP130'). A request that matches
        nothing is an error rather than a silently ignored line, and one that
        matches several residues - the same number in two chains - takes all
        of them and says which.

        :param topology:
            The OpenMM topology.
        :param coordinating_residues:
            The residues asked for, or None.

        :return:
            The residue indices, as a set.
        """

        if not coordinating_residues:
            return set()

        if isinstance(coordinating_residues, (str, int)):
            coordinating_residues = [coordinating_residues]

        residues = list(topology.residues())
        forced = set()

        for request in coordinating_residues:
            wanted = str(request).strip()
            matched = [
                residue for residue in residues
                if wanted in (str(residue.id), residue_label(residue))
            ]

            assert_msg_critical(
                len(matched) > 0, 'coordinating_residues asked for '
                f'{request}, which is not a residue of this structure')

            for residue in matched:
                forced.add(residue.index)

            found = ', '.join(f'{residue_label(residue)} '
                              f'(chain {residue.chain.id})' for residue in matched)
            self.ostream.print_info(
                f'Including {found} in the coordination sphere by request.')

        self.ostream.flush()

        return forced

    @on_master
    def _force_ligands(self, atoms, position_of, metal_indices, contacts,
                       notes, forced, metal_bond_cutoff, report_cutoff):
        """
        Makes the residues asked for ligands of their nearest metal.

        Only the closest sidechain donor of the residue is taken, and only
        when the cutoffs have not already made one of its atoms a ligand:
        a request is there to add coordination the distances missed, not to
        add a second contact to a residue that already has one.

        :param atoms:
            The candidate atoms.
        :param position_of:
            The callable returning the position of an atom index.
        :param metal_indices:
            The indices of the metal centers.
        :param contacts:
            The contacts found so far. Appended to, and promoted in, in place.
        :param notes:
            The list of review notes. Appended to in place.
        :param forced:
            The residue indices that are ligands whatever their distance.
        """

        if not forced:
            return

        already = {
            contact['res_index']
            for contact in contacts if contact['metals']
        }

        for res_index in sorted(forced):
            if res_index in already:
                continue

            donors = [
                atom for atom in atoms
                if atom.residue.index == res_index and atom.element is not None
                and atom.element.symbol in DONOR_ELEMENTS and atom.name not in
                BACKBONE_ATOM_NAMES and atom.index not in metal_indices
            ]

            assert_msg_critical(
                len(donors) > 0, 'coordinating_residues asked for '
                f'residue index {res_index}, whose sidechain has no '
                f'{list(DONOR_ELEMENTS)} atom to coordinate with')

            best = None
            for atom in donors:
                position = position_of(atom.index)
                for metal in metal_indices:
                    distance = float(np.linalg.norm(position - position_of(metal)))
                    if best is None or distance < best[0]:
                        best = (distance, atom, metal)

            distance, atom, metal = best
            label = residue_label(atom.residue)

            notes.append(
                f'{label} {atom.name} is {distance:.2f} A from a metal, '
                f'beyond the primary cutoff ({metal_bond_cutoff}), and '
                'was made a ligand because coordinating_residues asked for it')

            if distance > report_cutoff:
                self.ostream.print_warning(
                    f'{label} {atom.name} is {distance:.2f} A from its metal, '
                    'which is a long way for a bond. It is a ligand because '
                    'coordinating_residues asked for it.')
                self.ostream.flush()

            existing = [
                contact for contact in contacts if contact['index'] == atom.index
            ]

            if existing:
                # it was already reported as a near miss between the cutoffs
                existing[0]['metals'] = [metal]
                existing[0]['distances'] = [round(distance, 3)]
                continue

            contacts.append({
                'residue': label,
                'res_name': atom.residue.name,
                'res_index': atom.residue.index,
                'chain': atom.residue.chain.id,
                'atom': atom.name,
                'index': atom.index,
                'metals': [metal],
                'distances': [round(distance, 3)],
            })

    @on_master
    def _assign_binding_modes(self, ligands, notes, protected=None):
        """
        Classifies each ligand contact using the residue context.

        The modes a contact can end up with are monodentate, bidentate (two
        donor atoms of one residue on one metal), bridging_single (one donor
        atom on both metals), bridging_double (one donor atom on both metals
        while a second one reaches one of them), bridging_mu13 (a carboxylate
        reaching both metals through its two oxygens) and bridging (any other
        residue doing the same).

        Only the residues of BRIDGING_RESIDUES are allowed to reach two
        metals, through one donor atom or through two. A contact that would
        make any other residue bridge is trimmed back to its closest metal,
        and where the reach is spread over two atoms the whole far contact
        goes: carrying either further would put a metal bond into the force
        field that the chemistry does not support.

        :param ligands:
            The list of ligand contacts. Updated in place with a 'mode' entry,
            trimmed to one metal where the residue cannot bridge, and with any
            contact the rules discard removed from the list.
        :param notes:
            The list of review notes. Appended to in place.
        :param protected:
            The (residue index, atom name) pairs that a manual edit put there.
            The distance-driven demotions are heuristics reading an unrelaxed
            structure, and a contact asked for by hand outranks them, so a
            protected contact is kept and noted instead of being dropped. The
            rules that follow from what a residue can do chemically apply to
            it like to any other contact.
        """

        if protected is None:
            protected = set()

        by_residue = {}
        for ligand in ligands:
            by_residue.setdefault(ligand['res_index'], []).append(ligand)

        can_bridge = BRIDGING_RESIDUES
        carboxylates = CARBOXYLATE_RESIDUES

        # contacts the rules discard, taken out of the list at the end so that
        # the grouping is not disturbed halfway through
        discarded = set()

        for group in by_residue.values():
            res_name = group[0]['res_name']

            for ligand in group:
                if len(ligand['metals']) >= 2 and res_name not in can_bridge:
                    # the residue has no way of reaching two metals with one
                    # donor atom, so the longer contact is dropped rather than
                    # turned into a bond that the force field would then have
                    # to hold
                    closest = ligand['distances'].index(min(ligand['distances']))
                    dropped = [
                        f'{d:.2f} A' for i, d in enumerate(ligand['distances'])
                        if i != closest
                    ]
                    notes.append(
                        f'{ligand["residue"]} {ligand["atom"]} is within the '
                        f'primary cutoff of more than one metal, but '
                        f'{res_name} cannot bridge; keeping the closest '
                        f'contact and dropping {", ".join(dropped)}')
                    ligand['metals'] = [ligand['metals'][closest]]
                    ligand['distances'] = [ligand['distances'][closest]]

                if len(ligand['metals']) >= 2:
                    ligand['mode'] = 'bridging_single'
                else:
                    ligand['mode'] = 'monodentate'

            # One donor atom holding both metals while a second one of the
            # same residue reaches one of them is a single ligand binding
            # twice over, not a bridge with an unrelated contact next to it.
            # Only a residue of can_bridge gets this far: any other one was
            # trimmed to its closest metal just above.
            modes = {ligand['mode'] for ligand in group}
            if 'bridging_single' in modes and 'monodentate' in modes:
                for ligand in group:
                    ligand['mode'] = 'bridging_double'

            monodentate = all(ligand['mode'] == 'monodentate' for ligand in group)

            if len(group) >= 2 and monodentate:
                metals_hit = {ligand['metals'][0] for ligand in group}
                if len(metals_hit) >= 2:
                    if res_name not in can_bridge:
                        # two donor atoms on two metals is a bridge whatever
                        # it is called, so the residue is cut back to the one
                        # contact it can actually make. Keeping both would
                        # leave the residue labelled monodentate while still
                        # reporting, and bonding, two metals.
                        mode = None
                        nearest = min(group,
                                      key=lambda ligand: ligand['distances'][0])
                        far = [ligand for ligand in group if ligand is not nearest]
                        dropped = ', '.join(
                            f'{ligand["atom"]} at {ligand["distances"][0]:.2f} A'
                            for ligand in far)
                        notes.append(
                            f'{group[0]["residue"]} reaches two metals '
                            f'through {len(group)} atoms, but {res_name} '
                            f'cannot bridge; keeping {nearest["atom"]} at '
                            f'{nearest["distances"][0]:.2f} A and dropping '
                            f'{dropped}')
                        discarded.update(id(ligand) for ligand in far)
                        group[:] = [nearest]
                    elif res_name in carboxylates:
                        mode = 'bridging_mu13'
                    else:
                        mode = 'bridging'
                else:
                    mode = 'bidentate'

                if mode is not None:
                    for ligand in group:
                        ligand['mode'] = mode

            # A carboxylate whose two oxygens sit at very different distances
            # is not chelating: the far oxygen points away, and holding it at
            # the metal anyway would bend the group open.
            if (len(group) == 2 and res_name in carboxylates
                    and all(ligand['mode'] == 'bidentate' for ligand in group)):
                near, far = sorted(group, key=lambda ligand: ligand['distances'][0])
                separation = far['distances'][0] - near['distances'][0]
                asked_for = any((ligand['res_index'], ligand['atom']) in protected
                                for ligand in group)

                if separation > BIDENTATE_ASYMMETRY and asked_for:
                    notes.append(
                        f'{group[0]["residue"]} has {near["atom"]} at '
                        f'{near["distances"][0]:.2f} A and {far["atom"]} at '
                        f'{far["distances"][0]:.2f} A, a difference of '
                        f'{separation:.2f} A, which would normally read as '
                        f'monodentate through {near["atom"]}; both are kept '
                        'because a manual edit asked for them')
                elif separation > BIDENTATE_ASYMMETRY:
                    notes.append(
                        f'{group[0]["residue"]} has {near["atom"]} at '
                        f'{near["distances"][0]:.2f} A and {far["atom"]} at '
                        f'{far["distances"][0]:.2f} A, a difference of '
                        f'{separation:.2f} A; that is not a chelating '
                        f'bidentate, so it is read as monodentate through '
                        f'{near["atom"]} alone')
                    near['mode'] = 'monodentate'
                    discarded.add(id(far))
                    group[:] = [near]

            if res_name.startswith('HI') and len(group) >= 2:
                notes.append(
                    f'{group[0]["residue"]} appears to coordinate through '
                    'more than one ring nitrogen; imidazole geometry makes '
                    'this essentially impossible, so treat it as an artifact '
                    'and edit the binding modes')

        if discarded:
            ligands[:] = [
                ligand for ligand in ligands if id(ligand) not in discarded
            ]

    @on_master
    def _insert_contact(self, ligands, contact):
        """
        Puts a new ligand contact next to the others of its residue.

        A residue is read one row at a time, in the printed table and in the
        JSON alike, so a contact added by hand belongs with its siblings
        rather than at the end of the list.

        :param ligands:
            The ligand contacts. Updated in place.
        :param contact:
            The contact to insert.
        """

        position = None
        for index, ligand in enumerate(ligands):
            if ligand['res_index'] == contact['res_index']:
                position = index + 1

        if position is None:
            ligands.append(contact)
        else:
            ligands.insert(position, contact)

    @on_master
    def _add_contact_metal(self, contact, metal_index, distance):
        """
        Records one more metal on a ligand contact.

        The metals of a contact are kept sorted by index, the way
        _collect_ligands writes them, with the distances alongside.

        :param contact:
            The ligand contact. Updated in place.
        :param metal_index:
            The atom index of the metal.
        :param distance:
            The distance in Angstrom.
        """

        pairs = list(zip(contact['metals'], contact['distances']))
        pairs.append((metal_index, round(float(distance), 3)))
        pairs.sort(key=lambda pair: pair[0])

        contact['metals'] = [index for index, _ in pairs]
        contact['distances'] = [value for _, value in pairs]

    @on_master
    def _record_manual_bond(self,
                            records,
                            res_index,
                            atom_name,
                            metal_res_index,
                            action,
                            equilibrium=None):
        """
        Writes a manual edit down so that a re-detection can replay it.

        A bond is named by residue index, atom name and the residue index of
        the metal, none of which the renumbering of protonate touches, and all
        of which survive JSON. One edit per bond is kept: editing the same
        bond again overwrites what it did rather than stacking a second record
        behind it.

        :param records:
            The manual bond records. Appended to, in place.
        :param res_index:
            The residue index of the ligand.
        :param atom_name:
            The name of the donor atom.
        :param metal_res_index:
            The residue index of the metal center.
        :param action:
            Either 'add' or 'remove'.
        :param equilibrium:
            The distance in nanometers to pull this bond to in the crude pass,
            or None to measure it on the geometry like any other. Written onto
            the record only when there is one, so a record says exactly what
            was asked for. Read back by manual_bond_equilibria.
        """

        for record in records:
            if (record['res_index'] == res_index and record['atom'] == atom_name
                    and record['metal_res_index'] == metal_res_index):
                record['action'] = action
                # the key is only there when there is a distance to say, so
                # editing the bond again without one takes it back out
                record.pop('equilibrium', None)
                if equilibrium is not None:
                    record['equilibrium'] = equilibrium
                return

        new_record = {
            'res_index': res_index,
            'atom': atom_name,
            'metal_res_index': metal_res_index,
            'action': action,
        }

        if equilibrium is not None:
            new_record['equilibrium'] = equilibrium

        records.append(new_record)

    @on_master
    def _manual_protected(self, records):
        """
        The contacts that the classification rules must not take away again.

        :param records:
            The manual bond records.

        :return:
            The set of (residue index, atom name) pairs that were asked for.
        """

        return {(record['res_index'], record['atom'])
                for record in records if record['action'] == 'add'}

    @on_master
    def _apply_manual_bonds(self, ligands, records, metals, atoms, position_of, notes):
        """
        Replays the manual coordination edits onto a re-detected ligand list.

        The records name their atoms by residue index and atom name, so they
        outlive the renumbering of protonate and can be applied to any later
        detection of the same structure.

        :param ligands:
            The ligand contacts as detected. Updated in place.
        :param records:
            The manual bond records.
        :param metals:
            The metal entries of the binding modes.
        :param atoms:
            The candidate atoms the detection ran over.
        :param position_of:
            The callable returning the position in Angstrom of an atom index.
        :param notes:
            The list of review notes. Appended to in place.
        """

        if not records:
            return

        metal_by_res = {metal['res_index']: metal for metal in metals}
        atom_by_key = {(atom.residue.index, atom.name): atom for atom in atoms}

        for record in records:
            metal_entry = metal_by_res.get(record['metal_res_index'])
            atom = atom_by_key.get((record['res_index'], record['atom']))

            if metal_entry is None or atom is None:
                # a removal whose atom is not here has nothing left to undo,
                # which is the usual case: taking the last bond off a residue
                # takes the residue out of the active site as well
                if record['action'] == 'add':
                    notes.append(
                        f'the manual bond of residue index '
                        f'{record["res_index"]} {record["atom"]} could not be '
                        'applied because the atom or its metal is not part of '
                        'this active site')
                continue

            metal_index = metal_entry['index']
            metal_label = f'{metal_entry["element"]} (index {metal_index})'

            contact = None
            for ligand in ligands:
                if ligand['index'] == atom.index:
                    contact = ligand
                    break

            if record['action'] == 'remove':
                if contact is None or metal_index not in contact['metals']:
                    continue
                label = f'{contact["residue"]} {contact["atom"]}'
                position = contact['metals'].index(metal_index)
                contact['metals'].pop(position)
                contact['distances'].pop(position)
                if not contact['metals']:
                    ligands.remove(contact)
                notes.append(f'{label} is not bonded to {metal_label} '
                             'because remove_metal_bond asked for it')
                continue

            if contact is not None and metal_index in contact['metals']:
                continue

            try:
                distance = float(
                    np.linalg.norm(
                        position_of(atom.index) - position_of(metal_index)))
            except (KeyError, IndexError):
                notes.append(f'the manual bond of {atom.residue.name}'
                             f'{atom.residue.id} {atom.name} could not be applied '
                             'because the geometry holds no position for it')
                continue

            if contact is None:
                contact = {
                    'residue': residue_label(atom.residue),
                    'res_name': atom.residue.name,
                    'res_index': atom.residue.index,
                    'chain': atom.residue.chain.id,
                    'atom': atom.name,
                    'index': atom.index,
                    'metals': [],
                    'distances': [],
                }
                self._insert_contact(ligands, contact)

            self._add_contact_metal(contact, metal_index, distance)
            notes.append(f'{contact["residue"]} {contact["atom"]} is bonded '
                         f'to {metal_label} at {distance:.2f} A because '
                         'add_metal_bond asked for it')

    # ------------------------------------------------------------------
    # editing
    # ------------------------------------------------------------------

    @on_master
    def add_metal_bond(self,
                       request,
                       binding_modes,
                       topology,
                       positions,
                       resid,
                       metal,
                       metal_bond_cutoff,
                       report_cutoff,
                       atom=None,
                       chain=None,
                       equilibrium=None):
        """
        Bonds a residue to a metal center that the distances did not connect.

        The cutoffs read an unrelaxed structure, so a contact the design
        intends can sit just outside them, and reviewing the suggested
        coordination means being able to correct it. The edit is recorded in
        binding_modes['manual_bonds'] by residue index, atom name and metal
        residue index, none of which the renumbering of protonate touches, so
        that update_binding_modes puts the bond back after it has re-detected
        the coordination on a relaxed geometry.

        The modes of the whole residue are classified again afterwards, so a
        second oxygen added to a monodentate carboxylate turns the pair into a
        bidentate or a mu-1,3 bridge on its own. What the chemistry forbids is
        refused rather than quietly trimmed: only BRIDGING_RESIDUES can reach
        two metals, and the truncation is sidechain-only.

        :param binding_modes:
            The binding modes to edit. Not modified.
        :param topology:
            The OpenMM topology the binding modes index into.
        :param positions:
            The positions as an (N, 3) numpy array in Angstrom.
        :param resid:
            The residue, as an id ('130' or 130) or as a label ('ASP130').
        :param metal:
            The atom index of the metal center to bond to.
        :param metal_bond_cutoff:
            The bonding cutoff in Angstrom, which the distance of the new bond
            is reported against.
        :param report_cutoff:
            The reporting cutoff in Angstrom, beyond which the new bond is
            warned about as a long one.
        :param atom:
            The donor atom of the residue, as an atom name ('OE1') or as an
            atom index. Resolved automatically when left out, which is only
            possible where the choice is unambiguous: a sidechain with a
            single donor, or a carboxylate that coordinates nothing yet, whose
            two oxygens are interchangeable until one of them is bound.
        :param chain:
            The chain the residue is in. Only needed when the same residue
            number occurs in several chains.
        :param equilibrium:
            The distance in Angstrom the crude pre-QM pass should pull this
            bond to, instead of measuring it on the structure it was given. A
            bond is often added by hand because the residue is turned the wrong
            way round, and measuring the equilibrium on that geometry only pins
            the mistake in place; naming a distance swings the sidechain into
            position before any QM is paid for. It is recorded with the bond,
            so it survives every re-derivation, and only the seeded pass reads
            it -- once a Hessian exists the equilibria come from the geometry it
            was computed on. Note this is in Angstrom, while the
            metal_bond_equilibria setting is in nanometers.

        :return:
            The edited request, or the argument itself when the bond was
            already there.
        """

        positions = np.asarray(positions)

        metal_entry = self._resolve_metal(binding_modes, metal)
        residue = self._resolve_residue(topology, resid, chain)
        ligand_atom = self._resolve_ligand_atom(residue, atom, metal_entry, positions,
                                                binding_modes['ligands'])

        label = residue_label(residue)
        metal_label = f'{metal_entry["element"]} (index {metal_entry["index"]})'

        assert_msg_critical(
            ligand_atom.element is not None
            and ligand_atom.element.symbol in DONOR_ELEMENTS, 'add_metal_bond: '
            f'{label} {ligand_atom.name} is not one of '
            f'{list(DONOR_ELEMENTS)}, so it has no lone pair to donate '
            'to a metal')

        assert_msg_critical(
            ligand_atom.name not in BACKBONE_ATOM_NAMES, 'add_metal_bond: '
            f'{label} {ligand_atom.name} is a backbone atom, and the '
            'truncation scheme is sidechain-only, so it cannot be made a '
            'ligand')

        residue_contacts = [
            ligand for ligand in binding_modes['ligands']
            if ligand['res_index'] == residue.index
        ]

        already = [
            ligand for ligand in residue_contacts
            if ligand['index'] == ligand_atom.index
            and metal_entry['index'] in ligand['metals']
        ]

        if already:
            self.ostream.print_warning(
                f'{label} {ligand_atom.name} is already bonded to '
                f'{metal_label}; the request is left as it is')
            self.ostream.flush()
            return request

        reached = {metal_entry['index']}
        for ligand in residue_contacts:
            reached.update(ligand['metals'])

        assert_msg_critical(
            len(reached) < 2 or residue.name in BRIDGING_RESIDUES,
            'add_metal_bond: the bond would make '
            f'{label} reach {len(reached)} metals, and only '
            f'{list(BRIDGING_RESIDUES)} can bridge. An imidazole '
            'nitrogen has a single lone pair in the ring plane, so a second '
            'metal at bonding distance is a geometric artifact rather than a '
            'bond')

        distance = float(
            np.linalg.norm(positions[ligand_atom.index] -
                           positions[metal_entry['index']]))

        if distance > report_cutoff:
            self.ostream.print_warning(
                f'{label} {ligand_atom.name} is {distance:.2f} A from '
                f'{metal_label}, beyond even the secondary cutoff '
                f'({report_cutoff} A), which is a long way for a bond. '
                'It is a ligand because add_metal_bond asked for it.')
        elif distance > metal_bond_cutoff:
            self.ostream.print_warning(
                f'{label} {ligand_atom.name} is {distance:.2f} A from '
                f'{metal_label}, beyond the primary cutoff '
                f'({metal_bond_cutoff} A). It is a ligand because '
                'add_metal_bond asked for it.')

        if residue.name.startswith('HI') and any(
                ligand['index'] != ligand_atom.index
                for ligand in residue_contacts):
            self.ostream.print_warning(
                f'{label} would coordinate through more than one ring '
                'nitrogen; imidazole geometry makes that essentially '
                'impossible, so check the structure before running on it')

        if binding_modes.get('variants'):
            self.ostream.print_warning(
                'This structure has already been protonated. The edit changes '
                f'no hydrogen, so protonate again if the variant of {label} '
                'depends on this bond.')

        self.ostream.flush()

        assert_msg_critical(
            equilibrium is None or equilibrium > 0, 'add_metal_bond: the '
            f'equilibrium distance must be positive, got {equilibrium}')

        new_request = deepcopy(request)
        records = new_request.setdefault('manual_bonds', [])

        # the force field keeps its bond lengths in nanometers, while everything
        # a caller is shown here -- the cutoffs, the distances in these warnings
        # -- is in Angstrom, so the conversion happens once, at the boundary
        self._record_manual_bond(records,
                                 residue.index,
                                 ligand_atom.name,
                                 metal_entry['res_index'],
                                 'add',
                                 equilibrium=None if equilibrium is None else 0.1 *
                                 float(equilibrium))

        self.ostream.print_info(f'Bonded {label} {ligand_atom.name} to {metal_label} at '
                                f'{distance:.2f} A.')

        if equilibrium is not None:
            self.ostream.print_info(
                f'The crude pass will pull it to {equilibrium:.2f} A rather than '
                f'holding it at the {distance:.2f} A the structure has. The QM '
                'optimization and the fit that follows are not bound by it.')

        self.ostream.flush()

        return new_request

    @on_master
    def remove_metal_bond(self,
                          request,
                          binding_modes,
                          resid,
                          metal=None,
                          atom=None,
                          chain=None):
        """
        Takes a residue's bond to a metal center back out.

        The counterpart of add_metal_bond, and recorded the same way, so that
        a contact the cutoffs invented does not come back when
        update_binding_modes re-detects the coordination. The residue is
        looked up in the binding modes themselves, which carry its label and
        its chain, so no topology is needed.

        :param binding_modes:
            The binding modes to edit. Not modified.
        :param resid:
            The residue, as an id ('130' or 130) or as a label ('ASP130').
        :param metal:
            The atom index of the metal to unbind from. Only needed when the
            residue bridges two metals, since otherwise there is one bond to
            remove.
        :param atom:
            The donor atom to unbind, as an atom name or as an atom index.
            Left out, every contact the residue makes with that metal goes,
            which is what removing the bond of a bidentate means; naming an
            atom drops that arm alone and leaves the rest.
        :param chain:
            The chain the residue is in. Only needed when the same residue
            number occurs in several chains.

        :return:
            The edited binding modes.
        """

        wanted = str(resid).strip()
        matched = [
            ligand for ligand in binding_modes['ligands']
            if wanted in (ligand['residue'],
                          ligand['residue'][len(ligand['res_name']):]) and (
                              chain is None or str(ligand['chain']) == str(chain))
        ]

        assert_msg_critical(
            len(matched) > 0, 'remove_metal_bond: residue '
            f'{resid} coordinates no metal in these binding modes, so there '
            'is no bond to remove')

        residue_indices = {ligand['res_index'] for ligand in matched}
        found = ', '.join(
            sorted({
                f'{ligand["residue"]} (chain {ligand["chain"]})'
                for ligand in matched
            }))

        assert_msg_critical(
            len(residue_indices) == 1, 'remove_metal_bond: residue '
            f'{resid} matches {found}; pass chain= to say which one is meant')

        label = matched[0]['residue']
        res_index = matched[0]['res_index']
        reached = sorted(
            {index
             for ligand in matched
             for index in ligand['metals']})

        if metal is None:
            assert_msg_critical(
                len(reached) == 1, 'remove_metal_bond: '
                f'{label} bridges the metals at indices {reached}; pass '
                'metal= to say which bond to remove')
            metal_entry = self._resolve_metal(binding_modes, reached[0])
        else:
            metal_entry = self._resolve_metal(binding_modes, metal)
            assert_msg_critical(
                metal_entry['index'] in reached, 'remove_metal_bond: '
                f'{label} is not bonded to the metal at index '
                f'{metal_entry["index"]}; it reaches {reached}')

        metal_index = metal_entry['index']
        metal_label = f'{metal_entry["element"]} (index {metal_index})'

        targets = [ligand for ligand in matched if metal_index in ligand['metals']]

        if atom is not None:
            index = self._as_atom_index(atom)
            if index is None:
                name = str(atom).strip()
                targets = [ligand for ligand in targets if ligand['atom'] == name]
            else:
                targets = [ligand for ligand in targets if ligand['index'] == index]

            assert_msg_critical(
                len(targets) > 0, 'remove_metal_bond: '
                f'{label} {atom} is not bonded to {metal_label}; its bonds '
                'to it are ' +
                ', '.join(ligand['atom']
                          for ligand in matched if metal_index in ligand['metals']))

        if binding_modes.get('variants'):
            self.ostream.print_warning(
                'This structure has already been protonated. The edit changes '
                f'no hydrogen, so protonate again if the variant of {label} '
                'depends on this bond.')
            self.ostream.flush()

        new_request = deepcopy(request)
        records = new_request.setdefault('manual_bonds', [])

        removed = []
        for ligand in targets:
            removed.append(ligand['atom'])
            self._record_manual_bond(records, res_index, ligand['atom'],
                                     metal_entry['res_index'], 'remove')

        # a residue that is only a ligand because coordinating_residues named it
        # would have the bond put straight back on the next derivation, so the
        # request that would do that goes with the bond
        forced = new_request.get('coordinating_residues', [])
        still_bound = any(
            ligand['res_index'] == res_index and set(ligand['metals']) -
            {metal_index} for ligand in binding_modes['ligands'])

        if res_index in forced and not still_bound:
            new_request['coordinating_residues'] = [
                index for index in forced if index != res_index
            ]
            self.ostream.print_info(
                f'{label} no longer coordinates anything and was taken out of '
                'coordinating_residues, which would otherwise put the bond back on '
                'the next derivation.')

        self.ostream.print_info(
            f'Removed the bond(s) of {label} {", ".join(removed)} to '
            f'{metal_label}.')
        self.ostream.flush()

        return new_request

    @on_master
    def include_residue(self, request, binding_modes, topology, resid,
                        chain=None):
        """
        Puts a residue into the truncated active site.

        The cluster is otherwise exactly the residues that coordinate a
        metal. This adds one that does not have to -- a second-shell residue
        that hydrogen bonds to a ligand, or one whose sidechain the QM should
        see for any other reason. It is truncated and capped like every other
        residue, and it keeps whatever protonation the pH gives it unless
        update_protonation_state says otherwise.

        :param request:
            The record of what has been decided about the site. Not
            modified.
        :param binding_modes:
            The coordination as it stands, for what the site already holds.
        :param topology:
            The topology the residue is looked up in.
        :param resid:
            The residue, as an id ('58' or 58) or as a label ('TYR58').
        :param chain:
            The chain id, when the residue id occurs in more than one chain.

        :return:
            The edited request, or None when the residue was already in the
            site and there was nothing to do.
        """

        residue = self._resolve_residue(topology, resid, chain)
        self.check_truncatable(residue)

        label = residue_label(residue)
        already = residue.index in active_site_residues(binding_modes)
        excluded = residue.index in request.get('excluded_residues', [])

        if already and not excluded:
            self.ostream.print_warning(
                f'{label} is already part of the active site; nothing to '
                'include')
            self.ostream.flush()
            return None

        request = deepcopy(request)
        request['excluded_residues'] = [
            index for index in request.get('excluded_residues', [])
            if index != residue.index
        ]
        extra = set(request.get('extra_residues', []))
        extra.add(residue.index)
        request['extra_residues'] = sorted(extra)

        self.ostream.print_info(f'Added {label} to the active site.')
        self.ostream.flush()

        return request

    @on_master
    def remove_residue(self,
                       request,
                       binding_modes,
                       topology,
                       positions,
                       resid,
                       chain,
                       metal_elements,
                       metal_formal_charges,
                       metal_bond_cutoff,
                       report_cutoff):
        """
        Takes a residue out of the truncated active site.

        Any metal bonds it makes go with it, recorded the way
        remove_metal_bond records them, so that neither the residue nor its
        coordination comes back when the site is detected again on a relaxed
        geometry.

        A metal that would be left with no ligand at all, and the last
        residue of the site, are refused: what is left would not be an active
        site.

        :param request:
            The record of what has been decided about the site. Not
            modified.
        :param binding_modes:
            The coordination as it stands.
        :param topology:
            The topology the residue is looked up in, and the coordination
            is derived on again after each bond goes.
        :param positions:
            Its positions in Angstrom.
        :param resid:
            The residue, as an id ('130' or 130) or as a label ('ASP130').
        :param chain:
            The chain id, when the residue id occurs in more than one chain,
            else None.
        :param metal_elements:
            The elements treated as metal centers.
        :param metal_formal_charges:
            The formal charges assumed for the metal ions, by element.
        :param metal_bond_cutoff:
            The bonding cutoff in Angstrom of the re-detection.
        :param report_cutoff:
            The reporting cutoff in Angstrom of the re-detection.

        :return:
            The edited request.
        """

        residue = self._resolve_residue(topology, resid, chain)
        modes = binding_modes
        label = residue_label(residue)
        members = active_site_residues(modes)

        assert_msg_critical(
            residue.index in members, 'remove_residue: '
            f'{label} is not part of the active site, so there is '
            'nothing to remove')
        assert_msg_critical(
            len(members) > 1, 'remove_residue: '
            f'{label} is the only residue of the active site, and what '
            'would be left is not one')

        orphaned = sorted({
            f'{metal["element"]} (index {metal["index"]})'
            for metal in modes['metals']
            if not any(metal['index'] in ligand['metals']
                       for ligand in modes['ligands']
                       if ligand['res_index'] != residue.index)
        })
        assert_msg_critical(
            not orphaned, 'remove_residue: removing '
            f'{label} would leave {", ".join(orphaned)} with no ligand '
            'at all')

        detection = {
            'metal_elements': metal_elements,
            'metal_formal_charges': metal_formal_charges,
            'metal_bond_cutoff': metal_bond_cutoff,
            'report_cutoff': report_cutoff,
        }

        # the bonds have to go through remove_metal_bond, so that the
        # removals are recorded and a re-detection does not put back what
        # the geometry still looks like
        while any(ligand['res_index'] == residue.index
                  for ligand in modes['ligands']):
            bound = next(ligand for ligand in modes['ligands']
                         if ligand['res_index'] == residue.index)
            request = self.remove_metal_bond(request,
                                             modes,
                                             label,
                                             metal=bound['metals'][0])
            modes = self.derive_binding_modes(topology, positions, request,
                                              **detection)

        request = deepcopy(request)
        request['extra_residues'] = [
            index for index in request.get('extra_residues', [])
            if index != residue.index
        ]
        request['coordinating_residues'] = [
            index for index in request.get('coordinating_residues', [])
            if index != residue.index
        ]
        excluded = set(request.get('excluded_residues', []))
        excluded.add(residue.index)
        request['excluded_residues'] = sorted(excluded)

        self.ostream.print_info(f'Removed {label} from the active site.')
        self.ostream.flush()

        return request

    @on_master
    def _resolve_residue(self, topology, resid, chain):
        """
        Finds the one residue a manual edit is about.

        Matched the way _resolve_residues matches, against the residue id
        ('130' or 130) and against the label the binding modes report
        ('ASP130'), but an edit names a single residue, so a request that
        matches several of them is an error naming the chains rather than a
        set to work through.

        :param topology:
            The OpenMM topology.
        :param resid:
            The residue, as an id or as a label.
        :param chain:
            The chain the residue is in, or None.

        :return:
            The residue.
        """

        wanted = str(resid).strip()
        matched = [
            residue for residue in topology.residues()
            if wanted in (str(residue.id), residue_label(residue)) and (
                chain is None or str(residue.chain.id) == str(chain))
        ]

        in_chain = '' if chain is None else f' of chain {chain}'

        assert_msg_critical(
            len(matched) > 0, f'residue {resid}'
            f'{in_chain} is not a residue of this structure')

        found = ', '.join(f'{residue_label(residue)} '
                          f'(chain {residue.chain.id})' for residue in matched)

        assert_msg_critical(
            len(matched) == 1, f'residue {resid} matches {found}; '
            'pass chain= to say which one is meant')

        return matched[0]

    @on_master
    def _resolve_metal(self, binding_modes, metal):
        """
        Finds the metal entry an atom index names.

        The index is read with _as_atom_index rather than compared as given,
        because the number a caller has is the one off the active site's labels
        and those are strings: show_active_site draws '85' on the metal, and
        passing that straight back has to name the same atom as 85 does.

        :param binding_modes:
            The binding modes.
        :param metal:
            The atom index of a metal center, as an int or as the string the
            labels show.

        :return:
            The metal entry of the binding modes.
        """

        entries = {entry['index']: entry for entry in binding_modes['metals']}

        known = ', '.join(f'{entry["element"]} {index}'
                          for index, entry in sorted(entries.items()))

        wanted = self._as_atom_index(metal)

        assert_msg_critical(
            wanted is not None and wanted in entries,
            f'there is no metal center at atom index {metal}. '
            f'This site holds {known}, which is what its labels show')

        return entries[wanted]

    @on_master
    def _as_atom_index(self, atom):
        """
        Reads an atom argument as an atom index, or decides that it is a name.

        :param atom:
            The atom, as an index or as an atom name.

        :return:
            The index, or None when the argument is a name.
        """

        if isinstance(atom, (int, np.integer)) and not isinstance(atom, bool):
            return int(atom)

        text = str(atom).strip()

        return int(text) if text.isdigit() else None

    @on_master
    def _sidechain_donors(self, residue):
        """
        The atoms of a residue that could donate to a metal.

        :param residue:
            The residue.

        :return:
            The donor atoms, backbone excluded.
        """

        return [
            atom for atom in residue.atoms()
            if atom.element is not None and atom.element.symbol in DONOR_ELEMENTS
            and atom.name not in BACKBONE_ATOM_NAMES
        ]

    @on_master
    def _resolve_ligand_atom(self, residue, atom, metal_entry, positions, ligands):
        """
        Finds the donor atom of a residue that a manual bond is about.

        An atom given by name or by index is looked up and checked to belong
        to the residue. Left out, it is resolved only where the choice does
        not matter: a sidechain with a single donor has nothing to choose
        between, and the two oxygens of a carboxylate that coordinates nothing
        yet are interchangeable, so the nearer one is taken. Anything else -
        the two ring nitrogens of a histidine, the second oxygen of a
        carboxylate that already binds through the first - is a chemical
        choice, and is asked for rather than guessed.

        :param residue:
            The residue.
        :param atom:
            The atom, as a name, as an index, or None.
        :param metal_entry:
            The metal entry the bond is to.
        :param positions:
            The positions as an (N, 3) numpy array in Angstrom.
        :param ligands:
            The ligand contacts recorded so far, which say whether the residue
            already coordinates something.

        :return:
            The atom.
        """

        atoms = list(residue.atoms())
        donors = self._sidechain_donors(residue)
        label = residue_label(residue)
        names = [donor.name for donor in donors]

        if atom is not None:
            index = self._as_atom_index(atom)

            if index is not None:
                match = [entry for entry in atoms if entry.index == index]
                owner = None
                if not match:
                    for other in residue.chain.topology.atoms():
                        if other.index == index:
                            owner = (f'{other.residue.name}'
                                     f'{other.residue.id} {other.name}')
                            break
                assert_msg_critical(
                    len(match) == 1, 'add_metal_bond: atom index '
                    f'{index} is not part of {label}' +
                    (f', it is {owner}' if owner is not None else ''))
                return match[0]

            name = str(atom).strip()
            match = [entry for entry in atoms if entry.name == name]
            assert_msg_critical(
                len(match) == 1, 'add_metal_bond: '
                f'{label} has no atom named {name}; its sidechain donors are '
                f'{names}')
            return match[0]

        assert_msg_critical(
            len(donors) > 0, 'add_metal_bond: the '
            f'sidechain of {label} has no {list(DONOR_ELEMENTS)} atom to '
            'coordinate with')

        if len(donors) == 1:
            return donors[0]

        bound = [
            ligand for ligand in ligands
            if ligand['res_index'] == residue.index and ligand['metals']
        ]
        oxygens = [donor for donor in donors if donor.element.symbol == 'O']

        if (residue.name in CARBOXYLATE_RESIDUES and not bound
                and len(oxygens) == 2):
            # until one of them is bound the two oxygens are interchangeable,
            # so which one is picked does not matter; the nearer one keeps the
            # geometry closest to what it already is
            metal_position = positions[metal_entry['index']]
            return min(oxygens,
                       key=lambda donor: np.linalg.norm(positions[donor.index] -
                                                        metal_position))

        assert_msg_critical(
            False, 'add_metal_bond: which atom of '
            f'{label} binds the metal is a chemical choice rather than a '
            f'formality; pass atom= with one of {names}')

    # ------------------------------------------------------------------
    # protonation
    # ------------------------------------------------------------------

    @on_master
    def _histidine_variant(self, residue, positions, metal_positions):
        """
        Chooses the tautomer of a coordinating histidine.

        The choice is made only between the two ring nitrogens, on whichever
        of them sits closest to a metal, and the tautomer is then set so that
        this nitrogen carries no hydrogen. Ring carbons are ignored even when
        one of them is nearer to the metal than either nitrogen.

        :param residue:
            The histidine residue.
        :param positions:
            The positions as an (N, 3) numpy array in Angstrom.
        :param metal_positions:
            The positions of the metal centers.

        :return:
            The tuple of the variant name and a note, or None for no note.
        """

        def closest_metal_distance(atom):
            return min(
                float(np.linalg.norm(positions[atom.index] - metal_position))
                for metal_position in metal_positions)

        sidechain = [
            atom for atom in residue.atoms()
            if atom.name not in BACKBONE_ATOM_NAMES and atom.name != 'CA'
        ]
        ring_nitrogens = {
            atom.name: atom
            for atom in sidechain if atom.name in ('ND1', 'NE2')
        }

        if len(ring_nitrogens) < 2:
            return 'HID', (f'{residue_label(residue)} does not have both ring '
                           'nitrogens; defaulting to HID, please check')

        distances = {
            name: closest_metal_distance(atom)
            for name, atom in ring_nitrogens.items()
        }
        coordinating = min(distances, key=distances.get)

        # the coordinating nitrogen is the one that must not carry a hydrogen
        variant = 'HIE' if coordinating == 'ND1' else 'HID'

        note = None
        nearest = min(sidechain, key=closest_metal_distance)
        if nearest.name not in ('ND1', 'NE2'):
            note = (f'{residue_label(residue)} has {nearest.name} closer to a '
                    f'metal ({closest_metal_distance(nearest):.2f} A) than either '
                    f'ring nitrogen (ND1 {distances["ND1"]:.2f} A, NE2 '
                    f'{distances["NE2"]:.2f} A); the tautomer was still chosen '
                    'from the nitrogens, but the geometry looks distorted')

        return variant, note

    @on_master
    def suggest_variants(self, topology, positions, binding_modes,
                         protonation_overrides):
        """
        Chooses a protonation variant for each coordinating residue.

        coordinating carboxylates are deprotonated
        a coordinating cysteine is a thiolate
        and the histidine tautomer is set so that the coordinating nitrogen
        carries no hydrogen. Entries in protonation_overrides win over
        the rules.

        Every one of these rules is about a sidechain that coordinates a
        metal, so a residue that include_residue put in the cluster without
        one is deliberately left out of them: it keeps whatever Modeller
        picks for the pH, and update_protonation_state is how it is told
        otherwise. What Modeller picked is recorded by protonate.

        :param topology:
            The OpenMM topology.
        :param positions:
            The positions as an (N, 3) numpy array in Angstrom.
        :param binding_modes:
            The binding modes. Not modified; anything worth recording is
            returned as a note.
        :param protonation_overrides:
            The variants asked for by hand, keyed on residue label, id or
            index; None for none.

        :return:
            The tuple of the residue index to variant mapping and a list of
            notes for review.
        """

        residues = list(topology.residues())
        positions = np.asarray(positions)
        metal_positions = [
            positions[metal['index']] for metal in binding_modes['metals']
        ]

        by_residue = {}
        for ligand in binding_modes['ligands']:
            by_residue.setdefault(ligand['res_index'], []).append(ligand)

        notes = []
        variants = {}
        for res_index, group in by_residue.items():
            res_name = residues[res_index].name

            if res_name in ('ASP', 'ASH'):
                variant = 'ASP'
            elif res_name in ('GLU', 'GLH'):
                variant = 'GLU'
            elif res_name in ('CYS', 'CYX'):
                variant = 'CYX'
            elif res_name.startswith('HI'):
                variant, note = self._histidine_variant(residues[res_index], positions,
                                                        metal_positions)
                if note is not None:
                    notes.append(note)
            else:
                variant = None

            if variant is not None:
                variants[res_index] = variant

        overrides = protonation_overrides or {}
        for key, variant in overrides.items():
            by_id = [
                residue for residue in residues
                if str(residue.id) == str(key) or residue_label(residue) == str(key)
            ]
            matched = by_id or [
                residue for residue in residues if residue.index == key
            ]

            assert_msg_critical(
                len(matched) > 0, 'suggest_variants: override '
                f'residue {key} not found')

            named = ', '.join(f'{residue_label(residue)} '
                              f'(chain {residue.chain.id})' for residue in matched)

            assert_msg_critical(
                len(matched) == 1 or by_id, 'suggest_variants: override '
                f'{key} matches {named}; name the residue as a label such as '
                f'{matched[0].name}{matched[0].id} to say which is meant')

            for residue in matched:
                check_variant(residue, variant)
                variants[residue.index] = variant

        return variants, notes

    @on_master
    def protonate(self, topology, positions, binding_modes,
                  protonation_overrides):
        """
        Adds hydrogens with the protonation variants that the metal site
        requires.

        :param topology:
            The OpenMM topology.
        :param positions:
            The positions as an (N, 3) numpy array in Angstrom.
        :param binding_modes:
            The coordination, read for the variant rules. Not modified.
        :param protonation_overrides:
            The variants asked for by hand, keyed on residue label, id or
            index; None for none.

        :return:
            The tuple of the protonated topology, the positions in Angstrom,
            the residue index to variant mapping of what was actually built,
            and the notes the variant choice produced.
        """

        assert_msg_critical('openmm' in sys.modules,
                            'protonate: openmm is required')

        # Figure out the correct residue variants
        # based on the positions and the binding modes
        variants_by_index, notes = self.suggest_variants(
            topology, positions, binding_modes, protonation_overrides)

        variant_list = [None] * topology.getNumResidues()
        for res_index, variant in variants_by_index.items():
            variant_list[res_index] = variant

        modeller = mmapp.Modeller(topology, np.asarray(positions) * mmunit.angstrom)
        actual_variants = modeller.addHydrogens(variants=variant_list)

        topology_residues = list(topology.residues())
        for res_index in active_site_residues(binding_modes):
            if variants_by_index.get(res_index) is not None:
                continue
            chosen = actual_variants[res_index]
            variants_by_index[res_index] = (chosen if isinstance(chosen, str) else
                                            topology_residues[res_index].name)

        new_topology = modeller.topology
        new_positions = np.array(modeller.positions.value_in_unit(mmunit.angstrom))

        return new_topology, new_positions, variants_by_index, notes

    @on_master
    def update_protonation_state(self,
                                 protonation_overrides,
                                 topology,
                                 resid,
                                 variant,
                                 chain=None):
        """
        Records the protonation variant a residue is to be built with.

        Keyed by the residue's label rather than its index, because residue
        ids and residue indices overlap: a single chain numbered from one has
        index i for id i+1, so an index written here would also match its
        neighbour by id. The variant is checked against what OpenMM can build
        and what the charge table knows before it is recorded.

        :param protonation_overrides:
            The overrides as they stand, by residue label, or None. Not
            modified.
        :param topology:
            The topology the residue is looked up in.
        :param resid:
            The residue, as an id ('130' or 130) or as a label ('ASP130').
        :param variant:
            The variant to set, as OpenMM names it.
        :param chain:
            The chain id, when the residue id occurs in more than one chain.

        :return:
            The overrides with the residue's variant recorded.
        """

        residue = self._resolve_residue(topology, resid, chain)
        check_variant(residue, variant)

        overrides = dict(protonation_overrides or {})
        overrides[residue_label(residue)] = variant

        self.ostream.print_info(
            f'{residue_label(residue)} will be protonated as {variant}.')
        self.ostream.flush()

        return overrides

    # ------------------------------------------------------------------
    # extraction
    # ------------------------------------------------------------------

    @on_master
    def check_truncatable(self, residue):
        """
        Checks that a residue can be cut at its CA-CB bond.

        The truncation is sidechain-only and fixed, and a residue that only
        ever got here by coordinating a metal is always one the rule fits.
        include_residue takes any residue at all, so the two it does not fit
        are refused by name rather than left to fail later on a missing CB or,
        worse, to succeed with a dangling valence.

        :param residue:
            The residue to be truncated.
        """

        reason = UNTRUNCATABLE_RESIDUES.get(residue.name)

        assert_msg_critical(
            reason is None, 'check_truncatable: '
            f'{residue_label(residue)} {reason}, and the truncation of this '
            'module cuts every sidechain at CA-CB. It cannot be part of the '
            'active site')

    @on_master
    def extract_active_site(self, topology, positions, binding_modes,
                            cap_bond_length):
        """
        Builds the truncated QM active site.

        Sidechains are cut at the CA-CB bond and capped with a hydrogen placed
        along the CB to CA direction. No second-shell fragments and no
        backbone

        The connectivity is included in the returned data under 'connectivity_matrix'.
        Which atoms are bonded is a property of the site
        that was extracted, and nothing downstream should be able to pair the
        two up wrongly.

        It also comes with 'labels', one string per atom, holding what
        add_metal_bond and remove_metal_bond want to be told about that atom.
        Passing it to Molecule.show as atom_labels
        therefore draws the coordination edit straight onto the structure.
        These are not the element labels, which are read off the molecule with
        get_labels().

        To draw the site, hand Molecule.show the bonds that
        connectivity_bonds(active_site['connectivity_matrix']) returns rather
        than letting it perceive them by distance. That matters here for two
        reasons.

        :param topology:
            The protonated OpenMM topology.
        :param positions:
            The positions as an (N, 3) numpy array in Angstrom.
        :param binding_modes:
            The coordination, which says which residues the site holds.
        :param cap_bond_length:
            The C-H bond length in Angstrom of the capping hydrogens.

        :return:
            The active site dictionary. It records the active site indices of
            the capping hydrogens and of the beta carbons, the map back to
            the topology, the per-atom labels, the bonds and the charge.
        """

        self._check_supported_metals(binding_modes['metals'], 'extract_active_site')

        positions = np.asarray(positions)
        residues = list(topology.residues())

        res_indices = active_site_residues(binding_modes)

        # Check that all residues that coordinate a metal are included in the active site
        orphaned = sorted({
            ligand['residue']
            for ligand in binding_modes['ligands']
            if ligand['res_index'] not in res_indices
        })

        assert_msg_critical(
            not orphaned, 'extract_active_site: '
            f'{", ".join(orphaned)} coordinate(s) a metal but was excluded from '
            'the active site. Remove the metal bond before removing the residue')

        for res_index in res_indices:
            self.check_truncatable(residues[res_index])

        # Data for all active site elements
        labels = []
        atom_labels = []
        coords = []
        atom_map = {}
        cap_indices = []
        beta_carbon_indices = []
        metal_indices = []

        for metal in binding_modes['metals']:
            atom_map[len(coords)] = metal['index']
            metal_indices.append(len(coords))
            coords.append(positions[metal['index']])
            labels.append(metal['element'])
            atom_labels.append(str(metal['index']))

        for res_index in res_indices:
            residue = residues[res_index]
            res_atoms = list(residue.atoms())

            for atom in res_atoms:
                # Discard all backbone atoms
                if atom.name in BACKBONE_ATOM_NAMES:
                    continue

                # Discard the alpha carbon and replace it with a hydrogen at shorter bond length
                if atom.name == 'CA':
                    cb_atom = None
                    for other in res_atoms:
                        if other.name == 'CB':
                            cb_atom = other
                            break
                    assert_msg_critical(
                        cb_atom is not None, 'extract_active_site: '
                        f'residue {residue_label(residue)} has no CB to '
                        'cut at')
                    direction = positions[atom.index] - positions[cb_atom.index]
                    direction /= np.linalg.norm(direction)
                    # the cap is mapped to the CA it replaces, so that the
                    # CA-CB bond of the topology becomes the cap-CB bond of
                    # the active site
                    atom_map[len(coords)] = atom.index
                    cap_indices.append(len(coords))
                    coords.append(positions[cb_atom.index] +
                                  direction * cap_bond_length)
                    labels.append('H')
                    atom_labels.append('')

                # Include the rest of the atoms
                else:
                    if atom.name == 'CB':
                        beta_carbon_indices.append(len(coords))
                        # the beta carbon stands for its residue, so it is
                        # labelled with what add_metal_bond takes as resid
                        atom_labels.append(str(residue.id))
                    elif atom.element.symbol == 'H':
                        atom_labels.append('')
                    else:
                        atom_labels.append(atom.name)
                    atom_map[len(coords)] = atom.index
                    coords.append(positions[atom.index])
                    labels.append(atom.element.symbol)

        charge = sum(metal['formal_charge'] for metal in binding_modes['metals'])

        for res_index in res_indices:
            residue = residues[res_index]
            variant = binding_modes['variants'].get(res_index)
            if variant is None:
                variant = residue.name
                self.ostream.print_warning(
                    'No protonation variant recorded for '
                    f'{residue_label(residue)}; using the residue name for '
                    'the charge count')
            assert_msg_critical(
                variant in VARIANT_CHARGES, 'extract_active_site: no charge '
                f'known for variant {variant}')
            charge += VARIANT_CHARGES[variant]

        molecule = Molecule(labels, np.array(coords))
        molecule.set_charge(charge)
        molecule.set_multiplicity(1)

        # the charge and the multiplicity are on the molecule, which is
        # where they are read from; keeping a second copy here only gives
        # them somewhere to disagree
        active_site = {
            'molecule': molecule,
            'atom_map': atom_map,
            'cap_indices': cap_indices,
            'beta_carbon_indices': beta_carbon_indices,
            'metal_indices': metal_indices,
            'labels': atom_labels,
            'residues':
            [f'{residues[i].name}{residues[i].id}' for i in res_indices],
        }

        active_site['connectivity_matrix'] = self._build_connectivity(topology,
                                                                      active_site,
                                                                      binding_modes)

        return active_site

    @on_master
    def show_active_site(self, active_site, **kwargs):

        kwargs.setdefault('atom_labels', active_site['labels'])
        kwargs.setdefault('bonds',
                          connectivity_bonds(active_site['connectivity_matrix']))

        return active_site['molecule'].show(**kwargs)

    @on_master
    def _build_connectivity(self, topology, active_site, binding_modes):
        """
        Builds the connectivity matrix of the active site without perceiving any
        bonds.

        Covalent bonds come from the OpenMM topology, which is chemically
        correct by construction for standard residues, and the metal-ligand
        bonds come from the binding modes, which are explicit and reviewable.
        Molecule.get_connectivity_matrix is not used.

        Called by extract_active_site, which puts the result into the
        dictionary it returns.

        :param topology:
            The protonated OpenMM topology.
        :param active_site:
            The active site being extracted, for the atom map and the
            geometry.
        :param binding_modes:
            The coordination, for the metal-ligand bonds.

        :return:
            The connectivity matrix.
        """

        atom_map = active_site['atom_map']
        reverse_map = {
            top_index: site_index
            for site_index, top_index in atom_map.items()
        }

        n_atoms = len(atom_map)
        connectivity_matrix = np.zeros((n_atoms, n_atoms), dtype=int)

        for bond in topology.bonds():
            i, j = bond.atom1.index, bond.atom2.index
            if i in reverse_map and j in reverse_map:
                connectivity_matrix[reverse_map[i], reverse_map[j]] = 1
                connectivity_matrix[reverse_map[j], reverse_map[i]] = 1

        for ligand in binding_modes['ligands']:
            assert_msg_critical(
                ligand['index'] in reverse_map, '_build_connectivity: ligand atom '
                f'{ligand["residue"]} {ligand["atom"]} is not part of the '
                'extracted active site')
            lig_index = reverse_map[ligand['index']]
            for metal_index in ligand['metals']:
                metal = reverse_map[metal_index]
                connectivity_matrix[metal, lig_index] = 1
                connectivity_matrix[lig_index, metal] = 1

        # nothing but a metal bond should be long
        coords = active_site['molecule'].get_coordinates_in_angstrom()
        labels = active_site['molecule'].get_labels()
        metals = set(active_site['metal_indices'])

        for i in range(n_atoms):
            for j in range(i + 1, n_atoms):
                if not connectivity_matrix[i, j]:
                    continue
                if i in metals or j in metals:
                    continue
                # a covalent bond of a standard residue is never this long;
                # one that is says the structure is broken
                distance = np.linalg.norm(coords[i] - coords[j])
                if distance > 2.0:
                    self.ostream.print_warning(
                        f'Non-metal bond {i}-{j} ({labels[i]}-{labels[j]}) '
                        f'is {distance:.2f} A long')

        return connectivity_matrix

    @on_master
    def apply_metal_bonds(self, active_site, binding_modes):
        """
        Rewrites the metal-ligand bonds of an active site from the binding
        modes, leaving everything else alone.

        This is what a coordination edit made after the force field was built
        needs: the geometry, the charges and the Hessian all still belong to
        this site, and only which atoms the metals are bonded to has changed.
        Re-extracting would throw away an optimized geometry; re-detecting
        would read the distances the edit was made to overrule.

        The covalent bonds are untouched, and so is a metal-metal bond, which
        comes from the topology rather than from a ligand contact.

        :param active_site:
            The active site to rewrite. Not modified.
        :param binding_modes:
            The binding modes to read the coordination off, whose indices are
            those of the topology the active site was extracted from.

        :return:
            A new active site carrying the new connectivity, or the argument
            itself when nothing changed.
        """

        site_of = _site_index_map(active_site)
        metals = set(active_site['metal_indices'])

        matrix = np.array(active_site['connectivity_matrix'], copy=True)

        # only a metal-to-ligand entry is the binding modes' to say anything
        # about, so a metal-metal bond of the topology survives the rewrite
        for metal in metals:
            for other in range(matrix.shape[0]):
                if other in metals:
                    continue
                matrix[metal, other] = 0
                matrix[other, metal] = 0

        for ligand in binding_modes['ligands']:
            assert_msg_critical(
                ligand['index'] in site_of, 'apply_metal_bonds: '
                f'{ligand["residue"]} {ligand["atom"]} is not part of this '
                'active site, so it has no atom here to bond to a metal. Only a '
                'residue the cluster already holds can gain a metal bond once '
                'the force field is built; call build_active_site() again to '
                'extract the site afresh, which drops the fit')
            lig_index = site_of[ligand['index']]
            for metal_index in ligand['metals']:
                metal = site_of[metal_index]
                matrix[metal, lig_index] = 1
                matrix[lig_index, metal] = 1

        if np.array_equal(matrix, np.asarray(active_site['connectivity_matrix'])):
            return active_site

        new_active_site = dict(active_site)
        new_active_site['connectivity_matrix'] = matrix

        return new_active_site

    @on_master
    def connectivity_from_forcefield(self, active_site, forcefield):
        """
        Rewrites the metal-ligand bonds of an active site from a fitted force
        field, leaving everything else alone.

        This is what the weak bridge pruning needs: the fit is the only step
        that can decide a metal contact is not a bond after all, and it
        decides it on the force field, which it is handed and may edit. The
        active site it was built from is not, so without this the site goes
        on claiming a bond the force field no longer has -- and the site's
        matrix is what show_active_site draws, what the run summary counts,
        what a later Hessian walks its pairs out of, and what the manager
        describes a query site by, while a template is described by
        forcefield.bonds. That last disagreement is the expensive one: a
        structure would not reproduce the coarse topology of a template built
        from itself.

        The reverse of apply_metal_bonds and shaped the same way. The
        covalent bonds are untouched, and so is a metal-metal bond, which is
        not a ligand contact for the fit to have an opinion about.

        Note that the pruning cannot be undone from here: a rebuild derives
        the connectivity afresh and a rebuild is what drops the fit, which is
        the right way round. The one consequence worth knowing is that
        extract_pairs walked over a pruned matrix no longer covers the
        dropped pair, so a Hessian computed after a fit holds nothing for it.
        That is harmless -- the refit gives it no force constant and prunes
        it again -- but it is why the coverage is not a bug when it is
        missing.

        :param active_site:
            The active site to rewrite. Not modified.
        :param forcefield:
            The force field to read the metal bonds off.

        :return:
            A new active site carrying the new connectivity, or the argument
            itself when nothing changed.
        """

        metals = set(active_site['metal_indices'])

        matrix = np.array(active_site['connectivity_matrix'], copy=True)

        for metal in metals:
            for other in range(matrix.shape[0]):
                if other in metals:
                    continue
                matrix[metal, other] = 0
                matrix[other, metal] = 0

        for key in forcefield.bonds:
            first, second = key
            if not (metals & {first, second}):
                continue
            if first in metals and second in metals:
                continue
            matrix[first, second] = 1
            matrix[second, first] = 1

        if np.array_equal(matrix, np.asarray(active_site['connectivity_matrix'])):
            return active_site

        new_active_site = dict(active_site)
        new_active_site['connectivity_matrix'] = matrix

        return new_active_site

    @on_master
    def derive_site_coordination(self, topology, geometry, active_site,
                                 binding_modes, metal_bond_cutoff,
                                 report_cutoff):
        """
        Works out the coordination of an extracted active site from its own
        geometry.

        Restricted to the atoms the cluster holds, which is the difference
        between this and derive_binding_modes: once the site is extracted and
        a force field is keyed to it, a protein atom that drifts within the
        cutoff is not something the site can gain, so it is not looked at.

        The metals, the requests and the manual bond records are taken from
        the binding modes handed in, which is what makes deriving on the old
        geometry and on the new one comparable.

        :param topology:
            The protonated OpenMM topology, which the active site indexes into.
        :param geometry:
            The active site geometry, as a molecule or as an (N, 3) array in
            Angstrom. Ordered like the active site, not like the topology.
        :param active_site:
            The active site the geometry belongs to. Not modified.
        :param binding_modes:
            The modes to take the metals and the recorded decisions from.
        :param metal_bond_cutoff:
            The bonding cutoff in Angstrom.
        :param report_cutoff:
            The reporting cutoff in Angstrom.

        :return:
            The tuple of the ligand contacts and the notes deriving them
            produced.
        """

        if isinstance(geometry, Molecule):
            coordinates = geometry.get_coordinates_in_angstrom()
        else:
            coordinates = np.asarray(geometry, dtype=float)

        atom_map = active_site['atom_map']

        assert_msg_critical(
            len(coordinates) == len(atom_map),
            'derive_site_coordination: the geometry '
            f'has {len(coordinates)} atoms while the active site has '
            f'{len(atom_map)}')

        site_of = _site_index_map(active_site)

        atoms = list(topology.atoms())
        candidates = [atoms[top_index] for top_index in sorted(site_of)]

        def position_of(index):
            return coordinates[site_of[index]]

        notes = []
        ligands = self._collect_ligands(
            candidates, position_of, binding_modes['metals'], notes,
            set(binding_modes.get('coordinating_residues', [])),
            metal_bond_cutoff, report_cutoff)

        # a bond put there by hand is not something the distances can be
        # asked about again, so it is replayed onto every derivation
        records = binding_modes.get('manual_bonds', [])
        if records:
            self._apply_manual_bonds(ligands, records, binding_modes['metals'],
                                     candidates, position_of, notes)
            self._assign_binding_modes(ligands, notes,
                                       self._manual_protected(records))

        return ligands, notes

    @on_master
    def update_binding_modes(self, topology, geometry, active_site,
                             binding_modes, metal_bond_cutoff, report_cutoff):
        """
        Re-detects the coordination sphere on a new active site geometry.

        Relaxing the active site moves the metal-ligand distances: a contact
        that started just inside the primary cutoff can end up outside it, an
        asymmetric carboxylate can open into a monodentate one, and a second
        oxygen can rotate onto a metal. The rules of derive_binding_modes are
        applied again to the new coordinates so that what is fitted afterwards
        is the coordination the geometry actually has. Residues that were
        asked for through coordinating_residues stay ligands, since a distance is
        not what put them there.
        Nothing is modified in place, and the arguments themselves are returned
        when the coordination did not change, so the caller can overwrite what
        it holds unconditionally and still compare by identity.

        :param topology:
            The protonated OpenMM topology, which the active site indexes into.
        :param geometry:
            The new active site geometry, as a molecule or as an (N, 3) array
            in Angstrom. Ordered like the active site, not like the topology.
        :param active_site:
            The active site the geometry belongs to, whose connectivity is
            what gets checked. Not modified.
        :param binding_modes:
            The binding modes to check. Not modified.
        :param metal_bond_cutoff:
            The bonding cutoff in Angstrom.
        :param report_cutoff:
            The reporting cutoff in Angstrom.

        :return:
            The tuple of the binding modes, the active site and the flag that
            is True when the coordination changed. The active site is a new
            dictionary carrying the new connectivity when it did change, and
            the one that was passed in when it did not.
        """

        self._check_supported_metals(binding_modes['metals'], 'update_binding_modes')

        if isinstance(geometry, Molecule):
            coordinates = geometry.get_coordinates_in_angstrom()
        else:
            coordinates = np.asarray(geometry, dtype=float)

        site_of = _site_index_map(active_site)

        records = binding_modes.get('manual_bonds', [])

        new_ligands, notes = self.derive_site_coordination(
            topology, coordinates, active_site, binding_modes,
            metal_bond_cutoff, report_cutoff)

        def coordination(ligands):
            return {
                ligand['index']: (ligand['mode'], tuple(ligand['metals']))
                for ligand in ligands
            }

        old_coordination = coordination(binding_modes['ligands'])
        new_coordination = coordination(new_ligands)

        old_ligands = {
            ligand['index']: ligand
            for ligand in binding_modes['ligands']
        }
        new_by_index = {ligand['index']: ligand for ligand in new_ligands}

        # measured over the bonds the argument recorded, so that a contact
        # the update drops still contributes the distance it moved
        shifts = []
        for ligand in binding_modes['ligands']:
            assert_msg_critical(
                ligand['index'] in site_of, 'update_binding_modes: ligand atom '
                f'{ligand["residue"]} {ligand["atom"]} is not part of the '
                'active site the geometry belongs to')
            lig_index = site_of[ligand['index']]
            for metal_index, distance in zip(ligand['metals'], ligand['distances']):
                moved = float(
                    np.linalg.norm(coordinates[lig_index] -
                                   coordinates[site_of[metal_index]]))
                shifts.append(abs(moved - distance))

        largest_shift = max(shifts) if shifts else 0.0

        if new_coordination == old_coordination:
            self.print_binding_mode_update([], largest_shift, [])
            return binding_modes, active_site, False

        changes = []
        for index in sorted(set(old_coordination) | set(new_coordination)):
            old = old_ligands.get(index)
            new = new_by_index.get(index)

            if old is None:
                changes.append(('gained', f'{new["residue"]} {new["atom"]}',
                                f'{self._metal_contact_label(new, binding_modes)}, '
                                f'{new["mode"]}'))
            elif new is None:
                changes.append(('lost', f'{old["residue"]} {old["atom"]}', 'was '
                                f'{self._metal_contact_label(old, binding_modes)}, '
                                f'{old["mode"]}'))
            elif old_coordination[index] != new_coordination[index]:
                detail = []
                if old['mode'] != new['mode']:
                    detail.append(f'{old["mode"]} -> {new["mode"]}')
                if old['metals'] != new['metals']:
                    detail.append(self._metal_contact_label(new, binding_modes))
                changes.append(('changed', f'{new["residue"]} {new["atom"]}',
                                ', '.join(detail)))

        # a residue that lost every contact is still part of the truncated
        # active site, and only a new extraction can take it out
        dropped_residues = sorted({
            ligand['residue']
            for ligand in binding_modes['ligands']
            if ligand['index'] not in new_by_index
        } - {ligand['residue']
             for ligand in new_ligands})

        self.print_binding_mode_update(changes, largest_shift,
                                       dropped_residues)

        new_binding_modes = {
            'metals':
            deepcopy(binding_modes['metals']),
            'ligands':
            new_ligands,
            'variants':
            deepcopy(binding_modes.get('variants', {})),
            'coordinating_residues':
            sorted(binding_modes.get('coordinating_residues', [])),
            'extra_residues':
            sorted(binding_modes.get('extra_residues', [])),
            'excluded_residues':
            sorted(binding_modes.get('excluded_residues', [])),
            'manual_bonds':
            deepcopy(records),
            'notes':
            notes,
        }

        # only the metal-ligand bonds are read off the geometry, so the
        # covalent bonds of the matrix are left exactly as they were
        new_matrix = np.array(active_site['connectivity_matrix'], copy=True)

        for ligands, value in ((binding_modes['ligands'], 0), (new_ligands, 1)):
            for ligand in ligands:
                lig_index = site_of[ligand['index']]
                for metal_index in ligand['metals']:
                    metal = site_of[metal_index]
                    new_matrix[metal, lig_index] = value
                    new_matrix[lig_index, metal] = value

        self.print_binding_modes(new_binding_modes)

        new_active_site = dict(active_site)
        new_active_site['connectivity_matrix'] = new_matrix

        return new_binding_modes, new_active_site, True

    @on_master
    def _check_supported_metals(self, metals, method):
        """
        Rejects metal centers the builder is not validated for.

        The literature distances, the assumed formal charges and the
        coordination rules have only been checked against zinc, so a site
        built around any other metal would be produced with zinc's assumptions
        silently applied to it.

        :param metals:
            The list of metal entries of the binding modes.
        :param method:
            The name of the calling method, for the error message.
        """

        found = sorted({metal['element'] for metal in metals})
        unsupported = [
            element for element in found if element not in SUPPORTED_METAL_ELEMENTS
        ]

        assert_msg_critical(
            not unsupported, f'{method}: found {unsupported}, but '
            f'only {list(SUPPORTED_METAL_ELEMENTS)} is supported. The '
            'literature distances, formal charges and coordination rules have '
            'only been validated for zinc.')

    @on_master
    def _metal_contact_label(self, ligand, binding_modes):
        """
        Names the metals a contact reaches and how far away they are.

        :param ligand:
            The ligand contact.
        :param binding_modes:
            The binding modes, for the elements of the metals.

        :return:
            The formatted contact, such as 'Zn4 2.11 A'.
        """

        elements = {
            metal['index']: metal['element']
            for metal in binding_modes['metals']
        }

        return ', '.join(
            f'{elements.get(index, "metal")}{index} {distance:.2f} A'
            for index, distance in zip(ligand['metals'], ligand['distances']))

    # ------------------------------------------------------------------
    # the structural pass
    # ------------------------------------------------------------------

    @on_master
    def build_active_site(self,
                          topology,
                          positions,
                          request,
                          protonation_overrides,
                          report_detection,
                          cap_bond_length,
                          metal_elements,
                          metal_formal_charges,
                          metal_bond_cutoff,
                          report_cutoff,
                          coordinating_residues=None):
        """
        Runs the structural pass from a prepared structure to a truncated
        active site: detects the coordination, protonates the structure to
        match, detects again on what that returns and truncates the cluster.

        The detection runs twice on purpose. protonate renumbers the atoms
        and remaps nothing, so the coordination is simply derived again on
        the protonated topology and every atom index belongs to the topology
        it was derived on. That is safe because addHydrogens only adds
        hydrogens, which are not donors, so the same contacts are found
        either way.

        :param topology:
            The prepared topology.
        :param positions:
            Its positions in Angstrom.
        :param request:
            The record of what has been decided about the site. Not
            modified; the one returned carries the variants the protonation
            built and the residues coordinating_residues forced.
        :param protonation_overrides:
            Protonation variants by residue label, as
            update_protonation_state records them; None for none.
        :param report_detection:
            Whether the coordination found on the structure is printed. The
            truncated site always is.
        :param cap_bond_length:
            The C-H distance in Angstrom the capping hydrogens are placed at.
        :param metal_elements:
            The elements treated as metal centers.
        :param metal_formal_charges:
            The formal charges assumed for the metal ions, by element.
        :param metal_bond_cutoff:
            The bonding cutoff in Angstrom of the detection.
        :param report_cutoff:
            The reporting cutoff in Angstrom of the detection.
        :param coordinating_residues:
            Residues to make ligands whatever their distance, as ids or
            labels. Recorded into the returned request, so every later
            derivation honours them.

        :return:
            A dictionary holding the request ('request'), the protonated
            topology and positions ('protonated_topology',
            'protonated_positions'), the binding modes derived on them
            ('binding_modes') and the active site ('active_site').
        """

        detection = {
            'metal_elements': metal_elements,
            'metal_formal_charges': metal_formal_charges,
            'metal_bond_cutoff': metal_bond_cutoff,
            'report_cutoff': report_cutoff,
        }

        binding_modes = self.derive_binding_modes(
            topology,
            positions,
            request,
            coordinating_residues=coordinating_residues,
            **detection)

        if report_detection:
            self.print_binding_modes(binding_modes)

        # what was forced is a decision, so it is recorded like every other
        request = deepcopy(request)
        request['coordinating_residues'] = list(
            binding_modes['coordinating_residues'])

        protonated_topology, protonated_positions, variants, notes = (
            self.protonate(topology, positions, binding_modes,
                           protonation_overrides))

        request['variants'] = variants
        protonated_modes = self.derive_binding_modes(protonated_topology,
                                                     protonated_positions,
                                                     request, **detection)

        existing = protonated_modes.setdefault('notes', [])
        for note in notes:
            if note not in existing:
                existing.append(note)

        active_site = self.extract_active_site(protonated_topology,
                                               protonated_positions,
                                               protonated_modes,
                                               cap_bond_length)
        self.print_active_site(active_site, protonated_modes)

        return {
            'request': request,
            'protonated_topology': protonated_topology,
            'protonated_positions': protonated_positions,
            'binding_modes': protonated_modes,
            'active_site': active_site,
        }

    @on_master
    def site_coordination(self, topology, positions, geometry, active_site,
                          request, metal_elements, metal_formal_charges,
                          metal_bond_cutoff, report_cutoff):
        """
        The coordination of an extracted site under one of its geometries.

        Restricted to the atoms the cluster holds, unlike derive_binding_modes
        on the whole structure: once a force field is keyed to the site, a
        protein atom drifting into range is not something the site can gain.

        :param topology:
            The protonated topology the site was extracted from.
        :param positions:
            Its positions in Angstrom.
        :param geometry:
            The geometry to read, ordered like the active site.
        :param active_site:
            The active site.
        :param request:
            The record of what has been decided about the site.
        :param metal_elements:
            The elements treated as metal centers.
        :param metal_formal_charges:
            The formal charges assumed for the metal ions, by element.
        :param metal_bond_cutoff:
            The bonding cutoff in Angstrom.
        :param report_cutoff:
            The reporting cutoff in Angstrom.

        :return:
            The binding modes of the site under that geometry.
        """

        modes = self.derive_binding_modes(topology, positions, request,
                                          metal_elements, metal_formal_charges,
                                          metal_bond_cutoff, report_cutoff)
        ligands, notes = self.derive_site_coordination(
            topology, geometry, active_site, modes, metal_bond_cutoff,
            report_cutoff)

        modes = dict(modes)
        modes['ligands'] = ligands
        modes['notes'] = notes

        return modes

    @on_master
    def adopt_geometry(self, topology, positions, active_site, molecule,
                       request, metal_elements, metal_formal_charges,
                       metal_bond_cutoff, report_cutoff):
        """
        Puts a new geometry on an active site and detects the coordination
        again on it.

        :param topology:
            The protonated topology the site was extracted from.
        :param positions:
            Its positions in Angstrom.
        :param active_site:
            The active site. Not modified.
        :param molecule:
            The new geometry of the active site.
        :param request:
            The record of what has been decided about the site.
        :param metal_elements:
            The elements treated as metal centers.
        :param metal_formal_charges:
            The formal charges assumed for the metal ions, by element.
        :param metal_bond_cutoff:
            The bonding cutoff in Angstrom of the re-detection.
        :param report_cutoff:
            The reporting cutoff in Angstrom of the re-detection.

        :return:
            The active site carrying the new geometry, with its connectivity
            brought up to date.
        """

        before = self.site_coordination(topology, positions,
                                        active_site['molecule'], active_site,
                                        request, metal_elements,
                                        metal_formal_charges,
                                        metal_bond_cutoff, report_cutoff)

        moved = dict(active_site)
        moved['molecule'] = molecule

        _, moved, _ = self.update_binding_modes(topology, molecule, moved,
                                                before, metal_bond_cutoff,
                                                report_cutoff)

        return moved

    # ------------------------------------------------------------------
    # geometry
    # ------------------------------------------------------------------

    @on_master
    def _add_untemplated_impropers(self, system, forcefield):
        """
        Puts the impropers a residue template cannot express onto a system.

        OpenMM builds impropers out of a residue template by walking every atom
        with three or more bonds and taking combinations of what it is bonded
        to, then matching the template's first class against that central atom.
        An improper naming four atoms of any other shape is therefore written
        into the XML and then never matched -- no error, no term.

        The one a bidentate carboxylate gets is exactly that shape: its central
        carbon is bonded to its two oxygens and to the beta carbon, not to the
        metal the improper is meant to hold it planar with. Left to the
        template it is silently absent from the relaxation the site is meant to
        be doing it under, so it is added here by index, where nothing has to be
        bonded to anything.

        Which ones those are is read off forcefield.bonds, which is where the
        residue template's own bonds come from, so this asks the same question
        of the same graph OpenMM does rather than guessing at its answer.

        :param system:
            The OpenMM system built from the force field's own XML.
        :param forcefield:
            The force field generator the system was built from.

        :return:
            The improper keys that were added.
        """

        torsion_force = None
        for force in system.getForces():
            if isinstance(force, mm.PeriodicTorsionForce):
                torsion_force = force

        bonded = {}
        for first, second in forcefield.bonds:
            bonded.setdefault(first, set()).add(second)
            bonded.setdefault(second, set()).add(first)

        added = []

        for key in forcefield.impropers:
            central = key[0]
            if set(key[1:]) <= bonded.get(central, set()):
                continue

            assert_msg_critical(
                torsion_force is not None,
                '_add_untemplated_impropers: the force field carries an improper '
                'no residue template can express, but the system it built has no '
                'periodic torsion force to add it to')

            params = forcefield.impropers[key]
            # the central atom goes third, which is where OpenMM puts it when it
            # does build an improper itself
            torsion_force.addTorsion(key[1], key[2], central, key[3],
                                     params['periodicity'],
                                     params['phase'] * mmunit.degree,
                                     params['barrier'] * mmunit.kilojoule_per_mole)
            added.append(key)

        return added

    @on_master
    def minimize_active_site(self, active_site, forcefield, frozen_indices):
        """
        Minimizes the active site with its own force field.

        The caller says what is frozen, normally the beta carbons through
        constrained_indices. They are where the backbone holds the sidechains
        in place, so a free minimization of an isolated active site lets the
        truncated fragments drift apart and says nothing about the metal
        site. Note that the force field only carries electrostatics if
        build_forcefield was given partial charges.

        :param active_site:
            The extracted active site.
        :param forcefield:
            The force field to minimize on.
        :param frozen_indices:
            The active site indices to hold fixed; an empty list minimizes
            freely.

        :return:
            The minimized coordinates as an (N, 3) numpy array in Angstrom.
            The coordinates make a round trip through a PDB file, so they
            carry a rounding of 0.001 Angstrom.
        """

        assert_msg_critical('openmm' in sys.modules,
                            'minimize_active_site: openmm is '
                            'required')

        with tempfile.TemporaryDirectory() as temp_dir:
            stem = str(Path(temp_dir) / 'active_site')
            forcefield.write_openmm_files(stem)

            pdb = mmapp.PDBFile(f'{stem}.pdb')
            openmm_ff = mmapp.ForceField(f'{stem}.xml')
            system = openmm_ff.createSystem(pdb.topology,
                                            nonbondedMethod=mmapp.NoCutoff)

            # the atoms are written and read back in force field order, so a
            # force field index is a particle index here
            self._add_untemplated_impropers(system, forcefield)

            # a zero mass makes OpenMM hold the particle fixed
            for index in frozen_indices:
                system.setParticleMass(int(index), 0.0)

            integrator = mm.VerletIntegrator(0.001 * mmunit.picoseconds)
            simulation = mmapp.Simulation(pdb.topology, system, integrator)
            simulation.context.setPositions(pdb.positions)
            simulation.minimizeEnergy()

            state = simulation.context.getState(getPositions=True)
            coords = np.array(state.getPositions().value_in_unit(mmunit.angstrom))

        return coords

    @on_master
    def mm_optimize_active_site(self, active_site, forcefield,
                                constrain_metals, constrain_capping_hydrogens,
                                bond_change_warning):
        """
        Relaxes the active site on a crude force field of its own.

        A site comes out of a structure file with its ligands wherever the
        crystallographer or the design left them, and the constrained QM
        optimization then spends its first cycles cleaning that up. This pass
        does the same work at MM cost: the metal terms are the seeded ones
        build_forcefield puts on a force field built without a Hessian, every
        other term is what the generator assigns, and the beta carbons are
        frozen exactly as they are in the QM optimization.

        The force field is taken rather than built, so which one the pass runs
        on is the caller's decision: the crude seeded one for a run, or one
        seeded against the literature distances for a comparison. It carries no
        electrostatics unless it was given charges. That is what makes the pass
        crude, and why it stands in front of optimize_active_site rather than
        in place of it.

        :param active_site:
            The extracted active site.
        :param forcefield:
            The force field to relax on, as built by build_forcefield without a
            Hessian.
        :param constrain_metals:
            Whether to hold the metal centers as well as the beta carbons.
        :param constrain_capping_hydrogens:
            Whether the capping hydrogens are frozen along with the beta
            carbons.
        :param bond_change_warning:
            How far a metal-ligand bond may move, in Angstrom, before it is
            reported.

        :return:
            The relaxed molecule. The caller decides whether to put it back
            into the active site.
        """

        molecule = active_site['molecule']

        frozen_indices = constrained_indices(active_site,
                                             constrain_capping_hydrogens)

        if constrain_metals:
            frozen_indices = sorted(
                set(frozen_indices) | set(active_site['metal_indices']))

        coordinates = self.minimize_active_site(active_site, forcefield,
                                                frozen_indices)

        relaxed = Molecule(molecule.get_labels(), coordinates, 'angstrom')
        relaxed.set_charge(molecule.get_charge())
        relaxed.set_multiplicity(molecule.get_multiplicity())

        self.print_mm_optimization(active_site,
                                   forcefield,
                                   relaxed,
                                   frozen_indices,
                                   get_metal_keys(forcefield, active_site),
                                   SEEDED_EQUILIBRIUM_LABELS,
                                   bond_change_warning)

        return relaxed

    # ------------------------------------------------------------------
    # the seeded force field
    # ------------------------------------------------------------------

    @on_master
    def _manual_bond_records_by_key(self, active_site, topology, binding_modes):
        """
        Maps each hand-added metal bond onto the active site key it names.

        The records name their atoms the way they have to in order to survive
        the renumbering of protonate - by residue index, atom name and the
        residue index of the metal - so turning them into keys of the force
        field means going back through the atom map of the active site. Doing
        that once here is what keeps manual_bond_keys and
        manual_bond_equilibria reading the same bonds.

        :param active_site:
            The active site the keys are to index into.
        :param topology:
            The protonated topology the active site was extracted from.
        :param binding_modes:
            The binding modes carrying the records.

        :return:
            A dictionary of bond key to the record that asked for it. A record
            whose residue is no longer part of this active site is left out.
        """

        records = [
            record for record in binding_modes.get('manual_bonds', [])
            if record['action'] == 'add'
        ]

        if not records:
            return {}

        atoms = list(topology.atoms())
        site_of = _site_index_map(active_site)

        by_key = {}
        for top_index in site_of:
            atom = atoms[top_index]
            by_key[(atom.residue.index, atom.name)] = top_index

        metal_by_res = {
            atoms[metal['index']].residue.index: metal['index']
            for metal in binding_modes['metals']
        }

        found = {}
        for record in records:
            top_atom = by_key.get((record['res_index'], record['atom']))
            top_metal = metal_by_res.get(record['metal_res_index'])

            if top_atom is None or top_metal is None:
                # the residue is no longer part of this active site, so
                # there is no bond of this force field to speak about
                continue

            key = tuple(sorted((site_of[top_atom], site_of[top_metal])))
            found[key] = record

        return found

    @on_master
    def manual_bond_keys(self, active_site, topology, binding_modes):
        """
        The active site bond keys of the metal bonds that were asked for by
        hand.

        :param active_site:
            The active site the keys are to index into.
        :param topology:
            The protonated topology the active site was extracted from.
        :param binding_modes:
            The binding modes carrying the records.

        :return:
            The bond keys, as a set of sorted index pairs.
        """

        return set(self._manual_bond_records_by_key(active_site, topology,
                                                    binding_modes))

    @on_master
    def manual_bond_equilibria(self, active_site, topology, binding_modes):
        """
        The equilibrium distances asked for alongside a hand-added metal bond.

        A bond added by hand is often added because the geometry has the
        residue turned the wrong way round, and then measuring the equilibrium
        on that geometry only pins the mistake in place. Giving a distance
        makes the crude pass pull the contact to it instead, which is enough to
        swing a histidine round before any QM is paid for.

        Only the crude pass reads these. Once a Hessian exists the equilibria
        come from the optimized geometry through Seminario, which is the whole
        point of fitting.

        :param active_site:
            The active site the keys are to index into.
        :param topology:
            The protonated topology the active site was extracted from.
        :param binding_modes:
            The binding modes carrying the records.

        :return:
            A dictionary of bond key to equilibrium distance in nanometers,
            holding only the bonds that were given one.
        """

        return {
            key: record['equilibrium']
            for key, record in self._manual_bond_records_by_key(active_site, topology,
                                                                binding_modes).items()
            if record.get('equilibrium') is not None
        }

    @on_master
    def _add_metal_planarity_impropers(self, forcefield, active_site,
                                       force_constant):
        """
        Adds a weak improper nudging each metal into the plane of a
        coordinating histidine ring or a bidentate carboxylate.

        An imidazole nitrogen's coordinating lone pair, and a carboxylate's
        two, both sit in the plane of the group they belong to, so a metal
        they grip belongs there too. Read straight off
        forcefield.connectivity_matrix -- the final metal-ligand bonding,
        already pruned of any weak bridge arm when a Hessian was fitted -- so a
        mu-1,3 carboxylate (one oxygen per metal) never qualifies as bidentate,
        and a histidine reaching a metal through two ring nitrogens, which
        _assign_binding_modes already trims to one, is never seen with both.

        The central atom of each improper key is the donor heavy atom itself
        (the ring nitrogen, or the carboxylate carbon); the ordering of the
        other three does not matter here, since with periodicity 2 and phase
        180 degrees the energy is invariant to exchanging any two of them.

        A coordinating, deprotonated ring nitrogen has only two ring carbons
        as covalent neighbours of its own, so the metal bond is what fills the
        third slot a substituent hydrogen would otherwise occupy -- and GAFF's
        own generic sp2-planarity improper (populate_impropers) already picks
        it up as that slot's atom, over the same four atoms in some other
        order. Adding a second, differently-ordered term over the same set
        would double the effective restraint rather than set it, so any
        existing improper covering the same four atoms is replaced rather than
        added to.

        :param forcefield:
            The force field generator to add the impropers to. Modified in
            place.
        :param active_site:
            The active site, for the metal indices and the elements.
        :param force_constant:
            The improper barrier, in kJ/mol. Deliberately weak; see
            DEFAULT_METAL_PLANARITY_FORCE_CONSTANT.

        :return:
            The number of impropers added, for the caller to report.
        """

        matrix = forcefield.connectivity_matrix
        elements = active_site['molecule'].get_labels()
        metals = set(active_site['metal_indices'])

        def neighbors(index):
            return set(np.where(matrix[index])[0].tolist()) - {index}

        # the impropers indexed by the atom set they cover, so that install
        # does not walk the whole table once per restraint it places
        by_atoms = {}
        for existing in forcefield.impropers:
            if len(existing) == 4:
                by_atoms.setdefault(frozenset(existing), []).append(existing)

        def install(key):
            # a term already covering the same four atoms in some other order
            # -- GAFF's own generic sp2-planarity guess can produce exactly
            # this for a coordinating ring nitrogen -- is replaced rather than
            # layered under a second one, so the restraint this atom set gets
            # is the one force_constant names, not the sum of two
            target = frozenset(key)
            for existing in by_atoms.get(target, []):
                del forcefield.impropers[existing]
            by_atoms[target] = [key]
            forcefield.impropers[key] = {
                'type': 'Fourier',
                'barrier': force_constant,
                'phase': 180.0,
                'periodicity': 2,
                'comment': 'metal coordination planarity restraint',
            }

        added = 0

        for metal in sorted(metals):
            ligand_atoms = neighbors(metal) - metals

            # coordinating histidine ring nitrogen: bonded to the metal and to
            # exactly two other heavy atoms, its ring carbons
            for donor in sorted(ligand_atoms):
                if elements[donor] != 'N':
                    continue
                ring_neighbors = sorted(atom for atom in neighbors(donor) - metals
                                        if elements[atom] != 'H')
                if len(ring_neighbors) != 2:
                    continue
                install((donor, ring_neighbors[0], ring_neighbors[1], metal))
                added += 1

            # bidentate carboxylate: two oxygens on this metal sharing one
            # covalently bonded carbon
            oxygens = sorted(atom for atom in ligand_atoms if elements[atom] == 'O')
            for i, oxygen_1 in enumerate(oxygens):
                for oxygen_2 in oxygens[i + 1:]:
                    shared = (neighbors(oxygen_1) & neighbors(oxygen_2)) - metals
                    for carbon in sorted(shared):
                        if elements[carbon] != 'C':
                            continue
                        install((carbon, oxygen_1, oxygen_2, metal))
                        added += 1

        if added:
            self.ostream.print_info(
                f'Added {added} weak metal coordination planarity improper(s) '
                f'at {force_constant:.2f} kJ/mol.')
            self.ostream.flush()

        return added

    @on_master
    def build_forcefield(self,
                         active_site,
                         partial_charges,
                         metal_blind_typing,
                         reparameterize_metal_angles,
                         default_metal_angle_force_constant,
                         default_metal_bond_force_constant,
                         metal_angle_equilibria,
                         metal_bond_equilibria,
                         add_metal_planarity_impropers,
                         metal_planarity_force_constant,
                         mute_generator,
                         bond_equilibria=None):
        """
        Builds the active site force field with seeded metal terms.

        The metal terms are seeded by _seed_metal_terms rather than fitted:
        equilibria measured on the geometry unless a table or a request
        overrides them, and a flat default stiffness. That is the force field
        the crude pre-QM pass runs on, where getting the equilibrium geometry
        roughly right matters far more than the stiffness, and it is what
        QmParameterizer.fit_forcefield fits the metal terms of once there is a
        Hessian.

        :param active_site:
            The extracted active site.
        :param partial_charges:
            The partial charges of the active site. The charge of the capping
            hydrogens is redistributed over the remaining atoms before they are
            applied, since the caps do not exist in the protein. util.d4_charges
            is the cheap choice when none were fitted.
        :param metal_blind_typing:
            Whether the generator perceives the atom types as if the metal
            bonds were not there, so that a coordinating residue is typed as
            the amino acid it is. The covalent terms of the residues are what
            switching it off loses; a site whose ligands are genuinely not
            amino acids may want that.
        :param reparameterize_metal_angles:
            Whether the metal angles are seeded at all, or left at what the
            generator guessed.
        :param default_metal_angle_force_constant:
            The flat stiffness of every seeded metal angle, in kJ/mol/rad^2.
        :param default_metal_bond_force_constant:
            The flat stiffness of every seeded metal bond, in kJ/mol/nm^2.
        :param metal_angle_equilibria:
            Equilibrium angles by element triple, overriding the measured
            ones; None for none.
        :param metal_bond_equilibria:
            Equilibrium distances in nanometers by element pair, overriding
            the measured ones; None for none. LITERATURE_METAL_BONDS is ready
            to be assigned here.
        :param add_metal_planarity_impropers:
            Whether to add a weak improper nudging each metal into the plane
            of a coordinating histidine ring or a bidentate carboxylate; see
            _add_metal_planarity_impropers.
        :param metal_planarity_force_constant:
            The barrier of that improper, in kJ/mol.
        :param mute_generator:
            Whether the generator's own reporting is silenced. It names every
            parameter it looks up and every bond and angle it re-measures,
            which buries what this module has to say about the site -- and a
            shoehorning, which rebuilds the site after every edit, prints it
            all again each time. Set it False to see what GAFF did.
        :param bond_equilibria:
            Equilibrium distances in nanometers for individual metal bonds,
            from manual_bond_equilibria. Read by the crude pass alone.

        :return:
            The force field generator.
        """

        molecule = active_site['molecule']
        n_atoms = molecule.number_of_atoms()

        # The generator reports every parameter it looks up and every bond and
        # angle it re-measures, which for a site of this size is hundreds of
        # lines saying only that GAFF was used. It is given a stream of its own
        # rather than the shared one being muted, since OutputStream.mute is
        # reference counted and an unbalanced pair silences everything after it.
        forcefield = MMForceFieldGenerator(
            MPI.COMM_SELF, OutputStream(None) if mute_generator else self.ostream)
        # copied, not aliased: np.asarray hands back the caller's own array, and
        # the weak bridge pruning edits the generator's matrix, so sharing it
        # would have this function quietly rewriting the active site it was
        # given -- which no method here may do
        forcefield.connectivity_matrix = np.array(
            active_site['connectivity_matrix'], copy=True)
        forcefield.topology_update_flag = True
        # A metal bond is a bond like any other to the GAFF perception, so a
        # carboxylate oxygen gripping a metal is typed as an ether one and the
        # covalent terms around it lose their parameters. The bonds are still
        # made -- only the typing looks past them.
        forcefield.metal_blind_typing = metal_blind_typing

        forcefield.create_topology(molecule, resp=False)

        # The charges have to be applied after create_topology: setting
        # topology_update_flag, which is what makes the custom connectivity
        # take effect, also resets partial_charges to None, and resp=False
        # then fills them with zeros. Assigning them beforehand is silently
        # discarded.
        partial_charges = np.asarray(partial_charges)
        assert_msg_critical(
            partial_charges.shape == (n_atoms, ), 'build_forcefield: expected '
            f'{n_atoms} partial charges, got {partial_charges.shape}')
        # The capping hydrogens stand in for alpha carbons and do not
        # exist anywhere the force field is used, so their charge is
        # folded into the rest of the active site before anything is written.
        partial_charges = redistribute_cap_charges(active_site, partial_charges)
        forcefield.partial_charges = partial_charges
        for index in range(n_atoms):
            forcefield.atoms[index]['charge'] = partial_charges[index]

        self.annotate_atoms(forcefield, active_site)

        bonds, angles = get_metal_keys(forcefield, active_site)

        # switching the angles off leaves them at whatever the generator
        # guessed, here and in the fit
        if not reparameterize_metal_angles:
            angles = []

        self._seed_metal_terms(forcefield, active_site, bonds, angles,
                               default_metal_angle_force_constant,
                               default_metal_bond_force_constant,
                               metal_angle_equilibria, metal_bond_equilibria,
                               bond_equilibria)

        # A fit that later prunes a weak bridge arm takes every improper across
        # the dropped bond with it, so one added here cannot outlive its bond.
        if add_metal_planarity_impropers:
            self._add_metal_planarity_impropers(
                forcefield, active_site, metal_planarity_force_constant)

        # the typing that puts a term here is done by create_topology, and a
        # fit of the metal terms changes none of it
        self._print_check_atom_types(forcefield, active_site,
                                     metal_blind_typing)

        return forcefield

    @on_master
    def annotate_atoms(self, forcefield, active_site):
        """
        Writes what the truncation knows about an atom into its comment.

        A force field carries atom types, charges and bonds, and a geometry
        carries positions; neither of them says which carbon a sidechain was
        cut at, or which hydrogen stands in for an alpha carbon. That is known
        here, while the active site is still in hand, and the comment is the
        one field that survives into forcefield.json for anything downstream
        to read it back out of.

        The capping hydrogens do end up at exactly zero charge, but that is
        the charge correction doing its job rather than a marker, and it says
        nothing in a force field whose charges are all zero. The comment is
        the marker.

        Comments the generator already wrote are kept, and running this twice
        does not write the same note twice.

        :param forcefield:
            The force field generator to annotate.
        :param active_site:
            The active site, for the atoms the truncation created.
        """

        roles = (
            (BETA_CARBON_COMMENT, active_site['beta_carbon_indices']),
            (CAP_COMMENT, active_site['cap_indices']),
        )

        for note, indices in roles:
            for index in indices:
                atom = forcefield.atoms[index]
                comment = atom.get('comment', '') or ''

                if note in comment:
                    continue

                atom['comment'] = '; '.join(part for part in (comment, note)
                                            if part)

    @on_master
    def _lookup_equilibrium(self, table, elements):
        """
        Looks an element combination up in an equilibrium table.

        A bond and an angle read the same forwards and backwards, so a table
        only has to carry one of the two orders.

        :param table:
            The table, or None.
        :param elements:
            The element symbols of the term, in key order.

        :return:
            The equilibrium value, or None when the table does not hold it.
        """

        if not table:
            return None

        if elements in table:
            return table[elements]

        return table.get(elements[::-1])

    @on_master
    def _seed_metal_terms(self, forcefield, active_site, bonds, angles,
                          default_metal_angle_force_constant,
                          default_metal_bond_force_constant,
                          metal_angle_equilibria, metal_bond_equilibria,
                          bond_equilibria):
        """
        Seeds the metal terms for the crude MM pass.

        The equilibrium values are measured on the active site as the
        structure file gave it, which before any QM is run is the only
        description of the site there is. metal_bond_equilibria and
        metal_angle_equilibria override that per element combination;
        LITERATURE_METAL_BONDS is ready to be assigned to the first of them.
        The force constants are flat defaults, since nothing at this stage
        says anything about the stiffness of a metal term.

        :param bonds:
            The metal bond keys.
        :param angles:
            The metal angle keys. Empty when reparameterize_metal_angles is
            switched off.
        :param default_metal_angle_force_constant:
            The flat stiffness of every metal angle, in kJ/mol/rad^2.
        :param default_metal_bond_force_constant:
            The flat stiffness of every metal bond, in kJ/mol/nm^2.
        :param metal_angle_equilibria:
            Equilibrium angles by element triple, or None.
        :param metal_bond_equilibria:
            Equilibrium distances in nanometers by element pair, or None.
        :param bond_equilibria:
            Equilibrium distances in nanometers for individual bond keys. These
            beat both the measurement and the element table, since they name one
            bond rather than a kind of bond: manual_bond_equilibria turns what
            add_metal_bond was told into them.
        """

        labels = active_site['molecule'].get_labels()
        molecule = active_site['molecule']

        if bond_equilibria is None:
            bond_equilibria = {}

        for key in bonds:
            elements = tuple(labels[index] for index in key)
            equilibrium = bond_equilibria.get(tuple(sorted(key)))
            comment = SEEDED_FROM_REQUEST

            if equilibrium is None:
                equilibrium = self._lookup_equilibrium(metal_bond_equilibria, elements)
                comment = SEEDED_FROM_TABLE

            if equilibrium is None:
                # the getters index atoms from one, and the force field keeps
                # its bond lengths in nanometers
                equilibrium = 0.1 * molecule.get_distance_in_angstroms(
                    [index + 1 for index in key])
                comment = SEEDED_FROM_GEOMETRY
            forcefield.bonds[key]['equilibrium'] = equilibrium
            forcefield.bonds[key]['force_constant'] = (
                default_metal_bond_force_constant)
            forcefield.bonds[key]['comment'] = comment

        for key in angles:
            elements = tuple(labels[index] for index in key)
            equilibrium = self._lookup_equilibrium(metal_angle_equilibria, elements)
            if equilibrium is None:
                equilibrium = molecule.get_angle_in_degrees(
                    [index + 1 for index in key])
                comment = SEEDED_FROM_GEOMETRY
            else:
                comment = SEEDED_FROM_TABLE
            forcefield.angles[key]['equilibrium'] = equilibrium
            forcefield.angles[key]['force_constant'] = (
                default_metal_angle_force_constant)
            forcefield.angles[key]['comment'] = comment

        self.ostream.flush()

    @on_master
    def _print_check_atom_types(self, forcefield, active_site,
                                metal_blind_typing):
        """
        Warns about covalent terms left without parameters.

        GAFF types an atom from what it is bonded to, and a metal bond is a
        bond like any other to that perception. A carboxylate oxygen that
        grips a metal stops looking like a carbonyl oxygen and is typed as an
        ether one, which costs its carbon the carbonyl type in turn, and the
        combinations that leaves behind are not always in the parameter set.
        What is missing falls back to a flat constant -- 2.5e5 kJ/mol/nm^2 for
        a bond, 1000 kJ/mol/rad^2 for an angle -- several times the tabulated
        values it stands in for, on terms that carry hydrogen bending modes.

        metal_blind_typing is what stops that happening, and is on by
        default, so what this reports depends on it: with it on, a term in
        the ligand shell is a gap in the parameter set like any other, and
        with it off it is the coordination that put it there.

        Nothing downstream repairs these. reparameterize would fit a guessed
        term by default, but the fit here is restricted to the metal keys and
        the Hessian behind it to the pairs of the coordination, so a covalent
        term left flat stays flat.

        Only bonds and angles are looked at, which is the same line
        reparameterize draws when it picks its own terms: their fallbacks
        stand several times off the tabulated values they replace, while a
        guessed improper gets the generic 1.1 kcal/mol every force field uses
        and says nothing about the typing. Impropers also come out guessed in
        roughly the same number whatever the coordination does, so counting
        them would bury the part of the report that varies.

        The report is split by how far a term sits from the metal, because
        without the blind typing that is what separates the two cures. Inside
        the ligand shell -- a term holding a ligand atom or a neighbour of one
        -- the typing follows from the coordination, and it is the
        coordination that has to change. Outside it the gap is in the
        parameter set and has nothing to do with the metal.

        Terms are grouped by their atom types rather than listed one by one,
        since the types are what names the gap and a site can hold a dozen
        terms saying the same thing.

        :param forcefield:
            The force field generator, as it will be handed back.
        :param active_site:
            The active site, for the indices of the metal centers.
        :param metal_blind_typing:
            Whether the types were perceived with the metal bonds ignored,
            which decides what a term in the ligand shell is evidence of.
        """

        metals = set(active_site['metal_indices'])
        atoms = forcefield.atoms

        neighbours = {}
        for first, second in forcefield.bonds:
            neighbours.setdefault(first, set()).add(second)
            neighbours.setdefault(second, set()).add(first)

        ligands = {
            index
            for metal in metals
            for index in neighbours.get(metal, ())
        } - metals
        shell = (
            ligands
            | {index
               for ligand in ligands
               for index in neighbours.get(ligand, ())}) - metals

        def flat_keys(terms):
            return [
                key for key, term in terms.items()
                if term.get('comment') == UNPARAMETERIZED_COMMENT and not metals
                & set(key)
            ]

        flat = [(key, 'bond') for key in flat_keys(forcefield.bonds)]
        flat += [(key, 'angle') for key in flat_keys(forcefield.angles)]

        untyped = sorted(index for index, atom in atoms.items()
                         if UFF_TYPE_COMMENT in (atom.get('comment') or ''))

        if not flat and not untyped:
            return

        def describe(indices):
            return '-'.join(f'{atoms[index]["name"]}({atoms[index]["type"]})'
                            for index in indices)

        def report(terms):
            groups = {}
            for key, kind in terms:
                types = '-'.join(atoms[index]['type'] for index in key)
                groups.setdefault((kind, types), []).append(key)

            # the detail lines go through print_info: print_warning boxes
            # every line it is given, and a site can put a dozen groups here
            for (kind, types), keys in sorted(groups.items()):
                count = f' x{len(keys)}' if len(keys) > 1 else ''
                self.ostream.print_info(
                    f'  guessed {kind} {types}{count}, e.g. {describe(keys[0])}')

        if untyped:
            names = ', '.join(describe([index]) for index in untyped)
            self.ostream.print_warning(
                f'{len(untyped)} atom(s) fell back to UFF because GAFF could '
                f'not type them at all: {names}.')

        from_metal = [term for term in flat if shell & set(term[0])]
        elsewhere = [term for term in flat if term not in from_metal]

        if from_metal:
            if metal_blind_typing:
                cause = ('These were typed with the metal bonds ignored, so '
                         'the coordination is not what put them here: the '
                         'combination is missing from the parameter set.')
            else:
                cause = ('The metal bonds are part of what GAFF typed these '
                         'atoms from, so the coordination is what put them '
                         'here. Turning metal_blind_typing on types these '
                         'residues as the amino acids they are.')
            self.ostream.print_warning(
                f'{len(from_metal)} covalent term(s) in the ligand shell have '
                'no force field parameters and carry a flat guessed constant. '
                f'{cause} The fit does not repair them.')
            report(from_metal)

        if elsewhere:
            self.ostream.print_warning(
                f'{len(elsewhere)} covalent term(s) away from the metals also '
                'carry a guessed constant. Those are gaps in the parameter set '
                'rather than a consequence of the coordination.')
            report(elsewhere)

        self.ostream.flush()

    # ------------------------------------------------------------------
    # reporting
    # ------------------------------------------------------------------

    @on_master
    def print_binding_modes(self, binding_modes):
        """
        Prints the detected coordination sphere.
        """

        self.ostream.print_blank()
        self.ostream.print_header('Coordination sphere')
        self.ostream.print_header(19 * '-')

        for metal in binding_modes['metals']:
            self.ostream.print_header(
                param(f'metal {metal["element"]} (index {metal["index"]})',
                      f'charge {metal["formal_charge"]:+d}'))

        self.ostream.print_blank()
        valstr = '{:>10} {:>9} | {:>18} | {:>16}'.format('residue', 'atoms',
                                                         'distances (A)', 'mode')
        self.ostream.print_header(valstr)
        self.ostream.print_header(60 * '-')

        by_residue = {}
        for ligand in binding_modes['ligands']:
            by_residue.setdefault(ligand['res_index'], []).append(ligand)

        for group in by_residue.values():
            # a residue binding through several atoms is one ligand, so it gets
            # one row listing them side by side, the same way an atom bridging
            # two metals lists both of its distances. Merging is only
            # unambiguous while every atom contributes exactly one distance;
            # otherwise the distances could not be read back onto their atoms,
            # so that group stays one row per atom
            if any(len(ligand['distances']) != 1 for ligand in group):
                rows = [[ligand] for ligand in group]
            else:
                rows = [group]

            for row in rows:
                atoms = ', '.join(ligand['atom'] for ligand in row)
                distances = ', '.join(f'{d:.2f}' for ligand in row
                                      for d in ligand['distances'])
                modes = '/'.join(dict.fromkeys(ligand['mode'] for ligand in row))
                valstr = '{:>10} {:>9} | {:>18} | {:>16}'.format(
                    row[0]['residue'], atoms, distances, modes)
                self.ostream.print_header(valstr)

        for note in binding_modes['notes']:
            self.ostream.print_warning(note)

        self.ostream.print_blank()
        self.ostream.flush()

    @on_master
    def print_binding_mode_update(self, changes, largest_shift,
                                  dropped_residues):
        """
        Prints what re-detecting the coordination on a new geometry did.

        :param changes:
            The list of (kind, atom, detail) tuples describing the contacts
            that were gained, lost or reclassified. Empty when the coordination
            is unchanged.
        :param largest_shift:
            The largest change in Angstrom that the recorded metal-ligand
            bonds underwent.
        :param dropped_residues:
            The residues that no longer coordinate at all.
        """

        self.ostream.print_blank()
        self.ostream.print_header('Coordination update')
        self.ostream.print_header(19 * '-')
        self.ostream.print_header(param('largest bond change', f'{largest_shift:.2f} A'))

        if not changes:
            self.ostream.print_header(param('coordination', 'unchanged'))
            self.ostream.print_blank()
            self.ostream.print_info(
                'The new geometry gives the same coordination sphere; the '
                'binding modes and the connectivity matrix are kept as they '
                'are.')
            self.ostream.print_blank()
            self.ostream.flush()
            return

        self.ostream.print_header(param('contacts changed', len(changes)))
        self.ostream.print_blank()

        valstr = '{:>10} {:>12} | {:>46}'.format('change', 'atom', 'detail')
        self.ostream.print_header(valstr)
        self.ostream.print_header(72 * '-')

        for kind, atom, detail in changes:
            self.ostream.print_header('{:>10} {:>12} | {:>46}'.format(
                kind, atom, detail))

        self.ostream.print_blank()
        self.ostream.print_info(
            'The binding modes and the connectivity matrix were updated to '
            'the new geometry. Overwrite the ones you hold, or the fit will '
            'use a coordination the geometry no longer has.')

        for residue in dropped_residues:
            self.ostream.print_warning(
                f'{residue} no longer coordinates a metal, but it is still '
                'part of the truncated active site; extract the active site '
                'again to leave it out')

        self.ostream.print_blank()
        self.ostream.flush()

    @on_master
    def print_active_site(self, active_site, binding_modes):
        """
        Prints the composition of the truncated active site.
        """

        molecule = active_site['molecule']

        self.ostream.print_blank()
        self.ostream.print_header('Truncated active site')
        self.ostream.print_header(21 * '-')
        self.ostream.print_header(param('atoms', molecule.number_of_atoms()))
        self.ostream.print_header(param('charge', f'{int(molecule.get_charge()):+d}'))
        self.ostream.print_header(param('multiplicity',
                                        int(molecule.get_multiplicity())))
        self.ostream.print_header(
            param('capping hydrogens', len(active_site['cap_indices'])))
        self.ostream.print_header(
            param('bonds', int(active_site['connectivity_matrix'].sum() // 2)))
        print_param_list('residues', active_site['residues'], self.ostream)

        variants = sorted(binding_modes['variants'].values())
        print_param_list('protonation', variants, self.ostream)

        self.ostream.print_blank()
        self.ostream.flush()

    @on_master
    def print_mm_optimization(self,
                              active_site,
                              forcefield,
                              relaxed,
                              frozen_indices,
                              metal_keys,
                              equilibrium_labels,
                              bond_change_warning):
        """
        Prints what the crude MM relaxation did to the coordination sphere.

        The metal-ligand distances before and against after are the point of
        the table: the pass is there to clean up contacts and hydrogens, and
        a coordination sphere that moved more than a few hundredths of an
        Angstrom is the sign that it did something else instead.

        What the seeding did is read off the force field rather than off the
        settings that asked for it, so the table describes the terms the
        relaxation actually ran with.

        :param active_site:
            The active site, holding the geometry the pass started from.
        :param forcefield:
            The seeded force field the relaxation ran on.
        :param relaxed:
            The relaxed molecule.
        :param frozen_indices:
            The indices that were held fixed.
        :param metal_keys:
            The (bonds, angles) keys of the metal terms, from get_metal_keys.
        :param equilibrium_labels:
            The table from a seeded term's comment to the label its equilibrium
            source is printed as (SEEDED_EQUILIBRIUM_LABELS).
        :param bond_change_warning:
            How far a metal-ligand bond may move before it is reported.
        """

        def _seeded_constant(table, keys):
            """
            Returns the force constant the seeding put on a set of terms.

            :param table:
                The bond or angle table of the force field.
            :param keys:
                The keys of the metal terms.

            :return:
                The constant as a string, or 'varies' when they are not all the
                same, which the crude pass never makes them.
            """

            constants = {round(table[key]['force_constant'], 6) for key in keys}

            if len(constants) != 1:
                return 'varies'

            return f'{constants.pop():.0f}'

        def _seeded_equilibria(table, keys):
            """
            Returns where the seeding took its equilibrium values from.

            _seed_metal_terms writes that into the comment of every term it
            touches, so the force field says it without being asked again.

            :param table:
                The bond or angle table of the force field.
            :param keys:
                The keys of the metal terms.

            :return:
                'requested', 'given', 'measured', or 'mixed'.
            """

            sources = {table[key].get('comment') for key in keys}

            if len(sources) != 1:
                return 'mixed'

            return equilibrium_labels.get(sources.pop(), 'mixed')

        labels = active_site['molecule'].get_labels()
        molecule = active_site['molecule']
        metals = sorted(active_site['metal_indices'])
        bonds, angles = metal_keys

        before = molecule.get_coordinates_in_angstrom()
        after = relaxed.get_coordinates_in_angstrom()
        shift = np.linalg.norm(after - before, axis=1)

        self.ostream.print_blank()
        self.ostream.print_header('Crude MM relaxation')
        self.ostream.print_header(19 * '-')

        self.ostream.print_header(param('frozen atoms', len(frozen_indices)))
        self.ostream.print_header(
            param('metal centers',
                  'frozen' if set(metals) <= set(frozen_indices) else 'free'))
        self.ostream.print_header(
            param('metal bonds',
                  f'{len(bonds)}, k = {_seeded_constant(forcefield.bonds, bonds)}'))
        self.ostream.print_header(
            param(
                'metal angles', f'{len(angles)}, k = '
                f'{_seeded_constant(forcefield.angles, angles)}'
                if angles else 'left untouched'))
        self.ostream.print_header(
            param('bond equilibria', _seeded_equilibria(forcefield.bonds, bonds)))
        if angles:
            self.ostream.print_header(
                param('angle equilibria',
                      _seeded_equilibria(forcefield.angles, angles)))
        self.ostream.print_blank()

        valstr = '{:>12} {:>8} | {:>10} | {:>9} | {:>8}'.format(
            'atoms', 'elements', 'before (A)', 'after (A)', 'change')
        self.ostream.print_header(valstr)
        self.ostream.print_header(60 * '-')

        worst_bond = None
        for key in bonds:
            one_based = [index + 1 for index in key]
            was = molecule.get_distance_in_angstroms(one_based)
            now = relaxed.get_distance_in_angstroms(one_based)
            names = '-'.join(labels[index] for index in key)
            valstr = '{:>12} {:>8} | {:>10.2f} | {:>9.2f} | {:>+8.2f}'.format(
                str(key), names, was, now, now - was)
            self.ostream.print_header(valstr)
            if worst_bond is None or abs(now - was) > abs(worst_bond[1]):
                worst_bond = (names, now - was)

        for first in range(len(metals)):
            for second in range(first + 1, len(metals)):
                pair = (metals[first], metals[second])
                one_based = [index + 1 for index in pair]
                was = molecule.get_distance_in_angstroms(one_based)
                now = relaxed.get_distance_in_angstroms(one_based)
                names = '-'.join(labels[index] for index in pair)
                valstr = ('{:>12} {:>8} | {:>10.2f} | {:>9.2f} | '
                          '{:>+8.2f}').format(str(pair), names, was, now, now - was)
                self.ostream.print_header(valstr)

        self.ostream.print_blank()

        largest = int(np.argmax(shift))
        self.ostream.print_header(
            param('largest shift', f'{shift[largest]:.2f} A on '
                  f'{labels[largest]} {largest}'))
        self.ostream.print_header(param('mean shift', f'{shift.mean():.2f} A'))

        # A ligand swinging around its metal moves a long way in Cartesian
        # terms while the coordination sphere itself is untouched, so what
        # the pass has to be held to is the bond lengths, not the shifts.
        if worst_bond is not None:
            self.ostream.print_header(
                param('largest bond change',
                      f'{worst_bond[1]:+.2f} A on {worst_bond[0]}'))
            if abs(worst_bond[1]) > bond_change_warning:
                self.ostream.print_warning(
                    f'The crude relaxation changed a {worst_bond[0]} bond by '
                    f'{worst_bond[1]:+.2f} A. Check the metal terms it was '
                    'given before trusting the geometry it produced.')

        self.ostream.print_blank()
        self.ostream.flush()

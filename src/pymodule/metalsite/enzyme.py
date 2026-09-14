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
The enzyme phase of the metal site force field builder: what happens once
the metal terms are fitted and the whole protein is to carry them.
"""

import sys

from ..errorhandler import assert_msg_critical
from .qm import QmParameterizer
from .util import (Shell, on_master, param, get_metal_keys,
                   get_metal_impropers, redistribute_cap_charges,
                   backbone_charge_shift, _bond_separation)
from .openmmxml import (SITE_RESIDUE_NAME, protein_atom_parameters,
                        restructure_topology, build_templates,
                        forcefield_xml, metal_term_keys)

try:
    import openmm as mm
    import openmm.app as mmapp
    import openmm.unit as mmunit
except ImportError:
    pass


class EnzymeSystemBuilder(Shell):
    """
    Carries the fitted metal terms and charges of an active site into the
    whole enzyme: onto one OpenMM system by index, or into a force field
    XML that builds any system for it.

    Every method runs on the master rank and its result is broadcast. The
    three steps of create_enzyme_system that edit a system in place --
    redistribute_charges, redistribute_backbone_charges and
    rescale_14_exceptions -- are usable on their own on one rank; under
    MPI they edit the master's system alone, and it is the whole system
    create_enzyme_system returns that every rank receives.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.
    """

    @on_master
    def create_enzyme_system(self,
                             topology,
                             active_site,
                             forcefield,
                             partial_charges=None,
                             forcefield_files=('amber14-all.xml',
                                               'amber14/tip3pfb.xml')):
        """
        Injects the fitted metal terms into a force field system for the whole
        enzyme.

        The protein force field already covers everything except the metal, so
        only the metal bonds, angles and impropers are transferred -- every
        other term of the site belongs to a residue it parameterizes itself.
        The active site was built from the topology, so the atom map of the
        active site gives the correspondence directly and no graph matching is
        needed. Capping hydrogens are skipped, since they stand in for CA atoms
        that the protein force field parameterizes itself.

        The impropers are added to the periodic torsion force by index, which is
        what a metal improper needs: OpenMM builds an improper out of a residue
        template only for a central atom and three atoms bonded to it, so the
        one a bidentate carboxylate gets -- the carboxylate carbon against its
        two oxygens and the metal, which the carbon is not bonded to -- never
        reaches a system built that way. Added here it does.

        :param topology:
            The protonated OpenMM topology of the whole enzyme.
        :param active_site:
            The active site, for the map back to the topology.
        :param forcefield:
            The force field generator carrying the fitted metal parameters.
        :param partial_charges:
            The charges fitted on the active site, which replace the charges
            of the coordination sphere through redistribute_charges. D4
            charges are used when none are given, so that the system carries
            the same charges as the force field built beside it rather than
            the protein force field's own.
        :param forcefield_files:
            The OpenMM force field files for the protein.

        :return:
            The tuple of the OpenMM system and the list of added terms.
        """

        assert_msg_critical('openmm' in sys.modules,
                            'create_enzyme_system: openmm is '
                            'required')

        openmm_ff = mmapp.ForceField(*forcefield_files)
        system = openmm_ff.createSystem(topology, nonbondedMethod=mmapp.NoCutoff)

        if partial_charges is None:
            partial_charges = QmParameterizer(
                self.comm, self.ostream).d4_charges(active_site)

        self.redistribute_charges(system, topology, active_site,
                                  partial_charges)

        # the exceptions were worked out from the protein force field's own
        # charges, inside createSystem, and know nothing of the ones just
        # written over them
        for generator in openmm_ff.getGenerators():
            if isinstance(generator, mmapp.forcefield.NonbondedGenerator):
                rescaled = self.rescale_14_exceptions(
                    system, topology, generator.coulomb14scale)
                self.ostream.print_info(
                    f'Brought {rescaled} 1-4 exception(s) up to date with the '
                    'fitted charges.')

        atom_map = active_site['atom_map']
        caps = set(active_site['cap_indices'])
        bonds, angles = get_metal_keys(forcefield, active_site)
        impropers = get_metal_impropers(forcefield, active_site)

        bond_force = None
        angle_force = None
        torsion_force = None
        for force in system.getForces():
            if isinstance(force, mm.HarmonicBondForce):
                bond_force = force
            elif isinstance(force, mm.HarmonicAngleForce):
                angle_force = force
            elif isinstance(force, mm.PeriodicTorsionForce):
                torsion_force = force

        assert_msg_critical(
            bond_force is not None and angle_force is not None,
            'create_enzyme_system: the protein '
            'system has no harmonic bond or angle force to extend')

        assert_msg_critical(
            torsion_force is not None or not impropers,
            'create_enzyme_system: the force field carries metal impropers but '
            'the protein system has no periodic torsion force to extend')

        added = []

        for key in bonds:
            if caps & set(key):
                continue
            params = forcefield.bonds[key]
            bond_force.addBond(
                atom_map[key[0]], atom_map[key[1]],
                params['equilibrium'] * mmunit.nanometer, params['force_constant'] *
                mmunit.kilojoule_per_mole / mmunit.nanometer**2)
            added.append(('bond', key))

        for key in angles:
            if caps & set(key):
                continue
            params = forcefield.angles[key]
            angle_force.addAngle(
                atom_map[key[0]], atom_map[key[1]], atom_map[key[2]],
                params['equilibrium'] * mmunit.degree, params['force_constant'] *
                mmunit.kilojoule_per_mole / mmunit.radian**2)
            added.append(('angle', key))

        for key in impropers:
            if caps & set(key):
                continue
            params = forcefield.impropers[key]
            # an improper key is the central atom and its three substituents,
            # and OpenMM writes that as the torsion (first, second, central,
            # third) -- the order the same improper comes out in when the
            # generator's own XML is loaded through mmapp.ForceField, which is
            # what the active site is relaxed on. Adding it in the key's own
            # order instead would measure a different dihedral of the same four
            # atoms.
            torsion_force.addTorsion(atom_map[key[1]], atom_map[key[2]],
                                     atom_map[key[0]], atom_map[key[3]],
                                     params['periodicity'],
                                     params['phase'] * mmunit.degree,
                                     params['barrier'] * mmunit.kilojoule_per_mole)
            added.append(('improper', key))

        counts = {
            kind: sum(1 for term in added if term[0] == kind)
            for kind in ('bond', 'angle', 'improper')
        }
        self.ostream.print_info(
            f'Added {counts["bond"]} metal bond(s), {counts["angle"]} metal '
            f'angle(s) and {counts["improper"]} metal improper(s) to the enzyme '
            'system.')
        self.ostream.flush()

        return system, added

    @on_master
    def create_enzyme_forcefield(self,
                                 topology,
                                 positions,
                                 active_site,
                                 forcefield,
                                 partial_charges=None,
                                 forcefield_files=('amber14-all.xml',
                                                   'amber14/tip3pfb.xml'),
                                 site_residue_name=SITE_RESIDUE_NAME,
                                 drop_torsions_across_metal_bonds=True):
        """
        Writes the fitted metal site as an OpenMM force field XML.

        The counterpart of create_enzyme_system, and what it should be
        preferred to: that one puts the fitted terms onto one System and can
        describe nothing else, while this returns a force field and the
        topology it is for, which can be solvated, extended and rebuilt as
        often as the caller likes.

        The metal-ligand bonds become real bonds of the topology it returns,
        so the system built from it carries the 1-2 and 1-3 exclusions and the
        1-4 scaling a bonded metal model should have. The injected system has
        none of those, and that is the one place the two are meant to differ.

        :param topology:
            The protonated OpenMM topology of the whole enzyme.
        :param positions:
            Its positions as an (N, 3) array in Angstrom.
        :param active_site:
            The active site, for the map back to the topology.
        :param forcefield:
            The force field generator carrying the fitted metal parameters.
        :param partial_charges:
            The charges fitted on the active site. D4 charges are used when
            none are given, the way create_enzyme_system falls back, so that
            the file carries the same charges as the force field built beside
            it rather than the protein force field's own.
        :param forcefield_files:
            The OpenMM force field files for the protein. Every atom keeps the
            type these give it, so the file that comes back has to be loaded
            beside the same ones.
        :param site_residue_name:
            The name of the residue the site is moved into.
        :param drop_torsions_across_metal_bonds:
            Whether to zero the wildcard proper torsions the protein force
            field writes across a metal bond once it is a real one.

        :return:
            A dictionary holding the XML, the restructured topology and its
            positions, the residue templates and what restructure_topology
            returned.
        """

        assert_msg_critical('openmm.app' in sys.modules,
                            'create_enzyme_forcefield: openmm is required')

        if partial_charges is None:
            partial_charges = QmParameterizer(
                self.comm, self.ostream).d4_charges(active_site)

        protein_ff_not_used, protein_parameters = protein_atom_parameters(
            topology, forcefield_files)

        restructured = restructure_topology(topology,
                                            positions,
                                            active_site,
                                            forcefield,
                                            site_residue_name=site_residue_name)
        self.ostream.print_info(
            f'Moved {len(restructured["site_indices"])} atoms out of '
            f'{len(restructured["stub_residues"])} residues into residue '
            f'{site_residue_name}, and bonded '
            f'{len(restructured["metal_bonds"])} metal-ligand contacts.')
        self.ostream.flush()

        templates, correction = build_templates(topology, active_site,
                                                restructured, partial_charges,
                                                protein_parameters)
        shift = correction['shift']
        self.ostream.print_info(
            f'Built {len(templates)} residue templates; the coordination '
            f'region gives each of the {len(correction["uncovered"])} atoms '
            f'the active site does not cover {shift:+.4f} e.')
        self.ostream.flush()

        xml = forcefield_xml(
            active_site,
            forcefield,
            restructured,
            templates,
            forcefield_files=forcefield_files,
            drop_torsions_across_metal_bonds=drop_torsions_across_metal_bonds)
        bonds, angles, impropers = metal_term_keys(forcefield, active_site)
        self.ostream.print_info(
            f'Wrote {len(templates)} residue templates, {len(bonds)} metal '
            f'bonds, {len(angles)} metal angles and {len(impropers)} metal '
            'impropers.')
        self.ostream.flush()

        return {
            'xml': xml,
            'topology': restructured['topology'],
            'positions': restructured['positions'],
            'templates': templates,
            'backbone_shift': shift,
            'restructured': restructured,
        }

    @on_master
    def redistribute_charges(self, system, topology, active_site,
                             partial_charges):
        """
        Applies both charge corrections to a protein system.

        :param system:
            The OpenMM system to modify in place.
        :param topology:
            The protonated topology the system was built from.
        :param active_site:
            The active site, for the map back to the topology.
        :param partial_charges:
            The charges fitted on the active site, capping hydrogens included.

        :return:
            The tuple of the charge moved off the caps and the shift applied
            to each uncovered atom.
        """

        caps = active_site['cap_indices']
        cap_charge = float(sum(partial_charges[index] for index in caps))
        charges = redistribute_cap_charges(active_site, partial_charges)
        shift = self.redistribute_backbone_charges(system, topology,
                                                   active_site, charges)

        return cap_charge, shift

    @on_master
    def redistribute_backbone_charges(self, system, topology, active_site,
                                      partial_charges):
        """
        Restores the charge of the coordination region after the fitted
        charges are written into a protein system.

        :param system:
            The OpenMM system to modify in place.
        :param topology:
            The protonated topology the system was built from.
        :param active_site:
            The active site, for the map back to the topology.
        :param partial_charges:
            The active site charges, with the capping hydrogens already folded in.

        :return:
            The shift applied to each uncovered atom.
        """

        assert_msg_critical('openmm' in sys.modules,
                            'redistribute_backbone_charges: openmm '
                            'is required')

        nonbonded = None
        for force in system.getForces():
            if isinstance(force, mm.NonbondedForce):
                nonbonded = force
                break

        assert_msg_critical(
            nonbonded is not None, 'redistribute_backbone_charges: the '
            'system has no nonbonded force to write charges into')

        def get_charge(index):
            return nonbonded.getParticleParameters(index)[0].value_in_unit(
                mmunit.elementary_charge)

        correction = backbone_charge_shift(get_charge, topology, active_site,
                                           partial_charges)

        covered = correction['covered']
        residue_indices = correction['residue_indices']
        region = correction['region']
        uncovered = correction['uncovered']
        shift = correction['shift']
        total_before = correction['total_before']
        total_after = correction['total_after']
        difference = correction['difference']

        for atom in region:
            parameters = nonbonded.getParticleParameters(atom.index)
            if atom.index in covered:
                new_charge = covered[atom.index]
            else:
                new_charge = get_charge(atom.index) + shift
            nonbonded.setParticleParameters(atom.index, new_charge, parameters[1],
                                            parameters[2])

        residual = total_before - sum(get_charge(atom.index) for atom in region)
        assert_msg_critical(
            abs(residual) < 1.0e-6, 'redistribute_backbone_charges: the '
            f'charge moved by {residual:+.6f} e')

        self.ostream.print_blank()
        self.ostream.print_header('Charges written into the protein')
        self.ostream.print_header(31 * '-')
        self.ostream.print_header(
            param('coordination region',
                  f'{len(residue_indices)} residues, {len(region)} atoms'))
        self.ostream.print_header(
            param('covered by the active site', f'{len(covered)} atoms'))
        self.ostream.print_header(
            param('left to the protein', f'{len(uncovered)} atoms'))
        self.ostream.print_blank()
        self.ostream.print_header(
            param('region charge, protein', f'{total_before:+.4f} e'))
        self.ostream.print_header(
            param('region charge, fitted', f'{total_after:+.4f} e'))
        self.ostream.print_header(
            param('difference to recover', f'{difference:+.4f} e'))
        self.ostream.print_header(
            param('shift per uncovered atom', f'{shift:+.4f} e'))
        self.ostream.print_blank()
        # a shift and not a scaling: after the caps are folded in the active site
        # carries its own formal charge exactly, so the amount to recover is
        # always minus the sum of the backbone charges and the scale factor
        # that achieves it is identically zero
        self.ostream.print_info(
            'The region is corrected as a whole, so the ligand-to-metal '
            'donation the fit captured is kept; balancing each residue to its '
            'formal charge separately would undo it.')
        self.ostream.print_blank()
        self.ostream.flush()

        return shift

    @on_master
    def rescale_14_exceptions(self, system, topology, coulomb14scale):
        """
        Puts the fitted charges into the 1-4 exceptions as well.

        OpenMM works the exceptions out inside createSystem, from the charges
        the protein force field gave the atoms, and writing the fitted charges
        onto the particles afterwards leaves them behind: every 1-4
        electrostatic interaction of the coordination region goes on being
        scaled from the charge the fit replaced. Which pairs those are is
        worked out the way createSystem works it out, from the bond graph --
        a pair three bonds apart is scaled, anything closer is excluded
        outright and has nothing to rescale.

        :param system:
            The OpenMM system to modify in place.
        :param topology:
            The topology the system was built from.
        :param coulomb14scale:
            The 1-4 electrostatic scaling of the protein force field.

        :return:
            The number of exceptions that were brought up to date.
        """

        assert_msg_critical('openmm' in sys.modules,
                            'rescale_14_exceptions: openmm is required')

        nonbonded = None
        for force in system.getForces():
            if isinstance(force, mm.NonbondedForce):
                nonbonded = force
                break

        assert_msg_critical(
            nonbonded is not None, 'rescale_14_exceptions: the system has no '
            'nonbonded force to read exceptions from')

        bonded = {}
        for first, second in topology.bonds():
            bonded.setdefault(first.index, set()).add(second.index)
            bonded.setdefault(second.index, set()).add(first.index)

        def charge_of(index):
            return nonbonded.getParticleParameters(index)[0].value_in_unit(
                mmunit.elementary_charge)

        rescaled = 0
        for index in range(nonbonded.getNumExceptions()):
            first, second, charge_product, sigma, epsilon = (
                nonbonded.getExceptionParameters(index))
            if _bond_separation(bonded, first, second) != 3:
                continue
            updated = coulomb14scale * charge_of(first) * charge_of(second)
            if abs(updated -
                   charge_product.value_in_unit(mmunit.elementary_charge**2)
                   ) > 1.0e-10:
                rescaled += 1
            nonbonded.setExceptionParameters(index, first, second, updated, sigma,
                                             epsilon)

        return rescaled

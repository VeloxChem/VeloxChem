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
from enum import IntEnum
import numpy as np
import time
import sys

from ..veloxchemlib import mpi_master
from ..outputstream import OutputStream
from ..errorhandler import assert_msg_critical
from ..mmforcefieldgenerator import MMForceFieldGenerator
from . import util
from . import openmmxml
from .builder import ActiveSiteBuilder
from .qm import QmParameterizer
from .enzyme import EnzymeSystemBuilder

try:
    import openmm as mm
    import openmm.app as mmapp
    import openmm.unit as mmunit
except ImportError:
    pass

# ----------------------------------------------------------------------
# stage
# ----------------------------------------------------------------------


class Stage(IntEnum):
    """
    Staging state for this class
    - EMPTY: nothing built yet.
    - ACTIVE_SITE: a truncated cluster exists and is fully editable.
    - FITTED: the QM has been paid for and the metal terms are fitted. Metal
        bonds can be added or removed, but the site cannot be changed in
        any other way.
    - ENZYME: An enzyme system has been created with the fitted metal terms.
        Not terminal: a bond edit still refits, and drops the system it
        invalidates.
    """

    EMPTY = 0
    ACTIVE_SITE = 1
    FITTED = 2
    ENZYME = 3


# What each stage owns, and therefore what entering a lower one throws away.
#
# The Hessian and the charges are deliberately not here.
# They are not state, but a result of the expensive steps
STAGE_FIELDS = {
    Stage.ACTIVE_SITE: (),
    Stage.FITTED: ('_forcefield', ),
    Stage.ENZYME: ('_enzyme_system', '_enzyme_forcefield'),
}


class MetalSiteForceFieldBuilder:
    """
    Builds a bonded force field for the zinc center of a metallo-enzyme.

    Identifies the coordination sphere, fixes the protonation of the
    coordinating residues, truncates a QM active site at the CA-CB bonds,
    optimises the truncated active site, calculates a partial hessian and resp charges
    and fits the metal-ligand bond and angle parameters to a QM Hessian with the
    Seminario method. The fitted parameters can then be injected into a
    protein force field for the whole enzyme.

    The pipeline is three calls:

        molecule = builder.build_active_site('site.cif')
        forcefield = builder.build_forcefield()
        system, topology = builder.create_enzyme_system()

    The active site can be edited in between them, and how much of it can be
    edited is what separates the two gaps.

    The build_forcefield function calls all expensive QM steps, including by default
    a constrained optimization, a partial Hessian calculation, and a RESP charge calculation.

    Between the first and the second step, the
    site can be changed in the following ways: add_metal_bond ,remove_metal_bond,
    update_protonation_state, and include_residue, remove_residue.
    Between the second and the third step only bonds can be added or removed.

    The object holds the settings and the intermediates; the work itself is
    done by the phase classes it owns -- ActiveSiteBuilder, QmParameterizer
    and EnzymeSystemBuilder -- which keep no state and take everything they
    use. What a step produced is read back off a read-only property rather
    than juggled by the caller.

    Under MPI every method of a phase class is either run on the master rank
    and broadcast, or collective: every rank calls the three expensive QM
    steps at the same point, with the same molecule. This class runs on
    every rank and does nothing non-deterministic outside a shell call.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - metal_bond_cutoff: The distance in Angstrom within which a donor atom
          is taken to be bonded to a metal center. Generous on purpose: the
          scan reads an unrelaxed structure, where a stretched bridging
          contact is still a bond. A contact up to util.REPORT_CUTOFF_MARGIN
          further out is still reported for review.
        - prepare_protein: if the loaded pdb structure should be prepared with pdbfixer first
        - scf_drv: SCF driver instance that can be used to override the automatically generated one
        - xcfun: The exchange-correlation functional used when scf_drv is not
          provided
        - basis_set_label: The basis set label.
        - mute_scf: The flag for muting the output of the QM drivers.
        - do_qm_optimization: The flag for optimizing the active site before the
          Hessian is computed.
        - do_hessian: The flag for computing the Hessian. With this off, build_forcefield
          will use default force constants
        - do_resp: The flag for computing RESP charges. Falls back on D4 charges.
        - calculate_partial_hessian: The flag for restricting the Hessian to
          the atom pairs the metal terms are fitted from. False computes the
          whole thing, which is signifacntly more expensive but can useful if
          vibrational modes are of interest
        - partial_hessian_cutoff: The radius at which bonds are percieved for
          the partial Hessian. Intentionally slightly larger than the metal_bond_cutoff,
          so that a donor atom that is not bonded to a metal but is close enough
          to be covered by the partial Hessian can be added as a bond afterwards
          without paying for the Hessian again. None restricts it to the bonding as perceived.
        - constrain_capping_hydrogens: The flag for constraining the capping
          hydrogens in addition to the beta carbons.
        - average_metal_terms: The flag for averaging the fitted metal terms
          over equivalent atoms.
        - metal_hessian_fitting_method: The method the metal bonds and angles
          are fitted with: 'seminario' (default), 'improved-seminario' or
          'phf'/'phf(k)'.
        - metal_blind_typing: The flag for perceiving the GAFF atom types as
          if the metal bonds were not there, so that a coordinating residue
          is typed as the amino acid it is.
        - mute_forcefield_generator: The flag for keeping the force field
          generator's own commentary out of the output.
        - prune_weak_bridge_bonds: The flag for dropping the long arm of a
          bridging residue when the fit gave it no force constant.
        - add_metal_planarity_impropers: The flag for adding a weak improper
          nudging a metal into the plane of a coordinating histidine ring or
          a bidentate carboxylate.
        - reparameterize_metal_angles: The flag for touching the metal angles
          at all. False leaves every one of them at the value the generator
          guessed, in the crude pass and in the Hessian fit alike.
        - mm_constrain_metals: The flag for holding the metal centers fixed in
          the crude pass as well, on top of the beta carbons.
        - default_metal_bond_equilibria: Equilibrium metal-ligand distances in nm,
          keyed by element pair in either order, replacing the values measured
          on the input geometry. LITERATURE_METAL_BONDS is such a table.
        - default_metal_angle_equilibria: Equilibrium metal angles in degrees, keyed by
          element triple in either order, replacing the measured values.
        - protein_forcefield_files: The OpenMM force field files the enzyme
          system is built from and the enzyme force field XML is loaded
          beside.
        - folder: The folder every step writes its result to as soon as
          it has it, and reads back on a later run. Named after the creation
          time by default, so runs do not collide.
        - comm: The MPI communicator.
        - rank: The rank of the MPI process.
        - nodes: The number of MPI processes.
        - ostream: The output stream.
    """

    # ------------------------------------------------------------------
    # settings and state
    # ------------------------------------------------------------------

    def __init__(self, comm=None, ostream=None):

        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(sys.stdout)
            else:
                ostream = OutputStream(None)

        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()
        self.ostream = ostream

        # The settings. Each is read by the phase class _sites, _qm or
        # _enzyme builds when a phase runs, so a change made here at any
        # point reaches the next step.
        self.metal_bond_cutoff = 3.0
        self.prepare_protein = True

        self.scf_drv = None
        self.xcfun = 'PBE0'
        self.basis_set_label = 'def2-svp'
        self.mute_scf = True

        # workflow
        self.do_qm_optimization = True
        self.do_hessian = True
        self.do_resp = True
        self.calculate_partial_hessian = True
        # Deliberately more forgiving than metal_bond_cutoff: the perception
        # is a distance cutoff read on a single geometry, and the contact
        # that falls just outside it is exactly the one add_metal_bond is
        # reached for afterwards. Filling a block that was never computed
        # costs the whole Hessian again; widening what is computed costs a
        # few more blocks and changes nothing about what is fitted.
        self.partial_hessian_cutoff = 3.5
        self.constrain_capping_hydrogens = False
        self.average_metal_terms = False
        self.metal_hessian_fitting_method = 'seminario'

        self.metal_blind_typing = True
        self.mute_forcefield_generator = True
        self.prune_weak_bridge_bonds = True
        self.add_metal_planarity_impropers = True
        self.reparameterize_metal_angles = True
        self.default_metal_bond_equilibria = None
        self.default_metal_angle_equilibria = None
        self.mm_constrain_metals = False

        # the protein force field the fitted metal terms are added to
        self.protein_forcefield_files = ('amber14-all.xml',
                                         'amber14/tip3pfb.xml')

        self.folder = f'metal_site_{int(time.time())}'

        self._protonation_overrides = None  # todo what exactly is this, it should be private

        self._stage = Stage.EMPTY
        self._request = util.empty_request()
        self._mm_opt = True
        self._topology = None
        self._positions = None
        self._protonated_topology = None
        self._protonated_positions = None
        self._active_site = None
        self._forcefield = None
        # where the force field came from. A fit made here describes the
        # geometry its Hessian was computed on; one adopted from a template
        # describes the template, and a relaxation toward it therefore does
        # not invalidate it. Provenance is not derivable from the force
        # field, so it is recorded rather than worked out.
        self._adopted_forcefield = False
        self._hessian = None
        self._partial_charges = None
        self._enzyme_system = None
        self._enzyme_forcefield = None

    # The phase classes, built from the settings as they stand whenever a
    # phase runs. They hold nothing but those settings, the communicator and
    # the stream, so one made for a step is as good as one kept -- and a
    # setting changed, or a stream swapped for a silent one as the
    # shoehorning does, reaches the next step without anything keeping the
    # two in sync.

    def _sites(self):
        """
        An ActiveSiteBuilder on the current settings.
        """

        return ActiveSiteBuilder(
            self.comm,
            self.ostream,
            metal_bond_cutoff=self.metal_bond_cutoff,
            prepare_protein=self.prepare_protein,
            constrain_capping_hydrogens=self.constrain_capping_hydrogens,
            mm_constrain_metals=self.mm_constrain_metals,
            metal_blind_typing=self.metal_blind_typing,
            mute_forcefield_generator=self.mute_forcefield_generator,
            add_metal_planarity_impropers=self.add_metal_planarity_impropers,
            reparameterize_metal_angles=self.reparameterize_metal_angles,
            default_metal_bond_equilibria=self.default_metal_bond_equilibria,
            default_metal_angle_equilibria=self.default_metal_angle_equilibria)

    def _qm(self):
        """
        A QmParameterizer on the current settings.
        """

        return QmParameterizer(
            self.comm,
            self.ostream,
            scf_drv=self.scf_drv,
            xcfun=self.xcfun,
            basis_set_label=self.basis_set_label,
            mute_scf=self.mute_scf,
            constrain_capping_hydrogens=self.constrain_capping_hydrogens,
            partial_hessian_cutoff=self.partial_hessian_cutoff,
            average_metal_terms=self.average_metal_terms,
            metal_hessian_fitting_method=self.metal_hessian_fitting_method,
            prune_weak_bridge_bonds=self.prune_weak_bridge_bonds,
            reparameterize_metal_angles=self.reparameterize_metal_angles)

    def _enzyme(self):
        """
        An EnzymeSystemBuilder on the current settings.
        """

        return EnzymeSystemBuilder(
            self.comm,
            self.ostream,
            protein_forcefield_files=self.protein_forcefield_files)

    # ------------------------------------------------------------------
    # results
    # ------------------------------------------------------------------

    @property
    def stage(self):
        """
        How far along the workflow the builder is, as a lowercase name:
        'empty', 'active_site', 'fitted' or 'enzyme'.

        What it says is what decides which operations are legal, so this is
        worth reading when one is refused.
        """

        return self._stage.name.lower()

    @property
    def active_site_molecule(self):
        """
        The truncated active site, as the last step to touch it left it.
        """

        return None if self._active_site is None else (
            self._active_site['molecule'])

    @property
    def active_site(self):
        """
        The active site dictionary, as the last step to touch it left it:
        the molecule, the atom map back to the topology, the metal, cap and
        beta carbon indices, the labels and the connectivity matrix. This is
        what MetalForceFieldManager describes and matches against templates,
        and what adopt_forcefield takes back once a template has been
        transferred onto it.
        """

        return self._active_site

    @property
    def active_site_forcefield(self):
        """
        The force field carrying the fitted metal terms.
        """

        return self._forcefield

    @property
    def enzyme_system(self):
        """
        The OpenMM system of the whole enzyme.
        """

        return self._enzyme_system

    @property
    def enzyme_forcefield(self):
        """
        The metal site written as an OpenMM force field, as
        create_enzyme_forcefield returned it.
        """

        return self._enzyme_forcefield

    @property
    def enzyme_topology(self):
        """
        The protonated topology the enzyme system was built for. The same
        object as protonated_topology; read that one unless the enzyme
        system is what is being asked about.
        """

        return self._protonated_topology

    @property
    def enzyme_positions(self):
        """
        The positions of that topology, in Angstrom.
        """

        return self._protonated_positions

    def show_active_site(self, **kwargs):
        """
        Draws the extracted active site with all forcefield bonds explicitly visible.
        The atom labels are what the edit methods want to be told
        about an atom: the atom index on a metal, the resid on the beta
        carbon standing for its residue, the atom name on any other heavy
        atom.

        Every other keyword goes to Molecule.show untouched, so width, height,
        forming_bonds and the rest work as they always do. Passing bonds or
        atom_labels explicitly overrides what is worked out here.

        :param kwargs:
            Further keyword arguments for Molecule.show.

        :return:
            Whatever Molecule.show returns.
        """

        self._require('show_active_site', Stage.ACTIVE_SITE)

        return self._sites().show_active_site(self._active_site, **kwargs)

    def print_active_site(self):
        """
        Prints the composition of the active site as it now stands: atom and
        bond counts, the residues, and their protonation.
        """

        self._require('print_active_site', Stage.ACTIVE_SITE)

        self._sites().print_active_site(self._active_site,
                                        self.binding_modes)

    @property
    def optimization_constraints(self):
        """
        The constraint the constrained optimization runs under, as geomeTRIC
        reads it.

        What is frozen is what the fitted parameters end up describing, so
        this is here to be read before paying for the optimization that uses
        it. None when there is no active site yet, or when nothing is frozen.
        """

        if self._active_site is None:
            return None

        return util.freeze_constraints(
            util.constrained_indices(self._active_site,
                                     self.constrain_capping_hydrogens))

    @property
    def binding_modes(self):
        """
        Which residue coordinates which metal, at what distance.

        Derived when asked rather than kept, from the positions as they now
        stand and the decisions in the request, so it cannot fall behind a
        relaxation or an edit. The atom indices in it are those of the
        protonated topology, which is what the active site indexes into and
        the only numbering shown anywhere.
        """

        if self._protonated_topology is None:
            return None

        return self._sites().derive_binding_modes(self._protonated_topology,
                                                  self._protonated_positions,
                                                  self._request)

    @property
    def hessian(self):
        """
        The Hessian the metal terms were fitted from.
        """

        return self._hessian

    @property
    def partial_charges(self):
        """
        The charges of the active site, as they now stand.

        Before a fit these are the raw result of the charge calculation --
        RESP, or D4 when do_resp is off. Once the force field exists they are
        the ones it carries, with the capping hydrogens' share folded back
        into the site, which is the same array that reaches the enzyme system
        and partial_charges.txt.
        """

        return self._partial_charges

    # ------------------------------------------------------------------
    # the pipeline
    # ------------------------------------------------------------------

    def build_active_site(self,
                          cif_path=None,
                          mm_opt=True,
                          coordinating_residues=None):
        """
        Builds the truncated active site from a new pdb or rebuilds with edited binding
        modes from an earlier pdb

        With a path, the structure is read and repaired, the coordination
        sphere is detected and reported, and everything a previous call
        produced is dropped. Without one, the detection is not run again: the
        binding modes already on the builder are used.

        The coordinating residues are then protonated, the site is truncated
        at the CA-CB bonds, and unless mm_opt is switched off it is relaxed on
        a crude default force field.

        Calling this by hand is only needed to start from a different
        structure, to change mm_opt, or to reopen editing after
        build_forcefield.

        :param cif_path:
            The path to a .pdb, .cif or .pdbx file, or None to rebuild from
            the binding modes of the previous call.
        :param mm_opt:
            Whether to relax the extracted site on a crude force field.
        :param coordinating_residues:
            Coordinating residues to include, regardless of their distance from the active site
            as ids ('130') or labels ('ASP130'). Only used when a path is given.

        :return:
            The active site molecule.
        """

        assert_msg_critical('openmm' in sys.modules,
                            'MetalSiteForceFieldBuilder: openmm is required')

        if cif_path is not None:
            self._print_header(cif_path)
            topology, positions = self._sites().load_and_prepare_protein(
                cif_path)
        else:
            assert_msg_critical(
                self._topology is not None,
                'MetalSiteForceFieldBuilder.build_active_site: there is no '
                'structure to rebuild from. Call it with a cif_path first.')
            topology = self._topology
            positions = self._positions

        state = self._build_active_site(topology,
                                        positions,
                                        mm_opt,
                                        cif_path is not None,
                                        coordinating_residues)

        # an edit rebuilds the site the way this call built it, so what it
        # was told has to outlive the call
        self._mm_opt = mm_opt

        self._adopt(state)

        # everything downstream belonged to the site that was just replaced
        self._enter(Stage.ACTIVE_SITE)

        return self.active_site_molecule

    def _build_active_site(self,
                           topology,
                           positions,
                           mm_opt,
                           report_detection,
                           coordinating_residues=None):
        """
        The structural pass on a prepared structure, written to the folder
        and relaxed on the crude force field when asked.

        :param topology:
            The prepared topology.
        :param positions:
            Its positions in Angstrom.
        :param mm_opt:
            Whether to relax the extracted site on the crude force field.
        :param report_detection:
            Whether the coordination found on the structure is printed; a
            rebuild after an edit does not print it again.
        :param coordinating_residues:
            Residues to force as ligands, as build_active_site takes them.

        :return:
            The state _adopt takes on.
        """

        state = self._sites().build_active_site(
            topology,
            positions,
            self._request,
            self._protonation_overrides,
            report_detection,
            coordinating_residues=coordinating_residues)
        state['topology'] = topology
        state['positions'] = positions

        # PDBFile.writeFile wants a quantity: handed bare numbers it writes
        # nanometers as Angstrom and produces a structure shrunk tenfold
        self._save_intermediate(
            'protonated.pdb', lambda path: mmapp.PDBFile.writeFile(
                state['protonated_topology'],
                np.asarray(state['protonated_positions']) * mmunit.angstrom,
                str(path),
                keepIds=True))

        if mm_opt:
            state['active_site']['molecule'] = self._crude_relax(
                state['active_site'],
                self._sites().manual_bond_equilibria(state['active_site'],
                                                     state['protonated_topology'],
                                                     state['binding_modes']))

        return state

    # ------------------------------------------------------------------
    # reporting
    # ------------------------------------------------------------------

    def _functional_label(self):
        """
        Returns the functional actually in use.

        The xcfun setting only applies when no SCF driver is provided, so an
        assigned driver has to be asked for its own. Its xcfun stays a plain
        string until it runs, and becomes a functional object afterwards.

        :return:
            The name of the functional.
        """

        if self.scf_drv is not None:
            xcfun = self.scf_drv.xcfun
        else:
            xcfun = self.xcfun

        if xcfun is None:
            return 'Hartree-Fock'

        if hasattr(xcfun, 'get_func_label'):
            return xcfun.get_func_label()

        return str(xcfun)

    def _print_header(self, structure):
        """
        Prints the settings of the run.

        :param structure:
            The structure file being processed.
        """

        param = util.param

        self.ostream.print_blank()
        self.ostream.print_header('Metal Site Force Field Builder')
        self.ostream.print_header(32 * '=')
        self.ostream.print_blank()

        self.ostream.print_header(param('structure', Path(structure).name))
        self.ostream.print_header(
            param('primary cutoff', f'{self.metal_bond_cutoff:.2f} A'))
        self.ostream.print_header(
            param('secondary cutoff',
                  f'{self.metal_bond_cutoff + util.REPORT_CUTOFF_MARGIN:.2f} A'))
        self.ostream.print_header(param('basis set', self.basis_set_label))
        self.ostream.print_header(
            param('xc functional', self._functional_label()))
        self.ostream.print_header(
            param('SCF driver',
                  'given' if self.scf_drv is not None else 'default'))
        self.ostream.print_header(
            param('QM optimization', self.do_qm_optimization))
        if not self.do_hessian:
            hessian_line = 'skipped'
        elif self.calculate_partial_hessian:
            hessian_line = 'computed, partial'
        else:
            hessian_line = 'computed, full'
        self.ostream.print_header(param('Hessian', hessian_line))
        if (self.do_hessian and self.calculate_partial_hessian
                and self.partial_hessian_cutoff is not None):
            self.ostream.print_header(
                param('Hessian cutoff', f'{self.partial_hessian_cutoff:.2f} A'))
        self.ostream.print_header(
            param('partial charges', 'RESP' if self.do_resp else 'D4'))
        self.ostream.print_header(
            param(
                'constrained atoms', 'beta carbons + caps'
                if self.constrain_capping_hydrogens else 'beta carbons'))
        self.ostream.print_header(
            param(
                'weak bridge pruning', f'> {util.WEAK_BRIDGE_TOLERANCE:.2f} A'
                if self.prune_weak_bridge_bonds else 'off'))
        self.ostream.print_header(
            param('atom typing',
                  'metal-blind' if self.metal_blind_typing else 'as bonded'))
        self.ostream.print_header(param('MPI ranks', self.nodes))
        self.ostream.print_header(param('output folder', self.folder))
        self.ostream.print_blank()
        self.ostream.flush()

    # ------------------------------------------------------------------
    # edits
    #
    # Every edit brings the builder up to date itself, so an edit is never
    # a request waiting for a rebuild that could be forgotten. Which way it
    # does that is the stage: before the fit it rebuilds the site, after it
    # it refits the terms. The individual methods do not restate this.
    # ------------------------------------------------------------------

    def add_metal_bond(self,
                       resid,
                       metal,
                       atom=None,
                       chain=None,
                       equilibrium=None):
        """
        Bonds a residue to a metal center that the distances did not connect.

        The cutoffs read an unrelaxed structure, so a contact the design
        intends can sit just outside them. The edit is recorded by residue
        index, atom name and metal residue index, none of which the
        renumbering of the protonation touches, so a rebuild and a
        re-detection on a relaxed geometry both put the bond back.

        With remove_metal_bond, the only edit still allowed once the force
        field exists.

        :param resid:
            The residue, as an id ('130' or 130) or as a label ('ASP130').
        :param metal:
            The atom index of the metal to bind to, as the active site labels
            it.
        :param atom:
            The name of the coordinating atom ('OD1'). Worked out from the
            geometry when there is no ambiguity.
        :param chain:
            The chain id, when the residue id occurs in more than one chain.
        :param equilibrium:
            The distance in Angstrom the crude pass should pull this bond to,
            rather than holding it where the structure has it. A bond is often
            added by hand because the sidechain is turned the wrong way round,
            and an equilibrium measured on that geometry only pins the mistake
            in place. Nothing after the crude pass is bound by it.
        """

        self._require('add_metal_bond', Stage.ACTIVE_SITE)

        # an edit is a decision, so it goes into the request and nowhere else
        self._request = self._sites().add_metal_bond(self._request,
                                                     self.binding_modes,
                                                     self._protonated_topology,
                                                     self._protonated_positions,
                                                     resid,
                                                     metal,
                                                     atom=atom,
                                                     chain=chain,
                                                     equilibrium=equilibrium)
        self._reapply()

    def remove_metal_bond(self, resid, metal=None, atom=None, chain=None):
        """
        Takes a residue's bond to a metal center back out.

        The counterpart of add_metal_bond, recorded the same way, so that a
        contact the cutoffs invented does not come back when the coordination
        is detected again on a relaxed geometry.

        Taking the last bond off a residue takes the residue out of the site
        with it, unless include_residue asked for it. Use remove_residue to
        say so outright.

        :param resid:
            The residue, as an id ('130' or 130) or as a label ('ASP130').
        :param metal:
            The atom index of the metal to unbind from, as the active site
            labels it. Only needed when the residue reaches more than one.
        :param atom:
            The name of the coordinating atom. Only needed when the residue
            binds that metal through more than one atom.
        :param chain:
            The chain id, when the residue id occurs in more than one chain.
        """

        self._require('remove_metal_bond', Stage.ACTIVE_SITE)

        self._request = self._sites().remove_metal_bond(self._request,
                                                        self.binding_modes,
                                                        resid,
                                                        metal=metal,
                                                        atom=atom,
                                                        chain=chain)
        self._reapply()

    def update_protonation_state(self, resid, variant, chain=None):
        """
        Sets the protonation variant of one residue of the active site.

        The variants are chosen automatically from the coordination -- a
        coordinating carboxylate is deprotonated, a coordinating cysteine is a
        thiolate, and the histidine tautomer puts no hydrogen on the
        coordinating nitrogen -- and this is how one of them is overruled, or
        how a residue that coordinates nothing is told what it should be.

        The request goes into protonation_overrides -- the setting itself,
        not a store beside it -- so it survives every later rebuild.
        Hydrogens are added and removed by protonating again from before any
        were placed, which is why the whole site is rebuilt: re-protonating
        an already protonated topology is the bug that guards against.

        Refused after build_forcefield: the charges and parameters of a force
        field describe the protonation it was fitted to.

        :param resid:
            The residue, as an id ('130' or 130) or as a label ('ASP130').
        :param variant:
            The variant to set, as OpenMM names it: 'ASP'/'ASH', 'GLU'/'GLH',
            'CYS'/'CYX', 'HID'/'HIE'/'HIP'/'HIN' or 'LYS'/'LYN'. Note that
            CYX here is the cysteine without its HG, which for a metal-bound
            sidechain is the thiolate.
        :param chain:
            The chain id, when the residue id occurs in more than one chain.
        """

        self._require('update_protonation_state', Stage.ACTIVE_SITE,
                      Stage.ACTIVE_SITE)

        self._protonation_overrides = self._sites().update_protonation_state(
            self._protonation_overrides, self._topology, resid, variant, chain)
        self._reapply()

    def include_residue(self, resid, chain=None):
        """
        Puts a residue into the truncated active site.

        The cluster is otherwise exactly the residues that coordinate a
        metal. This adds one that does not have to -- a second-shell residue
        that hydrogen bonds to a ligand, or one whose sidechain the QM should
        see for any other reason. It is truncated and capped like every other
        residue, and it keeps whatever protonation the pH gives it unless
        update_protonation_state says otherwise.

        To make a residue coordinate a metal instead, use add_metal_bond, or
        build_active_site(coordinating_residues=...) to do it during the
        detection.

        Only possible before build_forcefield, since the geometry, the
        Hessian and the charges all describe the cluster they were computed
        for.

        :param resid:
            The residue, as an id ('58' or 58) or as a label ('TYR58').
        :param chain:
            The chain id, when the residue id occurs in more than one chain.
        """

        self._require('include_residue', Stage.ACTIVE_SITE, Stage.ACTIVE_SITE)

        request = self._sites().include_residue(self._request,
                                                self.binding_modes,
                                                self._topology, resid, chain)
        if request is None:
            return

        self._request = request
        self._reapply()

    def remove_residue(self, resid, chain=None):
        """
        Takes a residue out of the truncated active site.

        Any metal bonds it makes go with it, recorded the way
        remove_metal_bond records them, so that neither the residue nor its
        coordination comes back when the site is detected again on a relaxed
        geometry.

        A metal that would be left with no ligand at all, and the last
        residue of the site, are refused: what is left would not be an active
        site.

        Only possible before build_forcefield, for the same reason as
        include_residue.

        :param resid:
            The residue, as an id ('130' or 130) or as a label ('ASP130').
        :param chain:
            The chain id, when the residue id occurs in more than one chain.
        """

        self._require('remove_residue', Stage.ACTIVE_SITE, Stage.ACTIVE_SITE)

        self._request = self._sites().remove_residue(
            self._request, self.binding_modes, self._protonated_topology,
            self._protonated_positions, resid, chain)
        self._reapply()

    # ------------------------------------------------------------------
    # applying an edit
    # ------------------------------------------------------------------

    def _reapply(self):
        """
        Brings everything downstream of an edit up to date, whichever way
        this builder's stage says that is done. The one place that branch is
        written, and every edit method's last act.

        Before a fit the site is protonated and truncated again from the
        request as it stands, the way build_active_site did it without a
        path, so an edit does not silently change whether the crude pass
        happens. After a fit nothing expensive runs: the geometry, the
        Hessian and the charges all still describe this cluster, and only
        which atoms the metals are bonded to has changed, so the force field
        is built again from them with the new connectivity. Rebuilding it
        rather than patching it is what keeps the angles, torsions and
        impropers that cross an edited bond right, since the generator
        derives every one of them from the connectivity matrix. A term the
        Hessian holds nothing for is fitted to zero and reported by
        _check_force_constants; the bond is kept, and the warning says to
        recompute the Hessian rather than that the edit was refused.
        """

        if self._stage < Stage.FITTED:
            # not a new run, so the header is not printed again
            state = self._build_active_site(self._topology, self._positions,
                                            self._mm_opt, False)
            self._adopt(state)
            # the Hessian and the charges were computed for the site that was
            # just replaced; a file in the folder that still fits this one
            # is picked back up by build_forcefield, which validates it first
            self._enter(Stage.ACTIVE_SITE)
            return

        self.ostream.print_info(
            'Fitting the metal terms again on the edited coordination. '
            'Nothing is recomputed; a term the Hessian does not cover is '
            'reported below.')
        self.ostream.flush()

        coordination = self._sites().site_coordination(
            self._protonated_topology, self._protonated_positions,
            self._active_site['molecule'], self._active_site, self._request)
        self._active_site = self._sites().apply_metal_bonds(
            self._active_site, coordination)

        # entering FITTED drops the enzyme system built from the old terms
        self._fit_forcefield(self._hessian, self._partial_charges)

    # ------------------------------------------------------------------
    # steps
    #
    # The pipeline's stages, usable on their own. Each drops the fit above
    # it -- a geometry, a Hessian or a set of charges computed after a fit
    # describes something that fit does not.
    # ------------------------------------------------------------------

    def mm_optimize_active_site(self):
        """
        Relaxes the active site again, on the best force field there is.

        build_active_site does this once already unless it was told not to.
        This is the way to do it again, or after the fact on a site that was
        extracted without it. The relaxed geometry replaces the one on the
        builder, and which force field it runs on is _crude_relax's decision.

        The one step that does not always drop the fit above it. A fit made
        here describes the geometry its Hessian was computed on, so a new
        geometry drops it; a fit adopted from a template describes the
        template, and relaxing toward its parameters brings the two closer
        together, so that one is kept. The enzyme system goes either way --
        the positions it was built from have moved.

        :return:
            The relaxed active site molecule.
        """

        self._require('mm_optimize_active_site', Stage.ACTIVE_SITE)

        # None before a fit, which is what makes _crude_relax build the
        # seeded force field it is named for
        forcefield = self._forcefield
        adopted = self._adopted_forcefield

        if self._protonated_topology is None:
            manual_equilibria = None
        else:
            manual_equilibria = self._sites().manual_bond_equilibria(
                self._active_site, self._protonated_topology,
                self.binding_modes)

        molecule = self._crude_relax(self._active_site,
                                     manual_equilibria,
                                     forcefield=forcefield)
        self._active_site['molecule'] = molecule

        if adopted:
            # the parameters are the template's, not this geometry's
            self._enter(Stage.FITTED)
        else:
            # a force field fitted here describes the geometry it replaced
            self._enter(Stage.ACTIVE_SITE)

        return molecule

    def _crude_relax(self, active_site, bond_equilibria, forcefield=None):
        """
        Relaxes a site on a force field, and writes the result.

        The pre-QM pass, run once by build_active_site and again by
        mm_optimize_active_site. With no force field it builds the seeded one
        it is named for: equilibria off the geometry, a flat default
        stiffness, which is all that is known before a Hessian exists.

        A force field is taken rather than built once one exists, and the
        difference is not cosmetic: a seeded equilibrium is measured on the
        very geometry the pass is trying to improve, so a contact the cutoffs
        left open sits at its own minimum and cannot be moved. A fitted or
        transferred force field carries per-bond equilibria and force
        constants, and pulls the same contact in.

        :param active_site:
            The site to relax. Not modified.
        :param bond_equilibria:
            Distances in nanometers to pull individual metal bonds to, from
            manual_bond_equilibria. Passed in rather than read off the
            builder, because build_active_site relaxes the site it has just
            extracted, before that site is the one the builder holds. Read
            only when the seeded force field is built here.
        :param forcefield:
            The force field to relax on. Defaults to building the seeded
            one for this site.

        :return:
            The relaxed molecule.
        """

        sites = self._sites()

        if forcefield is None:
            forcefield = sites.build_forcefield(active_site,
                                                util.d4_charges(active_site),
                                                bond_equilibria=bond_equilibria)

        relaxed = sites.mm_optimize_active_site(active_site, forcefield)

        self._save_intermediate(util.MM_GEOMETRY_FILE,
                                lambda path: relaxed.write_xyz_file(str(path)))

        return relaxed

    def optimize_geometry(self):
        """
        Optimizes the active site with the beta carbons frozen.

        Freezing them keeps the spatial arrangement imposed by the protein
        backbone. Without it the site relaxes to a gas-phase geometry, and
        since the Seminario method takes the equilibrium values straight from
        the geometry, every fitted bond length and angle would then describe
        the wrong structure.

        The optimized geometry replaces the one on the builder, and the
        coordination is detected again on it, since a contact can close or
        open during the relaxation.

        Collective: every rank runs it, and the drivers parallelize inside.

        :return:
            The results of the optimization driver.
        """

        self._require('optimize_geometry', Stage.ACTIVE_SITE)
        active_site = self._active_site

        optimized, opt_results = self._qm().optimize_active_site(active_site)

        self._save_intermediate(
            util.GEOMETRY_FILE,
            lambda path: optimized.write_xyz_file(str(path)))
        self._adopt_geometry(optimized)

        return opt_results

    def calculate_hessian(self):
        """
        Computes the nuclear Hessian of the active site.

        Restricted to the atom pairs the metal terms are fitted from unless
        calculate_partial_hessian is switched off, since Seminario reads one
        block per metal bond and two per metal angle and the fit throws the
        rest away. Those pairs are taken with partial_hessian_cutoff's
        forgiveness, so a contact just outside the perceived coordination is
        covered as well and can be bonded afterwards without recomputing.

        Collective: every rank runs it, and the drivers parallelize inside.

        :return:
            The Hessian.
        """

        self._require('calculate_hessian', Stage.ACTIVE_SITE)
        active_site = self._active_site

        qm = self._qm()
        atom_pairs, atoms = qm.hessian_pairs(active_site)
        n_atoms = active_site['molecule'].number_of_atoms()

        if self.calculate_partial_hessian:
            restriction = (f'Hessian restricted to {len(atom_pairs)} atom '
                           f'pairs over {len(atoms)} of {n_atoms} atoms')
            if self.partial_hessian_cutoff is not None:
                restriction += (', covering every donor atom within '
                                f'{self.partial_hessian_cutoff:.2f} A of a '
                                'metal whether it is bonded to one or not')
            self.ostream.print_info(restriction + '.')
        else:
            self.ostream.print_info(
                f'Computing the full Hessian over all {n_atoms} atoms; the '
                f'metal terms read {len(atom_pairs)} of its blocks.')
        self.ostream.flush()

        hessian = qm.compute_hessian(
            active_site, atom_pairs if self.calculate_partial_hessian else None)

        self._save_intermediate(util.HESSIAN_FILE,
                                lambda path: np.savetxt(path, hessian))

        # a fit above this was made from a different Hessian
        self._enter(Stage.ACTIVE_SITE)
        self._hessian = hessian

        return hessian

    def calculate_partial_charges(self):
        """
        Computes the partial charges of the active site.

        RESP when do_resp is on, and D4 otherwise. Whichever it is, the same
        charges reach the force field, the enzyme system and the charges file,
        so a run without RESP cannot end up with a force field whose charges
        never reach the protein.

        The RESP fit is collective; D4 costs a fraction of a millisecond and
        runs on the master.

        :return:
            The charges.
        """

        self._require('calculate_partial_charges', Stage.ACTIVE_SITE)
        active_site = self._active_site

        if self.do_resp:
            charges = self._qm().compute_resp_charges(active_site)
        else:
            charges = self._qm().d4_charges(active_site)

        self._save_intermediate(util.CHARGES_FILE,
                                lambda path: np.savetxt(path, charges))

        # a fit above this carries different charges
        self._enter(Stage.ACTIVE_SITE)
        self._partial_charges = charges

        return charges

    # ------------------------------------------------------------------
    # the fit and the enzyme system
    # ------------------------------------------------------------------

    def build_forcefield(self,
                         hessian=None,
                         opt_geometry=None,
                         partial_charges=None):
        """
        Fits the metal terms and builds the force field of the active site.

        This is where the expensive work is triggered. For each of the
        geometry, the Hessian and the charges the precedence is the same: what
        is passed in here, else the file an earlier run left in folder,
        else computing it -- the constrained optimization when
        do_qm_optimization is on, the Hessian when do_hessian is on, and the
        charges as RESP or D4. Each of the three is validated against the
        extracted active site before it is used, and which of the three it
        came from is announced by the resolver rather than guessed at here.

        :param hessian:
            A matrix, or the path to a text file readable by numpy.loadtxt, to
            use instead of computing one.
        :param opt_geometry:
            A molecule, or the path to an xyz file, to use instead of running
            the constrained optimization.
        :param partial_charges:
            Charges, or the path to a text file, to use instead of computing
            them.

        :return:
            The force field generator carrying the fitted metal terms.
        """

        self._require('build_forcefield', Stage.ACTIVE_SITE)

        geometry = self._qm()._resolve_optimized_geometry(
            self._active_site['molecule'], self.folder, opt_geometry)

        if geometry is not None:
            self._adopt_geometry(geometry)
        elif self.do_qm_optimization:
            self.optimize_geometry()
        else:
            self.ostream.print_info(
                'Skipping the constrained optimization; using the geometry '
                'that build_active_site left on the builder.')
            self.ostream.flush()

        hessian = self._qm()._resolve_hessian(self._active_site, self.folder,
                                              hessian)

        if hessian is not None:
            self._hessian = hessian
        elif self.do_hessian:
            hessian = self.calculate_hessian()
        else:
            self.ostream.print_info(
                'Skipping the Hessian calculation; using default force '
                'constants for the metal terms.')
            self.ostream.flush()

        charges = self._qm()._resolve_partial_charges(
            self._active_site, self.folder, partial_charges)

        if charges is not None:
            self._partial_charges = charges
        else:
            charges = self.calculate_partial_charges()

        self._qm().print_partial_charges(
            self._protonated_topology, self._active_site, charges,
            util.redistribute_cap_charges(self._active_site, charges), {
                residue.index: util.residue_label(residue)
                for residue in self._protonated_topology.residues()
            })

        return self._fit_forcefield(hessian, charges)

    def adopt_forcefield(self, forcefield, active_site=None):
        """
        Adopts a force field that was fitted elsewhere as this site's force
        field, without running any QM.

        How MetalForceFieldManager.build_ff_from_template hands over a force
        field transferred from a matching template -- metal terms and charges
        an earlier run already paid for -- so that create_enzyme_system can be
        called here afterwards. The manager's one way into the Stage
        machinery, so that it never reaches into _forcefield or _stage by
        hand.

        :param forcefield:
            The force field to adopt, already carrying the fitted metal
            terms and charges.
        :param active_site:
            The active site the force field was built for, when it differs
            from the one already on the builder (a template's coordination
            can force metal bonds the site did not have). Defaults to
            leaving the active site as it is.

        :return:
            The adopted force field.
        """

        self._require('adopt_forcefield', Stage.ACTIVE_SITE)

        if active_site is not None:
            self._active_site = active_site

        self._forcefield = forcefield
        self._adopted_forcefield = True
        if forcefield.partial_charges is not None:
            self._partial_charges = np.asarray(forcefield.partial_charges)

        self._enter(Stage.FITTED)

        return self._forcefield

    def create_enzyme_system(self):
        """
        Injects the fitted metal terms into a force field system for the whole
        enzyme.

        The protein force field already covers everything except the metal, so
        only the metal bonds and angles are transferred, and the atom map of
        the active site gives the correspondence directly. The charges of the
        coordination sphere are replaced by the ones the force field carries,
        with the region charge restored over the local backbone.

        Nothing expensive is triggered: it is an error to call this before
        build_forcefield.

        :return:
            The tuple of the OpenMM system and the topology it was built for.
        """

        self._require('create_enzyme_system', Stage.FITTED)

        system, _ = self._enzyme().create_enzyme_system(
            self._protonated_topology, self._active_site, self._forcefield,
            self._partial_charges)

        self._save_intermediate(
            util.ENZYME_SYSTEM_FILE,
            lambda path: path.write_text(mm.XmlSerializer.serialize(system)))
        # written again beside the system, unchanged, so that the last step
        # to touch the folder leaves it complete whichever way the run got
        # here
        self._save_intermediate(
            util.FORCEFIELD_FILE,
            lambda path: MMForceFieldGenerator.save_forcefield_as_json(
                self._forcefield, str(path)))

        self._enzyme_system = system
        self._enter(Stage.ENZYME)

        return self._enzyme_system, self._protonated_topology

    def create_enzyme_forcefield(self):
        """
        Writes the fitted metal site as an OpenMM force field XML.

        The counterpart of create_enzyme_system and what a simulation
        should be built from: that one puts the fitted terms onto one
        System, which describes the topology it was built for and nothing
        else, while this writes a force field that can be loaded beside
        the protein force field to build whatever system is wanted --
        solvated, extended, or rebuilt after any change to the structure.

        The active site is moved into a residue of its own and the
        metal-ligand contacts become real bonds of the topology, so the
        system built from the file carries the 1-2 and 1-3 exclusions and
        the 1-4 scaling that a bonded metal model should have. Everything
        else -- every covalent term, every charge, every Lennard-Jones
        parameter -- comes out exactly as create_enzyme_system leaves it.

        Nothing expensive is triggered: it is an error to call this before
        build_forcefield.

        :return:
            A dictionary holding the XML, the restructured topology, its
            positions in Angstrom and the residue templates.
        """

        self._require('create_enzyme_forcefield', Stage.FITTED)

        result = self._enzyme().create_enzyme_forcefield(
            self._protonated_topology, self._protonated_positions,
            self._active_site, self._forcefield, self._partial_charges)

        # the XML and the topology it is for are written together, since
        # neither says anything without the other: the templates are keyed
        # on the chain and residue id of a topology whose active site has
        # been moved into one residue
        self._save_intermediate(openmmxml.SITE_XML_FILE,
                                lambda path: path.write_text(result['xml']))
        self._save_intermediate(
            openmmxml.SITE_TOPOLOGY_FILE,
            lambda path: mmapp.PDBFile.writeFile(
                result['topology'],
                # a bare list of numbers is written as Angstrom whatever
                # the values mean, so the unit goes on here
                result['positions'] * mmunit.angstrom,
                str(path),
                # the templates are keyed on chain and residue id, and
                # PDBFile numbers the residues from one unless told to keep
                # them, which would leave the file unable to load its own
                # topology back
                keepIds=True))

        # what restructure_topology returned is the phase's own business
        self._enzyme_forcefield = {
            key: result[key]
            for key in ('xml', 'topology', 'positions', 'templates',
                        'backbone_shift')
        }
        self._enter(Stage.ENZYME)

        return self._enzyme_forcefield

    # ------------------------------------------------------------------
    # the working folder
    # ------------------------------------------------------------------

    def _save_intermediate(self, name, writer):
        """
        Writes one intermediate as soon as the step that produced it finishes.

        Saving as the run goes rather than at the end means a run that dies
        part way through still leaves behind everything that did succeed, so a
        geometry optimization is not paid for twice because the Hessian after
        it failed.

        :param name:
            The file name, one of the file name constants of util.
        :param writer:
            A callable taking the path to write.
        """

        if self.rank != mpi_master():
            return

        folder = Path(self.folder)
        folder.mkdir(parents=True, exist_ok=True)
        writer(folder / name)

        self.ostream.print_info(f'Wrote {name} to {self.folder}')
        self.ostream.flush()

    # ------------------------------------------------------------------
    # stage machinery
    # ------------------------------------------------------------------

    def _enter(self, stage):
        """
        Records that the builder has reached a stage, and drops everything
        that belonged to a later one.

        The only place downstream state is invalidated. A step calls this
        once it has what the stage owns, and whatever described the site it
        replaced goes automatically, so no step has to remember what it
        made stale.

        Entering a stage the builder is already past is how a step that
        produces a new geometry, Hessian or set of charges says that the fit
        above it no longer describes them.

        :param stage:
            The stage now reached.
        """

        for level, names in STAGE_FIELDS.items():
            if level <= stage:
                continue
            for name in names:
                setattr(self, name, None)

        if self._forcefield is None:
            self._adopted_forcefield = False

        self._stage = stage

    def _require(self, method, minimum, maximum=None):
        """
        Checks that the workflow is far enough along for an operation, and
        not too far.

        One guard for the whole class: which operations are legal is a
        function of the stage, so there is one place that asks and one place
        that says why not.

        :param method:
            The name of the method that needs it.
        :param minimum:
            The earliest stage the operation makes sense at.
        :param maximum:
            The latest, for the operations the fit closes off.
        """

        # named from where the builder is rather than from what was asked
        # for, since what the caller needs is the next step to take, and
        # the step after that would fail for the same reason
        next_step = {
            Stage.EMPTY:
            'there is no active site yet. Call build_active_site first.',
            Stage.ACTIVE_SITE:
            'there is no force field yet. Call build_forcefield first.',
        }

        assert_msg_critical(
            self._stage >= minimum, f'MetalSiteForceFieldBuilder.{method}: '
            f'{next_step.get(self._stage, f"the workflow is at {self.stage}.")}'
        )

        assert_msg_critical(
            maximum is None or self._stage <= maximum,
            f'MetalSiteForceFieldBuilder.{method}: the force field has '
            'already been built, and its charges and parameters describe the '
            'protonation and the residues it was fitted to. Only '
            'add_metal_bond and remove_metal_bond can be used after '
            'build_forcefield. Call build_active_site() again to start over '
            'from the request, which drops the fit.')

    def _adopt(self, state):
        """
        Takes on what the structural pass produced.

        The one place the builder's structural state is written, so that
        build_active_site and a rebuild after an edit cannot come to write
        different sets of it.

        The Hessian and the charges go with it: they were computed for the
        site being replaced, and this is the one place a site is replaced.
        A file in the working folder that still fits the new one is picked
        back up by build_forcefield, which validates it before using it.

        :param state:
            What _build_active_site returned.
        """

        self._request = state['request']
        self._topology = state['topology']
        self._positions = state['positions']
        self._protonated_topology = state['protonated_topology']
        self._protonated_positions = state['protonated_positions']
        self._active_site = state['active_site']

        self._hessian = None
        self._partial_charges = None

    def _fit_forcefield(self, hessian, charges):
        """
        Builds the seeded force field and fits its metal terms.

        Shared by the first fit and by every refit after a bond edit, so
        the two cannot come to write different artifacts or print different
        tables. With no Hessian the seeded force field is the result: the
        metal terms keep their default force constants.

        :param hessian:
            The Hessian to fit from.
        :param charges:
            The charges the force field carries.

        :return:
            The force field, on every rank.
        """

        active_site = self._active_site

        # the seeded force field the fit starts from, or, with no Hessian,
        # the force field itself: default force constants on the metal terms
        forcefield = self._sites().build_forcefield(active_site, charges)

        if hessian is not None:
            # A bond asked for by hand is a decision, and the two things the
            # weak bridge pruning reads -- a zero force constant and a long
            # distance -- are exactly what such a bond looks like when the
            # Hessian does not cover it, so the pruning is told to leave it.
            if self._protonated_topology is None:
                protected = set()
            else:
                protected = self._sites().manual_bond_keys(
                    active_site, self._protonated_topology,
                    self.binding_modes)

            forcefield = self._qm().fit_forcefield(active_site,
                                                   forcefield,
                                                   hessian,
                                                   protected_bonds=protected)

        self._qm().print_metal_parameters(
            active_site, forcefield, util.get_metal_keys(forcefield,
                                                         active_site))

        # Everything the fit was made of and everything it produced, written
        # even when it was handed in or read back, so that the folder holds
        # the run rather than only the parts that happened to be computed.
        # The charges written are the ones the force field ended up carrying:
        # build_forcefield folds the capping hydrogens' charge into the rest
        # of the site, and writing the raw fit beside a force field carrying
        # the corrected one left two files disagreeing about what "the
        # charges" are. The force field goes out as JSON, which is what a
        # template reads.
        corrected = forcefield.partial_charges
        if corrected is None:
            corrected = charges
        self._save_intermediate(
            util.GEOMETRY_FILE,
            lambda path: active_site['molecule'].write_xyz_file(str(path)))
        if hessian is not None:
            self._save_intermediate(util.HESSIAN_FILE,
                                    lambda path: np.savetxt(path, hessian))
        self._save_intermediate(util.CHARGES_FILE,
                                lambda path: np.savetxt(path, corrected))
        self._save_intermediate(
            util.FORCEFIELD_FILE,
            lambda path: MMForceFieldGenerator.save_forcefield_as_json(
                forcefield, str(path)))

        self._forcefield = forcefield
        self._adopted_forcefield = False

        # The weak bridge pruning is the one step that can decide a metal
        # contact is not a bond after all, and it decides it on the force
        # field. Lifting that back onto the site is what stops the two
        # disagreeing about what the cluster is bonded like -- see
        # ActiveSiteBuilder.connectivity_from_forcefield.
        self._active_site = self._sites().connectivity_from_forcefield(
            self._active_site, self._forcefield)

        # The fit folds the capping hydrogens' charge into the rest of the
        # site, and from here on that is what "the charges" means: the same
        # array reaches the force field, the enzyme system, the file and the
        # partial_charges property, so no two of them can answer differently.
        # Before a fit they are the raw result of the charge calculation,
        # which is what there is to have.
        if self._forcefield.partial_charges is not None:
            self._partial_charges = np.asarray(self._forcefield.partial_charges)

        self._enter(Stage.FITTED)

        return self._forcefield

    # ------------------------------------------------------------------
    # state updates
    # ------------------------------------------------------------------

    def _adopt_geometry(self, molecule):
        """
        Puts a new geometry on the active site and detects the coordination
        again on it.

        :param molecule:
            The new geometry of the active site.
        """

        self._active_site = self._sites().adopt_geometry(
            self._protonated_topology, self._protonated_positions,
            self._active_site, molecule, self._request)

        # a force field fitted before this describes the geometry it replaced
        self._enter(Stage.ACTIVE_SITE)

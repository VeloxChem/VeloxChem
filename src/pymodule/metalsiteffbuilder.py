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
from copy import deepcopy
from enum import IntEnum
import numpy as np
import time
import sys

from .veloxchemlib import mpi_master
from .outputstream import OutputStream
from .errorhandler import assert_msg_critical
from . import metalsitecore as core

try:
    import openmm as mm
    import openmm.app as mmapp
    import openmm.unit as mmunit
except ImportError:
    pass


class Stage(IntEnum):
    """
    How far along the workflow a builder is, and therefore what it will do.

    Ordered, so that a requirement reads as a comparison. Which operations
    are legal is a function of this and nothing else: there is no second
    way to ask, and no inferring the answer from which field happens to be
    None.

    - EMPTY: nothing built yet.
    - ACTIVE_SITE: a truncated cluster exists and is fully editable.
    - FITTED: the QM has been paid for and the metal terms are fitted. Only
      the coordination is still open, because a bond edit can be refitted
      from the geometry, Hessian and charges that are already there, while
      the protonation and the residue membership are what those describe.
    - ENZYME: the fitted terms are in a system for the whole enzyme. Not
      terminal: a bond edit still refits, and drops the system it invalidates.
    """

    EMPTY = 0
    ACTIVE_SITE = 1
    FITTED = 2
    ENZYME = 3


# What each stage owns, and therefore what entering a lower one throws away.
# The single invalidation rule of the class: nothing is nulled by hand
# anywhere else, so a step cannot forget to drop what it invalidated.
#
# The Hessian and the charges are deliberately not here. They are inputs to
# the fit rather than the fit itself, and the public steps fill them one at a
# time, so clearing them whenever one of them is recomputed would wipe the
# other. What invalidates them is a different active site, and _adopt -- the
# one place a site is taken on -- is what drops them.
STAGE_FIELDS = {
    Stage.ACTIVE_SITE: (),
    Stage.FITTED: ('_forcefield',),
    Stage.ENZYME: ('_enzyme_system',),
}


class MetalSiteForceFieldBuilder:
    """
    Builds a bonded force field for the zinc center of a metallo-enzyme.

    Identifies the coordination sphere, fixes the protonation of the
    coordinating residues, truncates a QM active site at the CA-CB bonds, and
    fits the metal-ligand bond and angle parameters to a QM Hessian with the
    Seminario method. The fitted parameters can then be injected into a
    protein force field for the whole enzyme.

    The pipeline is three calls:

        molecule = builder.build_active_site('site.cif')
        forcefield = builder.build_forcefield()
        system, topology = builder.create_enzyme_system()

    The active site can be edited in between them, and how much of it can be
    edited is what separates the two gaps.

    Between the first and the second nothing expensive has run yet, so the
    site can be changed in every way it is described by: add_metal_bond and
    remove_metal_bond correct a coordination the distances got wrong,
    update_protonation_state overrules a variant, and include_residue and
    remove_residue change which residues the cluster is made of. Each of
    them extracts the site again on its own, so an edit is never a request
    waiting for a rebuild that might be forgotten.

    Between the second and the third only the coordination is still open.
    The charges and the fitted parameters describe the protonation and the
    residues they were fitted to, so those two are refused with an error
    naming the way back -- build_active_site again, which drops the fit.
    A bond edit is fitted again on the geometry, the Hessian and the charges
    that are already there, which costs nothing; a bond the Hessian holds no
    data for is fitted to zero and reported as such, and kept.

    The object holds the settings and the intermediates; the work itself is
    done by the functions of metalsitecore, which keep no state and take
    everything they use. What a step produced is read back off a read-only
    property rather than juggled by the caller.

    Under MPI everything cheap runs on the master rank and is broadcast, and
    the three expensive QM steps are collective: every rank calls them, at the
    same point, with the same molecule.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - metal_bond_cutoff: The distance in Angstrom within which a donor atom
          is taken to be bonded to a metal center.
        - report_cutoff: The distance in Angstrom out to which contacts are
          reported for review.
        - metal_elements: The elements treated as metal centers.
        - metal_formal_charges: The formal charges assumed for the metal ions.
        - bidentate_asymmetry: How much closer one carboxylate oxygen has to be
          than the other before the pair stops counting as chelating.
        - cap_bond_length: The C-H distance in Angstrom of the capping
          hydrogen that replaces CA.
        - protonation_overrides: The residue to variant mapping that overrides
          the automatic protonation rules. A residue is named by its id
          ('187'), by its label ('HIS187') or by its index, and a label is
          worth preferring: ids and indices overlap, so an index can name a
          residue by one reading and its neighbour by the other.
          update_protonation_state writes into this, by label, so a variant
          set by hand and one set as a setting are the same thing.
        - prepare_protein: The flag for repairing the structure with pdbfixer,
          which is what a protein force field needs in order to match
          templates.
        - scf_drv: The SCF driver, used as given. Created from xcfun if not
          provided. Set convergence thresholds and iteration limits on it.
        - xcfun: The exchange-correlation functional used when scf_drv is not
          provided. None runs Hartree-Fock.
        - basis_set_label: The basis set label.
        - mute_scf: The flag for muting the output of the QM drivers.
        - do_qm_optimization: The flag for optimizing the active site before the
          Hessian is computed.
        - do_resp: The flag for computing RESP charges. D4 charges are used
          when it is off.
        - calculate_partial_hessian: The flag for restricting the Hessian to
          the atom pairs the metal terms are fitted from. False computes the
          whole thing, which costs far more and changes no parameter.
        - constrain_capping_hydrogens: The flag for constraining the capping
          hydrogens in addition to the beta carbons.
        - average_metal_terms: The flag for averaging the fitted metal terms
          over equivalent atoms.
        - prune_weak_bridge_bonds: The flag for dropping the long arm of a
          bridging residue when the fit gave it no force constant.
        - weak_bridge_tolerance: How much longer than the residue's shortest
          metal bond, in Angstrom, the unfitted arm has to be before it is
          dropped.
        - reparameterize_metal_angles: The flag for touching the metal angles
          at all. False leaves every one of them at the value the generator
          guessed, in the crude pass and in the Hessian fit alike.
        - mm_max_iterations: The iteration limit of the crude relaxation. Zero
          runs until convergence.
        - mm_constrain_metals: The flag for holding the metal centers fixed in
          the crude pass as well, on top of the beta carbons.
        - mm_bond_change_warning: How far a metal-ligand bond may move in the
          crude pass, in Angstrom, before it is reported.
        - default_metal_bond_force_constant: The force constant given to every
          metal bond of the crude pass, in kJ/mol/nm^2.
        - default_metal_angle_force_constant: The force constant given to every
          metal angle of the crude pass, in kJ/mol/rad^2.
        - metal_bond_equilibria: Equilibrium metal-ligand distances in nm,
          keyed by element pair in either order, replacing the values measured
          on the input geometry. LITERATURE_METAL_BONDS is such a table.
        - metal_angle_equilibria: Equilibrium metal angles in degrees, keyed by
          element triple in either order, replacing the measured values.
        - protein_forcefield_files: The OpenMM force field files the enzyme
          system is built from.
        - output_folder: The folder every step writes its result to as soon as
          it has it, and reads back on a later run. Named after the creation
          time by default, so runs do not collide.
        - comm: The MPI communicator.
        - rank: The rank of the MPI process.
        - nodes: The number of MPI processes.
        - ostream: The output stream.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes the metal site force field builder.
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

        # coordination detection
        # the primary cutoff is generous on purpose: a stretched bridging
        # contact in an unrelaxed structure is still a bond
        self.metal_bond_cutoff = 3.0
        self.report_cutoff = 3.5
        self.metal_elements = tuple(core.METAL_ELEMENTS)
        self.metal_formal_charges = dict(core.METAL_FORMAL_CHARGES)

        # a carboxylate this much more lopsided than BIDENTATE_ASYMMETRY is
        # not chelating, whatever the second oxygen looks like on paper
        self.bidentate_asymmetry = core.BIDENTATE_ASYMMETRY

        # truncation and protonation
        self.cap_bond_length = 1.09
        self.protonation_overrides = None
        # the repair renumbers atoms and is what lets a protein force field
        # match its templates; a run that only wants the metal site can skip
        # it, and skip pdbfixer with it
        self.prepare_protein = True

        # QM settings. Anything beyond the functional and the basis set is
        # set on an SCF driver assigned to scf_drv, which is used as given.
        self.scf_drv = None
        self.xcfun = 'PBE0'
        self.basis_set_label = 'def2-svp'
        self.mute_scf = True

        # workflow
        self.do_qm_optimization = True
        self.do_resp = True
        # Only the metal terms are fitted, and Seminario reads one block of
        # the Hessian per bond and two per angle, so the analytical Hessian is
        # restricted to exactly those pairs. Everything else it would have
        # computed is thrown away by the fit, which is why this is on. Turn it
        # off for a Hessian that is worth something on its own -- a frequency
        # analysis, or a fit of terms the metals do not reach.
        self.calculate_partial_hessian = True
        # the beta carbon is where the backbone actually holds the sidechain;
        # the capping hydrogen only stands in for the alpha carbon
        self.constrain_capping_hydrogens = False
        # two topologically equivalent ligands of a strained site sit at
        # genuinely different distances, and that asymmetry is the signal
        self.average_metal_terms = False

        # A residue that bridges two metals, through one atom or through two,
        # can reach one of them far more weakly than the other. Where the fit
        # says so outright, by giving that arm no force constant at all, and
        # the geometry agrees by holding it this much further out than the
        # residue's shortest metal bond, it is dropped rather than kept as a
        # bond with no stiffness.
        self.prune_weak_bridge_bonds = True
        self.weak_bridge_tolerance = 0.25

        # Crude MM pass. The site is relaxed on a force field of its own
        # before any QM is run, so the optimization does not spend its first
        # cycles undoing whatever the structure file happened to hold. Nothing
        # at that point says anything about the stiffness of a metal term, so
        # the force constants are flat defaults and only the equilibrium
        # values carry information: they are measured on the input geometry
        # unless one of the tables below overrides them.
        self.mm_max_iterations = 0
        self.default_metal_bond_force_constant = (
            core.DEFAULT_METAL_BOND_FORCE_CONSTANT)
        self.default_metal_angle_force_constant = (
            core.DEFAULT_METAL_ANGLE_FORCE_CONSTANT)
        self.metal_bond_equilibria = None
        self.metal_angle_equilibria = None
        self.reparameterize_metal_angles = True
        # how far a metal-ligand bond may move in the crude pass before it is
        # worth saying so, in Angstrom
        self.mm_bond_change_warning = 0.25
        # Holding the metals as well keeps the metal-metal distance and the
        # shape of the site exactly as the structure had them, and lets the
        # crude pass tidy up the ligands around a fixed frame.
        self.mm_constrain_metals = False

        # the protein force field the fitted metal terms are added to
        self.protein_forcefield_files = ('amber14-all.xml',
                                         'amber14/tip3pfb.xml')

        # Output. The folder is named after the time the builder was created,
        # so two runs never overwrite one another; assign output_folder to put
        # the files somewhere of your own choosing, or to point at an earlier
        # run and reuse whatever it managed to finish.
        self.output_folder = f'metal_site_{int(time.time())}'

        # State. Two kinds of thing, and the difference is the whole of it.
        #
        # The request is what has been decided: the metal bonds set by hand,
        # the protonation variants, and which residues were added to or taken
        # out of the cluster. It is keyed on residue indices, labels and atom
        # names, none of which adding hydrogens disturbs, so it outlives every
        # renumbering and is replayed onto every derivation.
        #
        # The coordination is not kept at all. It is a function of a geometry
        # and the request, so it is derived again whenever the answer could
        # have changed, and the atom indices in it always belong to the
        # topology it was derived on. That is what removes the old pair of
        # pre- and post-protonation copies and the translation between them.
        #
        # The rest is geometry and the results of the expensive steps. The
        # prepared topology is kept apart from the protonated one because
        # build_active_site rebuilds from the prepared one, and protonating an
        # already protonated topology is the bug that guards against.
        self._stage = Stage.EMPTY
        self._request = core.site_request()
        self._mm_opt = True
        self._topology = None
        self._positions = None
        self._protonated_topology = None
        self._protonated_positions = None
        self._active_site = None
        self._forcefield = None
        self._hessian = None
        self._partial_charges = None
        self._enzyme_system = None

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
    def enzyme_topology(self):
        """
        The protonated topology the enzyme system was built for.
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
        Draws the active site as it now stands.

        The bonds and the atom labels are the site's own rather than
        perceived from the geometry: a metal-ligand bond is a decision that
        nothing perceives, and the labels are what add_metal_bond and
        remove_metal_bond want to be told about an atom, so a coordination
        edit can be read straight off the picture.

        Every other keyword reaches Molecule.show untouched.

        :param kwargs:
            Further keyword arguments for Molecule.show.

        :return:
            Whatever Molecule.show returns.
        """

        self._require('show_active_site', Stage.ACTIVE_SITE)

        return core.show_active_site(self._active_site, **kwargs)

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

        return core.freeze_constraints(
            self._active_site,
            constrain_capping_hydrogens=self.constrain_capping_hydrogens)

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

        return self._derive(self._protonated_topology,
                            self._protonated_positions)

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
        Builds the truncated active site, from a structure or from the
        binding modes as they now stand.

        With a path, the structure is read and repaired, the coordination
        sphere is detected and reported, and everything a previous call
        produced is dropped. Without one, the detection is not run again: the
        binding modes already on the builder are used, which is what makes a
        correction with add_metal_bond or remove_metal_bond survive a rebuild.

        The coordinating residues are then protonated, the site is truncated
        at the CA-CB bonds, and unless mm_opt is switched off it is relaxed on
        a crude force field of its own so that a later QM optimization does not
        spend its first cycles cleaning up the structure file.

        The edit methods rebuild the site themselves, so calling this again
        by hand is only needed to start from a different structure, to change
        mm_opt, or to reopen editing after build_forcefield.

        :param cif_path:
            The path to a .pdb, .cif or .pdbx file, or None to rebuild from
            the binding modes of the previous call.
        :param mm_opt:
            Whether to relax the extracted site on a crude force field.
        :param coordinating_residues:
            Residues that coordinate a metal however far out they sit, as ids
            ('130') or labels ('ASP130'). Only used when a path is given.

        :return:
            The active site molecule.
        """

        assert_msg_critical('openmm' in sys.modules,
                            'MetalSiteForceFieldBuilder: openmm is required')

        state = self._on_master(lambda: self._build_active_site(
            cif_path, mm_opt, coordinating_residues))

        # an edit rebuilds the site the way this call built it, so what it
        # was told has to outlive the call
        self._mm_opt = mm_opt

        self._adopt(state)

        # everything downstream belonged to the site that was just replaced
        self._enter(Stage.ACTIVE_SITE)

        return self.active_site_molecule

    def _build_active_site(self, cif_path, mm_opt, coordinating_residues):
        """
        The work of build_active_site, on one rank.

        :return:
            The state build_active_site adopts and broadcasts.
        """

        if cif_path is not None:
            self._print_header(cif_path)
            topology, positions = core.load_and_prepare_protein(
                cif_path, prepare=self.prepare_protein)
        else:
            assert_msg_critical(
                self._topology is not None,
                'MetalSiteForceFieldBuilder.build_active_site: there is no '
                'structure to rebuild from. Call it with a cif_path first.')
            topology = self._topology
            positions = self._positions

        binding_modes = self._derive(
            topology, positions, coordinating_residues=coordinating_residues)

        if cif_path is not None:
            core._print_binding_modes(binding_modes, ostream=self.ostream)

        protonated_topology, protonated_positions, variants, notes = (
            core.protonate(topology,
                           positions,
                           binding_modes,
                           protonation_overrides=self.protonation_overrides))

        # what the protonation actually built is a decision, so it joins the
        # request; the coordination is then derived again on the protonated
        # topology, where the atom indices are the ones everything downstream
        # uses. Nothing is remapped: a derivation is always about the
        # topology it was handed.
        self._request['variants'] = variants
        protonated_modes = self._derive(protonated_topology,
                                        protonated_positions)
        core._merge_notes(protonated_modes, notes)

        active_site = core.extract_active_site(
            protonated_topology,
            protonated_positions,
            protonated_modes,
            cap_bond_length=self.cap_bond_length,
            ostream=self.ostream)
        core._print_active_site(active_site,
                                protonated_modes,
                                ostream=self.ostream)

        self._save_intermediate(
            'protonated.pdb', lambda path: self._write_pdb(
                path, protonated_topology, protonated_positions))

        if mm_opt:
            active_site['molecule'] = self._crude_relax(
                active_site,
                core.manual_bond_equilibria(active_site, protonated_topology,
                                            protonated_modes))

        return {
            'request': self._request,
            'topology': topology,
            'positions': positions,
            'protonated_topology': protonated_topology,
            'protonated_positions': protonated_positions,
            'active_site': active_site,
        }

    def add_metal_bond(self,
                       resid,
                       metal,
                       atom=None,
                       chain=None,
                       equilibrium=None):
        """
        Bonds a residue to a metal center that the distances did not connect.

        The cutoffs read an unrelaxed structure, so a contact the design
        intends can sit just outside them. The edit is recorded in the binding
        modes by residue index, atom name and metal residue index, none of
        which the renumbering of the protonation touches, so that a rebuild
        and a re-detection on a relaxed geometry both put the bond back.

        The active site is brought up to date on its own: before
        build_forcefield it is extracted again from the corrected
        coordination, and afterwards the metal terms are fitted again on the
        geometry, the Hessian and the charges that are already there. It is
        the one edit that is still allowed once the force field exists.

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
            The distance in Angstrom the crude pass should pull this bond
            to, rather than holding it where the structure has it. A bond is
            often added by hand because the sidechain is turned the wrong
            way round, and an equilibrium measured on that geometry only
            pins the mistake in place; naming a distance swings a histidine
            into place before any QM is paid for. The QM optimization and
            the Hessian fit that follow are not bound by it.
        """

        self._require('add_metal_bond', Stage.ACTIVE_SITE)

        def edit(request, modes):
            return core.add_metal_bond(request,
                                       modes,
                                       self._protonated_topology,
                                       self._protonated_positions,
                                       resid,
                                       metal,
                                       atom=atom,
                                       chain=chain,
                                       equilibrium=equilibrium,
                                       ostream=self.ostream,
                                       **self.detection_settings())

        self._edit_coordination(edit)

    def remove_metal_bond(self, resid, metal=None, atom=None, chain=None):
        """
        Takes a residue's bond to a metal center back out.

        The counterpart of add_metal_bond, and recorded the same way, so that
        a contact the cutoffs invented does not come back when the
        coordination is detected again on a relaxed geometry. The active site
        is brought up to date on its own, and this too is still allowed once
        the force field exists.

        Note that taking the last bond off a residue takes the residue out of
        the active site with it, unless include_residue asked for it. Use
        remove_residue to say so outright.

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

        def edit(request, modes):
            return core.remove_metal_bond(
                request,
                modes,
                resid,
                metal=metal,
                atom=atom,
                chain=chain,
                bidentate_asymmetry=self.bidentate_asymmetry,
                ostream=self.ostream)

        self._edit_coordination(edit)

    def update_protonation_state(self, resid, variant, chain=None):
        """
        Sets the protonation variant of one residue of the active site.

        The variants are chosen automatically from the coordination -- a
        coordinating carboxylate is deprotonated, a coordinating cysteine is a
        thiolate, and the histidine tautomer puts no hydrogen on the
        coordinating nitrogen -- and this is how one of them is overruled, or
        how a residue that coordinates nothing is told what it should be.

        The request goes into protonation_overrides, so it is the same thing
        as the setting and survives every later rebuild. Hydrogens are added
        and removed by protonating the structure again from before any
        hydrogens were placed, which is also why the whole active site is
        rebuilt: nothing expensive has run at this point, and re-protonating
        an already protonated topology is not something to reach for.

        Only possible before build_forcefield, since the charges and the
        parameters of the force field describe the protonation it was fitted
        to.

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

        def work():
            residue = core._resolve_residue(self._topology, resid, chain)
            core.check_variant(residue, variant)

            # keyed by the label rather than the index, because residue ids
            # and residue indices overlap: a single chain numbered from one
            # has index i for id i+1, so an index written here would also
            # match its neighbour by id
            overrides = dict(self.protonation_overrides or {})
            overrides[f'{residue.name}{residue.id}'] = variant

            self.ostream.print_info(
                f'{residue.name}{residue.id} will be protonated as {variant}.')
            self.ostream.flush()

            return overrides

        self.protonation_overrides = self._on_master(work)
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

        self._require('include_residue', Stage.ACTIVE_SITE,
                      Stage.ACTIVE_SITE)

        def work():
            residue = core._resolve_residue(self._topology, resid, chain)
            core.check_truncatable(residue)

            request = self._request
            label = f'{residue.name}{residue.id}'

            already = residue.index in core.active_site_residues(
                self.binding_modes)
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

        request = self._on_master(work)

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

        self._require('remove_residue', Stage.ACTIVE_SITE,
                      Stage.ACTIVE_SITE)

        def work():
            residue = core._resolve_residue(self._topology, resid, chain)

            request = self._request
            modes = self.binding_modes
            label = f'{residue.name}{residue.id}'
            members = core.active_site_residues(modes)

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

            # the bonds have to go through remove_metal_bond, so that the
            # removals are recorded and a re-detection does not put back what
            # the geometry still looks like
            while any(ligand['res_index'] == residue.index
                      for ligand in modes['ligands']):
                bound = next(ligand for ligand in modes['ligands']
                             if ligand['res_index'] == residue.index)
                request = core.remove_metal_bond(
                    request,
                    modes,
                    label,
                    metal=bound['metals'][0],
                    bidentate_asymmetry=self.bidentate_asymmetry,
                    ostream=self.ostream)
                modes = self._derive(self._protonated_topology,
                                     self._protonated_positions,
                                     request=request)

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

        self._request = self._on_master(work)
        self._reapply()

    def mm_optimize_active_site(self):
        """
        Relaxes the active site again on a crude force field of its own.

        build_active_site does this once already unless it was told not to.
        This is the way to do it again, or to do it after the fact on a site
        that was extracted without it. The relaxed geometry replaces the one
        on the builder.

        :return:
            The relaxed active site molecule.
        """

        self._require('mm_optimize_active_site', Stage.ACTIVE_SITE)

        molecule = self._on_master(lambda: self._crude_relax(
            self._active_site, self._manual_equilibria()))
        self._active_site['molecule'] = molecule

        # a force field fitted before this describes the geometry it replaced
        self._enter(Stage.ACTIVE_SITE)

        return molecule

    def _crude_relax(self, active_site, bond_equilibria=None):
        """
        Relaxes a site on a force field of its own, and writes the result.

        The pre-QM pass, run once by build_active_site and again by
        mm_optimize_active_site. The force field is built here rather than
        taken, because this is the one caller that wants the seeded one: the
        equilibria come off the geometry and the stiffness is a flat default,
        which is all that is known before a Hessian exists.

        Runs on one rank; the caller broadcasts.

        :param active_site:
            The site to relax. Not modified.
        :param bond_equilibria:
            Distances in nanometers to pull individual metal bonds to, from
            manual_bond_equilibria. Passed in rather than read off the
            builder, because build_active_site relaxes the site it has just
            extracted, before that site is the one the builder holds.

        :return:
            The relaxed molecule.
        """

        forcefield = core.build_forcefield(
            active_site,
            comm=MPI.COMM_SELF,
            ostream=self.ostream,
            bond_equilibria=bond_equilibria,
            **self.fit_settings())
        relaxed = core.mm_optimize_active_site(
            active_site,
            forcefield,
            constrain_metals=self.mm_constrain_metals,
            constrain_capping_hydrogens=self.constrain_capping_hydrogens,
            max_iterations=self.mm_max_iterations,
            bond_change_warning=self.mm_bond_change_warning,
            ostream=self.ostream)

        self._save_intermediate(core.MM_GEOMETRY_FILE,
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

        optimized, opt_results = core.optimize_active_site(
            active_site,
            constrain_capping_hydrogens=self.constrain_capping_hydrogens,
            comm=self.comm,
            ostream=self.ostream,
            **self._qm_kwargs())

        self._save_intermediate(
            core.GEOMETRY_FILE,
            lambda path: optimized.write_xyz_file(str(path)))
        self._adopt_geometry(optimized)

        return opt_results

    def calculate_hessian(self):
        """
        Computes the nuclear Hessian of the active site.

        Restricted to the atom pairs the metal terms are fitted from unless
        calculate_partial_hessian is switched off, since Seminario reads one
        block per metal bond and two per metal angle and the fit throws the
        rest away.

        Collective: every rank runs it, and the drivers parallelize inside.

        :return:
            The Hessian.
        """

        self._require('calculate_hessian', Stage.ACTIVE_SITE)
        active_site = self._active_site

        atom_pairs, atoms = core.extract_pairs(
            active_site['connectivity_matrix'],
            active_site['metal_indices'],
            bond_count=2)
        n_atoms = active_site['molecule'].number_of_atoms()

        if self.calculate_partial_hessian:
            self.ostream.print_info(
                f'Hessian restricted to {len(atom_pairs)} atom pairs over '
                f'{len(atoms)} of {n_atoms} atoms.')
        else:
            self.ostream.print_info(
                f'Computing the full Hessian over all {n_atoms} atoms; the '
                f'metal terms read {len(atom_pairs)} of its blocks.')
        self.ostream.flush()

        hessian = core.compute_hessian(
            active_site,
            atom_pairs=atom_pairs if self.calculate_partial_hessian else None,
            comm=self.comm,
            ostream=self.ostream,
            **self._qm_kwargs())

        self._save_intermediate(core.HESSIAN_FILE,
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
            charges = core.compute_resp_charges(active_site,
                                                mute_scf=self.mute_scf,
                                                comm=self.comm,
                                                ostream=self.ostream)
        else:
            charges = self._on_master(
                lambda: core.d4_charges(active_site, ostream=self.ostream))

        self._save_intermediate(core.CHARGES_FILE,
                                lambda path: np.savetxt(path, charges))

        # a fit above this carries different charges
        self._enter(Stage.ACTIVE_SITE)
        self._partial_charges = charges

        return charges

    def build_forcefield(self,
                         hessian=None,
                         opt_geometry=None,
                         partial_charges=None):
        """
        Fits the metal terms and builds the force field of the active site.

        This is where the expensive work is triggered. For each of the
        geometry, the Hessian and the charges the precedence is the same: what
        is passed in here, else the file an earlier run left in output_folder,
        else computing it -- the constrained optimization when
        do_qm_optimization is on, the Hessian always, and the charges as RESP
        or D4. Each of the three is validated against the extracted active
        site before it is used.

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

        geometry = self._on_master(lambda: core._resolve_optimized_geometry(
            self._active_site,
            folder=self.output_folder,
            optimized_geometry=opt_geometry,
            ostream=self.ostream))

        if geometry is not None:
            self.ostream.print_info(
                'Using the geometry given to build_forcefield; skipping the '
                'constrained optimization.')
            self.ostream.flush()
            self._adopt_geometry(geometry)
        elif self.do_qm_optimization:
            self.optimize_geometry()

        hessian = self._on_master(
            lambda: core._resolve_hessian(self._active_site,
                                          folder=self.output_folder,
                                          hessian=hessian,
                                          ostream=self.ostream))

        if hessian is not None:
            self.ostream.print_info(
                'Using the Hessian given to build_forcefield; skipping the '
                'QM Hessian.')
            self.ostream.flush()
            self._hessian = hessian
        else:
            hessian = self.calculate_hessian()

        charges = self._on_master(lambda: core._resolve_partial_charges(
            self._active_site,
            folder=self.output_folder,
            partial_charges=partial_charges,
            ostream=self.ostream))

        if charges is not None:
            self.ostream.print_info(
                'Using the charges given to build_forcefield; skipping the '
                'charge calculation.')
            self.ostream.flush()
            self._partial_charges = charges
        else:
            charges = self.calculate_partial_charges()

        self._on_master(
            lambda: core._print_partial_charges(self._protonated_topology,
                                                self._active_site,
                                                charges,
                                                ostream=self.ostream))

        return self._fit_and_broadcast(hessian, charges)

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
        assert_msg_critical(
            self._forcefield is not None,
            'MetalSiteForceFieldBuilder.create_enzyme_system: there is no '
            'force field yet. Call build_forcefield first.')

        self._enzyme_system = self._on_master(self._create_enzyme_system)
        self._enter(Stage.ENZYME)

        return self._enzyme_system, self._protonated_topology

    def _create_enzyme_system(self):
        """
        The work of create_enzyme_system, on one rank.

        The force field is written out again alongside the system. Building
        the enzyme reads it rather than changing it, so the file it replaces
        holds the same parameters -- but writing it here means the last step
        to touch the folder leaves it complete, whichever way the run
        reached this point.
        """

        system, _ = core.create_enzyme_system(
            self._protonated_topology,
            self._active_site,
            self._forcefield,
            partial_charges=self._partial_charges,
            forcefield_files=self.protein_forcefield_files,
            ostream=self.ostream)

        self._save_intermediate(
            core.ENZYME_SYSTEM_FILE,
            lambda path: path.write_text(mm.XmlSerializer.serialize(system)))
        self._save_intermediate(
            core.FORCEFIELD_FILE,
            lambda path: core.save_forcefield(path, self._forcefield))

        return system

    # ------------------------------------------------------------------
    # the working folder
    # ------------------------------------------------------------------

    def _working_folder(self):
        """
        Returns the working folder, creating it if it does not exist.

        Called from _save_intermediate alone, which has already returned on
        every rank but the master, so there is no rank guard here: guarding
        in two places is how the two guards drift apart.

        :return:
            The folder as a Path.
        """

        folder = Path(self.output_folder)
        folder.mkdir(parents=True, exist_ok=True)

        return folder

    def _save_intermediate(self, name, writer):
        """
        Writes one intermediate as soon as the step that produced it finishes.

        Saving as the run goes rather than at the end means a run that dies
        part way through still leaves behind everything that did succeed, so a
        geometry optimization is not paid for twice because the Hessian after
        it failed.

        :param name:
            The file name, one of the file name constants of metalsitecore.
        :param writer:
            A callable taking the path to write.
        """

        if self.rank != mpi_master():
            return

        path = self._working_folder() / name
        writer(path)

        self.ostream.print_info(f'Wrote {name} to {self.output_folder}')
        self.ostream.flush()

    def _write_run_artifacts(self, forcefield, hessian, charges):
        """
        Writes everything the fit was made of and everything it produced.

        The geometry, the Hessian and the charges are written again here even
        when they were handed in or read back, so that the folder holds the
        run rather than only the parts of it that happened to be computed.
        The force field goes out as JSON, which is what a template reads.

        The charges written are the ones the force field ended up carrying,
        not the ones handed in: build_forcefield folds the capping hydrogens'
        charge into the rest of the site, and writing the raw fit beside a
        force field carrying the corrected one left two files disagreeing
        about what "the charges" are. Reading the file back and correcting it
        again is harmless, since redistribute_cap_charges is idempotent.

        :param forcefield:
            The fitted force field.
        :param hessian:
            The Hessian it was fitted from.
        :param charges:
            The charges it was given, before the cap correction.
        """

        corrected = forcefield.partial_charges
        if corrected is None:
            corrected = charges

        self._save_intermediate(
            core.GEOMETRY_FILE, lambda path: self._active_site['molecule'].
            write_xyz_file(str(path)))
        self._save_intermediate(core.HESSIAN_FILE,
                                lambda path: np.savetxt(path, hessian))
        self._save_intermediate(core.CHARGES_FILE,
                                lambda path: np.savetxt(path, corrected))
        self._save_intermediate(
            core.FORCEFIELD_FILE,
            lambda path: core.save_forcefield(path, forcefield))

    @staticmethod
    def _write_pdb(path, topology, positions):
        """
        Writes a topology and its positions as a PDB file.

        PDBFile.writeFile wants a quantity rather than bare numbers, and
        handing it a bare list writes nanometers as Angstrom and produces a
        structure shrunk tenfold.

        :param path:
            The file to write.
        :param topology:
            The OpenMM topology.
        :param positions:
            The positions in Angstrom.
        """

        with Path(path).open('w') as handle:
            mmapp.PDBFile.writeFile(topology,
                                    np.asarray(positions) * mmunit.angstrom,
                                    handle,
                                    keepIds=True)

    # ------------------------------------------------------------------
    # plumbing
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
            self._stage >= minimum,
            f'MetalSiteForceFieldBuilder.{method}: '
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
        Takes on what a structural step produced, on every rank.

        The step itself runs on the master and its result is broadcast, so
        this is the one place the builder's structural state is written and
        the two callers cannot come to write different sets of it.

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

    def _derive(self,
                topology,
                positions,
                coordinating_residues=None,
                request=None):
        """
        Works out the coordination of a structure from its geometry and the
        request.

        The one place the detection is called from, so that every derivation
        reads the same settings and replays the same decisions. Cheap enough
        -- it is a distance scan -- that keeping a copy instead would buy
        nothing and cost the chance of it going stale.

        :param topology:
            The topology to derive on.
        :param positions:
            Its positions in Angstrom.
        :param coordinating_residues:
            Residues to force as ligands, merged into the request's own list.
        :param request:
            The request to derive under. Defaults to the builder's own, and
            is given explicitly only while an edit is being worked out.

        :return:
            The binding modes.
        """

        return core.suggest_binding_modes(
            topology,
            positions,
            coordinating_residues=coordinating_residues,
            metal_elements=self.metal_elements,
            metal_formal_charges=self.metal_formal_charges,
            ostream=self.ostream,
            request=self._request if request is None else request,
            **self.detection_settings())

    def _on_master(self, work):
        """
        Runs work on the master rank and hands the result to every rank.

        Everything but the three QM steps is cheap enough that running it on
        one rank and broadcasting the result costs nothing, and it keeps the
        force field construction away from the collectives inside the drivers
        it uses. A failure on the master is broadcast as well and raised
        everywhere, so a rank cannot be left waiting for a result that is
        never coming.

        :param work:
            A callable taking no arguments.

        :return:
            What work returned, on every rank.
        """

        outcome = None

        if self.rank == mpi_master():
            try:
                outcome = ('value', work())
            except Exception as error:
                outcome = ('error', error)

        if self.nodes > 1:
            outcome = self.comm.bcast(outcome, root=mpi_master())

        kind, payload = outcome

        if kind == 'error':
            raise payload

        return payload

    def _adopt_geometry(self, molecule):
        """
        Puts a new geometry on the active site and detects the coordination
        again on it.

        A contact can leave the cutoff or a carboxylate can open up during a
        relaxation, and the connectivity the Hessian and the fit read comes
        out of that detection.

        :param molecule:
            The new geometry of the active site.
        """

        self._active_site = self._on_master(
            lambda: self._update_binding_modes(molecule))

        # a force field fitted before this describes the geometry it replaced
        self._enter(Stage.ACTIVE_SITE)

    def _update_binding_modes(self, molecule):
        """
        The work of _adopt_geometry, on one rank.

        The coordination before the move is derived rather than remembered,
        which is what the comparison needs and what nothing has to keep in
        step.
        """

        before = self._site_coordination(self._active_site['molecule'])

        self._active_site['molecule'] = molecule
        _, active_site, _ = core.update_binding_modes(
            self._protonated_topology,
            molecule,
            self._active_site,
            before,
            ostream=self.ostream,
            **self.detection_settings())

        return active_site

    def _site_coordination(self, geometry):
        """
        The coordination of the extracted site under one of its geometries.

        Restricted to the atoms the cluster holds, unlike the whole-structure
        derivation behind the binding_modes property: once a force field is
        keyed to the site, a protein atom drifting into range is not
        something the site can gain.

        :param geometry:
            The geometry to read, ordered like the active site.

        :return:
            The binding modes of the site under that geometry.
        """

        modes = self._derive(self._protonated_topology,
                             self._protonated_positions)
        ligands, notes = core.derive_site_coordination(
            self._protonated_topology,
            geometry,
            self._active_site,
            modes,
            ostream=self.ostream,
            **self.detection_settings())

        modes = dict(modes)
        modes['ligands'] = ligands
        modes['notes'] = notes

        return modes

    def _edit_coordination(self, edit):
        """
        Applies one coordination edit and brings the builder up to date.

        An edit is a decision, so it goes into the request, which is the
        same thing at either stage: there is only one numbering now, and the
        request survives every derivation. What differs between the stages
        is only what is done about it afterwards -- the site is extracted
        again before the fit, and the metal terms are fitted again after it.

        :param edit:
            A callable taking the request and the coordination it is being
            made against, and returning the edited request.
        """

        self._request = self._on_master(
            lambda: edit(self._request, self.binding_modes))

        self._reapply()

    def _reapply(self):
        """
        Brings everything downstream of an edit up to date.

        An edit is not a request to be remembered for later: the builder that
        made it describes the edited site immediately, at whichever stage it
        is in. Before the force field exists that means extracting the site
        again from the corrected binding modes; afterwards it means fitting
        the metal terms again on what is already there.
        """

        if self._stage < Stage.FITTED:
            self._rebuild_active_site()
        else:
            self._refit_forcefield()

    def _rebuild_active_site(self):
        """
        Protonates and truncates again from the binding modes as they stand.

        The same work build_active_site does without a path, and run the same
        way it was run then, so that an edit does not silently change whether
        the crude pass happens. The run header is not printed again: this is
        not a new run.
        """

        state = self._on_master(
            lambda: self._build_active_site(None, self._mm_opt, None))

        self._adopt(state)

        # the site they were computed for has just been replaced; a file in
        # the folder that still fits this one is picked back up by
        # build_forcefield, which validates it before it uses it
        self._enter(Stage.ACTIVE_SITE)

    def _refit_forcefield(self):
        """
        Fits the metal terms again after the coordination was edited.

        Nothing expensive runs. The geometry, the Hessian and the charges all
        still describe this cluster, and only which atoms the metals are
        bonded to has changed, so the force field is built again from them
        with the new connectivity. Rebuilding it rather than patching it is
        what keeps the angles, torsions and impropers that cross an edited
        bond right: the generator derives every one of them from the
        connectivity matrix.

        A term the Hessian holds nothing for is fitted to zero and reported by
        _check_force_constants, which is what a bond added by hand looks like
        when the Hessian was computed for a coordination without it. The bond
        is kept: the warning says to recompute the Hessian, not that the edit
        was refused.
        """

        self.ostream.print_info(
            'Fitting the metal terms again on the edited coordination. '
            'Nothing is recomputed; a term the Hessian does not cover is '
            'reported below.')
        self.ostream.flush()

        self._active_site = self._on_master(lambda: core.apply_metal_bonds(
            self._active_site,
            self._site_coordination(self._active_site['molecule'])))

        # entering FITTED drops the enzyme system built from the old terms
        self._fit_and_broadcast(self._hessian, self._partial_charges)

    def _fit_and_broadcast(self, hessian, charges):
        """
        Fits the metal terms on one rank and hands the force field to every
        rank.

        The fit itself is cheap and is kept off the collectives inside the
        generator, so it runs on the master alone. MMForceFieldGenerator does
        not pickle on its own -- it owns an output stream -- which is why the
        crossing goes through core.broadcast_forcefield rather than a plain
        bcast.

        Shared by the first fit and by every refit after a bond edit, so the
        two cannot come to write different artifacts or print different
        tables.

        :param hessian:
            The Hessian to fit from.
        :param charges:
            The charges the force field carries.

        :return:
            The force field, on every rank.
        """

        forcefield = None

        if self.rank == mpi_master():
            forcefield = core.build_forcefield(
                self._active_site,
                hessian=hessian,
                partial_charges=charges,
                comm=MPI.COMM_SELF,
                ostream=self.ostream,
                protected_bonds=self._protected_bonds(),
                **self.fit_settings())
            core._print_metal_parameters(self._active_site,
                                         forcefield,
                                         ostream=self.ostream)
            self._write_run_artifacts(forcefield, hessian, charges)

        self._forcefield = core.broadcast_forcefield(forcefield,
                                                     comm=self.comm,
                                                     ostream=self.ostream)

        # The fit folds the capping hydrogens' charge into the rest of the
        # site, and from here on that is what "the charges" means: the same
        # array reaches the force field, the enzyme system, the file and the
        # partial_charges property, so no two of them can answer differently.
        # Before a fit they are the raw result of the charge calculation,
        # which is what there is to have.
        if self._forcefield.partial_charges is not None:
            self._partial_charges = np.asarray(
                self._forcefield.partial_charges)

        self._enter(Stage.FITTED)

        return self._forcefield

    def _manual_equilibria(self):
        """
        The equilibrium distances add_metal_bond was told to pull bonds to.

        Read by the crude pass alone. A hand-added bond often exists because
        the sidechain is turned the wrong way round, and measuring its
        equilibrium on that geometry would hold it there; this is what lets
        the relaxation move it instead.

        :return:
            The bond keys mapped to distances in nanometers.
        """

        if self._protonated_topology is None:
            return {}

        return core.manual_bond_equilibria(self._active_site,
                                           self._protonated_topology,
                                           self.binding_modes)

    def _protected_bonds(self):
        """
        The metal bonds the weak bridge pruning must leave alone.

        A bond that was asked for by hand is a decision, and the two things
        the pruning reads -- a zero force constant and a long distance -- are
        exactly what such a bond looks like when the Hessian does not cover
        it. Without this the fit would take it straight back out.

        :return:
            The bond keys, as a set.
        """

        if self._protonated_topology is None:
            return set()

        return core.manual_bond_keys(self._active_site,
                                     self._protonated_topology,
                                     self.binding_modes)

    def detection_settings(self):
        """
        The settings the coordination detection reads, as keyword arguments.

        Public because MetalForceFieldManager measures a structure with the
        settings of the builder it carries, and splatting this into its own
        core call is what stops the two drifting apart on what a cutoff
        means. Note that metal_elements and metal_formal_charges are not in
        here -- they are passed alongside it, since the detection takes them
        but the re-detection on a new geometry does not.

        :return:
            The keyword arguments.
        """

        return {
            'metal_bond_cutoff': self.metal_bond_cutoff,
            'report_cutoff': self.report_cutoff,
            'bidentate_asymmetry': self.bidentate_asymmetry,
        }

    def _qm_kwargs(self):
        """
        The settings the QM drivers read.

        :return:
            The keyword arguments.
        """

        return {
            'scf_drv': self.scf_drv,
            'xcfun': self.xcfun,
            'basis_set_label': self.basis_set_label,
            'mute_scf': self.mute_scf,
        }

    def fit_settings(self):
        """
        The settings the force field construction reads, as keyword
        arguments, both for the fit and for the seeding a force field built
        without a Hessian falls back on.

        Public for the same reason as detection_settings: the manager builds
        force fields of its own, and reading them off the builder is what
        keeps a transferred force field and a fitted one made the same way.

        :return:
            The keyword arguments.
        """

        return {
            'average_metal_terms': self.average_metal_terms,
            'prune_weak_bridge_bonds': self.prune_weak_bridge_bonds,
            'weak_bridge_tolerance': self.weak_bridge_tolerance,
            'reparameterize_metal_angles': self.reparameterize_metal_angles,
            'default_metal_bond_force_constant':
                self.default_metal_bond_force_constant,
            'default_metal_angle_force_constant':
                self.default_metal_angle_force_constant,
            'metal_bond_equilibria': self.metal_bond_equilibria,
            'metal_angle_equilibria': self.metal_angle_equilibria,
        }

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

        param = core._param

        self.ostream.print_blank()
        self.ostream.print_header('Metal Site Force Field Builder')
        self.ostream.print_header(32 * '=')
        self.ostream.print_blank()

        self.ostream.print_header(param('structure', Path(structure).name))
        self.ostream.print_header(
            param('primary cutoff', f'{self.metal_bond_cutoff:.2f} A'))
        self.ostream.print_header(
            param('secondary cutoff', f'{self.report_cutoff:.2f} A'))
        self.ostream.print_header(param('basis set', self.basis_set_label))
        self.ostream.print_header(
            param('xc functional', self._functional_label()))
        self.ostream.print_header(
            param('SCF driver',
                  'given' if self.scf_drv is not None else 'default'))
        self.ostream.print_header(
            param('QM optimization', self.do_qm_optimization))
        self.ostream.print_header(
            param(
                'Hessian', 'computed, partial'
                if self.calculate_partial_hessian else 'computed, full'))
        self.ostream.print_header(
            param('partial charges', 'RESP' if self.do_resp else 'D4'))
        self.ostream.print_header(
            param(
                'constrained atoms', 'beta carbons + caps'
                if self.constrain_capping_hydrogens else 'beta carbons'))
        self.ostream.print_header(
            param(
                'weak bridge pruning', f'> {self.weak_bridge_tolerance:.2f} A'
                if self.prune_weak_bridge_bonds else 'off'))
        self.ostream.print_header(param('MPI ranks', self.nodes))
        self.ostream.print_header(param('output folder', self.output_folder))
        self.ostream.print_blank()
        self.ostream.flush()

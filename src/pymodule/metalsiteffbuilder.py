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

        # State. Everything a step produces is kept here and handed out
        # through a property, so that a caller carries a builder rather than
        # six dictionaries. The prepared topology is kept apart from the
        # protonated one: build_active_site rebuilds from the prepared one, and
        # protonating an already protonated topology is the bug that guards
        # against.
        self._structure = None
        self._mm_opt = True
        self._topology = None
        self._positions = None
        self._binding_modes = None
        self._protonated_topology = None
        self._protonated_positions = None
        self._protonated_binding_modes = None
        self._active_site = None
        self._forcefield = None
        self._hessian = None
        self._partial_charges = None
        self._enzyme_system = None

    # ------------------------------------------------------------------
    # results
    # ------------------------------------------------------------------

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

        The coordination as it now stands, which after the protonation is the
        one whose indices match the active site. The uncorrected set the edits
        are applied to is the one build_active_site rebuilds from.
        """

        if self._protonated_binding_modes is not None:
            return self._protonated_binding_modes

        return self._binding_modes

    @property
    def hessian(self):
        """
        The Hessian the metal terms were fitted from.
        """

        return self._hessian

    @property
    def partial_charges(self):
        """
        The charges of the active site, as fitted or as computed with D4.
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

        self._topology = state['topology']
        self._positions = state['positions']
        self._binding_modes = state['binding_modes']
        self._protonated_topology = state['protonated_topology']
        self._protonated_positions = state['protonated_positions']
        self._protonated_binding_modes = state['protonated_binding_modes']
        self._active_site = state['active_site']
        if cif_path is not None:
            self._structure = cif_path

        # everything downstream belonged to the site that was just replaced
        self._forcefield = None
        self._hessian = None
        self._partial_charges = None
        self._enzyme_system = None

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
            binding_modes = core.suggest_binding_modes(
                topology,
                positions,
                coordinating_residues=coordinating_residues,
                metal_elements=self.metal_elements,
                metal_formal_charges=self.metal_formal_charges,
                ostream=self.ostream,
                **self._detection_kwargs())
            core._print_binding_modes(binding_modes, ostream=self.ostream)
        else:
            assert_msg_critical(
                self._topology is not None,
                'MetalSiteForceFieldBuilder.build_active_site: there is no '
                'structure to rebuild from. Call it with a cif_path first.')
            topology = self._topology
            positions = self._positions
            binding_modes = self._binding_modes

        protonated_topology, protonated_positions, protonated_modes = (
            core.protonate(topology,
                           positions,
                           binding_modes,
                           protonation_overrides=self.protonation_overrides))

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
            'binding_modes.json',
            lambda path: core.save_binding_modes(path, binding_modes))
        self._save_intermediate(
            'protonated.pdb', lambda path: self._write_pdb(
                path, protonated_topology, protonated_positions))

        if mm_opt:
            forcefield = core.build_forcefield(active_site,
                                               comm=MPI.COMM_SELF,
                                               ostream=self.ostream,
                                               **self._fit_kwargs())
            relaxed = core.mm_optimize_active_site(
                active_site,
                forcefield,
                constrain_metals=self.mm_constrain_metals,
                constrain_capping_hydrogens=self.constrain_capping_hydrogens,
                max_iterations=self.mm_max_iterations,
                bond_change_warning=self.mm_bond_change_warning,
                ostream=self.ostream)
            active_site['molecule'] = relaxed
            self._save_intermediate(
                core.MM_GEOMETRY_FILE,
                lambda path: relaxed.write_xyz_file(str(path)))

        return {
            'topology': topology,
            'positions': positions,
            'binding_modes': binding_modes,
            'protonated_topology': protonated_topology,
            'protonated_positions': protonated_positions,
            'protonated_binding_modes': protonated_modes,
            'active_site': active_site,
        }

    def add_metal_bond(self, resid, metal, atom=None, chain=None):
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
        """

        self._require_binding_modes('add_metal_bond')

        def edit(modes, topology, positions):
            return core.add_metal_bond(modes,
                                       topology,
                                       positions,
                                       resid,
                                       self._metal_index(metal, modes),
                                       atom=atom,
                                       chain=chain,
                                       ostream=self.ostream,
                                       **self._detection_kwargs())

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

        self._require_binding_modes('remove_metal_bond')

        def edit(modes, topology, positions):
            return core.remove_metal_bond(
                modes,
                resid,
                metal=None if metal is None else self._metal_index(
                    metal, modes),
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

        self._require_binding_modes('update_protonation_state')
        self._require_editable('update_protonation_state')

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

        self._require_binding_modes('include_residue')
        self._require_editable('include_residue')

        def work():
            residue = core._resolve_residue(self._topology, resid, chain)
            core.check_truncatable(residue)

            modes = self._binding_modes
            label = f'{residue.name}{residue.id}'

            already = residue.index in core.active_site_residues(modes)
            excluded = residue.index in modes.get('excluded_residues', [])

            if already and not excluded:
                self.ostream.print_warning(
                    f'{label} is already part of the active site; nothing to '
                    'include')
                self.ostream.flush()
                return None

            modes = deepcopy(modes)
            modes['excluded_residues'] = [
                index for index in modes.get('excluded_residues', [])
                if index != residue.index
            ]
            extra = set(modes.get('extra_residues', []))
            extra.add(residue.index)
            modes['extra_residues'] = sorted(extra)

            core._merge_notes(modes, [
                f'{label} is part of the active site because include_residue '
                'asked for it, whether or not it coordinates a metal'
            ])

            self.ostream.print_info(f'Added {label} to the active site.')
            self.ostream.flush()

            return modes

        modes = self._on_master(work)

        if modes is None:
            return

        self._binding_modes = modes
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

        self._require_binding_modes('remove_residue')
        self._require_editable('remove_residue')

        def work():
            residue = core._resolve_residue(self._topology, resid, chain)

            modes = self._binding_modes
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
                modes = core.remove_metal_bond(
                    modes,
                    label,
                    metal=bound['metals'][0],
                    bidentate_asymmetry=self.bidentate_asymmetry,
                    ostream=self.ostream)

            modes = deepcopy(modes)
            modes['extra_residues'] = [
                index for index in modes.get('extra_residues', [])
                if index != residue.index
            ]
            modes['coordinating_residues'] = [
                index for index in modes.get('coordinating_residues', [])
                if index != residue.index
            ]
            excluded = set(modes.get('excluded_residues', []))
            excluded.add(residue.index)
            modes['excluded_residues'] = sorted(excluded)

            core._merge_notes(modes, [
                f'{label} is not part of the active site because '
                'remove_residue asked for it'
            ])

            self.ostream.print_info(f'Removed {label} from the active site.')
            self.ostream.flush()

            return modes

        self._binding_modes = self._on_master(work)
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

        self._require_active_site('mm_optimize_active_site')

        molecule = self._on_master(self._mm_optimize_active_site)
        self._active_site['molecule'] = molecule

        return molecule

    def _mm_optimize_active_site(self):
        """
        The work of mm_optimize_active_site, on one rank.
        """

        forcefield = core.build_forcefield(self._active_site,
                                           comm=MPI.COMM_SELF,
                                           ostream=self.ostream,
                                           **self._fit_kwargs())
        relaxed = core.mm_optimize_active_site(
            self._active_site,
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

        active_site = self._require_active_site('optimize_geometry')

        optimized, opt_results = core.optimize_active_site(
            active_site,
            constrain_capping_hydrogens=self.constrain_capping_hydrogens,
            comm=self.comm,
            ostream=self.ostream,
            **self._qm_kwargs())

        self._save_intermediate(core.GEOMETRY_FILE,
                                lambda path: optimized.write_xyz_file(
                                    str(path)))
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

        active_site = self._require_active_site('calculate_hessian')

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

        active_site = self._require_active_site('calculate_partial_charges')

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

        self._require_active_site('build_forcefield')

        geometry = self._on_master(
            lambda: core._resolve_optimized_geometry(
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

        hessian = self._on_master(lambda: core._resolve_hessian(
            self._active_site,
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

        self._on_master(lambda: core._print_partial_charges(
            self._protonated_topology,
            self._active_site,
            charges,
            ostream=self.ostream))

        forcefield = None
        if self.rank == mpi_master():
            forcefield = core.build_forcefield(
                self._active_site,
                hessian=hessian,
                partial_charges=charges,
                comm=MPI.COMM_SELF,
                ostream=self.ostream,
                protected_bonds=self._protected_bonds(),
                **self._fit_kwargs())
            core._print_metal_parameters(self._active_site,
                                         forcefield,
                                         ostream=self.ostream)
            self._write_run_artifacts(forcefield, hessian, charges)

        self._forcefield = core.broadcast_forcefield(forcefield,
                                                     comm=self.comm,
                                                     ostream=self.ostream)

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

        self._require_active_site('create_enzyme_system')
        assert_msg_critical(
            self._forcefield is not None,
            'MetalSiteForceFieldBuilder.create_enzyme_system: there is no '
            'force field yet. Call build_forcefield first.')

        self._enzyme_system = self._on_master(self._create_enzyme_system)

        return self._enzyme_system, self._protonated_topology

    def _create_enzyme_system(self):
        """
        The work of create_enzyme_system, on one rank.
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

        return system

    # ------------------------------------------------------------------
    # the working folder
    # ------------------------------------------------------------------

    def _working_folder(self):
        """
        Returns the working folder, creating it if it does not exist.

        :return:
            The folder as a Path.
        """

        folder = Path(self.output_folder)

        if self.rank == mpi_master():
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
        The force field goes out as JSON for what comes after the run, and as
        OpenMM files for anything that wants to simulate the site on its own.

        :param forcefield:
            The fitted force field.
        :param hessian:
            The Hessian it was fitted from.
        :param charges:
            The charges it carries.
        """

        self._save_intermediate(
            core.GEOMETRY_FILE, lambda path: self._active_site['molecule'].
            write_xyz_file(str(path)))
        self._save_intermediate(core.HESSIAN_FILE,
                                lambda path: np.savetxt(path, hessian))
        self._save_intermediate(core.CHARGES_FILE,
                                lambda path: np.savetxt(path, charges))
        self._save_intermediate(
            core.FORCEFIELD_FILE,
            lambda path: core.save_forcefield(path, forcefield))
        self._save_intermediate(
            'forcefield.xml', lambda path: forcefield.write_openmm_files(
                str(path.with_suffix(''))))

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

        state = self._on_master(lambda: self._update_binding_modes(molecule))

        self._protonated_binding_modes = state['binding_modes']
        self._active_site = state['active_site']

    def _update_binding_modes(self, molecule):
        """
        The work of _adopt_geometry, on one rank.
        """

        self._active_site['molecule'] = molecule
        binding_modes, active_site, _ = core.update_binding_modes(
            self._protonated_topology,
            molecule,
            self._active_site,
            self._protonated_binding_modes,
            bidentate_asymmetry=self.bidentate_asymmetry,
            ostream=self.ostream)

        return {'binding_modes': binding_modes, 'active_site': active_site}

    def _edit_coordination(self, edit):
        """
        Applies one coordination edit and brings the builder up to date.

        Which binding modes the edit is applied to is the whole difference
        between the two stages. Before the force field exists the edit belongs
        on the modes the site is rebuilt from, which are the ones from before
        the hydrogens were placed. Afterwards there is a geometry, a Hessian
        and a set of charges that all belong to this cluster, and rebuilding
        would throw them away, so the edit is applied to the protonated modes
        the active site indexes into and the metal terms are fitted again.

        :param edit:
            A callable taking the binding modes, the topology and the
            positions they index into, and returning the edited modes.
        """

        if self._forcefield is None:
            self._binding_modes = self._on_master(lambda: edit(
                self._binding_modes, self._topology, self._positions))
        else:
            self._protonated_binding_modes = self._on_master(
                lambda: edit(self._protonated_binding_modes, self.
                             _protonated_topology, self._protonated_positions))

        self._reapply()

    def _metal_index(self, metal, binding_modes):
        """
        Translates a metal atom index into the numbering of one set of
        binding modes.

        A metal is named by the number the active site labels it with, which
        is an index into the protonated topology, since that is the only one
        anything draws or prints for the user to read. The modes an edit is
        applied to before the force field exists are the ones from before the
        hydrogens were placed, where the same atom has a different index, so
        the two are matched up by position: protonate remaps the entries in
        place and keeps their order.

        :param metal:
            The atom index of the metal, as the active site labels it.
        :param binding_modes:
            The binding modes whose numbering is wanted.

        :return:
            The atom index of that metal in those binding modes.
        """

        reference = self._protonated_binding_modes or self._binding_modes
        wanted = int(metal)

        position = None
        for index, entry in enumerate(reference['metals']):
            if entry['index'] == wanted:
                position = index
                break

        known = ', '.join(
            f'{entry["element"]} {entry["index"]}'
            for entry in reference['metals'])

        assert_msg_critical(
            position is not None, 'MetalSiteForceFieldBuilder: there is no '
            f'metal center at atom index {wanted}. The active site holds '
            f'{known}, which is what its labels show')

        return binding_modes['metals'][position]['index']

    def _reapply(self):
        """
        Brings everything downstream of an edit up to date.

        An edit is not a request to be remembered for later: the builder that
        made it describes the edited site immediately, at whichever stage it
        is in. Before the force field exists that means extracting the site
        again from the corrected binding modes; afterwards it means fitting
        the metal terms again on what is already there.
        """

        if self._forcefield is None:
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

        self._topology = state['topology']
        self._positions = state['positions']
        self._binding_modes = state['binding_modes']
        self._protonated_topology = state['protonated_topology']
        self._protonated_positions = state['protonated_positions']
        self._protonated_binding_modes = state['protonated_binding_modes']
        self._active_site = state['active_site']

        # the site they were computed for has just been replaced; a file in
        # the folder that still fits this one is picked back up by
        # build_forcefield, which validates it before it uses it
        self._hessian = None
        self._partial_charges = None
        self._enzyme_system = None

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
            self._active_site, self._protonated_binding_modes))

        forcefield = None
        if self.rank == mpi_master():
            forcefield = core.build_forcefield(
                self._active_site,
                hessian=self._hessian,
                partial_charges=self._partial_charges,
                comm=MPI.COMM_SELF,
                ostream=self.ostream,
                protected_bonds=self._protected_bonds(),
                **self._fit_kwargs())
            core._print_metal_parameters(self._active_site,
                                         forcefield,
                                         ostream=self.ostream)
            self._write_run_artifacts(forcefield, self._hessian,
                                      self._partial_charges)

        self._forcefield = core.broadcast_forcefield(forcefield,
                                                     comm=self.comm,
                                                     ostream=self.ostream)

        # it was built from the metal terms that have just been refitted
        self._enzyme_system = None

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

        if self._protonated_binding_modes is None:
            return set()

        return core.manual_bond_keys(self._active_site,
                                     self._protonated_topology,
                                     self._protonated_binding_modes)

    def _require_editable(self, method):
        """
        Checks that the force field has not been built yet.

        The protonation and the residues of the active site are what the
        charges and the fitted parameters describe, so changing either after
        the fit would leave a force field describing a cluster that no longer
        exists. Coordination is the exception: only which atoms the metals
        bond to changes, and that can be fitted again from what is already
        there.

        :param method:
            The name of the method that needs it.
        """

        assert_msg_critical(
            self._forcefield is None,
            f'MetalSiteForceFieldBuilder.{method}: the force field has '
            'already been built, and its charges and parameters describe the '
            'protonation and the residues it was fitted to. Only '
            'add_metal_bond and remove_metal_bond can be used after '
            'build_forcefield. Call build_active_site() again to start over '
            'from the binding modes, which drops the fit.')

    def _require_active_site(self, method):
        """
        Returns the active site, or says which call is missing.

        :param method:
            The name of the method that needs it.

        :return:
            The active site.
        """

        assert_msg_critical(
            self._active_site is not None,
            f'MetalSiteForceFieldBuilder.{method}: there is no active site '
            'yet. Call build_active_site first.')

        return self._active_site

    def _require_binding_modes(self, method):
        """
        Checks that there are binding modes to edit.

        :param method:
            The name of the method that needs them.
        """

        assert_msg_critical(
            self._binding_modes is not None,
            f'MetalSiteForceFieldBuilder.{method}: there is no coordination '
            'to correct yet. Call build_active_site first.')

    def _detection_kwargs(self):
        """
        The settings the coordination detection reads.

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

    def _fit_kwargs(self):
        """
        The settings the force field construction reads, both for the fit and
        for the seeding a force field built without a Hessian falls back on.

        :return:
            The keyword arguments.
        """

        return {
            'average_metal_terms':
                self.average_metal_terms,
            'prune_weak_bridge_bonds':
                self.prune_weak_bridge_bonds,
            'weak_bridge_tolerance':
                self.weak_bridge_tolerance,
            'reparameterize_metal_angles':
                self.reparameterize_metal_angles,
            'default_metal_bond_force_constant':
                self.default_metal_bond_force_constant,
            'default_metal_angle_force_constant':
                self.default_metal_angle_force_constant,
            'metal_bond_equilibria':
                self.metal_bond_equilibria,
            'metal_angle_equilibria':
                self.metal_angle_equilibria,
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

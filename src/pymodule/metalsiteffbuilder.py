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
import tempfile
import json
import time
import sys

from .veloxchemlib import mpi_master
from .molecule import Molecule
from .molecularbasis import MolecularBasis
from .outputstream import OutputStream
from .mmforcefieldgenerator import MMForceFieldGenerator
from .respchargesdriver import RespChargesDriver
from .scfrestdriver import ScfRestrictedDriver
from .scfunrestdriver import ScfUnrestrictedDriver
from .scfgradientdriver import ScfGradientDriver
from .scfhessiandriver import ScfHessianDriver
from .optimizationdriver import OptimizationDriver
from .errorhandler import assert_msg_critical

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


class MetalSiteForceFieldBuilder:
    """
    Builds a bonded force field for the zinc center of a metallo-enzyme.

    Identifies the coordination sphere, fixes the protonation of the
    coordinating residues, truncates a QM active site at the CA-CB bonds, and fits
    the metal-ligand bond and angle parameters to a QM Hessian with the
    Seminario method. The fitted parameters can then be injected into a
    protein force field for the whole enzyme.

    Only settings are kept on the object. Everything a step produces - the
    binding modes, the truncated active site, the connectivity matrix and the
    force field - is returned by the method that builds it and passed to the
    methods that need it, so a call site shows exactly what each step consumes.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - metal_bond_cutoff: The distance in Angstrom within which a donor atom is
          taken to be bonded to a metal center.
        - report_cutoff: The distance in Angstrom out to which contacts are
          reported for review.
        - metal_elements: The elements treated as metal centers.
        - metal_formal_charges: The formal charges assumed for the metal ions.
        - cap_bond_length: The C-H distance in Angstrom of the capping
          hydrogen that replaces CA.
        - protonation_overrides: The residue id to variant mapping that
          overrides the automatic protonation rules.
        - scf_drv: The SCF driver, used as given. Created from xcfun if not
          provided. Set convergence thresholds and iteration limits on it.
        - xcfun: The exchange-correlation functional used when scf_drv is not
          provided. None runs Hartree-Fock.
        - basis_set_label: The basis set label.
        - mute_scf: The flag for muting the output of the QM drivers.
        - optimized_geometry: A molecule, or the path to an xyz file, to use
          instead of running the constrained optimization. Falls back to
          opt_active_site.xyz in the working folder.
        - hessian: A matrix, or the path to a text file readable by
          numpy.loadtxt, to use instead of computing a Hessian. Falls back to
          hessian.txt in the working folder.
        - partial_charges: Charges, or the path to a text file, to use instead
          of computing RESP charges. Falls back to partial_charges.txt in the
          working folder.
        - do_qm_optimization: The flag for optimizing the active site before the
          Hessian is computed.
        - do_mm_optimization: The flag for relaxing the active site on a crude
          force field of its own before any QM is run.
        - mm_max_iterations: The iteration limit of that relaxation. Zero runs
          until convergence.
        - default_metal_bond_force_constant: The force constant given to every
          metal bond of the crude pass, in kJ/mol/nm^2.
        - default_metal_angle_force_constant: The force constant given to every
          metal angle of the crude pass, in kJ/mol/rad^2.
        - metal_bond_equilibria: Equilibrium metal-ligand distances in nm,
          keyed by element pair in either order, replacing the values measured
          on the input geometry. LITERATURE_METAL_BONDS is such a table.
        - metal_angle_equilibria: Equilibrium metal angles in degrees, keyed by
          element triple in either order, replacing the measured values.
        - reparameterize_metal_angles: The flag for touching the metal angles
          at all. False leaves every one of them at the value the generator
          guessed, in the crude pass and in the Hessian fit alike.
        - mm_bond_change_warning: How far a metal-ligand bond may move in the
          crude pass, in Angstrom, before it is reported.
        - mm_constrain_metals: The flag for holding the metal centers fixed in
          the crude pass as well, on top of the beta carbons.
        - do_resp: The flag for computing RESP charges.
        - constrain_capping_hydrogens: The flag for constraining the capping
          hydrogens in addition to the beta carbons.
        - average_metal_terms: The flag for averaging the fitted metal terms
          over equivalent atoms.
        - folder: The working folder. Every expensive step writes its result
          here as soon as it has it, and reads it back on a later run. Named
          after the creation time by default, so runs do not collide.
        - build_enzyme_system: The flag for building a force field system
          for the whole enzyme at the end of a run, which is what applies the
          fitted charges to the protein.
        - save_output: The flag for writing the results at the end of a run.
        - comm: The MPI communicator.
        - rank: The rank of the MPI process.
        - nodes: The number of MPI processes.
        - ostream: The output stream.
    """

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

    # Residues whose sidechain can bridge two metal centers. A carboxylate has
    # two donor oxygens and a thiolate sulfur has several lone pairs, so both
    # can reach two metals at once. An imidazole nitrogen has a single lone
    # pair in the ring plane, so a histidine bound to two metals is a geometric
    # artifact rather than a bridge, however short the second contact looks.
    BRIDGING_RESIDUES = ('ASP', 'ASH', 'GLU', 'GLH', 'CYS', 'CYM', 'CYX')

    # Atoms dropped when truncating a sidechain at the CA-CB bond. CA itself
    # is replaced by a capping hydrogen rather than dropped.
    BACKBONE_ATOM_NAMES = ('N', 'C', 'O', 'OXT', 'H', 'H2', 'H3', 'HA', 'HA2',
                           'HA3', 'HXT')

    # Net charge of each protonation variant, for the active site charge
    # bookkeeping.
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

    # Literature equilibrium metal-ligand distances in nm. The crude pre-QM
    # pass measures its equilibrium values on the input geometry instead;
    # assign this to metal_bond_equilibria to impose these in their place.
    LITERATURE_METAL_BONDS = {
        ('Zn', 'N'): 0.205,
        ('Zn', 'O'): 0.200,
        ('Zn', 'S'): 0.230,
    }

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
        self.metal_elements = tuple(self.METAL_ELEMENTS)
        self.metal_formal_charges = dict(self.METAL_FORMAL_CHARGES)

        # a carboxylate this much more lopsided than BIDENTATE_ASYMMETRY is
        # not chelating, whatever the second oxygen looks like on paper
        self.bidentate_asymmetry = self.BIDENTATE_ASYMMETRY

        # truncation and protonation
        self.cap_bond_length = 1.09
        self.protonation_overrides = None

        # QM settings. Anything beyond the functional and the basis set is
        # set on an SCF driver assigned to scf_drv, which is used as given.
        self.scf_drv = None
        self.xcfun = 'PBE0'
        self.basis_set_label = 'def2-svp'
        self.mute_scf = True

        # Precomputed input. Each may be given as an object or as a path.
        # Whatever is provided is validated against the extracted active site and
        # the step that would have produced it is skipped.
        self.optimized_geometry = None
        self.hessian = None
        self.partial_charges = None

        # workflow
        self.do_qm_optimization = True
        self.do_resp = True
        # the beta carbon is where the backbone actually holds the sidechain;
        # the capping hydrogen only stands in for the alpha carbon
        self.constrain_capping_hydrogens = False
        # two topologically equivalent ligands of a strained site sit at
        # genuinely different distances, and that asymmetry is the signal
        self.average_metal_terms = False

        # Crude MM pass. The site is relaxed on a force field of its own
        # before any QM is run, so the optimization does not spend its first
        # cycles undoing whatever the structure file happened to hold. Nothing
        # at that point says anything about the stiffness of a metal term, so
        # the force constants are flat defaults and only the equilibrium
        # values carry information: they are measured on the input geometry
        # unless one of the tables below overrides them.
        self.do_mm_optimization = True
        self.mm_max_iterations = 0
        self.default_metal_bond_force_constant = 100000.0
        self.default_metal_angle_force_constant = 200.0
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

        # Whether compute() goes on to build a force field system for the
        # whole enzyme. That is where the fitted charges reach the protein and
        # where redistribute_backbone_charges restores the region charge, so
        # it also pulls in prepare_protein, which the protein force field needs
        # in order to match templates.
        self.build_enzyme_system = True

        # Output. The folder is named after the time the builder was created,
        # so two runs never overwrite one another; assign folder to put the
        # files somewhere of your own choosing, or to point at an earlier run
        # and reuse whatever it managed to finish.
        self.folder = f'metal_site_{int(time.time())}'
        self.save_output = True

    def compute(self, structure):
        """
        Runs the complete pipeline on a structure.

        Detects the coordination sphere, protonates the coordinating
        residues, truncates the QM active site, optionally optimizes it with the
        beta carbons frozen, computes a pair-restricted Hessian and RESP
        charges, and fits the metal-ligand terms.

        A geometry, Hessian or set of partial charges supplied through
        optimized_geometry, hessian or partial_charges is validated against
        the extracted active site and used in place of the step that would
        have produced it. The force field is always fitted afresh from
        whatever geometry and Hessian the run ends up with.

        :param structure:
            The path to a PDB or mmCIF file.

        :return:
            The results dictionary.
        """

        assert_msg_critical('openmm' in sys.modules,
                            'MetalSiteForceFieldBuilder: openmm is required')

        self._print_header(structure)

        topology, positions = self.load_structure(structure)

        if self.build_enzyme_system:
            # every index downstream refers to this topology, so the repair
            # has to happen before anything is extracted from it
            topology, positions = self.prepare_protein(topology, positions)

        binding_modes = self.suggest_binding_modes(topology, positions)
        self._print_binding_modes(binding_modes)

        topology, positions, binding_modes = self.protonate(
            topology, positions, binding_modes)
        active_site = self.extract_active_site(topology, positions,
                                               binding_modes)
        connectivity_matrix = self.build_connectivity(topology, active_site,
                                                      binding_modes)
        self._print_active_site(active_site, binding_modes, connectivity_matrix)

        # anything supplied through the settings is validated against the
        # active site and replaces the step that would have produced it
        geometry = self._resolve_optimized_geometry(active_site)
        hessian = self._resolve_hessian(active_site)
        charges = self._resolve_partial_charges(active_site)

        # a supplied geometry is the geometry to work with, so there is
        # nothing for the crude pass to improve on
        if geometry is None and self.do_mm_optimization:
            active_site['molecule'] = self.mm_optimize_active_site(
                active_site, connectivity_matrix)

        if geometry is not None:
            self.ostream.print_info(
                'Using the supplied optimized geometry; skipping the '
                'constrained optimization.')
            self.ostream.flush()
        elif self.do_qm_optimization:
            geometry = self.optimize_active_site(active_site)

        if geometry is not None:
            active_site['molecule'] = geometry
            binding_modes, connectivity_matrix, _ = self.update_binding_modes(
                topology, geometry, active_site, binding_modes,
                connectivity_matrix)

        molecule = active_site['molecule']

        atom_pairs, atoms = self.extract_pairs(connectivity_matrix,
                                               active_site['metal_indices'],
                                               bond_count=2)

        if hessian is not None:
            self.ostream.print_info(
                'Using the supplied Hessian; skipping the QM Hessian.')
            self.ostream.flush()
        else:
            self.ostream.print_info(
                f'Hessian restricted to {len(atom_pairs)} atom pairs over '
                f'{len(atoms)} of {molecule.number_of_atoms()} atoms.')
            self.ostream.flush()
            hessian = self.compute_hessian(active_site, atom_pairs)

        if charges is not None:
            self.ostream.print_info(
                'Using the supplied partial charges; skipping the RESP '
                'calculation.')
            self.ostream.flush()
        elif self.do_resp:
            charges = self.compute_resp_charges(active_site)

        if charges is not None:
            self._print_partial_charges(topology, active_site, charges)

        forcefield = self.build_forcefield(active_site, connectivity_matrix,
                                           hessian, charges)
        # written here rather than inside build_forcefield, which the
        # crude MM pass calls as well: only this one is the fit
        self._save_intermediate(
            self.FORCEFIELD_FILE,
            lambda path: self.save_forcefield(path, forcefield))

        self._print_metal_parameters(active_site, forcefield)

        enzyme_system = None
        if self.build_enzyme_system:
            enzyme_system, added = self.create_enzyme_system(
                topology, active_site, forcefield, charges)

        results = {
            'topology': topology,
            'positions': positions,
            'enzyme_system': enzyme_system,
            'binding_modes': binding_modes,
            'active_site': active_site,
            'connectivity_matrix': connectivity_matrix,
            'atom_pairs': atom_pairs,
            'hessian': hessian,
            'partial_charges': charges,
            'forcefield': forcefield,
        }

        if self.save_output:
            self.save_results(results)

        return results

    def working_folder(self):
        """
        Returns the working folder, creating it if it does not exist.

        :return:
            The folder as a Path.
        """

        folder = Path(self.folder)

        if self.rank == mpi_master():
            folder.mkdir(parents=True, exist_ok=True)

        return folder

    def _save_intermediate(self, name, writer):
        """
        Writes one intermediate as soon as the step that produced it finishes.

        Saving inside the steps rather than only at the end of compute() means
        a run that dies part way through still leaves behind everything that
        did succeed, so a geometry optimization is not paid for twice because
        the Hessian after it failed.

        :param name:
            The file name, one of the class-level file name attributes.
        :param writer:
            A callable taking the path to write.
        """

        if not self.save_output or self.rank != mpi_master():
            return

        path = self.working_folder() / name
        writer(path)

        self.ostream.print_info(f'Wrote {name} to {self.folder}')
        self.ostream.flush()

    def _folder_file(self, name):
        """
        Returns the path of an intermediate in the working folder, or None
        when it is not there.

        :param name:
            The file name, one of the class-level file name attributes.

        :return:
            The path, or None.
        """

        path = Path(self.folder) / name

        return path if path.is_file() else None

    def save_results(self, results, folder=None):
        """
        Writes everything a run produced into one folder.

        Entries that are absent are skipped rather than treated as an error, so
        a structural run that stopped before any QM was paid for still leaves
        its binding modes behind to be reviewed. Every file a later run can
        read back is named by the class-level file name attributes, so that
        saving and reloading cannot drift apart.

        :param results:
            The results dictionary of compute(), or any subset of it.
        :param folder:
            The folder to write into. Defaults to the working folder.

        :return:
            The folder as a Path.
        """

        if folder is None:
            folder = self.working_folder()
        else:
            folder = Path(folder)
            if self.rank == mpi_master():
                folder.mkdir(parents=True, exist_ok=True)

        if self.rank != mpi_master():
            return folder

        written = []

        binding_modes = results.get('binding_modes')
        if binding_modes is not None:
            self.save_binding_modes(folder / 'binding_modes.json',
                                    binding_modes)
            written.append('binding_modes.json')

        topology = results.get('topology')
        positions = results.get('positions')
        if topology is not None and positions is not None:
            with (folder / 'protonated.pdb').open('w') as fh:
                mmapp.PDBFile.writeFile(topology,
                                        np.asarray(positions) * mmunit.angstrom,
                                        fh,
                                        keepIds=True)
            written.append('protonated.pdb')

        active_site = results.get('active_site')
        if active_site is not None:
            active_site['molecule'].write_xyz_file(
                str(folder / self.GEOMETRY_FILE))
            written.append(self.GEOMETRY_FILE)

        hessian = results.get('hessian')
        if hessian is not None:
            np.savetxt(folder / self.HESSIAN_FILE, hessian)
            written.append(self.HESSIAN_FILE)

        partial_charges = results.get('partial_charges')
        if partial_charges is not None:
            np.savetxt(folder / self.CHARGES_FILE, partial_charges)
            written.append(self.CHARGES_FILE)

        forcefield = results.get('forcefield')
        if forcefield is not None:
            forcefield.write_openmm_files(str(folder / 'forcefield'))
            written.extend(['forcefield.xml', 'forcefield.pdb'])
            self.save_forcefield(folder / self.FORCEFIELD_FILE, forcefield)
            written.append(self.FORCEFIELD_FILE)

        enzyme_system = results.get('enzyme_system')
        if enzyme_system is not None:
            (folder / self.ENZYME_SYSTEM_FILE).write_text(
                mm.XmlSerializer.serialize(enzyme_system))
            written.append(self.ENZYME_SYSTEM_FILE)

        if written:
            self.ostream.print_info(f'Wrote {", ".join(written)} to {folder}')
            self.ostream.flush()

        return folder

    def load_structure(self, structure):
        """
        Reads a PDB or mmCIF structure.

        :param structure:
            The path to a .pdb, .cif or .pdbx file.

        :return:
            The tuple of the OpenMM topology and the positions as an (N, 3)
            numpy array in Angstrom.
        """

        assert_msg_critical(
            'openmm' in sys.modules,
            'MetalSiteForceFieldBuilder.load_structure: openmm is required')

        path = Path(structure)

        assert_msg_critical(
            path.is_file(),
            f'MetalSiteForceFieldBuilder.load_structure: {path} not found')

        if path.suffix.lower() in ('.cif', '.pdbx', '.mmcif'):
            pdb = mmapp.PDBxFile(str(path))
        else:
            pdb = mmapp.PDBFile(str(path))

        positions = np.array(pdb.positions.value_in_unit(mmunit.angstrom))

        return pdb.topology, positions

    def prepare_protein(self, topology, positions):
        """
        Adds missing heavy atoms so that a protein force field can match
        templates.  Missing residues are deliberately not built.
        Only needed before create_enzyme_system.

        :param topology:
            The OpenMM topology.
        :param positions:
            The positions as an (N, 3) numpy array in Angstrom.

        :return:
            The tuple of the repaired topology and positions in Angstrom.
        """

        assert_msg_critical(
            'pdbfixer' in sys.modules or 'PDBFixer' in globals(),
            'MetalSiteForceFieldBuilder.prepare_protein: pdbfixer is required')

        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / 'input.pdb'
            with path.open('w') as fh:
                mmapp.PDBFile.writeFile(topology,
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
            element for element in found
            if element not in self.SUPPORTED_METAL_ELEMENTS
        ]

        assert_msg_critical(
            not unsupported,
            f'MetalSiteForceFieldBuilder.{method}: found {unsupported}, but '
            f'only {list(self.SUPPORTED_METAL_ELEMENTS)} is supported. The '
            'literature distances, formal charges and coordination rules have '
            'only been validated for zinc.')

    def suggest_binding_modes(self, topology, positions, include_residue=None):
        """
        Proposes the coordination topology of the metal centers from geometry.

        Unrelaxed designed structures are unreliable enough that the intended
        ligand assignment has to win over raw distances, so the dictionary can
        be written to and read back from JSON and is meant to be reviewed before
        use. Everything downstream reads the dictionary, never the geometry.

        include_residue is the same principle applied before the fact: a
        residue named there is made a ligand of its nearest metal however far
        away it sits, through its closest sidechain donor. That is how a
        design whose ligand has not been placed properly yet, or a contact the
        cutoffs are too tight for, gets into the active site. The request is
        recorded in the dictionary, so update_binding_modes honours it too and
        a relaxation cannot drop what was asked for.

        :param topology:
            The OpenMM topology.
        :param positions:
            The positions as an (N, 3) numpy array in Angstrom.
        :param include_residue:
            Residues that must be ligands whatever their distance, each given
            as a residue id ('130' or 130) or as a residue label ('ASP130').

        :return:
            The binding modes dictionary.
        """

        positions = np.asarray(positions)

        metals = []
        for atom in topology.atoms():
            if atom.element is None:
                continue
            if atom.element.symbol not in self.metal_elements:
                continue
            symbol = atom.element.symbol
            metals.append({
                'index':
                atom.index,
                'element':
                symbol,
                'chain':
                atom.residue.chain.id,
                'resid':
                atom.residue.id,
                'res_index':
                atom.residue.index,
                'formal_charge':
                self.metal_formal_charges.get(symbol, 2),
            })

        assert_msg_critical(
            len(metals) > 0,
            'MetalSiteForceFieldBuilder.suggest_binding_modes: no metal atom '
            f'found. Recognized elements: {self.metal_elements}')

        self._check_supported_metals(metals, 'suggest_binding_modes')

        forced = self._resolve_residues(topology, include_residue)

        notes = []
        ligands = self._collect_ligands(topology.atoms(),
                                        lambda index: positions[index], metals,
                                        notes, forced)

        binding_modes = {
            'metals': metals,
            'ligands': ligands,
            'variants': {},
            'cutoffs': {
                'primary': self.metal_bond_cutoff,
                'secondary': self.report_cutoff,
            },
            'include_residue': sorted(forced),
            'notes': notes,
        }

        return binding_modes

    def _collect_ligands(self, atoms, position_of, metals, notes, forced=None):
        """
        Builds the classified ligand contact list of a set of candidate atoms.

        Shared by suggest_binding_modes, which offers it every atom of the
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

        :return:
            The list of ligand contacts, each carrying its mode.
        """

        metal_indices = [metal['index'] for metal in metals]
        atoms = list(atoms)
        contacts = []

        for atom in atoms:
            if atom.element is None or atom.index in metal_indices:
                continue
            if atom.element.symbol not in self.DONOR_ELEMENTS:
                continue

            position = position_of(atom.index)
            distances = {
                index: float(np.linalg.norm(position - position_of(index)))
                for index in metal_indices
            }
            closest = min(distances.values())
            if closest > self.report_cutoff:
                continue

            label = f'{atom.residue.name}{atom.residue.id}'

            if atom.name in self.BACKBONE_ATOM_NAMES:
                notes.append(
                    f'backbone atom {label} {atom.name} is {closest:.2f} A '
                    'from a metal; the truncation scheme is sidechain-only, '
                    'so it is not treated as a ligand')
                continue

            bonded_to = sorted([
                index for index, dist in distances.items()
                if dist <= self.metal_bond_cutoff
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
                    f'primary ({self.metal_bond_cutoff}) and secondary '
                    f'({self.report_cutoff}) cutoffs; review whether it '
                    'should be a ligand')

        self._force_ligands(atoms, position_of, metal_indices, contacts, notes,
                            forced)

        ligands = [contact for contact in contacts if contact['metals']]

        self._assign_binding_modes(ligands, notes, self.bidentate_asymmetry)

        for metal in metals:
            n_ligands = sum(1 for ligand in ligands
                            if metal['index'] in ligand['metals'])
            if n_ligands < 3:
                notes.append(
                    f'metal {metal["element"]} (index {metal["index"]}) has '
                    f'only {n_ligands} ligand(s) within the primary cutoff; '
                    'check for a missing bridging ligand or water')

        return ligands

    def _resolve_residues(self, topology, include_residue):
        """
        Turns the residues asked for into residue indices of the topology.

        A request is matched against the residue id ('130' or 130) and against
        the label the binding modes report ('ASP130'). A request that matches
        nothing is an error rather than a silently ignored line, and one that
        matches several residues - the same number in two chains - takes all
        of them and says which.

        :param topology:
            The OpenMM topology.
        :param include_residue:
            The residues asked for, or None.

        :return:
            The residue indices, as a set.
        """

        if not include_residue:
            return set()

        if isinstance(include_residue, (str, int)):
            include_residue = [include_residue]

        residues = list(topology.residues())
        forced = set()

        for request in include_residue:
            wanted = str(request).strip()
            matched = [
                residue for residue in residues
                if wanted in (str(residue.id), f'{residue.name}{residue.id}')
            ]

            assert_msg_critical(
                len(matched) > 0,
                'MetalSiteForceFieldBuilder: include_residue asked for '
                f'{request}, which is not a residue of this structure')

            for residue in matched:
                forced.add(residue.index)

            found = ', '.join(f'{residue.name}{residue.id} '
                              f'(chain {residue.chain.id})'
                              for residue in matched)
            self.ostream.print_info(
                f'Including {found} in the coordination sphere by request.')

        self.ostream.flush()

        return forced

    def _force_ligands(self, atoms, position_of, metal_indices, contacts, notes,
                       forced):
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
                and atom.element.symbol in self.DONOR_ELEMENTS
                and atom.name not in self.BACKBONE_ATOM_NAMES
                and atom.index not in metal_indices
            ]

            assert_msg_critical(
                len(donors) > 0,
                'MetalSiteForceFieldBuilder: include_residue asked for '
                f'residue index {res_index}, whose sidechain has no '
                f'{list(self.DONOR_ELEMENTS)} atom to coordinate with')

            best = None
            for atom in donors:
                position = position_of(atom.index)
                for metal in metal_indices:
                    distance = float(
                        np.linalg.norm(position - position_of(metal)))
                    if best is None or distance < best[0]:
                        best = (distance, atom, metal)

            distance, atom, metal = best
            label = f'{atom.residue.name}{atom.residue.id}'

            notes.append(
                f'{label} {atom.name} is {distance:.2f} A from a metal, '
                f'beyond the primary cutoff ({self.metal_bond_cutoff}), and '
                'was made a ligand because include_residue asked for it')

            if distance > self.report_cutoff:
                self.ostream.print_warning(
                    f'{label} {atom.name} is {distance:.2f} A from its metal, '
                    'which is a long way for a bond. It is a ligand because '
                    'include_residue asked for it.')
                self.ostream.flush()

            existing = [
                contact for contact in contacts
                if contact['index'] == atom.index
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

    @staticmethod
    def _assign_binding_modes(ligands, notes, bidentate_asymmetry=None):
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
        :param bidentate_asymmetry:
            How much closer one carboxylate oxygen has to be than the other
            before the pair is read as monodentate. Defaults to
            BIDENTATE_ASYMMETRY.
        """

        if bidentate_asymmetry is None:
            bidentate_asymmetry = (
                MetalSiteForceFieldBuilder.BIDENTATE_ASYMMETRY)

        by_residue = {}
        for ligand in ligands:
            by_residue.setdefault(ligand['res_index'], []).append(ligand)

        can_bridge = MetalSiteForceFieldBuilder.BRIDGING_RESIDUES
        carboxylates = ('ASP', 'GLU', 'ASH', 'GLH')

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
                    closest = ligand['distances'].index(min(
                        ligand['distances']))
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

            monodentate = all(ligand['mode'] == 'monodentate'
                              for ligand in group)

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
                        far = [
                            ligand for ligand in group if ligand is not nearest
                        ]
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
                near, far = sorted(group,
                                   key=lambda ligand: ligand['distances'][0])
                separation = far['distances'][0] - near['distances'][0]
                if separation > bidentate_asymmetry:
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

    def save_binding_modes(self, filename, binding_modes):
        """
        Writes the binding modes to a JSON file for review.

        :param filename:
            The name of the JSON file.
        """

        if self.rank == mpi_master():
            Path(filename).write_text(json.dumps(binding_modes, indent=2))

    def load_binding_modes(self, filename):
        """
        Reads reviewed binding modes back from a JSON file.

        :param filename:
            The name of the JSON file.

        :return:
            The binding modes dictionary.
        """

        binding_modes = json.loads(Path(filename).read_text())
        # JSON keys are always strings
        binding_modes['variants'] = {
            int(key): value
            for key, value in binding_modes.get('variants', {}).items()
        }
        self._check_supported_metals(binding_modes.get('metals', []),
                                     'load_binding_modes')

        return binding_modes

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
            if atom.name not in self.BACKBONE_ATOM_NAMES and atom.name != 'CA'
        ]
        ring_nitrogens = {
            atom.name: atom
            for atom in sidechain if atom.name in ('ND1', 'NE2')
        }

        if len(ring_nitrogens) < 2:
            return 'HID', (
                f'{residue.name}{residue.id} does not have both ring '
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
            note = (
                f'{residue.name}{residue.id} has {nearest.name} closer to a '
                f'metal ({closest_metal_distance(nearest):.2f} A) than either '
                f'ring nitrogen (ND1 {distances["ND1"]:.2f} A, NE2 '
                f'{distances["NE2"]:.2f} A); the tautomer was still chosen '
                'from the nitrogens, but the geometry looks distorted')

        return variant, note

    def suggest_variants(self, topology, positions, binding_modes):
        """
        Chooses a protonation variant for each coordinating residue.

        coordinating carboxylates are deprotonated
        a coordinating cysteine is a thiolate
        and the histidine tautomer is set so that the coordinating nitrogen
        carries no hydrogen. Entries in self.protonation_overrides win over
        the rules.

        :param topology:
            The OpenMM topology.
        :param positions:
            The positions as an (N, 3) numpy array in Angstrom.
        :param binding_modes:
            The binding modes. Not modified; anything worth recording is
            returned as a note.

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
                variant, note = self._histidine_variant(residues[res_index],
                                                        positions,
                                                        metal_positions)
                if note is not None:
                    notes.append(note)
            else:
                variant = None

            if variant is not None:
                variants[res_index] = variant

        overrides = self.protonation_overrides or {}
        for key, variant in overrides.items():
            matched = False
            for residue in residues:
                if str(residue.id) == str(key) or residue.index == key:
                    variants[residue.index] = variant
                    matched = True
            assert_msg_critical(
                matched,
                'MetalSiteForceFieldBuilder.suggest_variants: override '
                f'residue {key} not found')

        return variants, notes

    def protonate(self, topology, positions, binding_modes):
        """
        Adds hydrogens with the protonation variants that the metal site
        requires.

        Adding hydrogens renumbers the atoms, so every index stored in the
        binding modes is invalidated. A corrected copy is returned rather than
        the argument being rewritten in place, which keeps the renumbering
        visible at the call site.

        :param topology:
            The OpenMM topology.
        :param positions:
            The positions as an (N, 3) numpy array in Angstrom.
        :param binding_modes:
            The binding modes. Not modified.

        :return:
            The tuple of the protonated topology, the positions in Angstrom
            and the binding modes with their indices remapped.
        """

        assert_msg_critical(
            'openmm' in sys.modules,
            'MetalSiteForceFieldBuilder.protonate: openmm is required')

        variants_by_index, notes = self.suggest_variants(
            topology, positions, binding_modes)

        # the atom indices of the binding modes are about to be invalidated,
        # so the caller gets a corrected copy rather than a silently rewritten
        # version of what it passed in
        binding_modes = deepcopy(binding_modes)
        binding_modes['notes'].extend(notes)

        variant_list = [None] * topology.getNumResidues()
        for res_index, variant in variants_by_index.items():
            variant_list[res_index] = variant

        modeller = mmapp.Modeller(topology,
                                  np.asarray(positions) * mmunit.angstrom)
        modeller.addHydrogens(variants=variant_list)

        new_topology = modeller.topology
        new_positions = np.array(
            modeller.positions.value_in_unit(mmunit.angstrom))

        # residue indices survive adding hydrogens, atom indices do not
        lookup = {
            (atom.residue.index, atom.name): atom.index
            for atom in new_topology.atoms()
        }
        old_atoms = list(topology.atoms())

        def remap(old_index):
            old_atom = old_atoms[old_index]
            key = (old_atom.residue.index, old_atom.name)
            assert_msg_critical(
                key in lookup, 'MetalSiteForceFieldBuilder.protonate: atom '
                f'{old_atom.name} of residue {old_atom.residue.name}'
                f'{old_atom.residue.id} vanished while adding hydrogens')
            return lookup[key]

        # remap takes an index into the old topology, so the ligands' stored
        # metal indices can be remapped independently of the metal entries
        for ligand in binding_modes['ligands']:
            ligand['index'] = remap(ligand['index'])
            ligand['metals'] = [remap(index) for index in ligand['metals']]

        for metal in binding_modes['metals']:
            metal['index'] = remap(metal['index'])

        binding_modes['variants'] = variants_by_index

        return new_topology, new_positions, binding_modes

    def extract_active_site(self, topology, positions, binding_modes):
        """
        Builds the truncated QM active site.

        Sidechains are cut at the CA-CB bond and capped with a hydrogen placed
        along the CB to CA direction. No second-shell fragments and no
        backbone: the scheme is fixed, which also guarantees that the RESP and
        Hessian calculations see the same truncation.

        :param topology:
            The protonated OpenMM topology.
        :param positions:
            The positions as an (N, 3) numpy array in Angstrom.

        :return:
            The active site dictionary. It records the active site indices of
            the capping hydrogens and of the beta carbons, the map back to
            the topology, and the charge.
        """

        self._check_supported_metals(binding_modes['metals'],
                                     'extract_active_site')

        positions = np.asarray(positions)
        residues = list(topology.residues())

        res_indices = sorted(
            {ligand['res_index']
             for ligand in binding_modes['ligands']})

        labels = []
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

        for res_index in res_indices:
            residue = residues[res_index]
            res_atoms = list(residue.atoms())

            for atom in res_atoms:
                if atom.name in self.BACKBONE_ATOM_NAMES:
                    continue

                if atom.name == 'CA':
                    cb_atom = None
                    for other in res_atoms:
                        if other.name == 'CB':
                            cb_atom = other
                            break
                    assert_msg_critical(
                        cb_atom is not None,
                        'MetalSiteForceFieldBuilder.extract_active_site: '
                        f'residue {residue.name}{residue.id} has no CB to '
                        'cut at')
                    direction = positions[atom.index] - positions[cb_atom.index]
                    direction /= np.linalg.norm(direction)
                    # the cap is mapped to the CA it replaces, so that the
                    # CA-CB bond of the topology becomes the cap-CB bond of
                    # the active site
                    atom_map[len(coords)] = atom.index
                    cap_indices.append(len(coords))
                    coords.append(positions[cb_atom.index] +
                                  direction * self.cap_bond_length)
                    labels.append('H')
                else:
                    if atom.name == 'CB':
                        beta_carbon_indices.append(len(coords))
                    atom_map[len(coords)] = atom.index
                    coords.append(positions[atom.index])
                    labels.append(atom.element.symbol)

        charge = sum(metal['formal_charge']
                     for metal in binding_modes['metals'])

        for res_index in res_indices:
            residue = residues[res_index]
            variant = binding_modes['variants'].get(res_index)
            if variant is None:
                variant = residue.name
                self.ostream.print_warning(
                    'No protonation variant recorded for '
                    f'{residue.name}{residue.id}; using the residue name for '
                    'the charge count')
            assert_msg_critical(
                variant in self.VARIANT_CHARGES,
                'MetalSiteForceFieldBuilder.extract_active_site: no charge '
                f'known for variant {variant}')
            charge += self.VARIANT_CHARGES[variant]

        molecule = Molecule(labels, np.array(coords))
        molecule.set_charge(charge)
        molecule.set_multiplicity(1)

        active_site = {
            'molecule': molecule,
            'labels': labels,
            'atom_map': atom_map,
            'cap_indices': cap_indices,
            'beta_carbon_indices': beta_carbon_indices,
            'metal_indices': metal_indices,
            'charge': charge,
            'multiplicity': 1,
            'residues':
            [f'{residues[i].name}{residues[i].id}' for i in res_indices],
            'res_indices': res_indices,
        }

        return active_site

    def build_connectivity(self,
                           topology,
                           active_site,
                           binding_modes,
                           warn_above=2.0):
        """
        Builds the connectivity matrix of the active site without perceiving any
        bonds.

        Covalent bonds come from the OpenMM topology, which is chemically
        correct by construction for standard residues, and the metal-ligand
        bonds come from the binding modes, which are explicit and reviewable.
        Molecule.get_connectivity_matrix is not used.

        :param topology:
            The protonated OpenMM topology.
        :param warn_above:
            The length in Angstrom above which a non-metal bond is reported.

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
                ligand['index'] in reverse_map,
                'MetalSiteForceFieldBuilder.build_connectivity: ligand atom '
                f'{ligand["residue"]} {ligand["atom"]} is not part of the '
                'extracted active site')
            lig_index = reverse_map[ligand['index']]
            for metal_index in ligand['metals']:
                metal = reverse_map[metal_index]
                connectivity_matrix[metal, lig_index] = 1
                connectivity_matrix[lig_index, metal] = 1

        # nothing but a metal bond should be long
        coords = active_site['molecule'].get_coordinates_in_angstrom()
        labels = active_site['labels']
        metals = set(active_site['metal_indices'])

        for i in range(n_atoms):
            for j in range(i + 1, n_atoms):
                if not connectivity_matrix[i, j]:
                    continue
                if i in metals or j in metals:
                    continue
                distance = np.linalg.norm(coords[i] - coords[j])
                if distance > warn_above:
                    self.ostream.print_warning(
                        f'Non-metal bond {i}-{j} ({labels[i]}-{labels[j]}) '
                        f'is {distance:.2f} A long')

        self.ostream.flush()

        return connectivity_matrix

    def update_binding_modes(self, topology, geometry, active_site,
                             binding_modes, connectivity_matrix):
        """
        Re-detects the coordination sphere on a new active site geometry.

        Relaxing the active site moves the metal-ligand distances: a contact
        that started just inside the primary cutoff can end up outside it, an
        asymmetric carboxylate can open into a monodentate one, and a second
        oxygen can rotate onto a metal. The rules of suggest_binding_modes are
        applied again to the new coordinates so that what is fitted afterwards
        is the coordination the geometry actually has. Residues that were
        asked for through include_residue stay ligands, since a distance is
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
            The active site the geometry belongs to. Not modified.
        :param binding_modes:
            The binding modes to check. Not modified.
        :param connectivity_matrix:
            The connectivity matrix to check. Not modified.

        :return:
            The tuple of the binding modes, the connectivity matrix and the
            flag that is True when the coordination changed.
        """

        self._check_supported_metals(binding_modes['metals'],
                                     'update_binding_modes')

        if isinstance(geometry, Molecule):
            coordinates = geometry.get_coordinates_in_angstrom()
        else:
            coordinates = np.asarray(geometry, dtype=float)

        atom_map = active_site['atom_map']

        assert_msg_critical(
            len(coordinates) == len(atom_map),
            'MetalSiteForceFieldBuilder.update_binding_modes: the geometry '
            f'has {len(coordinates)} atoms while the active site has '
            f'{len(atom_map)}')

        # a capping hydrogen is mapped to the CA it replaces, so its position
        # must not be read back as that CA's
        cap_atoms = {
            atom_map[site_index]
            for site_index in active_site['cap_indices']
        }
        site_of = {
            top_index: site_index
            for site_index, top_index in atom_map.items()
            if top_index not in cap_atoms
        }

        atoms = list(topology.atoms())
        candidates = [atoms[top_index] for top_index in sorted(site_of)]

        notes = []
        new_ligands = self._collect_ligands(
            candidates, lambda index: coordinates[site_of[index]],
            binding_modes['metals'], notes,
            set(binding_modes.get('include_residue', [])))

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
                ligand['index'] in site_of,
                'MetalSiteForceFieldBuilder.update_binding_modes: ligand atom '
                f'{ligand["residue"]} {ligand["atom"]} is not part of the '
                'active site the geometry belongs to')
            lig_index = site_of[ligand['index']]
            for metal_index, distance in zip(ligand['metals'],
                                             ligand['distances']):
                moved = float(
                    np.linalg.norm(coordinates[lig_index] -
                                   coordinates[site_of[metal_index]]))
                shifts.append(abs(moved - distance))

        largest_shift = max(shifts) if shifts else 0.0

        if new_coordination == old_coordination:
            self._print_binding_mode_update([], largest_shift)
            return binding_modes, connectivity_matrix, False

        changes = []
        for index in sorted(set(old_coordination) | set(new_coordination)):
            old = old_ligands.get(index)
            new = new_by_index.get(index)

            if old is None:
                changes.append(
                    ('gained', f'{new["residue"]} {new["atom"]}',
                     f'{self._metal_contact_label(new, binding_modes)}, '
                     f'{new["mode"]}'))
            elif new is None:
                changes.append(
                    ('lost', f'{old["residue"]} {old["atom"]}', 'was '
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

        self._print_binding_mode_update(changes, largest_shift,
                                        dropped_residues)

        new_binding_modes = {
            'metals': deepcopy(binding_modes['metals']),
            'ligands': new_ligands,
            'variants': deepcopy(binding_modes.get('variants', {})),
            'cutoffs': {
                'primary': self.metal_bond_cutoff,
                'secondary': self.report_cutoff,
            },
            'notes': notes,
        }

        # only the metal-ligand bonds are read off the geometry, so the
        # covalent bonds of the matrix are left exactly as they were
        new_matrix = np.array(connectivity_matrix, copy=True)

        for ligands, value in ((binding_modes['ligands'], 0), (new_ligands, 1)):
            for ligand in ligands:
                lig_index = site_of[ligand['index']]
                for metal_index in ligand['metals']:
                    metal = site_of[metal_index]
                    new_matrix[metal, lig_index] = value
                    new_matrix[lig_index, metal] = value

        self._print_binding_modes(new_binding_modes)

        return new_binding_modes, new_matrix, True

    @staticmethod
    def _metal_contact_label(ligand, binding_modes):
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

    @staticmethod
    def extract_pairs(connectivity_matrix,
                      source_atoms,
                      bond_count=2,
                      initial_bond_range=None,
                      coordinates=None):
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
        :param initial_bond_range:
            The distance in Angstrom within which the first shell is taken.
            Lets the function run on a topology in which the metal has no bonds
            yet. Requires coordinates.
        :param coordinates:
            The coordinates in Angstrom, only needed with initial_bond_range.

        :return:
            The tuple of the sorted pair list and the sorted atom list.
        """

        connectivity_matrix = np.asarray(connectivity_matrix)
        source_atoms = list(source_atoms)

        assert_msg_critical(
            initial_bond_range is None or coordinates is not None,
            'MetalSiteForceFieldBuilder.extract_pairs: initial_bond_range '
            'requires coordinates')

        def neighbors(index, depth):
            if depth == 0 and initial_bond_range is not None:
                coords = np.asarray(coordinates)
                distances = np.linalg.norm(coords - coords[index], axis=1)
                found = set(
                    np.where(distances <= initial_bond_range)[0].tolist())
                found.discard(index)
                return found
            return set(np.where(connectivity_matrix[index])[0].tolist())

        visited = set(source_atoms)
        pairs = set()
        frontier = list(source_atoms)

        for depth in range(bond_count):
            next_frontier = []
            for index in frontier:
                for neighbor in neighbors(index, depth):
                    pairs.add((min(index, neighbor), max(index, neighbor)))
                    if neighbor not in visited:
                        visited.add(neighbor)
                        next_frontier.append(neighbor)
            frontier = next_frontier

        pairs = sorted(pairs)
        atoms = sorted({index for pair in pairs for index in pair})

        return pairs, atoms

    def _resolve_optimized_geometry(self, active_site):
        """
        Validates a geometry supplied through optimized_geometry, or left
        behind in the working folder by an earlier run. The element sequence is checked against the extracted active site
        :param active_site:
            The active site, to validate against.

        :return:
            The molecule, or None if there is nothing to use.
        """

        source = self.optimized_geometry

        if source is None:
            source = self._folder_file(self.GEOMETRY_FILE)
            if source is None:
                return None
            self.ostream.print_info(
                f'Reusing {self.GEOMETRY_FILE} from {self.folder}.')
            self.ostream.flush()

        if isinstance(source, Molecule):
            molecule = Molecule(source)
        else:
            path = Path(source)
            assert_msg_critical(
                path.is_file(),
                'MetalSiteForceFieldBuilder: optimized_geometry file '
                f'{path} not found')
            molecule = Molecule.read_xyz_file(str(path))

        active_site = active_site['molecule']

        assert_msg_critical(
            molecule.number_of_atoms() == active_site.number_of_atoms(),
            'MetalSiteForceFieldBuilder: optimized_geometry has '
            f'{molecule.number_of_atoms()} atoms but the extracted active site '
            f'has {active_site.number_of_atoms()}')

        assert_msg_critical(
            list(molecule.get_labels()) == list(active_site.get_labels()),
            'MetalSiteForceFieldBuilder: the elements of '
            'optimized_geometry do not match the extracted active site, so it '
            'describes a different structure')

        molecule.set_charge(active_site.get_charge())
        molecule.set_multiplicity(active_site.get_multiplicity())

        return molecule

    def _resolve_hessian(self, active_site):
        """
        Validates a Hessian supplied through the hessian setting, or left
        behind in the working folder by an earlier run.

        :param active_site:
            The active site, to validate the shape against.

        :return:
            The Hessian, or None if there is nothing to use.
        """

        source = self.hessian

        if source is None:
            source = self._folder_file(self.HESSIAN_FILE)
            if source is None:
                return None
            self.ostream.print_info(
                f'Reusing {self.HESSIAN_FILE} from {self.folder}.')
            self.ostream.flush()

        assert_msg_critical(
            not isinstance(source, (str, Path)) or Path(source).is_file(),
            f'MetalSiteForceFieldBuilder: hessian file {source} not found')

        if isinstance(source, (str, Path)):
            hessian = np.loadtxt(source)
        else:
            hessian = np.asarray(source)

        n_atoms = active_site['molecule'].number_of_atoms()
        expected = (3 * n_atoms, 3 * n_atoms)

        assert_msg_critical(
            hessian.shape == expected,
            f'MetalSiteForceFieldBuilder: hessian has shape {hessian.shape} '
            f'but the extracted active site needs {expected}')

        return hessian

    def _resolve_partial_charges(self, active_site):
        """
        Validates charges supplied through the partial_charges setting, or
        left behind in the working folder by an earlier run.

        :param active_site:
            The active site, to validate the count against.

        :return:
            The partial charges, or None if there is nothing to use.
        """

        source = self.partial_charges

        if source is None:
            source = self._folder_file(self.CHARGES_FILE)
            if source is None:
                return None
            self.ostream.print_info(
                f'Reusing {self.CHARGES_FILE} from {self.folder}.')
            self.ostream.flush()

        if isinstance(source, (str, Path)):
            assert_msg_critical(
                Path(source).is_file(),
                'MetalSiteForceFieldBuilder: partial_charges file '
                f'{source} not found')
            charges = np.loadtxt(source)
        else:
            charges = np.asarray(source)

        n_atoms = active_site['molecule'].number_of_atoms()

        assert_msg_critical(
            charges.shape == (n_atoms, ),
            f'MetalSiteForceFieldBuilder: partial_charges has shape '
            f'{charges.shape} but the extracted active site has {n_atoms} atoms'
        )

        total = float(np.sum(charges))
        expected = active_site['charge']

        if abs(total - expected) > 1.0e-3:
            self.ostream.print_warning(
                f'The supplied partial charges sum to {total:+.3f}, but the '
                f'active site charge is {expected:+d}')
            self.ostream.flush()

        return charges

    def _get_scf_driver(self, molecule):
        """
        Returns the SCF driver and basis set for the active site

        The driver assigned to scf_drv is used exactly as given, so it carries
        every QM setting beyond the functional and the basis set.

        :param molecule:
            The active site molecule.

        :return:
            The tuple of the SCF driver and the basis set.
        """

        if self.scf_drv is None:
            if molecule.get_multiplicity() != 1:
                scf_drv = ScfUnrestrictedDriver(self.comm, self.ostream)
            else:
                scf_drv = ScfRestrictedDriver(self.comm, self.ostream)
            if self.xcfun is not None:
                scf_drv.xcfun = self.xcfun
            self.scf_drv = scf_drv

        basis = MolecularBasis.read(molecule, self.basis_set_label)

        return self.scf_drv, basis

    def _print_muted_notice(self, step):
        """
        Announces a long calculation whose output is being suppressed.

        :param step:
            A description of the step about to run.
        """

        if self.mute_scf:
            self.ostream.print_info(
                f'Running {step} with muted QM output. Set mute_scf to False '
                'to follow it.')
        else:
            self.ostream.print_info(f'Running {step}.')
        self.ostream.flush()

    def constrained_indices(self, active_site):
        """
        Returns the active site indices held fixed during optimization.

        The beta carbons are constrained by default: that is where the
        backbone actually holds the sidechain in place. The capping hydrogens
        merely stand in for the alpha carbons, and constraining them as well
        makes the truncated fragment rigid, which is available through
        constrain_capping_hydrogens.

        :return:
            The sorted list of active site indices.
        """

        indices = set(active_site['beta_carbon_indices'])

        if self.constrain_capping_hydrogens:
            indices |= set(active_site['cap_indices'])

        return sorted(indices)

    def optimize_active_site(self, active_site, frozen_indices=None):
        """
        Optimizes the active site with the beta carbons frozen.

        Freezing them keeps the spatial arrangement imposed by the protein
        backbone. Without it the site relaxes to a gas-phase geometry, and
        since the Seminario method takes the equilibrium values straight from
        the geometry, every fitted bond length and angle would then describe
        the wrong structure.

        :param frozen_indices:
            The zero-based active site indices to freeze. Defaults to
            constrained_indices().

        :return:
            The optimized molecule. The caller decides whether to put it back
            into the active site.
        """

        if frozen_indices is None:
            frozen_indices = self.constrained_indices(active_site)

        molecule = active_site['molecule']

        constraints = None
        if len(frozen_indices) > 0:
            # geomeTRIC indexes atoms from one
            selection = ','.join(
                str(index + 1) for index in sorted(frozen_indices))
            constraints = [f'freeze xyz {selection}']

        self._print_muted_notice(
            f'the constrained optimization with {len(frozen_indices)} '
            'atom(s) frozen')

        scf_drv, basis = self._get_scf_driver(molecule)

        grad_drv = ScfGradientDriver(scf_drv)
        opt_drv = OptimizationDriver(grad_drv)
        opt_drv.constraints = constraints

        if self.mute_scf:
            grad_drv.ostream.mute()
            opt_drv.ostream.mute()

        opt_results = opt_drv.compute(molecule, basis)

        if self.mute_scf:
            grad_drv.ostream.unmute()
            opt_drv.ostream.unmute()

        optimized = Molecule.read_xyz_string(opt_results['final_geometry'])
        optimized.set_charge(molecule.get_charge())
        optimized.set_multiplicity(molecule.get_multiplicity())

        self._save_intermediate(
            self.GEOMETRY_FILE,
            lambda path: optimized.write_xyz_file(str(path)))

        return optimized

    def compute_hessian(self, active_site, atom_pairs=None):
        """
        Computes the nuclear Hessian of the active site.

        With atom_pairs the analytical Hessian is restricted to those blocks
        and everything else is left at zero, which is all the Seminario method
        needs when only the metal terms are being fitted. The diagonal blocks
        are added by the Hessian driver itself, so only the off-diagonal pairs
        need to be given.

        :param atom_pairs:
            The list of zero-based (i, j) tuples, typically from
            extract_pairs. None computes the full Hessian.

        :return:
            The Hessian as a (3N, 3N) numpy array in Hartree per Bohr squared.
        """

        molecule = active_site['molecule']

        scf_drv, basis = self._get_scf_driver(molecule)

        assert_msg_critical(
            scf_drv.solvation_model is None,
            'MetalSiteForceFieldBuilder.compute_hessian: ScfHessianDriver '
            'does not support a solvation model')

        # The Hessian driver reuses scf_drv.scf_results without checking
        # which geometry they belong to, so the SCF is run here for the
        # current one. This costs nothing: the driver would otherwise run the
        # same SCF itself.
        self._print_muted_notice('the SCF for the Hessian')

        self._print_muted_notice('the Hessian')

        hessian_drv = ScfHessianDriver(scf_drv)
        # the numerical path ignores atom_pairs entirely and would silently
        # compute the full Hessian instead
        hessian_drv.numerical = False
        if atom_pairs is None:
            hessian_drv.atom_pairs = None
        else:
            hessian_drv.atom_pairs = [tuple(pair) for pair in atom_pairs]

        if self.mute_scf:
            hessian_drv.ostream.mute()

        hessian_drv.compute(molecule, basis)

        if self.mute_scf:
            hessian_drv.ostream.unmute()

        hessian = np.copy(hessian_drv.hessian)

        self._save_intermediate(self.HESSIAN_FILE,
                                lambda path: np.savetxt(path, hessian))

        return hessian

    def compute_resp_charges(self, active_site):
        """
        Computes RESP charges for the active site.

        :return:
            The partial charges as an (N,) numpy array.
        """

        molecule = active_site['molecule']

        self._print_muted_notice('the RESP charge fit at Hartree-Fock/6-31G*')

        resp_drv = RespChargesDriver(self.comm, self.ostream)

        if self.mute_scf:
            resp_drv.ostream.mute()

        # Neither a basis nor SCF results are passed: the driver then defaults
        # to Hartree-Fock with 6-31G*, which is what RESP charges are meant to
        # be fitted to, and runs its own SCF. Handing it the active site's own
        # functional and basis would silently fit the charges at a level the
        # RESP parameters were never derived for.
        charges = resp_drv.compute(molecule)

        if self.mute_scf:
            resp_drv.ostream.unmute()

        charges = self.comm.bcast(charges, root=mpi_master())
        charges = np.array(charges)

        self._save_intermediate(self.CHARGES_FILE,
                                lambda path: np.savetxt(path, charges))

        return charges

    def get_metal_keys(self, forcefield, active_site):
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

    def build_forcefield(self,
                         active_site,
                         connectivity_matrix,
                         hessian=None,
                         partial_charges=None):
        """
        Builds the active site force field and fits the metal terms.

        Without a Hessian the metal terms are seeded by _seed_metal_terms
        instead of fitted, which is the force field the crude pre-QM pass
        runs on: getting the equilibrium geometry roughly right matters far
        more than the stiffness at that stage.

        :param hessian:
            The Hessian as a (3N, 3N) numpy array. The string 'xtb' is
            rejected: it would make reparameterize re-optimize the active site
            without constraints, silently replacing every equilibrium value
            with a gas-phase one.
        :param partial_charges:
            The partial charges fitted on the active site. The charge of the
            capping hydrogens is redistributed over the remaining atoms before
            they are applied, since the caps do not exist in the protein.

        :return:
            The force field generator.
        """

        assert_msg_critical(
            not isinstance(hessian, str),
            'MetalSiteForceFieldBuilder.build_forcefield: the Hessian must be '
            'a numpy array.')

        molecule = active_site['molecule']
        n_atoms = molecule.number_of_atoms()

        forcefield = MMForceFieldGenerator(self.comm, self.ostream)
        forcefield.connectivity_matrix = np.asarray(connectivity_matrix)
        forcefield.topology_update_flag = True

        # the generator shares the output stream of the builder, so muting it
        # has to be balanced before anything of our own is printed
        forcefield.ostream.mute()
        forcefield.create_topology(molecule, resp=False)
        forcefield.ostream.unmute()

        # The charges have to be applied after create_topology: setting
        # topology_update_flag, which is what makes the custom connectivity
        # take effect, also resets partial_charges to None, and resp=False
        # then fills them with zeros. Assigning them beforehand is silently
        # discarded.
        if partial_charges is not None:
            partial_charges = np.asarray(partial_charges)
            assert_msg_critical(
                partial_charges.shape == (n_atoms, ),
                'MetalSiteForceFieldBuilder.build_forcefield: expected '
                f'{n_atoms} partial charges, got {partial_charges.shape}')
            # The capping hydrogens stand in for alpha carbons and do not
            # exist anywhere the force field is used, so their charge is
            # folded into the rest of the active site before anything is written.
            partial_charges = self.redistribute_cap_charges(
                active_site, partial_charges)
            forcefield.partial_charges = partial_charges
            for index in range(n_atoms):
                forcefield.atoms[index]['charge'] = partial_charges[index]

        self.annotate_atoms(forcefield, active_site)

        bonds, angles = self.get_metal_keys(forcefield, active_site)

        # switching the angles off has to leave them at whatever the generator
        # guessed in both passes, so it is gated once, here
        if not self.reparameterize_metal_angles:
            angles = []

        if hessian is None:
            self._seed_metal_terms(forcefield, active_site, bonds, angles)
        else:
            hessian = np.asarray(hessian)
            assert_msg_critical(
                hessian.shape == (3 * n_atoms, 3 * n_atoms),
                'MetalSiteForceFieldBuilder.build_forcefield: Hessian shape '
                f'{hessian.shape} does not match {(3 * n_atoms, 3 * n_atoms)}')
            forcefield.ostream.mute()
            forcefield.reparameterize(
                hessian,
                reparameterize_keys=bonds + angles,
                average_metal_terms=self.average_metal_terms)
            forcefield.ostream.unmute()
            self._check_force_constants(forcefield, active_site, bonds, angles)

        return forcefield

    def save_forcefield(self, filename, forcefield):
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

        if self.rank == mpi_master():
            MMForceFieldGenerator.save_forcefield_as_json(
                forcefield, str(filename))

    def load_forcefield(self, filename, active_site=None):
        """
        Reads a force field back from a JSON file.

        The file carries the parameters alone, so the molecule of the active
        site is put back onto the generator: without it the force field cannot
        be written as OpenMM files or minimized. partial_charges is restored
        from the charges of the atoms, which are the redistributed ones that
        build_forcefield wrote.

        :param filename:
            The name of the JSON file.
        :param active_site:
            The active site the force field belongs to. When given, the atom
            count and the elements are checked against it and its molecule is
            attached to the generator.

        :return:
            The force field generator.
        """

        assert_msg_critical(
            Path(filename).is_file(),
            f'MetalSiteForceFieldBuilder: forcefield file {filename} not found')

        forcefield = MMForceFieldGenerator.load_forcefield_from_json_file(
            str(filename))

        if active_site is not None:
            self._check_forcefield(forcefield, active_site)
            forcefield.molecule = active_site['molecule']

        forcefield.partial_charges = np.array(
            [atom['charge'] for atom in forcefield.atoms.values()])

        return forcefield

    @staticmethod
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

    def _check_forcefield(self, forcefield, active_site):
        """
        Checks that a force field describes the extracted active site.

        Its bonds and angles are keyed by plain atom indices, which fit any
        cluster of the same size, so a force field belonging to another site
        would otherwise be applied to this one without a word.

        :param forcefield:
            The force field generator.
        :param active_site:
            The active site, to validate against.
        """

        labels = list(active_site['labels'])
        elements = self._forcefield_elements(forcefield)

        assert_msg_critical(
            len(elements) == len(labels),
            f'MetalSiteForceFieldBuilder: the force field has {len(elements)} '
            f'atoms but the extracted active site has {len(labels)}')

        assert_msg_critical(
            elements == labels,
            'MetalSiteForceFieldBuilder: the elements of the force field do '
            'not match the extracted active site, so it describes a different '
            'structure')

    def annotate_atoms(self, forcefield, active_site):
        """
        Writes what the truncation knows about an atom into its comment.

        A force field carries atom types, charges and bonds, and a geometry
        carries positions; neither of them says which carbon a sidechain was
        cut at. That is known here, while the active site is still in hand,
        and the comment is the one field that survives into forcefield.json
        for anything downstream to read it back out of.

        Comments the generator already wrote are kept, and running this twice
        does not write the same note twice.

        :param forcefield:
            The force field generator to annotate.
        :param active_site:
            The active site, for the atoms the truncation created.
        """

        for index in active_site['beta_carbon_indices']:
            atom = forcefield.atoms[index]
            comment = atom.get('comment', '') or ''

            if self.BETA_CARBON_COMMENT in comment:
                continue

            atom['comment'] = '; '.join(part
                                        for part in (comment,
                                                     self.BETA_CARBON_COMMENT)
                                        if part)

    @staticmethod
    def _lookup_equilibrium(table, elements):
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

    def _seed_metal_terms(self, forcefield, active_site, bonds, angles):
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
        """

        labels = active_site['labels']
        molecule = active_site['molecule']

        for key in bonds:
            elements = tuple(labels[index] for index in key)
            equilibrium = self._lookup_equilibrium(self.metal_bond_equilibria,
                                                   elements)
            if equilibrium is None:
                # the getters index atoms from one, and the force field keeps
                # its bond lengths in nanometers
                equilibrium = 0.1 * molecule.get_distance_in_angstroms(
                    [index + 1 for index in key])
                comment = 'measured on the input geometry'
            else:
                comment = 'given equilibrium'
            forcefield.bonds[key]['equilibrium'] = equilibrium
            forcefield.bonds[key]['force_constant'] = (
                self.default_metal_bond_force_constant)
            forcefield.bonds[key]['comment'] = comment

        for key in angles:
            elements = tuple(labels[index] for index in key)
            equilibrium = self._lookup_equilibrium(self.metal_angle_equilibria,
                                                   elements)
            if equilibrium is None:
                equilibrium = molecule.get_angle_in_degrees(
                    [index + 1 for index in key])
                comment = 'measured on the input geometry'
            else:
                comment = 'given equilibrium'
            forcefield.angles[key]['equilibrium'] = equilibrium
            forcefield.angles[key]['force_constant'] = (
                self.default_metal_angle_force_constant)
            forcefield.angles[key]['comment'] = comment

        self.ostream.flush()

    def _check_force_constants(self, forcefield, active_site, bonds, angles):
        """
        Warns about metal terms whose force constant was fitted to zero.

        The Seminario method clamps negative Hessian projections to zero, and
        a negative projection means the geometry is not a stationary point
        along that coordinate. On an unrelaxed structure this typically wipes
        out the long, strained metal-ligand bonds, which are exactly the ones
        that matter for a bridged binuclear site.

        :param bonds:
            The metal bond keys.
        :param angles:
            The metal angle keys.
        """

        labels = active_site['labels']

        zero_bonds = [
            key for key in bonds
            if forcefield.bonds[key]['force_constant'] == 0.0
        ]
        zero_angles = [
            key for key in angles
            if forcefield.angles[key]['force_constant'] == 0.0
        ]

        if not zero_bonds and not zero_angles:
            return

        self.ostream.print_warning(
            f'{len(zero_bonds)} of {len(bonds)} metal bond(s) and '
            f'{len(zero_angles)} of {len(angles)} metal angle(s) got a zero '
            'force constant. The Hessian is most likely not evaluated at a '
            'stationary point; run optimize_active_site first.')

        for key in zero_bonds:
            names = '-'.join(labels[index] for index in key)
            self.ostream.print_warning(f'  zero force constant: {key} {names}')

        self.ostream.flush()

    def minimize_active_site(self,
                             active_site,
                             forcefield,
                             frozen_indices=None,
                             max_iterations=0):
        """
        Minimizes the active site with its own force field.

        The beta carbons are frozen by default. They are where the backbone
        holds the sidechains in place, so a free minimization of an isolated
        active site lets the truncated fragments drift apart and says nothing
        about the metal site. Note that the force field only carries
        electrostatics if build_forcefield was given partial charges.

        :param frozen_indices:
            The active site indices to hold fixed. Defaults to
            constrained_indices(); an empty list minimizes freely.
        :param max_iterations:
            The maximum number of minimization steps. Zero runs until
            convergence.

        :return:
            The minimized coordinates as an (N, 3) numpy array in Angstrom.
            The coordinates make a round trip through a PDB file, so they
            carry a rounding of 0.001 Angstrom.
        """

        assert_msg_critical(
            'openmm' in sys.modules,
            'MetalSiteForceFieldBuilder.minimize_active_site: openmm is '
            'required')

        if frozen_indices is None:
            frozen_indices = self.constrained_indices(active_site)

        with tempfile.TemporaryDirectory() as temp_dir:
            stem = str(Path(temp_dir) / 'active_site')
            forcefield.write_openmm_files(stem)

            pdb = mmapp.PDBFile(f'{stem}.pdb')
            openmm_ff = mmapp.ForceField(f'{stem}.xml')
            system = openmm_ff.createSystem(pdb.topology,
                                            nonbondedMethod=mmapp.NoCutoff)

            # a zero mass makes OpenMM hold the particle fixed
            for index in frozen_indices:
                system.setParticleMass(int(index), 0.0)

            integrator = mm.VerletIntegrator(0.001 * mmunit.picoseconds)
            simulation = mmapp.Simulation(pdb.topology, system, integrator)
            simulation.context.setPositions(pdb.positions)
            simulation.minimizeEnergy(maxIterations=max_iterations)

            state = simulation.context.getState(getPositions=True)
            coords = np.array(state.getPositions().value_in_unit(
                mmunit.angstrom))

        return coords

    def mm_optimize_active_site(self,
                                active_site,
                                connectivity_matrix,
                                frozen_indices=None,
                                constrain_metals=None):
        """
        Relaxes the active site on a crude force field of its own.

        A site comes out of a structure file with its ligands wherever the
        crystallographer or the design left them, and the constrained QM
        optimization then spends its first cycles cleaning that up. This pass
        does the same work at MM cost: the metal terms are seeded by
        _seed_metal_terms, every other term is what the generator assigns, and
        the beta carbons are frozen exactly as they are in the QM
        optimization.

        It runs before the RESP charges exist, so the force field it uses
        carries no electrostatics at all. That is what makes it crude, and why
        it stands in front of optimize_active_site rather than in place of it.

        :param active_site:
            The extracted active site.
        :param connectivity_matrix:
            The connectivity of the active site.
        :param frozen_indices:
            The active site indices to hold fixed. Defaults to
            constrained_indices().
        :param constrain_metals:
            Whether to hold the metal centers as well, added to whatever
            frozen_indices amounts to. Defaults to mm_constrain_metals.

        :return:
            The relaxed molecule. The caller decides whether to put it back
            into the active site.
        """

        molecule = active_site['molecule']

        if constrain_metals is None:
            constrain_metals = self.mm_constrain_metals

        if frozen_indices is None:
            frozen_indices = self.constrained_indices(active_site)

        if constrain_metals:
            frozen_indices = sorted(
                set(frozen_indices) | set(active_site['metal_indices']))

        forcefield = self.build_forcefield(active_site, connectivity_matrix)
        coordinates = self.minimize_active_site(
            active_site,
            forcefield,
            frozen_indices=frozen_indices,
            max_iterations=self.mm_max_iterations)

        relaxed = Molecule(active_site['labels'], coordinates, 'angstrom')
        relaxed.set_charge(molecule.get_charge())
        relaxed.set_multiplicity(molecule.get_multiplicity())

        self._print_mm_optimization(active_site, forcefield, relaxed,
                                    frozen_indices)

        self._save_intermediate(self.MM_GEOMETRY_FILE,
                                lambda path: relaxed.write_xyz_file(str(path)))

        return relaxed

    @staticmethod
    def redistribute_cap_charges(active_site, partial_charges):
        """
        Moves the charge of the capping hydrogens onto the rest of the active site.

        The operation is idempotent: the caps end up at zero, so a second pass
        finds nothing left to move. Both build_forcefield and
        redistribute_charges apply it, and either may be handed charges that
        have already been through it.

        :param active_site:
            The active site, for the indices of the capping hydrogens.
        :param partial_charges:
            The charges fitted on the whole active site.

        :return:
            A copy of the charges with the caps emptied and their charge spread
            over the other atoms.
        """

        charges = np.array(partial_charges, dtype=float)
        caps = set(active_site['cap_indices'])
        rest = [index for index in range(charges.size) if index not in caps]

        assert_msg_critical(
            len(rest) > 0,
            'MetalSiteForceFieldBuilder.redistribute_cap_charges: the active site '
            'is nothing but capping hydrogens')

        cap_charge = sum(charges[index] for index in caps)
        charges[list(caps)] = 0.0
        charges[rest] += cap_charge / len(rest)

        return charges

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
        charges = self.redistribute_cap_charges(active_site, partial_charges)
        shift = self.redistribute_backbone_charges(system, topology,
                                                   active_site, charges)

        return cap_charge, shift

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

        assert_msg_critical(
            'openmm' in sys.modules,
            'MetalSiteForceFieldBuilder.redistribute_backbone_charges: openmm '
            'is required')

        charges = np.asarray(partial_charges)
        caps = set(active_site['cap_indices'])
        atom_map = active_site['atom_map']

        # the caps map to the alpha carbons they replaced, so they have to be
        # left out or CA would be counted as covered by the active site
        covered = {
            atom_map[index]: charges[index]
            for index in range(len(charges)) if index not in caps
        }

        # The atom map indexes the topology the active site was extracted
        # from. Handing over a different one - most easily by running
        # prepare_protein after the extraction rather than before it - silently
        # writes the active site charges onto whatever atoms happen to hold those
        # indices, so the elements are checked before anything is modified.
        atoms = list(topology.atoms())
        labels = active_site['labels']
        mismatched = [
            index for index in covered
            if index >= len(atoms) or atoms[index].element is None
            or atoms[index].element.symbol != labels[next(
                c for c, t in atom_map.items() if t == index and c not in caps)]
        ]

        assert_msg_critical(
            not mismatched,
            'MetalSiteForceFieldBuilder.redistribute_backbone_charges: the '
            'atom map does not match this topology. Extract the active site '
            'from the same topology the system is built from; '
            'prepare_protein renumbers the atoms, so it must run first.')

        nonbonded = None
        for force in system.getForces():
            if isinstance(force, mm.NonbondedForce):
                nonbonded = force
                break

        assert_msg_critical(
            nonbonded is not None,
            'MetalSiteForceFieldBuilder.redistribute_backbone_charges: the '
            'system has no nonbonded force to write charges into')

        def get_charge(index):
            return nonbonded.getParticleParameters(index)[0].value_in_unit(
                mmunit.elementary_charge)

        atoms = list(topology.atoms())
        residue_indices = {atoms[index].residue.index for index in covered}
        region = [
            atom for residue in topology.residues()
            if residue.index in residue_indices for atom in residue.atoms()
        ]
        uncovered = [atom for atom in region if atom.index not in covered]

        total_before = sum(get_charge(atom.index) for atom in region)
        total_after = (sum(covered.values()) +
                       sum(get_charge(atom.index) for atom in uncovered))
        difference = total_before - total_after

        assert_msg_critical(
            len(uncovered) > 0 or abs(difference) < 1.0e-6,
            'MetalSiteForceFieldBuilder.redistribute_backbone_charges: the '
            f'coordination region is off by {difference:+.4f} e with no '
            'uncovered atom to absorb it')

        shift = difference / len(uncovered) if uncovered else 0.0

        for atom in region:
            parameters = nonbonded.getParticleParameters(atom.index)
            if atom.index in covered:
                new_charge = covered[atom.index]
            else:
                new_charge = get_charge(atom.index) + shift
            nonbonded.setParticleParameters(atom.index, new_charge,
                                            parameters[1], parameters[2])

        residual = total_before - sum(get_charge(atom.index) for atom in region)
        assert_msg_critical(
            abs(residual) < 1.0e-6,
            'MetalSiteForceFieldBuilder.redistribute_backbone_charges: the '
            f'charge moved by {residual:+.6f} e')

        self.ostream.print_blank()
        self.ostream.print_header('Charges written into the protein')
        self.ostream.print_header(31 * '-')
        self.ostream.print_header(
            self._param(
                'coordination region',
                f'{len(residue_indices)} residues, {len(region)} atoms'))
        self.ostream.print_header(
            self._param('covered by the active site', f'{len(covered)} atoms'))
        self.ostream.print_header(
            self._param('left to the protein', f'{len(uncovered)} atoms'))
        self.ostream.print_blank()
        self.ostream.print_header(
            self._param('region charge, protein', f'{total_before:+.4f} e'))
        self.ostream.print_header(
            self._param('region charge, fitted', f'{total_after:+.4f} e'))
        self.ostream.print_header(
            self._param('difference to recover', f'{difference:+.4f} e'))
        self.ostream.print_header(
            self._param('shift per uncovered atom', f'{shift:+.4f} e'))
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
        only the metal bonds and angles are transferred. The active site was built
        from the topology, so the atom map of the active site gives the
        correspondence directly and no graph matching is needed. Capping
        hydrogens are skipped, since they stand in for CA atoms that the
        protein force field parameterizes itself.

        :param topology:
            The protonated OpenMM topology of the whole enzyme.
        :param active_site:
            The active site, for the map back to the topology.
        :param forcefield:
            The force field generator carrying the fitted metal parameters.
        :param partial_charges:
            The charges fitted on the active site. When given they replace the
            charges of the coordination sphere through redistribute_charges;
            otherwise the protein force field's own charges are left alone.
        :param forcefield_files:
            The OpenMM force field files for the protein.

        :return:
            The tuple of the OpenMM system and the list of added terms.
        """

        assert_msg_critical(
            'openmm' in sys.modules,
            'MetalSiteForceFieldBuilder.create_enzyme_system: openmm is '
            'required')

        openmm_ff = mmapp.ForceField(*forcefield_files)
        system = openmm_ff.createSystem(topology,
                                        nonbondedMethod=mmapp.NoCutoff)

        if partial_charges is not None:
            self.redistribute_charges(system, topology, active_site,
                                      partial_charges)

        atom_map = active_site['atom_map']
        caps = set(active_site['cap_indices'])
        bonds, angles = self.get_metal_keys(forcefield, active_site)

        bond_force = None
        angle_force = None
        for force in system.getForces():
            if isinstance(force, mm.HarmonicBondForce):
                bond_force = force
            elif isinstance(force, mm.HarmonicAngleForce):
                angle_force = force

        assert_msg_critical(
            bond_force is not None and angle_force is not None,
            'MetalSiteForceFieldBuilder.create_enzyme_system: the protein '
            'system has no harmonic bond or angle force to extend')

        added = []

        for key in bonds:
            if caps & set(key):
                continue
            params = forcefield.bonds[key]
            bond_force.addBond(
                atom_map[key[0]], atom_map[key[1]],
                params['equilibrium'] * mmunit.nanometer,
                params['force_constant'] * mmunit.kilojoule_per_mole /
                mmunit.nanometer**2)
            added.append(('bond', key))

        for key in angles:
            if caps & set(key):
                continue
            params = forcefield.angles[key]
            angle_force.addAngle(
                atom_map[key[0]], atom_map[key[1]], atom_map[key[2]],
                params['equilibrium'] * mmunit.degree,
                params['force_constant'] * mmunit.kilojoule_per_mole /
                mmunit.radian**2)
            added.append(('angle', key))

        self.ostream.print_info(
            f'Added {sum(1 for term in added if term[0] == "bond")} metal '
            f'bond(s) and {sum(1 for term in added if term[0] == "angle")} '
            'metal angle(s) to the enzyme system.')
        self.ostream.flush()

        self._save_intermediate(
            self.ENZYME_SYSTEM_FILE,
            lambda path: path.write_text(mm.XmlSerializer.serialize(system)))

        return system, added

    @staticmethod
    def _param(label, value, label_width=26, value_width=20):
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

    def _print_param_list(self, label, items, value_width=20):
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
            self.ostream.print_header(
                self._param(label if i == 0 else '', chunk))

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

        self.ostream.print_blank()
        self.ostream.print_header('Metal Site Force Field Builder')
        self.ostream.print_header(32 * '=')
        self.ostream.print_blank()

        self.ostream.print_header(self._param('structure',
                                              Path(structure).name))
        self.ostream.print_header(
            self._param('primary cutoff', f'{self.metal_bond_cutoff:.2f} A'))
        self.ostream.print_header(
            self._param('secondary cutoff', f'{self.report_cutoff:.2f} A'))
        self.ostream.print_header(self._param('basis set',
                                              self.basis_set_label))
        self.ostream.print_header(
            self._param('xc functional', self._functional_label()))
        self.ostream.print_header(
            self._param('SCF driver',
                        'given' if self.scf_drv is not None else 'default'))
        self.ostream.print_header(
            self._param(
                'MM pre-optimization', 'skipped' if self.optimized_geometry
                is not None else self.do_mm_optimization))
        self.ostream.print_header(
            self._param(
                'QM optimization', 'supplied' if self.optimized_geometry
                is not None else self.do_qm_optimization))
        self.ostream.print_header(
            self._param('Hessian',
                        'supplied' if self.hessian is not None else 'computed'))
        self.ostream.print_header(
            self._param(
                'RESP charges', 'supplied'
                if self.partial_charges is not None else self.do_resp))
        self.ostream.print_header(
            self._param(
                'constrained atoms', 'beta carbons + caps'
                if self.constrain_capping_hydrogens else 'beta carbons'))
        self.ostream.print_header(
            self._param('enzyme system', self.build_enzyme_system))
        self.ostream.print_header(
            self._param('output folder',
                        self.folder if self.save_output else 'none'))
        self.ostream.print_blank()
        self.ostream.flush()

    def _print_binding_modes(self, binding_modes):
        """
        Prints the detected coordination sphere.
        """

        self.ostream.print_blank()
        self.ostream.print_header('Coordination sphere')
        self.ostream.print_header(19 * '-')

        for metal in binding_modes['metals']:
            self.ostream.print_header(
                self._param(
                    f'metal {metal["element"]} (index {metal["index"]})',
                    f'charge {metal["formal_charge"]:+d}'))

        self.ostream.print_blank()
        valstr = '{:>10} {:>9} | {:>18} | {:>16}'.format(
            'residue', 'atoms', 'distances (A)', 'mode')
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
                modes = '/'.join(dict.fromkeys(ligand['mode']
                                               for ligand in row))
                valstr = '{:>10} {:>9} | {:>18} | {:>16}'.format(
                    row[0]['residue'], atoms, distances, modes)
                self.ostream.print_header(valstr)

        for note in binding_modes['notes']:
            self.ostream.print_warning(note)

        self.ostream.print_blank()
        self.ostream.flush()

    def _print_binding_mode_update(self,
                                   changes,
                                   largest_shift,
                                   dropped_residues=None):
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
        self.ostream.print_header(
            self._param('largest bond change', f'{largest_shift:.2f} A'))

        if not changes:
            self.ostream.print_header(self._param('coordination', 'unchanged'))
            self.ostream.print_blank()
            self.ostream.print_info(
                'The new geometry gives the same coordination sphere; the '
                'binding modes and the connectivity matrix are kept as they '
                'are.')
            self.ostream.print_blank()
            self.ostream.flush()
            return

        self.ostream.print_header(self._param('contacts changed', len(changes)))
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

        for residue in dropped_residues or []:
            self.ostream.print_warning(
                f'{residue} no longer coordinates a metal, but it is still '
                'part of the truncated active site; extract the active site '
                'again to leave it out')

        self.ostream.print_blank()
        self.ostream.flush()

    def _print_active_site(self, active_site, binding_modes,
                           connectivity_matrix):
        """
        Prints the composition of the truncated active site.
        """

        molecule = active_site['molecule']

        self.ostream.print_blank()
        self.ostream.print_header('Truncated active site')
        self.ostream.print_header(21 * '-')
        self.ostream.print_header(
            self._param('atoms', molecule.number_of_atoms()))
        self.ostream.print_header(
            self._param('charge', f'{active_site["charge"]:+d}'))
        self.ostream.print_header(
            self._param('multiplicity', active_site['multiplicity']))
        self.ostream.print_header(
            self._param('capping hydrogens', len(active_site['cap_indices'])))
        self.ostream.print_header(
            self._param('bonds', int(connectivity_matrix.sum() // 2)))
        self._print_param_list('residues', active_site['residues'])

        variants = sorted(binding_modes['variants'].values())
        self._print_param_list('protonation', variants)

        self.ostream.print_blank()
        self.ostream.flush()

    def _print_partial_charges(self, topology, active_site, partial_charges):
        """
        Prints the fitted charges and what the capping correction did to them.

        :param topology:
            The protonated topology, for the residue each active site atom belongs
            to.
        :param active_site:
            The active site.
        :param partial_charges:
            The charges as fitted, capping hydrogens included.
        """

        charges = np.asarray(partial_charges)
        caps = sorted(active_site['cap_indices'])
        metals = active_site['metal_indices']
        labels = active_site['labels']
        n_atoms = len(charges)
        rest = [index for index in range(n_atoms) if index not in caps]
        cap_charge = float(sum(charges[index] for index in caps))

        self.ostream.print_blank()
        self.ostream.print_header('Partial charges')
        self.ostream.print_header(15 * '-')
        self.ostream.print_header(
            self._param('active site charge', f'{active_site["charge"]:+d}'))
        self.ostream.print_header(
            self._param('fitted total', f'{charges.sum():+.4f} e'))
        self.ostream.print_header(
            self._param('on capping hydrogens', f'{cap_charge:+.4f} e'))
        self.ostream.print_header(
            self._param(f'spread over {len(rest)} atoms',
                        f'{cap_charge / len(rest):+.4f} e each'))
        self.ostream.print_blank()

        # group what is left by the residue each atom came from
        atoms = list(topology.atoms())
        by_residue = {}
        for index in rest:
            residue = atoms[active_site['atom_map'][index]].residue
            by_residue.setdefault(residue, []).append(index)

        corrected = self.redistribute_cap_charges(active_site, charges)

        valstr = '{:>16} {:>7} | {:>12}'.format('fragment', 'atoms', 'charge')
        self.ostream.print_header(valstr)
        self.ostream.print_header(45 * '-')

        for residue, indices in by_residue.items():
            total = sum(corrected[index] for index in indices)
            if len(indices) == 1 and indices[0] in metals:
                name = f'{labels[indices[0]]} (metal)'
            else:
                name = f'{residue.name}{residue.id}'
            valstr = '{:>16} {:>7} | {:>12.4f}'.format(name, len(indices),
                                                       total)
            self.ostream.print_header(valstr)

        self.ostream.print_blank()
        self.ostream.flush()

    def _print_mm_optimization(self, active_site, forcefield, relaxed,
                               frozen_indices):
        """
        Prints what the crude MM relaxation did to the coordination sphere.

        The metal-ligand distances before and against after are the point of
        the table: the pass is there to clean up contacts and hydrogens, and
        a coordination sphere that moved more than a few hundredths of an
        Angstrom is the sign that it did something else instead.
        """

        labels = active_site['labels']
        molecule = active_site['molecule']
        metals = sorted(active_site['metal_indices'])
        bonds, angles = self.get_metal_keys(forcefield, active_site)

        before = molecule.get_coordinates_in_angstrom()
        after = relaxed.get_coordinates_in_angstrom()
        shift = np.linalg.norm(after - before, axis=1)

        self.ostream.print_blank()
        self.ostream.print_header('Crude MM relaxation')
        self.ostream.print_header(19 * '-')

        self.ostream.print_header(
            self._param('frozen atoms', len(frozen_indices)))
        self.ostream.print_header(
            self._param(
                'metal centers',
                'frozen' if set(metals) <= set(frozen_indices) else 'free'))
        self.ostream.print_header(
            self._param(
                'iteration limit', self.mm_max_iterations
                if self.mm_max_iterations > 0 else 'convergence'))
        self.ostream.print_header(
            self._param(
                'metal bonds', f'{len(bonds)}, k = '
                f'{self.default_metal_bond_force_constant:.0f}'))
        self.ostream.print_header(
            self._param(
                'metal angles', f'{len(angles)}, k = '
                f'{self.default_metal_angle_force_constant:.0f}'
                if self.reparameterize_metal_angles else 'left untouched'))
        self.ostream.print_header(
            self._param('bond equilibria',
                        'given' if self.metal_bond_equilibria else 'measured'))
        if self.reparameterize_metal_angles:
            self.ostream.print_header(
                self._param(
                    'angle equilibria',
                    'given' if self.metal_angle_equilibria else 'measured'))
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
                          '{:>+8.2f}').format(str(pair), names, was, now,
                                              now - was)
                self.ostream.print_header(valstr)

        self.ostream.print_blank()

        largest = int(np.argmax(shift))
        self.ostream.print_header(
            self._param(
                'largest shift', f'{shift[largest]:.2f} A on '
                f'{labels[largest]} {largest}'))
        self.ostream.print_header(
            self._param('mean shift', f'{shift.mean():.2f} A'))

        # A ligand swinging around its metal moves a long way in Cartesian
        # terms while the coordination sphere itself is untouched, so what
        # the pass has to be held to is the bond lengths, not the shifts.
        if worst_bond is not None:
            self.ostream.print_header(
                self._param('largest bond change',
                            f'{worst_bond[1]:+.2f} A on {worst_bond[0]}'))
            if abs(worst_bond[1]) > self.mm_bond_change_warning:
                self.ostream.print_warning(
                    f'The crude relaxation changed a {worst_bond[0]} bond by '
                    f'{worst_bond[1]:+.2f} A. Check the metal terms it was '
                    'given before trusting the geometry it produced.')

        self.ostream.print_blank()
        self.ostream.flush()

    def _print_metal_parameters(self, active_site, forcefield):
        """
        Prints the fitted metal bonds and angles.
        """

        labels = active_site['labels']
        metals = set(active_site['metal_indices'])
        coords = active_site['molecule'].get_coordinates_in_angstrom()
        bonds, angles = self.get_metal_keys(forcefield, active_site)

        self.ostream.print_blank()
        self.ostream.print_header('Metal bonds')
        self.ostream.print_header(11 * '-')
        valstr = '{:>12} {:>7} | {:>9} | {:>21}'.format('atoms', 'elements',
                                                        'r0 (A)',
                                                        'k (kcal/mol/A^2)')
        self.ostream.print_header(valstr)
        self.ostream.print_header(60 * '-')

        for key in bonds:
            params = forcefield.bonds[key]
            names = '-'.join(labels[index] for index in key)
            # kJ/mol/nm^2 to kcal/mol/A^2
            force_constant = params['force_constant'] / 100.0 / 4.184
            valstr = '{:>12} {:>7} | {:>9.3f} | {:>21.1f}'.format(
                str(key), names, params['equilibrium'] * 10.0, force_constant)
            self.ostream.print_header(valstr)

        self.ostream.print_blank()
        self.ostream.print_header('Metal angles')
        self.ostream.print_header(12 * '-')
        valstr = '{:>14} {:>9} | {:>12} | {:>19}'.format(
            'atoms', 'elements', 'theta0 (deg)', 'k (kJ/mol/rad^2)')
        self.ostream.print_header(valstr)
        self.ostream.print_header(60 * '-')

        for key in angles:
            params = forcefield.angles[key]
            names = '-'.join(labels[index] for index in key)
            bridging = key[0] in metals and key[2] in metals
            # the marker gets a fixed-width field of its own, otherwise the
            # centering of print_header would shift the marked line
            valstr = '{:>14} {:>9} | {:>12.1f} | {:>19.1f} {:<9}'.format(
                str(key), names, params['equilibrium'],
                params['force_constant'], 'bridging' if bridging else '')
            self.ostream.print_header(valstr)

        self.ostream.print_blank()

        for metal_a in sorted(metals):
            for metal_b in sorted(metals):
                if metal_a >= metal_b:
                    continue
                distance = np.linalg.norm(coords[metal_a] - coords[metal_b])
                self.ostream.print_header(
                    self._param(f'{labels[metal_a]}-{labels[metal_b]} distance',
                                f'{distance:.3f} A'))

        self.ostream.print_blank()
        self.ostream.flush()

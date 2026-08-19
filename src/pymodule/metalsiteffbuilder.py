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
import tempfile
import json
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
    coordinating residues, truncates a QM cluster at the CA-CB bonds, and fits
    the metal-ligand bond and angle parameters to a QM Hessian with the
    Seminario method. The fitted parameters can then be injected into a
    protein force field for the whole enzyme.

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
        - optimized_geometry: A molecule or xyz file to use instead of running
          the constrained optimization.
        - hessian: A Hessian matrix or text file to use instead of computing
          one.
        - partial_charges: Partial charges or a text file to use instead of
          computing RESP charges.
        - do_qm_optimization: The flag for optimizing the cluster before the
          Hessian is computed.
        - do_resp: The flag for computing RESP charges.
        - constrain_capping_hydrogens: The flag for constraining the capping
          hydrogens in addition to the beta carbons.
        - average_metal_terms: The flag for averaging the fitted metal terms
          over equivalent atoms.
        - binding_modes: The coordination topology of the metal centers.
        - active_site: The truncated QM cluster and its mapping back to the
          structure.
        - connectivity_matrix: The connectivity matrix of the cluster.
        - forcefield: The force field generator carrying the fitted terms.
        - results: The results of the last complete run.
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
    SUPPORTED_METAL_ELEMENTS = ('Zn',)

    # Formal charges assumed for bare metal ions. Only used for the cluster
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

    # Atoms dropped when truncating a sidechain at the CA-CB bond. CA itself
    # is replaced by a capping hydrogen rather than dropped.
    BACKBONE_ATOM_NAMES = ('N', 'C', 'O', 'OXT', 'H', 'H2', 'H3', 'HA', 'HA2',
                           'HA3', 'HXT')

    # Net charge of each protonation variant, for the cluster charge
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

    # Literature equilibrium metal-ligand distances in nm, used for the crude
    # pre-QM pass. Deliberately paired with weak force constants: at that
    # stage the equilibrium geometry matters far more than the stiffness.
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

        # truncation and protonation
        self.cap_bond_length = 1.09
        self.protonation_overrides = None

        # QM settings. Anything beyond the functional and the basis set is
        # set on an SCF driver assigned to scf_drv, which is used as given.
        self.scf_drv = None
        self.xcfun = 'PBE0'
        self.basis_set_label = 'def2-svp'
        self.mute_scf = True

        # Precomputed input. Whatever is provided here is validated against
        # the extracted cluster and the corresponding step is skipped.
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

        # results
        self.binding_modes = None
        self.active_site = None
        self.connectivity_matrix = None
        self.forcefield = None
        self.results = None

    def compute(self, structure):
        """
        Runs the complete pipeline on a structure.

        Detects the coordination sphere, protonates the coordinating
        residues, truncates the QM cluster, optionally optimizes it with the
        beta carbons frozen, computes a pair-restricted Hessian and RESP
        charges, and fits the metal-ligand terms.

        A geometry, Hessian or set of partial charges supplied through
        optimized_geometry, hessian or partial_charges is validated against
        the extracted cluster and used in place of the step that would have
        produced it.

        :param structure:
            The path to a PDB or mmCIF file.

        :return:
            The results dictionary, also stored in self.results.
        """

        assert_msg_critical('openmm' in sys.modules,
                            'MetalSiteForceFieldBuilder: openmm is required')

        self._print_header(structure)

        topology, positions = self.load_structure(structure)
        self.suggest_binding_modes(topology, positions)
        self._print_binding_modes()

        topology, positions = self.protonate(topology, positions)
        self.extract_active_site(topology, positions)
        self.build_connectivity(topology)
        self._print_active_site()

        # anything supplied through the settings is validated against the
        # cluster and replaces the step that would have produced it
        geometry = self._resolve_optimized_geometry()
        hessian = self._resolve_hessian()
        charges = self._resolve_partial_charges()

        if geometry is not None:
            self.active_site['molecule'] = geometry
            self.ostream.print_info(
                'Using the supplied optimized geometry; skipping the '
                'constrained optimization.')
            self.ostream.flush()
        elif self.do_qm_optimization:
            self.optimize_active_site()

        molecule = self.active_site['molecule']

        atom_pairs, atoms = self.extract_pairs(
            self.connectivity_matrix,
            self.active_site['metal_indices'],
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
            hessian = self.compute_hessian(atom_pairs)

        if charges is not None:
            self.ostream.print_info(
                'Using the supplied partial charges; skipping the RESP '
                'calculation.')
            self.ostream.flush()
        elif self.do_resp:
            charges = self.compute_resp_charges()

        self.build_forcefield(hessian, charges)
        self._print_metal_parameters()

        self.results = {
            'topology': topology,
            'positions': positions,
            'binding_modes': self.binding_modes,
            'active_site': self.active_site,
            'connectivity_matrix': self.connectivity_matrix,
            'atom_pairs': atom_pairs,
            'hessian': hessian,
            'partial_charges': charges,
            'forcefield': self.forcefield,
        }

        return self.results

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

    def suggest_binding_modes(self, topology, positions):
        """
        Proposes the coordination topology of the metal centers from geometry.

        Unrelaxed designed structures are unreliable enough that the intended
        ligand assignment has to win over raw distances, so the dictionary can
        be written to and read back from JSON and is meant to be reviewed before
        use. Everything downstream reads the dictionary, never the geometry.

        :param topology:
            The OpenMM topology.
        :param positions:
            The positions as an (N, 3) numpy array in Angstrom.

        :return:
            The binding modes dictionary, also stored in self.binding_modes.
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
                'index': atom.index,
                'element': symbol,
                'chain': atom.residue.chain.id,
                'resid': atom.residue.id,
                'res_index': atom.residue.index,
                'formal_charge': self.metal_formal_charges.get(symbol, 2),
            })

        assert_msg_critical(
            len(metals) > 0,
            'MetalSiteForceFieldBuilder.suggest_binding_modes: no metal atom '
            f'found. Recognized elements: {self.metal_elements}')

        self._check_supported_metals(metals, 'suggest_binding_modes')

        metal_indices = [metal['index'] for metal in metals]
        notes = []
        contacts = []

        for atom in topology.atoms():
            if atom.element is None or atom.index in metal_indices:
                continue
            if atom.element.symbol not in self.DONOR_ELEMENTS:
                continue

            distances = {
                metal['index']: float(
                    np.linalg.norm(positions[atom.index] -
                                   positions[metal['index']])
                ) for metal in metals
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
                'residue': label,
                'res_name': atom.residue.name,
                'res_index': atom.residue.index,
                'chain': atom.residue.chain.id,
                'atom': atom.name,
                'index': atom.index,
                'metals': bonded_to,
                'distances': [round(distances[i], 3) for i in bonded_to],
            })

            if not bonded_to:
                notes.append(
                    f'{label} {atom.name} at {closest:.2f} A is between the '
                    f'primary ({self.metal_bond_cutoff}) and secondary '
                    f'({self.report_cutoff}) cutoffs; review whether it '
                    'should be a ligand')

        ligands = [contact for contact in contacts if contact['metals']]

        self._assign_binding_modes(ligands, notes)

        for metal in metals:
            n_ligands = sum(
                1 for ligand in ligands if metal['index'] in ligand['metals'])
            if n_ligands < 3:
                notes.append(
                    f'metal {metal["element"]} (index {metal["index"]}) has '
                    f'only {n_ligands} ligand(s) within the primary cutoff; '
                    'check for a missing bridging ligand or water')

        self.binding_modes = {
            'metals': metals,
            'ligands': ligands,
            'variants': {},
            'cutoffs': {
                'primary': self.metal_bond_cutoff,
                'secondary': self.report_cutoff,
            },
            'notes': notes,
        }

        return self.binding_modes

    @staticmethod
    def _assign_binding_modes(ligands, notes):
        """
        Classifies each ligand contact using the residue context.

        :param ligands:
            The list of ligand contacts. Updated in place with a 'mode' entry.
        :param notes:
            The list of review notes. Appended to in place.
        """

        by_residue = {}
        for ligand in ligands:
            by_residue.setdefault(ligand['res_index'], []).append(ligand)

        for group in by_residue.values():
            res_name = group[0]['res_name']

            for ligand in group:
                if len(ligand['metals']) >= 2:
                    ligand['mode'] = 'bridging_single'
                else:
                    ligand['mode'] = 'monodentate'

            monodentate = all(
                ligand['mode'] == 'monodentate' for ligand in group)

            if len(group) >= 2 and monodentate:
                metals_hit = {ligand['metals'][0] for ligand in group}
                if len(metals_hit) >= 2:
                    if res_name in ('ASP', 'GLU', 'ASH', 'GLH'):
                        mode = 'bridging_mu13'
                    else:
                        mode = 'bridging'
                else:
                    mode = 'bidentate'
                for ligand in group:
                    ligand['mode'] = mode

            if res_name.startswith('HI') and len(group) >= 2:
                notes.append(
                    f'{group[0]["residue"]} appears to coordinate through '
                    'more than one ring nitrogen; imidazole geometry makes '
                    'this essentially impossible, so treat it as an artifact '
                    'and edit the binding modes')

    def save_binding_modes(self, filename):
        """
        Writes the binding modes to a JSON file for review.

        :param filename:
            The name of the JSON file.
        """

        assert_msg_critical(
            self.binding_modes is not None,
            'MetalSiteForceFieldBuilder.save_binding_modes: no binding modes '
            'available. Run suggest_binding_modes first.')

        if self.rank == mpi_master():
            Path(filename).write_text(json.dumps(self.binding_modes, indent=2))

    def load_binding_modes(self, filename):
        """
        Reads reviewed binding modes back from a JSON file.

        :param filename:
            The name of the JSON file.

        :return:
            The binding modes dictionary, also stored in self.binding_modes.
        """

        binding_modes = json.loads(Path(filename).read_text())
        # JSON keys are always strings
        binding_modes['variants'] = {
            int(key): value
            for key, value in binding_modes.get('variants', {}).items()
        }
        self._check_supported_metals(binding_modes.get('metals', []),
                                     'load_binding_modes')

        self.binding_modes = binding_modes

        return self.binding_modes

    def suggest_variants(self, topology):
        """
        Chooses a protonation variant for each coordinating residue.

        Automated pKa predictors are not parameterized on metal-coordinating
        residues, so the rules used here are deterministic: coordinating
        carboxylates are deprotonated, a coordinating cysteine is a thiolate,
        and the histidine tautomer is set so that the coordinating nitrogen
        carries no hydrogen. Entries in self.protonation_overrides win over
        the rules.

        :param topology:
            The OpenMM topology.

        :return:
            The dictionary mapping residue index to variant name.
        """

        assert_msg_critical(
            self.binding_modes is not None,
            'MetalSiteForceFieldBuilder.suggest_variants: no binding modes '
            'available. Run suggest_binding_modes first.')

        residues = list(topology.residues())

        by_residue = {}
        for ligand in self.binding_modes['ligands']:
            by_residue.setdefault(ligand['res_index'], []).append(ligand)

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
                closest = min(group, key=lambda x: min(x['distances']))
                if closest['atom'] == 'ND1':
                    variant = 'HIE'
                elif closest['atom'] == 'NE2':
                    variant = 'HID'
                else:
                    variant = 'HID'
                    self.binding_modes['notes'].append(
                        f'{closest["residue"]} coordinates through '
                        f'{closest["atom"]}, which is not a ring nitrogen; '
                        'defaulting to HID, please check')
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

        return variants

    def protonate(self, topology, positions):
        """
        Adds hydrogens with the protonation variants that the metal site
        requires.

        Adding hydrogens renumbers the atoms, so the indices in
        self.binding_modes are remapped to the new topology and the applied
        variants are recorded there.

        :param topology:
            The OpenMM topology.
        :param positions:
            The positions as an (N, 3) numpy array in Angstrom.

        :return:
            The tuple of the protonated topology and positions in Angstrom.
        """

        assert_msg_critical(
            'openmm' in sys.modules,
            'MetalSiteForceFieldBuilder.protonate: openmm is required')

        variants_by_index = self.suggest_variants(topology)

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
        for ligand in self.binding_modes['ligands']:
            ligand['index'] = remap(ligand['index'])
            ligand['metals'] = [remap(index) for index in ligand['metals']]

        for metal in self.binding_modes['metals']:
            metal['index'] = remap(metal['index'])

        self.binding_modes['variants'] = variants_by_index

        return new_topology, new_positions

    def extract_active_site(self, topology, positions):
        """
        Builds the truncated QM cluster.

        Sidechains are cut at the CA-CB bond and capped with a hydrogen placed
        along the CB to CA direction. No second-shell fragments and no
        backbone: the scheme is fixed, which also guarantees that the RESP and
        Hessian calculations see the same truncation.

        :param topology:
            The protonated OpenMM topology.
        :param positions:
            The positions as an (N, 3) numpy array in Angstrom.

        :return:
            The active site dictionary, also stored in self.active_site. It
            records the cluster indices of the capping hydrogens and of the
            beta carbons, the map back to the topology, and the charge.
        """

        assert_msg_critical(
            self.binding_modes is not None,
            'MetalSiteForceFieldBuilder.extract_active_site: no binding modes '
            'available. Run suggest_binding_modes first.')

        self._check_supported_metals(self.binding_modes['metals'],
                                     'extract_active_site')

        positions = np.asarray(positions)
        residues = list(topology.residues())

        res_indices = sorted(
            {ligand['res_index'] for ligand in self.binding_modes['ligands']})

        labels = []
        coords = []
        atom_map = {}
        cap_indices = []
        beta_carbon_indices = []
        metal_indices = []

        for metal in self.binding_modes['metals']:
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
                    # the cluster
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

        charge = sum(
            metal['formal_charge'] for metal in self.binding_modes['metals'])

        for res_index in res_indices:
            residue = residues[res_index]
            variant = self.binding_modes['variants'].get(res_index)
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

        self.active_site = {
            'molecule': molecule,
            'labels': labels,
            'atom_map': atom_map,
            'cap_indices': cap_indices,
            'beta_carbon_indices': beta_carbon_indices,
            'metal_indices': metal_indices,
            'charge': charge,
            'multiplicity': 1,
            'residues': [
                f'{residues[i].name}{residues[i].id}' for i in res_indices
            ],
            'res_indices': res_indices,
        }

        return self.active_site

    def build_connectivity(self, topology, warn_above=2.0):
        """
        Builds the connectivity matrix of the cluster without perceiving any
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
            The connectivity matrix, also stored in self.connectivity_matrix.
        """

        assert_msg_critical(
            self.active_site is not None,
            'MetalSiteForceFieldBuilder.build_connectivity: no active site '
            'available. Run extract_active_site first.')

        atom_map = self.active_site['atom_map']
        reverse_map = {
            top_index: cluster_index
            for cluster_index, top_index in atom_map.items()
        }

        n_atoms = len(atom_map)
        connectivity_matrix = np.zeros((n_atoms, n_atoms), dtype=int)

        for bond in topology.bonds():
            i, j = bond.atom1.index, bond.atom2.index
            if i in reverse_map and j in reverse_map:
                connectivity_matrix[reverse_map[i], reverse_map[j]] = 1
                connectivity_matrix[reverse_map[j], reverse_map[i]] = 1

        for ligand in self.binding_modes['ligands']:
            assert_msg_critical(
                ligand['index'] in reverse_map,
                'MetalSiteForceFieldBuilder.build_connectivity: ligand atom '
                f'{ligand["residue"]} {ligand["atom"]} is not part of the '
                'extracted cluster')
            lig_index = reverse_map[ligand['index']]
            for metal_index in ligand['metals']:
                metal = reverse_map[metal_index]
                connectivity_matrix[metal, lig_index] = 1
                connectivity_matrix[lig_index, metal] = 1

        # nothing but a metal bond should be long
        coords = self.active_site['molecule'].get_coordinates_in_angstrom()
        labels = self.active_site['labels']
        metals = set(self.active_site['metal_indices'])

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
        self.connectivity_matrix = connectivity_matrix

        return self.connectivity_matrix

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

    def _resolve_optimized_geometry(self):
        """
        Validates a geometry supplied through optimized_geometry.

        :return:
            The molecule, or None if nothing was supplied.
        """

        if self.optimized_geometry is None:
            return None

        if isinstance(self.optimized_geometry, Molecule):
            molecule = Molecule(self.optimized_geometry)
        else:
            path = Path(self.optimized_geometry)
            assert_msg_critical(
                path.is_file(),
                'MetalSiteForceFieldBuilder: optimized_geometry file '
                f'{path} not found')
            molecule = Molecule.read_xyz_file(str(path))

        cluster = self.active_site['molecule']

        assert_msg_critical(
            molecule.number_of_atoms() == cluster.number_of_atoms(),
            'MetalSiteForceFieldBuilder: optimized_geometry has '
            f'{molecule.number_of_atoms()} atoms but the extracted cluster '
            f'has {cluster.number_of_atoms()}')

        assert_msg_critical(
            list(molecule.get_labels()) == list(cluster.get_labels()),
            'MetalSiteForceFieldBuilder: the elements of '
            'optimized_geometry do not match the extracted cluster, so it '
            'describes a different structure')

        molecule.set_charge(cluster.get_charge())
        molecule.set_multiplicity(cluster.get_multiplicity())

        return molecule

    def _resolve_hessian(self):
        """
        Validates a Hessian supplied through the hessian setting.

        :return:
            The Hessian, or None if nothing was supplied.
        """

        if self.hessian is None:
            return None

        assert_msg_critical(
            not isinstance(self.hessian, str) or Path(self.hessian).is_file(),
            f'MetalSiteForceFieldBuilder: hessian file {self.hessian} not '
            'found')

        if isinstance(self.hessian, str):
            hessian = np.loadtxt(self.hessian)
        else:
            hessian = np.asarray(self.hessian)

        n_atoms = self.active_site['molecule'].number_of_atoms()
        expected = (3 * n_atoms, 3 * n_atoms)

        assert_msg_critical(
            hessian.shape == expected,
            f'MetalSiteForceFieldBuilder: hessian has shape {hessian.shape} '
            f'but the extracted cluster needs {expected}')

        return hessian

    def _resolve_partial_charges(self):
        """
        Validates charges supplied through the partial_charges setting.

        :return:
            The partial charges, or None if nothing was supplied.
        """

        if self.partial_charges is None:
            return None

        if isinstance(self.partial_charges, str):
            assert_msg_critical(
                Path(self.partial_charges).is_file(),
                'MetalSiteForceFieldBuilder: partial_charges file '
                f'{self.partial_charges} not found')
            charges = np.loadtxt(self.partial_charges)
        else:
            charges = np.asarray(self.partial_charges)

        n_atoms = self.active_site['molecule'].number_of_atoms()

        assert_msg_critical(
            charges.shape == (n_atoms,),
            f'MetalSiteForceFieldBuilder: partial_charges has shape '
            f'{charges.shape} but the extracted cluster has {n_atoms} atoms')

        total = float(np.sum(charges))
        expected = self.active_site['charge']

        if abs(total - expected) > 1.0e-3:
            self.ostream.print_warning(
                f'The supplied partial charges sum to {total:+.3f}, but the '
                f'cluster charge is {expected:+d}')
            self.ostream.flush()

        return charges

    def _get_scf_driver(self, molecule):
        """
        Returns the SCF driver and basis set for the cluster, without running
        anything.

        The driver assigned to scf_drv is used exactly as given, so it carries
        every QM setting beyond the functional and the basis set. Nothing is
        computed here: the gradient driver runs its own SCF at every geometry
        of an optimization, so an SCF at setup time would be thrown away.

        :param molecule:
            The cluster molecule.

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

    def _run_scf(self, molecule, basis):
        """
        Runs the SCF for a given geometry.

        Both the gradient and the Hessian driver reuse whatever is already in
        scf_drv.scf_results and only fall back to running an SCF when it is
        empty. They do not check that those results belong to the geometry
        they were handed, so the SCF has to be run explicitly here whenever
        the geometry may have moved since the driver was last used.

        :param molecule:
            The cluster molecule.
        :param basis:
            The basis set.

        :return:
            The SCF results.
        """

        if self.mute_scf:
            self.scf_drv.ostream.mute()

        scf_results = self.scf_drv.compute(molecule, basis)

        if self.mute_scf:
            self.scf_drv.ostream.unmute()

        assert_msg_critical(
            self.scf_drv.is_converged,
            'MetalSiteForceFieldBuilder: SCF did not converge in '
            f'{self.scf_drv.max_iter} iterations. Binuclear metal clusters '
            'converge slowly; raise max_iter or loosen conv_thresh on the '
            'SCF driver assigned to scf_drv.')

        return scf_results

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

    def constrained_indices(self):
        """
        Returns the cluster indices held fixed during optimization.

        The beta carbons are constrained by default: that is where the
        backbone actually holds the sidechain in place. The capping hydrogens
        merely stand in for the alpha carbons, and constraining them as well
        makes the truncated fragment rigid, which is available through
        constrain_capping_hydrogens.

        :return:
            The sorted list of cluster indices.
        """

        assert_msg_critical(
            self.active_site is not None,
            'MetalSiteForceFieldBuilder.constrained_indices: no active site '
            'available. Run extract_active_site first.')

        indices = set(self.active_site['beta_carbon_indices'])

        if self.constrain_capping_hydrogens:
            indices |= set(self.active_site['cap_indices'])

        return sorted(indices)

    def optimize_active_site(self, frozen_indices=None):
        """
        Optimizes the cluster with the beta carbons frozen.

        Freezing them keeps the spatial arrangement imposed by the protein
        backbone. Without it the site relaxes to a gas-phase geometry, and
        since the Seminario method takes the equilibrium values straight from
        the geometry, every fitted bond length and angle would then describe
        the wrong structure.

        :param frozen_indices:
            The zero-based cluster indices to freeze. Defaults to
            constrained_indices().

        :return:
            The optimized molecule, also stored in the active site.
        """

        assert_msg_critical(
            self.active_site is not None,
            'MetalSiteForceFieldBuilder.optimize_active_site: no active site '
            'available. Run extract_active_site first.')

        if frozen_indices is None:
            frozen_indices = self.constrained_indices()

        molecule = self.active_site['molecule']

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

        self.active_site['molecule'] = optimized

        return optimized

    def compute_hessian(self, atom_pairs=None):
        """
        Computes the nuclear Hessian of the cluster.

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

        assert_msg_critical(
            self.active_site is not None,
            'MetalSiteForceFieldBuilder.compute_hessian: no active site '
            'available. Run extract_active_site first.')

        molecule = self.active_site['molecule']

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
        self._run_scf(molecule, basis)

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

        return np.copy(hessian_drv.hessian)

    def compute_resp_charges(self):
        """
        Computes RESP charges for the cluster.

        :return:
            The partial charges as an (N,) numpy array.
        """

        assert_msg_critical(
            self.active_site is not None,
            'MetalSiteForceFieldBuilder.compute_resp_charges: no active site '
            'available. Run extract_active_site first.')

        molecule = self.active_site['molecule']

        self._print_muted_notice('the RESP charge fit at Hartree-Fock/6-31G*')

        resp_drv = RespChargesDriver(self.comm, self.ostream)

        if self.mute_scf:
            resp_drv.ostream.mute()

        # Neither a basis nor SCF results are passed: the driver then defaults
        # to Hartree-Fock with 6-31G*, which is what RESP charges are meant to
        # be fitted to, and runs its own SCF. Handing it the cluster's own
        # functional and basis would silently fit the charges at a level the
        # RESP parameters were never derived for.
        charges = resp_drv.compute(molecule)

        if self.mute_scf:
            resp_drv.ostream.unmute()

        charges = self.comm.bcast(charges, root=mpi_master())

        return np.array(charges)

    def get_metal_keys(self, forcefield=None):
        """
        Returns the bond and angle keys that involve a metal center.

        :param forcefield:
            The force field generator. Defaults to self.forcefield.

        :return:
            The tuple of the bond key list and the angle key list.
        """

        if forcefield is None:
            forcefield = self.forcefield

        assert_msg_critical(
            forcefield is not None,
            'MetalSiteForceFieldBuilder.get_metal_keys: no force field '
            'available. Run build_forcefield first.')

        metals = set(self.active_site['metal_indices'])

        bonds = [key for key in forcefield.bonds if metals & set(key)]
        angles = [key for key in forcefield.angles if metals & set(key)]

        return bonds, angles

    def build_forcefield(self, hessian=None, partial_charges=None):
        """
        Builds the cluster force field and fits the metal terms.

        Without a Hessian the metal terms are seeded from literature distances
        with deliberately weak force constants, which is the crude pre-QM
        pass: getting the equilibrium geometry roughly right matters far more
        than the stiffness at that stage.

        :param hessian:
            The Hessian as a (3N, 3N) numpy array. The string 'xtb' is
            rejected: it would make reparameterize re-optimize the cluster
            without constraints, silently replacing every equilibrium value
            with a gas-phase one.
        :param partial_charges:
            The partial charges. Defaults to zeros.

        :return:
            The force field generator, also stored in self.forcefield.
        """

        assert_msg_critical(
            self.active_site is not None,
            'MetalSiteForceFieldBuilder.build_forcefield: no active site '
            'available. Run extract_active_site first.')

        assert_msg_critical(
            self.connectivity_matrix is not None,
            'MetalSiteForceFieldBuilder.build_forcefield: no connectivity '
            'available. Run build_connectivity first.')

        assert_msg_critical(
            not isinstance(hessian, str),
            'MetalSiteForceFieldBuilder.build_forcefield: the Hessian must be '
            "a numpy array. Passing 'xtb' makes reparameterize re-optimize "
            'the cluster without constraints, which silently replaces every '
            'equilibrium value with a gas-phase one.')

        molecule = self.active_site['molecule']
        n_atoms = molecule.number_of_atoms()

        forcefield = MMForceFieldGenerator(self.comm, self.ostream)
        forcefield.connectivity_matrix = np.asarray(self.connectivity_matrix)
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
                partial_charges.shape == (n_atoms,),
                'MetalSiteForceFieldBuilder.build_forcefield: expected '
                f'{n_atoms} partial charges, got {partial_charges.shape}')
            forcefield.partial_charges = partial_charges
            for index in range(n_atoms):
                forcefield.atoms[index]['charge'] = partial_charges[index]

        self.forcefield = forcefield
        bonds, angles = self.get_metal_keys(forcefield)

        if hessian is None:
            self._seed_metal_terms(bonds, angles)
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
            self._check_force_constants(bonds, angles)

        return self.forcefield

    def _seed_metal_terms(self, bonds, angles):
        """
        Seeds the metal terms from literature values for the crude pre-QM
        pass.

        :param bonds:
            The metal bond keys.
        :param angles:
            The metal angle keys.
        """

        labels = self.active_site['labels']

        for key in bonds:
            elements = (labels[key[0]], labels[key[1]])
            equilibrium = self.LITERATURE_METAL_BONDS.get(elements)
            if equilibrium is None:
                equilibrium = self.LITERATURE_METAL_BONDS.get(elements[::-1])
            if equilibrium is None:
                self.ostream.print_warning(
                    f'No literature distance for {elements[0]}-{elements[1]}, '
                    'keeping the guessed value')
                continue
            self.forcefield.bonds[key]['equilibrium'] = equilibrium
            self.forcefield.bonds[key]['force_constant'] = 100000.0
            self.forcefield.bonds[key]['comment'] = 'literature estimate'

        for key in angles:
            self.forcefield.angles[key]['equilibrium'] = 109.5
            self.forcefield.angles[key]['force_constant'] = 200.0
            self.forcefield.angles[key]['comment'] = 'literature estimate'

        self.ostream.flush()

    def _check_force_constants(self, bonds, angles):
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

        labels = self.active_site['labels']

        zero_bonds = [
            key for key in bonds
            if self.forcefield.bonds[key]['force_constant'] == 0.0
        ]
        zero_angles = [
            key for key in angles
            if self.forcefield.angles[key]['force_constant'] == 0.0
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

    def minimize_active_site(self, frozen_indices=None, max_iterations=0):
        """
        Minimizes the cluster with its own force field.

        The beta carbons are frozen by default. They are where the backbone
        holds the sidechains in place, so a free minimization of an isolated
        cluster lets the truncated fragments drift apart and says nothing
        about the metal site. Note that the force field only carries
        electrostatics if build_forcefield was given partial charges.

        :param frozen_indices:
            The cluster indices to hold fixed. Defaults to
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

        assert_msg_critical(
            self.forcefield is not None,
            'MetalSiteForceFieldBuilder.minimize_active_site: no force field '
            'available. Run build_forcefield first.')

        if frozen_indices is None:
            frozen_indices = self.constrained_indices()

        with tempfile.TemporaryDirectory() as temp_dir:
            stem = str(Path(temp_dir) / 'cluster')
            self.forcefield.write_openmm_files(stem)

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

    def create_enzyme_system(self,
                             topology,
                             forcefield_files=('amber14-all.xml',
                                               'amber14/tip3pfb.xml')):
        """
        Injects the fitted metal terms into a force field system for the whole
        enzyme.

        The protein force field already covers everything except the metal, so
        only the metal bonds and angles are transferred. The cluster was built
        from the topology, so the atom map of the active site gives the
        correspondence directly and no graph matching is needed. Capping
        hydrogens are skipped, since they stand in for CA atoms that the
        protein force field parameterizes itself.

        :param topology:
            The protonated OpenMM topology of the whole enzyme.
        :param forcefield_files:
            The OpenMM force field files for the protein.

        :return:
            The tuple of the OpenMM system and the list of added terms.
        """

        assert_msg_critical(
            'openmm' in sys.modules,
            'MetalSiteForceFieldBuilder.create_enzyme_system: openmm is '
            'required')

        assert_msg_critical(
            self.forcefield is not None,
            'MetalSiteForceFieldBuilder.create_enzyme_system: no force field '
            'available. Run build_forcefield first.')

        openmm_ff = mmapp.ForceField(*forcefield_files)
        system = openmm_ff.createSystem(topology,
                                        nonbondedMethod=mmapp.NoCutoff)

        atom_map = self.active_site['atom_map']
        caps = set(self.active_site['cap_indices'])
        bonds, angles = self.get_metal_keys()

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
            params = self.forcefield.bonds[key]
            bond_force.addBond(
                atom_map[key[0]], atom_map[key[1]],
                params['equilibrium'] * mmunit.nanometer,
                params['force_constant'] * mmunit.kilojoule_per_mole /
                mmunit.nanometer**2)
            added.append(('bond', key))

        for key in angles:
            if caps & set(key):
                continue
            params = self.forcefield.angles[key]
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
        self.ostream.print_blank()
        self.ostream.flush()

    def _print_binding_modes(self):
        """
        Prints the detected coordination sphere.
        """

        self.ostream.print_blank()
        self.ostream.print_header('Coordination sphere')
        self.ostream.print_header(19 * '-')

        for metal in self.binding_modes['metals']:
            self.ostream.print_header(
                self._param(
                    f'metal {metal["element"]} (index {metal["index"]})',
                    f'charge {metal["formal_charge"]:+d}'))

        self.ostream.print_blank()
        valstr = '{:>10} {:>5} | {:>18} | {:>16}'.format(
            'residue', 'atom', 'distances (A)', 'mode')
        self.ostream.print_header(valstr)
        self.ostream.print_header(60 * '-')

        for ligand in self.binding_modes['ligands']:
            distances = ', '.join(f'{d:.2f}' for d in ligand['distances'])
            valstr = '{:>10} {:>5} | {:>18} | {:>16}'.format(
                ligand['residue'], ligand['atom'], distances, ligand['mode'])
            self.ostream.print_header(valstr)

        for note in self.binding_modes['notes']:
            self.ostream.print_warning(note)

        self.ostream.print_blank()
        self.ostream.flush()

    def _print_active_site(self):
        """
        Prints the composition of the truncated cluster.
        """

        molecule = self.active_site['molecule']

        self.ostream.print_blank()
        self.ostream.print_header('Truncated active site')
        self.ostream.print_header(21 * '-')
        self.ostream.print_header(
            self._param('atoms', molecule.number_of_atoms()))
        self.ostream.print_header(
            self._param('charge', f'{self.active_site["charge"]:+d}'))
        self.ostream.print_header(
            self._param('multiplicity', self.active_site['multiplicity']))
        self.ostream.print_header(
            self._param('capping hydrogens',
                        len(self.active_site['cap_indices'])))
        self.ostream.print_header(
            self._param('bonds', int(self.connectivity_matrix.sum() // 2)))
        self._print_param_list('residues', self.active_site['residues'])

        variants = sorted(self.binding_modes['variants'].values())
        self._print_param_list('protonation', variants)

        self.ostream.print_blank()
        self.ostream.flush()

    def _print_metal_parameters(self):
        """
        Prints the fitted metal bonds and angles.
        """

        labels = self.active_site['labels']
        metals = set(self.active_site['metal_indices'])
        coords = self.active_site['molecule'].get_coordinates_in_angstrom()
        bonds, angles = self.get_metal_keys()

        self.ostream.print_blank()
        self.ostream.print_header('Metal bonds')
        self.ostream.print_header(11 * '-')
        valstr = '{:>12} {:>7} | {:>9} | {:>21}'.format('atoms', 'elements',
                                                        'r0 (A)',
                                                        'k (kcal/mol/A^2)')
        self.ostream.print_header(valstr)
        self.ostream.print_header(60 * '-')

        for key in bonds:
            params = self.forcefield.bonds[key]
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
            params = self.forcefield.angles[key]
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

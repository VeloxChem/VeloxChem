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
The functional core of the metal site force field builder.

Every function here takes what it uses and returns what it produces: no state
is kept between calls and nothing is read off an object.  MetalSiteForceFieldBuilder
orchestrates these and holds the intermediates; MetalForceFieldManager calls
the same functions directly.
"""

from mpi4py import MPI
from pathlib import Path
from copy import deepcopy
import numpy as np
import tempfile
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

# Elements recognized when scanning a structure for metal centers. A
# center that is found but not supported is reported rather than ignored.
METAL_ELEMENTS = ('Zn', 'Fe', 'Cu', 'Mg', 'Mn', 'Co', 'Ni', 'Ca', 'Cd')

# Elements the builder is validated for. The literature distances, the
# formal charges and the coordination rules have only been checked against
# zinc sites, so anything else is rejected instead of silently treated as
# if it behaved the same way.
SUPPORTED_METAL_ELEMENTS = ('Zn',)

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

# Residues whose sidechain ends in a carboxylate, whose two oxygens are
# interchangeable for as long as neither of them is bound.
CARBOXYLATE_RESIDUES = ('ASP', 'ASH', 'GLU', 'GLH')

# Atoms dropped when truncating a sidechain at the CA-CB bond. CA itself
# is replaced by a capping hydrogen rather than dropped.
BACKBONE_ATOM_NAMES = ('N', 'C', 'O', 'OXT', 'H', 'H2', 'H3', 'HA', 'HA2',
                       'HA3', 'HXT')

# Residues the CA-CB truncation cannot cut, and why. A residue that only
# coordinates a metal is always one of ASP/GLU/CYS/HIS and never lands
# here, but include_residue takes any residue at all, so the rule the
# truncation relies on is written down where it is enforced.
UNTRUNCATABLE_RESIDUES = {
    'GLY': 'has no CB to cut at',
    'PRO': 'has a sidechain that closes back onto the backbone nitrogen, '
           'so cutting at CA-CB leaves CD with a dangling valence and no cap',
}

# Net charge of each protonation variant, for the active site charge
# bookkeeping. Note that CYX here is Modeller's variant of that name - a
# cysteine with no HG, which for a metal-bound sidechain is a thiolate -
# and not Amber's disulfide-bridged CYX, which is neutral.
# Every residue the truncation can cut needs an entry, since
# include_residue can put any of them in the cluster.
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
    'ALA': 0,
    'GLY': 0,
    'ILE': 0,
    'LEU': 0,
    'PHE': 0,
    'PRO': 0,
    'VAL': 0,
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
CAP_COMMENT = 'capping hydrogen'

# Literature equilibrium metal-ligand distances in nm. The crude pre-QM
# pass measures its equilibrium values on the input geometry instead;
# assign this to metal_bond_equilibria to impose these in their place.
LITERATURE_METAL_BONDS = {
    ('Zn', 'N'): 0.205,
    ('Zn', 'O'): 0.200,
    ('Zn', 'S'): 0.230,
}

# Force constants the crude pre-QM pass gives every metal term, in
# kJ/mol/nm^2 and kJ/mol/rad^2. Nothing at that stage says anything about the
# stiffness of a metal term, so they are flat: only the equilibrium values
# carry information until the Hessian is fitted.
DEFAULT_METAL_BOND_FORCE_CONSTANT = 100000.0
DEFAULT_METAL_ANGLE_FORCE_CONSTANT = 200.0

# ----------------------------------------------------------------------
# utilities
# ----------------------------------------------------------------------


def _stream(ostream):
    """
    Returns the stream to report through.

    A core function reports only when it is given somewhere to report to, so a
    caller that wants the numbers and not the commentary simply leaves the
    argument out.

    :param ostream:
        The output stream, or None.

    :return:
        The given stream, or a silent one.
    """

    return OutputStream(None) if ostream is None else ostream


def _folder_file(name, folder=None):
    """
    Returns the path of an intermediate in the working folder, or None
    when it is not there.

    :param name:
        The file name, one of the class-level file name attributes.

    :return:
        The path, or None.
    """

    if folder is None:
        return None

    path = Path(folder) / name

    return path if path.is_file() else None


def broadcast_forcefield(forcefield, comm=None, ostream=None):
    """
    Hands a force field generator from the master rank to every other one.

    A generator cannot be pickled: it owns an output stream, and a stream
    around sys.stdout does not survive the crossing. The stream is
    therefore set aside for the trip and put back on both sides, which
    keeps everything the JSON on disk leaves out - the pairs, the
    connectivity matrix and the atom type tables.

    :param forcefield:
        The force field on the master rank, ignored elsewhere.
    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream to attach to the generator on arrival.

    :return:
        The force field, on every rank.
    """

    comm = MPI.COMM_WORLD if comm is None else comm

    if comm.Get_size() == 1:
        return forcefield

    stream = None
    if comm.Get_rank() == mpi_master():
        stream = forcefield.ostream
        forcefield.ostream = None

    forcefield = comm.bcast(forcefield, root=mpi_master())

    forcefield.ostream = stream if stream is not None else _stream(ostream)

    return forcefield


# ----------------------------------------------------------------------
# loading
# ----------------------------------------------------------------------


def load_and_prepare_protein(structure, prepare=True):
    """
    Reads a structure and repairs it for a protein force field.

    The two are one step because the repair renumbers atoms: every index
    downstream refers to the topology this returns, so a structure that is
    read and then prepared separately is a chance to extract from the
    wrong one.

    :param structure:
        The path to a .pdb, .cif or .pdbx file.
    :param prepare:
        Whether to run the repair. It is what gives a protein force field a
        topology it can template, so a run that ends in an enzyme system
        needs it; a run that only wants the metal site can skip it, and skip
        pdbfixer with it.

    :return:
        The tuple of the OpenMM topology and the positions as an (N, 3)
        numpy array in Angstrom.
    """

    topology, positions = _load_structure(structure)

    if prepare:
        topology, positions = _prepare_protein(topology, positions)

    return topology, positions


def _load_structure(structure):
    """
    Reads a PDB or mmCIF structure.

    :param structure:
        The path to a .pdb, .cif or .pdbx file.

    :return:
        The tuple of the OpenMM topology and the positions as an (N, 3)
        numpy array in Angstrom.
    """

    assert_msg_critical('openmm' in sys.modules,
                        'load_and_prepare_protein: openmm is '
                        'required')

    path = Path(structure)

    assert_msg_critical(path.is_file(), f'load_and_prepare_protein: {path} not '
                        'found')

    if path.suffix.lower() in ('.cif', '.pdbx', '.mmcif'):
        pdb = mmapp.PDBxFile(str(path))
    else:
        pdb = mmapp.PDBFile(str(path))

    positions = np.array(pdb.positions.value_in_unit(mmunit.angstrom))

    return pdb.topology, positions


def _prepare_protein(topology, positions):
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

    assert_msg_critical('pdbfixer' in sys.modules or 'PDBFixer' in globals(),
                        'prepare_protein: pdbfixer is required')

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


def _check_supported_metals(metals, method):
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


# ----------------------------------------------------------------------
# detection
# ----------------------------------------------------------------------


def site_request():
    """
    A new, empty record of what the user has decided about a metal site.

    This is the only thing worth keeping between derivations. The
    coordination itself is not: it is a function of the geometry and of
    this record, so suggest_binding_modes is called again whenever the
    answer could have changed rather than a copy being carried forward and
    kept in step.

    Everything in here is keyed on residue indices, residue labels and
    atom names. None of those is disturbed by adding hydrogens, which is
    what lets the record outlive the renumbering that invalidates every
    atom index around it.

    - manual_bonds: the metal bonds added or removed by hand, as records
      of (residue index, atom name, metal residue index, action).
    - variants: the protonation variant of each residue, by residue index,
      as the protonation actually built it.
    - coordinating_residues: residues that are ligands whatever their
      distance.
    - extra_residues / excluded_residues: residues put into or taken out
      of the truncated cluster by hand.

    :return:
        The empty request.
    """

    return {
        'manual_bonds': [],
        'variants': {},
        'coordinating_residues': [],
        'extra_residues': [],
        'excluded_residues': [],
    }


def suggest_binding_modes(topology,
                          positions,
                          coordinating_residues=None,
                          metal_elements=METAL_ELEMENTS,
                          metal_formal_charges=None,
                          ostream=None,
                          bidentate_asymmetry=BIDENTATE_ASYMMETRY,
                          metal_bond_cutoff=3.0,
                          report_cutoff=3.5,
                          request=None):
    """
    Derives the coordination topology of the metal centers from geometry.

    This is a query, not a record: it is called again whenever the answer
    could have changed -- after the hydrogens are placed, after a
    relaxation, after an edit -- and nothing keeps what it returned. The
    atom indices it writes therefore always belong to the topology it was
    handed, which is what makes them safe to use.

    What survives instead is the request: the metal bonds decided by hand
    and the protonation variants, both keyed on residue indices and atom
    names, neither of which adding hydrogens disturbs. They are replayed
    here, so a derivation reproduces every edit that was ever made.

    :param topology:
        The OpenMM topology.
    :param positions:
        The positions as an (N, 3) numpy array in Angstrom.
    :param coordinating_residues:
        Residues that must be ligands whatever their distance, each given
        as a residue id ('130' or 130) or as a residue label ('ASP130').
        Merged into the request's own list.
    :param metal_elements:
        The elements treated as metal centers. Note that detecting one the
        builder is not validated for is a hard failure, so a zinc structure
        that also holds a calcium or magnesium ion has to narrow this to
        ('Zn',) by hand.
    :param metal_formal_charges:
        The formal charges assumed for the metal ions. Defaults to
        METAL_FORMAL_CHARGES.
    :param request:
        The stored request -- manual metal bonds, protonation variants and
        the residue membership overrides. Replayed onto the detection.

    :return:
        The binding modes dictionary, freshly derived.
    """

    ostream = _stream(ostream)

    if metal_formal_charges is None:
        metal_formal_charges = METAL_FORMAL_CHARGES

    if request is None:
        request = site_request()

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
        len(metals) > 0, 'suggest_binding_modes: no metal atom '
        f'found. Recognized elements: {metal_elements}')

    _check_supported_metals(metals, 'suggest_binding_modes')

    forced = _resolve_residues(topology, coordinating_residues, ostream=ostream)
    forced |= set(request.get('coordinating_residues', []))

    atoms = list(topology.atoms())

    def position_of(index):
        return positions[index]

    notes = []
    ligands = _collect_ligands(atoms,
                               position_of,
                               metals,
                               notes,
                               forced,
                               bidentate_asymmetry=bidentate_asymmetry,
                               metal_bond_cutoff=metal_bond_cutoff,
                               report_cutoff=report_cutoff,
                               ostream=ostream)

    # A bond decided by hand is not something the distances can be asked
    # about again, so the records are replayed onto every derivation. This
    # is what makes deriving twice give the same answer.
    records = request.get('manual_bonds', [])
    if records:
        _apply_manual_bonds(ligands, records, metals, atoms, position_of, notes)
        _assign_binding_modes(ligands,
                              notes,
                              bidentate_asymmetry=bidentate_asymmetry,
                              protected=_manual_protected(records))

    binding_modes = {
        'metals': metals,
        'ligands': ligands,
        'variants': dict(request.get('variants', {})),
        'coordinating_residues': sorted(forced),
        # what the caller asked for on top of the coordination, rather
        # than the membership those requests work out to: which residues
        # the cluster holds is derived by active_site_residues, so a
        # re-detection cannot leave the two disagreeing
        'extra_residues': sorted(request.get('extra_residues', [])),
        'excluded_residues': sorted(request.get('excluded_residues', [])),
        'manual_bonds': deepcopy(records),
        'notes': notes,
    }

    return binding_modes


def _collect_ligands(atoms,
                     position_of,
                     metals,
                     notes,
                     forced=None,
                     bidentate_asymmetry=BIDENTATE_ASYMMETRY,
                     metal_bond_cutoff=3.0,
                     report_cutoff=3.5,
                     ostream=None):
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

        label = f'{atom.residue.name}{atom.residue.id}'

        if atom.name in BACKBONE_ATOM_NAMES:
            notes.append(
                f'backbone atom {label} {atom.name} is {closest:.2f} A '
                'from a metal; the truncation scheme is sidechain-only, '
                'so it is not treated as a ligand')
            continue

        bonded_to = sorted([
            index for index, dist in distances.items()
            if dist <= metal_bond_cutoff
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
                f'primary ({metal_bond_cutoff}) and secondary '
                f'({report_cutoff}) cutoffs; review whether it '
                'should be a ligand')

    _force_ligands(atoms,
                   position_of,
                   metal_indices,
                   contacts,
                   notes,
                   forced,
                   metal_bond_cutoff=metal_bond_cutoff,
                   report_cutoff=report_cutoff,
                   ostream=ostream)

    ligands = [contact for contact in contacts if contact['metals']]

    _assign_binding_modes(ligands,
                          notes,
                          bidentate_asymmetry=bidentate_asymmetry)

    for metal in metals:
        n_ligands = sum(
            1 for ligand in ligands if metal['index'] in ligand['metals'])
        if n_ligands < 3:
            notes.append(
                f'metal {metal["element"]} (index {metal["index"]}) has '
                f'only {n_ligands} ligand(s) within the primary cutoff; '
                'check for a missing bridging ligand or water')

    return ligands


def _resolve_residues(topology, coordinating_residues, ostream=None):
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

    ostream = _stream(ostream)

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
            if wanted in (str(residue.id), f'{residue.name}{residue.id}')
        ]

        assert_msg_critical(
            len(matched) > 0, 'coordinating_residues asked for '
            f'{request}, which is not a residue of this structure')

        for residue in matched:
            forced.add(residue.index)

        found = ', '.join(f'{residue.name}{residue.id} '
                          f'(chain {residue.chain.id})' for residue in matched)
        ostream.print_info(
            f'Including {found} in the coordination sphere by request.')

    ostream.flush()

    return forced


def _force_ligands(atoms,
                   position_of,
                   metal_indices,
                   contacts,
                   notes,
                   forced,
                   metal_bond_cutoff=3.0,
                   report_cutoff=3.5,
                   ostream=None):
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

    ostream = _stream(ostream)

    if not forced:
        return

    already = {
        contact['res_index'] for contact in contacts if contact['metals']
    }

    for res_index in sorted(forced):
        if res_index in already:
            continue

        donors = [
            atom for atom in atoms
            if atom.residue.index == res_index and atom.element is not None and
            atom.element.symbol in DONOR_ELEMENTS and atom.name not in
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
        label = f'{atom.residue.name}{atom.residue.id}'

        notes.append(
            f'{label} {atom.name} is {distance:.2f} A from a metal, '
            f'beyond the primary cutoff ({metal_bond_cutoff}), and '
            'was made a ligand because coordinating_residues asked for it')

        if distance > report_cutoff:
            ostream.print_warning(
                f'{label} {atom.name} is {distance:.2f} A from its metal, '
                'which is a long way for a bond. It is a ligand because '
                'coordinating_residues asked for it.')
            ostream.flush()

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


def _assign_binding_modes(ligands,
                          notes,
                          bidentate_asymmetry=None,
                          protected=None):
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

    if bidentate_asymmetry is None:
        bidentate_asymmetry = BIDENTATE_ASYMMETRY

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
        if (len(group) == 2 and res_name in carboxylates and
                all(ligand['mode'] == 'bidentate' for ligand in group)):
            near, far = sorted(group, key=lambda ligand: ligand['distances'][0])
            separation = far['distances'][0] - near['distances'][0]
            asked_for = any((ligand['res_index'], ligand['atom']) in protected
                            for ligand in group)

            if separation > bidentate_asymmetry and asked_for:
                notes.append(
                    f'{group[0]["residue"]} has {near["atom"]} at '
                    f'{near["distances"][0]:.2f} A and {far["atom"]} at '
                    f'{far["distances"][0]:.2f} A, a difference of '
                    f'{separation:.2f} A, which would normally read as '
                    f'monodentate through {near["atom"]}; both are kept '
                    'because a manual edit asked for them')
            elif separation > bidentate_asymmetry:
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


def _insert_contact(ligands, contact):
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


def _add_contact_metal(contact, metal_index, distance):
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


def _record_manual_bond(records,
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
        if (record['res_index'] == res_index and record['atom'] == atom_name and
                record['metal_res_index'] == metal_res_index):
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


def _manual_protected(records):
    """
    The contacts that the classification rules must not take away again.

    :param records:
        The manual bond records.

    :return:
        The set of (residue index, atom name) pairs that were asked for.
    """

    return {(record['res_index'], record['atom'])
            for record in records
            if record['action'] == 'add'}


def _merge_notes(binding_modes, notes):
    """
    Adds review notes without repeating the ones already there.

    Editing the coordination classifies it again, and most of what that
    produces is what the last pass produced, so the notes are merged by
    their text rather than appended to.

    :param binding_modes:
        The binding modes. Its notes are updated in place.
    :param notes:
        The notes to merge in.
    """

    existing = binding_modes.setdefault('notes', [])

    for note in notes:
        if note not in existing:
            existing.append(note)


def _apply_manual_bonds(ligands, records, metals, atoms, position_of, notes):
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
                'residue': f'{atom.residue.name}{atom.residue.id}',
                'res_name': atom.residue.name,
                'res_index': atom.residue.index,
                'chain': atom.residue.chain.id,
                'atom': atom.name,
                'index': atom.index,
                'metals': [],
                'distances': [],
            }
            _insert_contact(ligands, contact)

        _add_contact_metal(contact, metal_index, distance)
        notes.append(f'{contact["residue"]} {contact["atom"]} is bonded '
                     f'to {metal_label} at {distance:.2f} A because '
                     'add_metal_bond asked for it')


# ----------------------------------------------------------------------
# editing
# ----------------------------------------------------------------------


def add_metal_bond(request,
                   binding_modes,
                   topology,
                   positions,
                   resid,
                   metal,
                   atom=None,
                   chain=None,
                   bidentate_asymmetry=BIDENTATE_ASYMMETRY,
                   metal_bond_cutoff=3.0,
                   report_cutoff=3.5,
                   equilibrium=None,
                   ostream=None):
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

    ostream = _stream(ostream)

    positions = np.asarray(positions)

    metal_entry = _resolve_metal(binding_modes, metal)
    residue = _resolve_residue(topology, resid, chain)
    ligand_atom = _resolve_ligand_atom(residue, atom, metal_entry, positions,
                                       binding_modes['ligands'])

    label = f'{residue.name}{residue.id}'
    metal_label = f'{metal_entry["element"]} (index {metal_entry["index"]})'

    assert_msg_critical(
        ligand_atom.element is not None and
        ligand_atom.element.symbol in DONOR_ELEMENTS, 'add_metal_bond: '
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
        if ligand['index'] == ligand_atom.index and
        metal_entry['index'] in ligand['metals']
    ]

    if already:
        ostream.print_warning(
            f'{label} {ligand_atom.name} is already bonded to '
            f'{metal_label}; the request is left as it is')
        ostream.flush()
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
        ostream.print_warning(
            f'{label} {ligand_atom.name} is {distance:.2f} A from '
            f'{metal_label}, beyond even the secondary cutoff '
            f'({report_cutoff} A), which is a long way for a bond. '
            'It is a ligand because add_metal_bond asked for it.')
    elif distance > metal_bond_cutoff:
        ostream.print_warning(
            f'{label} {ligand_atom.name} is {distance:.2f} A from '
            f'{metal_label}, beyond the primary cutoff '
            f'({metal_bond_cutoff} A). It is a ligand because '
            'add_metal_bond asked for it.')

    if residue.name.startswith('HI') and any(
            ligand['index'] != ligand_atom.index
            for ligand in residue_contacts):
        ostream.print_warning(
            f'{label} would coordinate through more than one ring '
            'nitrogen; imidazole geometry makes that essentially '
            'impossible, so check the structure before running on it')

    if binding_modes.get('variants'):
        ostream.print_warning(
            'This structure has already been protonated. The edit changes '
            f'no hydrogen, so protonate again if the variant of {label} '
            'depends on this bond.')

    ostream.flush()

    assert_msg_critical(
        equilibrium is None or equilibrium > 0, 'add_metal_bond: the '
        f'equilibrium distance must be positive, got {equilibrium}')

    new_request = deepcopy(request)
    records = new_request.setdefault('manual_bonds', [])

    # the force field keeps its bond lengths in nanometers, while everything
    # a caller is shown here -- the cutoffs, the distances in these warnings
    # -- is in Angstrom, so the conversion happens once, at the boundary
    _record_manual_bond(records,
                        residue.index,
                        ligand_atom.name,
                        metal_entry['res_index'],
                        'add',
                        equilibrium=None
                        if equilibrium is None else 0.1 * float(equilibrium))

    ostream.print_info(f'Bonded {label} {ligand_atom.name} to {metal_label} at '
                       f'{distance:.2f} A.')

    if equilibrium is not None:
        ostream.print_info(
            f'The crude pass will pull it to {equilibrium:.2f} A rather than '
            f'holding it at the {distance:.2f} A the structure has. The QM '
            'optimization and the fit that follows are not bound by it.')

    ostream.flush()

    return new_request


def remove_metal_bond(request,
                      binding_modes,
                      resid,
                      metal=None,
                      atom=None,
                      chain=None,
                      bidentate_asymmetry=BIDENTATE_ASYMMETRY,
                      ostream=None):
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

    ostream = _stream(ostream)

    wanted = str(resid).strip()
    matched = [
        ligand for ligand in binding_modes['ligands']
        if wanted in (ligand['residue'],
                      ligand['residue'][len(ligand['res_name']):]) and
        (chain is None or str(ligand['chain']) == str(chain))
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
        {index for ligand in matched for index in ligand['metals']})

    if metal is None:
        assert_msg_critical(
            len(reached) == 1, 'remove_metal_bond: '
            f'{label} bridges the metals at indices {reached}; pass '
            'metal= to say which bond to remove')
        metal_entry = _resolve_metal(binding_modes, reached[0])
    else:
        metal_entry = _resolve_metal(binding_modes, metal)
        assert_msg_critical(
            metal_entry['index'] in reached, 'remove_metal_bond: '
            f'{label} is not bonded to the metal at index '
            f'{metal_entry["index"]}; it reaches {reached}')

    metal_index = metal_entry['index']
    metal_label = f'{metal_entry["element"]} (index {metal_index})'

    targets = [ligand for ligand in matched if metal_index in ligand['metals']]

    if atom is not None:
        index = _as_atom_index(atom)
        if index is None:
            name = str(atom).strip()
            targets = [ligand for ligand in targets if ligand['atom'] == name]
        else:
            targets = [ligand for ligand in targets if ligand['index'] == index]

        assert_msg_critical(
            len(targets) > 0, 'remove_metal_bond: '
            f'{label} {atom} is not bonded to {metal_label}; its bonds '
            'to it are ' + ', '.join(ligand['atom']
                                     for ligand in matched
                                     if metal_index in ligand['metals']))

    if binding_modes.get('variants'):
        ostream.print_warning(
            'This structure has already been protonated. The edit changes '
            f'no hydrogen, so protonate again if the variant of {label} '
            'depends on this bond.')
        ostream.flush()

    new_request = deepcopy(request)
    records = new_request.setdefault('manual_bonds', [])

    removed = []
    for ligand in targets:
        removed.append(ligand['atom'])
        _record_manual_bond(records, res_index, ligand['atom'],
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
        ostream.print_info(
            f'{label} no longer coordinates anything and was taken out of '
            'coordinating_residues, which would otherwise put the bond back on '
            'the next derivation.')

    ostream.print_info(
        f'Removed the bond(s) of {label} {", ".join(removed)} to '
        f'{metal_label}.')
    ostream.flush()

    return new_request


def _resolve_residue(topology, resid, chain=None):
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
        if wanted in (str(residue.id), f'{residue.name}{residue.id}') and
        (chain is None or str(residue.chain.id) == str(chain))
    ]

    in_chain = '' if chain is None else f' of chain {chain}'

    assert_msg_critical(
        len(matched) > 0, f'residue {resid}'
        f'{in_chain} is not a residue of this structure')

    found = ', '.join(f'{residue.name}{residue.id} '
                      f'(chain {residue.chain.id})' for residue in matched)

    assert_msg_critical(
        len(matched) == 1, f'residue {resid} matches {found}; '
        'pass chain= to say which one is meant')

    return matched[0]


def _resolve_metal(binding_modes, metal):
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

    wanted = _as_atom_index(metal)

    assert_msg_critical(
        wanted is not None and wanted in entries,
        f'there is no metal center at atom index {metal}. '
        f'This site holds {known}, which is what its labels show')

    return entries[wanted]


def _as_atom_index(atom):
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


def _sidechain_donors(residue):
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


def _resolve_ligand_atom(residue, atom, metal_entry, positions, ligands):
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
    donors = _sidechain_donors(residue)
    label = f'{residue.name}{residue.id}'
    names = [donor.name for donor in donors]

    if atom is not None:
        index = _as_atom_index(atom)

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

    if (residue.name in CARBOXYLATE_RESIDUES and not bound and
            len(oxygens) == 2):
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


# ----------------------------------------------------------------------
# protonation
# ----------------------------------------------------------------------


def _histidine_variant(residue, positions, metal_positions):
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
        atom.name: atom for atom in sidechain if atom.name in ('ND1', 'NE2')
    }

    if len(ring_nitrogens) < 2:
        return 'HID', (f'{residue.name}{residue.id} does not have both ring '
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
        note = (f'{residue.name}{residue.id} has {nearest.name} closer to a '
                f'metal ({closest_metal_distance(nearest):.2f} A) than either '
                f'ring nitrogen (ND1 {distances["ND1"]:.2f} A, NE2 '
                f'{distances["NE2"]:.2f} A); the tautomer was still chosen '
                'from the nitrogens, but the geometry looks distorted')

    return variant, note


def known_variants(res_name):
    """
    The protonation variants OpenMM will accept for a residue.

    Modeller keeps them in hydrogens.xml and refuses anything else from
    deep inside addHydrogens, with a bare ValueError. Asking it what it
    knows before it is called turns that into an error naming the legal
    set, and keeps the check and the thing it guards from drifting apart.

    :param res_name:
        The residue name, as the topology has it.

    :return:
        The variant names, as a tuple. Empty for a residue whose
        protonation Modeller does not offer a choice about.
    """

    assert_msg_critical('openmm' in sys.modules,
                        'known_variants: openmm is required')

    mmapp.Modeller._loadStandardHydrogenDefinitions()
    spec = mmapp.Modeller._residueHydrogens.get(res_name)

    if spec is None:
        return ()

    return tuple(spec.variants)


def check_variant(residue, variant):
    """
    Checks that a protonation variant can be asked for and paid for.

    Two things have to hold, and they come from different places: OpenMM
    has to know how to build the variant, and the charge bookkeeping of
    the active site has to know what it is worth. A variant that passes
    one and not the other fails much later and much less clearly - either
    inside addHydrogens or on the cluster charge - so both are checked
    here.

    :param residue:
        The residue the variant is for.
    :param variant:
        The variant name.
    """

    label = f'{residue.name}{residue.id}'
    legal = known_variants(residue.name)

    assert_msg_critical(
        len(legal) > 0, 'check_variant: '
        f'{label} has no protonation variants to choose between; OpenMM '
        'builds it one way only')

    assert_msg_critical(
        variant in legal, 'check_variant: '
        f'{variant} is not a protonation variant of {label}. The variants '
        f'OpenMM knows for {residue.name} are {list(legal)}')

    assert_msg_critical(
        variant in VARIANT_CHARGES, 'check_variant: no charge known for '
        f'variant {variant}, so the charge of the active site could not be '
        'counted. Add it to VARIANT_CHARGES')


def suggest_variants(topology,
                     positions,
                     binding_modes,
                     protonation_overrides=None):
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
            variant, note = _histidine_variant(residues[res_index], positions,
                                               metal_positions)
            if note is not None:
                notes.append(note)
        else:
            variant = None

        if variant is not None:
            variants[res_index] = variant

    # An override is matched against the residue id, the ASP130 label and,
    # for a caller that has a residue in hand, the residue index. Note that
    # ids and indices overlap almost completely -- a single chain numbered
    # from one has index i for id i+1 -- so a key that hits both is a
    # collision rather than two ways of saying the same thing, and naming it
    # beats protonating a residue nobody asked about.
    overrides = protonation_overrides or {}
    for key, variant in overrides.items():
        by_id = [
            residue for residue in residues if str(residue.id) == str(key) or
            f'{residue.name}{residue.id}' == str(key)
        ]
        matched = by_id or [
            residue for residue in residues if residue.index == key
        ]

        assert_msg_critical(
            len(matched) > 0, 'suggest_variants: override '
            f'residue {key} not found')

        named = ', '.join(f'{residue.name}{residue.id} '
                          f'(chain {residue.chain.id})' for residue in matched)

        assert_msg_critical(
            len(matched) == 1 or by_id, 'suggest_variants: override '
            f'{key} matches {named}; name the residue as a label such as '
            f'{matched[0].name}{matched[0].id} to say which is meant')

        for residue in matched:
            check_variant(residue, variant)
            variants[residue.index] = variant

    return variants, notes


def protonate(topology, positions, binding_modes, protonation_overrides=None):
    """
    Adds hydrogens with the protonation variants that the metal site
    requires.

    Adding hydrogens renumbers the atoms, which invalidates every atom
    index the coordination holds. Nothing is remapped: the coordination is
    derived from geometry, so the caller runs suggest_binding_modes again
    on what comes back and gets indices that are correct by construction.
    That is safe because addHydrogens only adds hydrogens -- they are not
    donors and they do not move a heavy atom -- so the same contacts are
    found either way.

    Residue indices do survive, which is why the variants are keyed by one.

    :param topology:
        The OpenMM topology.
    :param positions:
        The positions as an (N, 3) numpy array in Angstrom.
    :param binding_modes:
        The coordination, read for the variant rules. Not modified.

    :return:
        The tuple of the protonated topology, the positions in Angstrom,
        the residue index to variant mapping of what was actually built,
        and the notes the variant choice produced.
    """

    assert_msg_critical('openmm' in sys.modules,
                        'protonate: openmm is required')

    variants_by_index, notes = suggest_variants(
        topology,
        positions,
        binding_modes,
        protonation_overrides=protonation_overrides)

    variant_list = [None] * topology.getNumResidues()
    for res_index, variant in variants_by_index.items():
        variant_list[res_index] = variant

    modeller = mmapp.Modeller(topology, np.asarray(positions) * mmunit.angstrom)
    actual_variants = modeller.addHydrogens(variants=variant_list)

    # A residue put in the cluster without a metal bond is left to the pH
    # by suggest_variants, so what it ended up as is only known here.
    # Recording it means the charge count reads what was built rather than
    # falling back on the residue name and warning about it.
    # addHydrogens only names a variant it chose between - the cysteines and
    # the histidines - and reports None for a residue it left at the default
    # for the pH. That default is the variant the residue is named after:
    # every alternative in hydrogens.xml is gated behind a maxph the default
    # pH of 7 does not reach.
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


# ----------------------------------------------------------------------
# extraction
# ----------------------------------------------------------------------


def check_truncatable(residue):
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
        f'{residue.name}{residue.id} {reason}, and the truncation of this '
        'module cuts every sidechain at CA-CB. It cannot be part of the '
        'active site')


def active_site_residues(binding_modes):
    """
    The residues the truncated active site holds.

    A residue is in the cluster when it coordinates a metal, or when
    include_residue asked for it, and out of it when remove_residue said
    so. Only those two requests are stored; the membership itself is
    worked out here every time it is needed, so that a re-detection which
    gains or loses a contact cannot leave a stored set behind.

    :param binding_modes:
        The binding modes.

    :return:
        The residue indices, sorted.
    """

    residues = {ligand['res_index'] for ligand in binding_modes['ligands']}
    residues |= set(binding_modes.get('extra_residues', []))
    residues -= set(binding_modes.get('excluded_residues', []))

    return sorted(residues)


def extract_active_site(topology,
                        positions,
                        binding_modes,
                        cap_bond_length=1.09,
                        ostream=None):
    """
    Builds the truncated QM active site.

    Sidechains are cut at the CA-CB bond and capped with a hydrogen placed
    along the CB to CA direction. No second-shell fragments and no
    backbone: the scheme is fixed, which also guarantees that the RESP and
    Hessian calculations see the same truncation.

    The connectivity comes with it, under 'connectivity_matrix'. It is not
    a separate object: which atoms are bonded is a property of the site
    that was extracted, and nothing downstream should be able to pair the
    two up wrongly.

    It also comes with 'labels', one string per atom, holding what
    add_metal_bond and remove_metal_bond want to be told about that atom:
    the atom index for a metal, the resid on the beta carbon that stands
    for its residue, the atom name on any other heavy atom, and an empty
    string on every hydrogen. Passing it to Molecule.show as atom_labels
    therefore draws the coordination edit straight onto the structure.
    These are not the element labels, which are read off the molecule with
    get_labels().

    To draw the site, hand Molecule.show the bonds that
    connectivity_bonds(active_site['connectivity_matrix']) returns rather
    than letting it perceive them by distance. That matters here for two
    reasons: the metal-ligand bonds are a decision rather than a distance
    and nothing perceives them, and OpenMM places hydrogens at about 1.19
    Angstrom, past the C-H threshold of Molecule.get_connectivity_matrix.
    The bonds are derived where they are drawn rather than stored beside
    the matrix, so the two cannot drift apart.

    :param topology:
        The protonated OpenMM topology.
    :param positions:
        The positions as an (N, 3) numpy array in Angstrom.

    :return:
        The active site dictionary. It records the active site indices of
        the capping hydrogens and of the beta carbons, the map back to
        the topology, the per-atom labels, the bonds and the charge.
    """

    ostream = _stream(ostream)

    _check_supported_metals(binding_modes['metals'], 'extract_active_site')

    positions = np.asarray(positions)
    residues = list(topology.residues())

    res_indices = active_site_residues(binding_modes)

    # a ligand whose residue was excluded would be left with a metal bond
    # to an atom the cluster does not hold, which _build_connectivity
    # would only find out about several steps later
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
        check_truncatable(residues[res_index])

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
            if atom.name in BACKBONE_ATOM_NAMES:
                continue

            if atom.name == 'CA':
                cb_atom = None
                for other in res_atoms:
                    if other.name == 'CB':
                        cb_atom = other
                        break
                assert_msg_critical(
                    cb_atom is not None, 'extract_active_site: '
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
                              direction * cap_bond_length)
                labels.append('H')
                atom_labels.append('')
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
            ostream.print_warning(
                'No protonation variant recorded for '
                f'{residue.name}{residue.id}; using the residue name for '
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
        'residues': [
            f'{residues[i].name}{residues[i].id}' for i in res_indices
        ],
    }

    active_site['connectivity_matrix'] = _build_connectivity(topology,
                                                             active_site,
                                                             binding_modes,
                                                             ostream=ostream)

    return active_site


def connectivity_bonds(connectivity_matrix):
    """
    Reads a connectivity matrix as the list of bonds Molecule.show draws.

    The pairs are zero-indexed and plain ints, since they are handed to
    RDKit, which does not take numpy integers.

    :param connectivity_matrix:
        The connectivity of an active site.

    :return:
        The bonds, as index pairs.
    """

    matrix = np.asarray(connectivity_matrix)

    return [(int(i), int(j))
            for i, j in zip(*np.triu_indices_from(matrix, k=1))
            if matrix[i, j]]


def show_active_site(active_site, **kwargs):
    """
    Draws the extracted active site.

    Two things about a truncated metal site are not recoverable from the
    geometry, so they are handed to Molecule.show rather than left to it:

    - the bonds. Perceiving them by distance loses every metal-ligand bond,
      which is a decision rather than a distance and which nothing
      perceives, and invents C-H ones, since OpenMM places hydrogens at
      about 1.19 Angstrom -- past the threshold of
      Molecule.get_connectivity_matrix.
    - the atom labels. These are what the edit methods want to be told
      about an atom: the atom index on a metal, the resid on the beta
      carbon standing for its residue, the atom name on any other heavy
      atom, and nothing on a hydrogen. Drawing them puts the argument for
      add_metal_bond or remove_metal_bond straight onto the structure.

    Every other keyword goes to Molecule.show untouched, so width, height,
    forming_bonds and the rest work as they always do. Passing bonds or
    atom_labels explicitly overrides what is worked out here.

    :param active_site:
        The active site to draw.
    :param kwargs:
        Further keyword arguments for Molecule.show.

    :return:
        Whatever Molecule.show returns.
    """

    kwargs.setdefault('atom_labels', active_site['labels'])
    kwargs.setdefault(
        'bonds', connectivity_bonds(active_site['connectivity_matrix']))

    return active_site['molecule'].show(**kwargs)


def _build_connectivity(topology,
                        active_site,
                        binding_modes,
                        warn_above=2.0,
                        ostream=None):
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
    :param warn_above:
        The length in Angstrom above which a non-metal bond is reported.

    :return:
        The connectivity matrix.
    """

    ostream = _stream(ostream)

    atom_map = active_site['atom_map']
    reverse_map = {
        top_index: site_index for site_index, top_index in atom_map.items()
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
            distance = np.linalg.norm(coords[i] - coords[j])
            if distance > warn_above:
                ostream.print_warning(
                    f'Non-metal bond {i}-{j} ({labels[i]}-{labels[j]}) '
                    f'is {distance:.2f} A long')

    ostream.flush()

    return connectivity_matrix


def apply_metal_bonds(active_site, binding_modes):
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

    caps = set(active_site['cap_indices'])
    site_of = {
        top_index: site_index
        for site_index, top_index in active_site['atom_map'].items()
        if site_index not in caps
    }
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


def derive_site_coordination(topology,
                             geometry,
                             active_site,
                             binding_modes,
                             bidentate_asymmetry=BIDENTATE_ASYMMETRY,
                             metal_bond_cutoff=3.0,
                             report_cutoff=3.5,
                             ostream=None):
    """
    Works out the coordination of an extracted active site from its own
    geometry.

    Restricted to the atoms the cluster holds, which is the difference
    between this and suggest_binding_modes: once the site is extracted and
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

    # a capping hydrogen is mapped to the CA it replaces, so its position
    # must not be read back as that CA's
    cap_atoms = {
        atom_map[site_index] for site_index in active_site['cap_indices']
    }
    site_of = {
        top_index: site_index
        for site_index, top_index in atom_map.items()
        if top_index not in cap_atoms
    }

    atoms = list(topology.atoms())
    candidates = [atoms[top_index] for top_index in sorted(site_of)]

    def position_of(index):
        return coordinates[site_of[index]]

    notes = []
    ligands = _collect_ligands(candidates,
                               position_of,
                               binding_modes['metals'],
                               notes,
                               set(
                                   binding_modes.get('coordinating_residues',
                                                     [])),
                               bidentate_asymmetry=bidentate_asymmetry,
                               metal_bond_cutoff=metal_bond_cutoff,
                               report_cutoff=report_cutoff,
                               ostream=ostream)

    # a bond put there by hand is not something the distances can be
    # asked about again, so it is replayed onto every derivation
    records = binding_modes.get('manual_bonds', [])
    if records:
        _apply_manual_bonds(ligands, records, binding_modes['metals'],
                            candidates, position_of, notes)
        _assign_binding_modes(ligands,
                              notes,
                              bidentate_asymmetry=bidentate_asymmetry,
                              protected=_manual_protected(records))

    return ligands, notes


def update_binding_modes(topology,
                         geometry,
                         active_site,
                         binding_modes,
                         bidentate_asymmetry=BIDENTATE_ASYMMETRY,
                         metal_bond_cutoff=3.0,
                         report_cutoff=3.5,
                         ostream=None):
    """
    Re-detects the coordination sphere on a new active site geometry.

    Relaxing the active site moves the metal-ligand distances: a contact
    that started just inside the primary cutoff can end up outside it, an
    asymmetric carboxylate can open into a monodentate one, and a second
    oxygen can rotate onto a metal. The rules of suggest_binding_modes are
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

    :return:
        The tuple of the binding modes, the active site and the flag that
        is True when the coordination changed. The active site is a new
        dictionary carrying the new connectivity when it did change, and
        the one that was passed in when it did not.
    """

    _check_supported_metals(binding_modes['metals'], 'update_binding_modes')

    if isinstance(geometry, Molecule):
        coordinates = geometry.get_coordinates_in_angstrom()
    else:
        coordinates = np.asarray(geometry, dtype=float)

    atom_map = active_site['atom_map']

    cap_atoms = {
        atom_map[site_index] for site_index in active_site['cap_indices']
    }
    site_of = {
        top_index: site_index
        for site_index, top_index in atom_map.items()
        if top_index not in cap_atoms
    }

    records = binding_modes.get('manual_bonds', [])

    new_ligands, notes = derive_site_coordination(
        topology,
        coordinates,
        active_site,
        binding_modes,
        bidentate_asymmetry=bidentate_asymmetry,
        metal_bond_cutoff=metal_bond_cutoff,
        report_cutoff=report_cutoff,
        ostream=ostream)

    def coordination(ligands):
        return {
            ligand['index']: (ligand['mode'], tuple(ligand['metals']))
            for ligand in ligands
        }

    old_coordination = coordination(binding_modes['ligands'])
    new_coordination = coordination(new_ligands)

    old_ligands = {
        ligand['index']: ligand for ligand in binding_modes['ligands']
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
        _print_binding_mode_update([], largest_shift, ostream=ostream)
        return binding_modes, active_site, False

    changes = []
    for index in sorted(set(old_coordination) | set(new_coordination)):
        old = old_ligands.get(index)
        new = new_by_index.get(index)

        if old is None:
            changes.append(('gained', f'{new["residue"]} {new["atom"]}',
                            f'{_metal_contact_label(new, binding_modes)}, '
                            f'{new["mode"]}'))
        elif new is None:
            changes.append(('lost', f'{old["residue"]} {old["atom"]}', 'was '
                            f'{_metal_contact_label(old, binding_modes)}, '
                            f'{old["mode"]}'))
        elif old_coordination[index] != new_coordination[index]:
            detail = []
            if old['mode'] != new['mode']:
                detail.append(f'{old["mode"]} -> {new["mode"]}')
            if old['metals'] != new['metals']:
                detail.append(_metal_contact_label(new, binding_modes))
            changes.append(('changed', f'{new["residue"]} {new["atom"]}',
                            ', '.join(detail)))

    # a residue that lost every contact is still part of the truncated
    # active site, and only a new extraction can take it out
    dropped_residues = sorted({
        ligand['residue']
        for ligand in binding_modes['ligands']
        if ligand['index'] not in new_by_index
    } - {ligand['residue'] for ligand in new_ligands})

    _print_binding_mode_update(changes,
                               largest_shift,
                               dropped_residues,
                               ostream=ostream)

    new_binding_modes = {
        'metals': deepcopy(binding_modes['metals']),
        'ligands': new_ligands,
        'variants': deepcopy(binding_modes.get('variants', {})),
        'coordinating_residues': sorted(
            binding_modes.get('coordinating_residues', [])),
        'extra_residues': sorted(binding_modes.get('extra_residues', [])),
        'excluded_residues': sorted(binding_modes.get('excluded_residues', [])),
        'manual_bonds': deepcopy(records),
        'notes': notes,
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

    _print_binding_modes(new_binding_modes, ostream=ostream)

    new_active_site = dict(active_site)
    new_active_site['connectivity_matrix'] = new_matrix

    return new_binding_modes, new_active_site, True


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
        metal['index']: metal['element'] for metal in binding_modes['metals']
    }

    return ', '.join(
        f'{elements.get(index, "metal")}{index} {distance:.2f} A'
        for index, distance in zip(ligand['metals'], ligand['distances']))


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
        'extract_pairs: initial_bond_range '
        'requires coordinates')

    def neighbors(index, depth):
        if depth == 0 and initial_bond_range is not None:
            coords = np.asarray(coordinates)
            distances = np.linalg.norm(coords - coords[index], axis=1)
            found = set(np.where(distances <= initial_bond_range)[0].tolist())
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


# ----------------------------------------------------------------------
# geometry
# ----------------------------------------------------------------------


def constrained_indices(active_site, constrain_capping_hydrogens=False):
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

    if constrain_capping_hydrogens:
        indices |= set(active_site['cap_indices'])

    return sorted(indices)


def freeze_constraints(active_site,
                       frozen_indices=None,
                       constrain_capping_hydrogens=False):
    """
    Returns the geomeTRIC constraint the optimization is run under.

    What is frozen is what the fitted parameters end up describing, so it is
    worth being able to read it before paying for the optimization that uses
    it.

    :param active_site:
        The extracted active site.
    :param frozen_indices:
        The zero-based active site indices to freeze. Defaults to
        constrained_indices().
    :param constrain_capping_hydrogens:
        Whether the capping hydrogens are frozen along with the beta carbons,
        when frozen_indices is left to the default.

    :return:
        The constraint as geomeTRIC reads it, or None when nothing is frozen.
    """

    if frozen_indices is None:
        frozen_indices = constrained_indices(
            active_site,
            constrain_capping_hydrogens=constrain_capping_hydrogens)

    if len(frozen_indices) == 0:
        return None

    # geomeTRIC indexes atoms from one
    selection = ','.join(str(index + 1) for index in sorted(frozen_indices))

    return f'freeze xyz {selection}'


def minimize_active_site(active_site,
                         forcefield,
                         frozen_indices=None,
                         max_iterations=0,
                         constrain_capping_hydrogens=False):
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

    assert_msg_critical('openmm' in sys.modules,
                        'minimize_active_site: openmm is '
                        'required')

    if frozen_indices is None:
        frozen_indices = constrained_indices(
            active_site,
            constrain_capping_hydrogens=constrain_capping_hydrogens)

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
        coords = np.array(state.getPositions().value_in_unit(mmunit.angstrom))

    return coords


def mm_optimize_active_site(active_site,
                            forcefield,
                            frozen_indices=None,
                            constrain_metals=False,
                            constrain_capping_hydrogens=False,
                            max_iterations=0,
                            bond_change_warning=0.25,
                            ostream=None):
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
    :param frozen_indices:
        The active site indices to hold fixed. Defaults to
        constrained_indices().
    :param constrain_metals:
        Whether to hold the metal centers as well, added to whatever
        frozen_indices amounts to.
    :param constrain_capping_hydrogens:
        Whether the capping hydrogens are frozen along with the beta
        carbons, when frozen_indices is left to the default.
    :param max_iterations:
        The iteration limit. Zero runs until convergence.
    :param bond_change_warning:
        How far a metal-ligand bond may move before it is reported.

    :return:
        The relaxed molecule. The caller decides whether to put it back
        into the active site.
    """

    ostream = _stream(ostream)

    molecule = active_site['molecule']

    if frozen_indices is None:
        frozen_indices = constrained_indices(
            active_site,
            constrain_capping_hydrogens=constrain_capping_hydrogens)

    if constrain_metals:
        frozen_indices = sorted(
            set(frozen_indices) | set(active_site['metal_indices']))

    coordinates = minimize_active_site(active_site,
                                       forcefield,
                                       frozen_indices=frozen_indices,
                                       max_iterations=max_iterations)

    relaxed = Molecule(molecule.get_labels(), coordinates, 'angstrom')
    relaxed.set_charge(molecule.get_charge())
    relaxed.set_multiplicity(molecule.get_multiplicity())

    _print_mm_optimization(active_site,
                           forcefield,
                           relaxed,
                           frozen_indices,
                           max_iterations=max_iterations,
                           bond_change_warning=bond_change_warning,
                           ostream=ostream)

    return relaxed


# ----------------------------------------------------------------------
# qm
# ----------------------------------------------------------------------


def _get_scf_driver(molecule,
                    basis_set_label='def2-svp',
                    scf_drv=None,
                    xcfun='PBE0',
                    comm=None,
                    ostream=None):
    """
    Returns the SCF driver and basis set for the active site

    The driver assigned to scf_drv is used exactly as given, so it carries
    every QM setting beyond the functional and the basis set.

    :param molecule:
        The active site molecule.

    :return:
        The tuple of the SCF driver and the basis set.
    """

    comm = MPI.COMM_WORLD if comm is None else comm
    ostream = _stream(ostream)

    if scf_drv is None:
        if molecule.get_multiplicity() != 1:
            scf_drv = ScfUnrestrictedDriver(comm, ostream)
        else:
            scf_drv = ScfRestrictedDriver(comm, ostream)
        if xcfun is not None:
            scf_drv.xcfun = xcfun
        scf_drv = scf_drv

    basis = MolecularBasis.read(molecule, basis_set_label)

    return scf_drv, basis


def _print_muted_notice(step, mute_scf=True, ostream=None):
    """
    Announces a long calculation whose output is being suppressed.

    :param step:
        A description of the step about to run.
    """

    ostream = _stream(ostream)

    if mute_scf:
        ostream.print_info(
            f'Running {step} with muted QM output. Set mute_scf to False '
            'to follow it.')
    else:
        ostream.print_info(f'Running {step}.')
    ostream.flush()


def _run_scf(scf_drv, molecule, basis, mute_scf=True):
    """
    Runs the SCF for a geometry, whatever the driver already holds.

    Both the gradient and the Hessian driver reuse whatever is already in
    scf_drv.scf_results and only fall back to running an SCF when it is
    empty. They do not check that those results belong to the geometry they
    were handed, so the SCF is run explicitly whenever the geometry may have
    moved since the driver was last used.

    :param scf_drv:
        The SCF driver.
    :param molecule:
        The active site molecule.
    :param basis:
        The basis set.
    :param mute_scf:
        Whether to mute the driver while it runs.
    """

    if mute_scf:
        scf_drv.ostream.mute()

    scf_drv.compute(molecule, basis)

    if mute_scf:
        scf_drv.ostream.unmute()


def optimize_active_site(active_site,
                         frozen_indices=None,
                         mute_scf=True,
                         ostream=None,
                         basis_set_label='def2-svp',
                         scf_drv=None,
                         xcfun='PBE0',
                         constrain_capping_hydrogens=False,
                         comm=None):
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
        The tuple of the optimized molecule and the results of the
        optimization driver. The caller decides whether to put the molecule
        back into the active site.
    """

    if frozen_indices is None:
        frozen_indices = constrained_indices(
            active_site,
            constrain_capping_hydrogens=constrain_capping_hydrogens)

    molecule = active_site['molecule']

    constraint = freeze_constraints(active_site, frozen_indices=frozen_indices)
    constraints = None if constraint is None else [constraint]

    _print_muted_notice(
        f'the constrained optimization with {len(frozen_indices)} '
        'atom(s) frozen',
        mute_scf=mute_scf,
        ostream=ostream)

    scf_drv, basis = _get_scf_driver(molecule,
                                     basis_set_label=basis_set_label,
                                     scf_drv=scf_drv,
                                     xcfun=xcfun,
                                     comm=comm,
                                     ostream=ostream)

    grad_drv = ScfGradientDriver(scf_drv)
    opt_drv = OptimizationDriver(grad_drv)
    opt_drv.constraints = constraints

    if mute_scf:
        grad_drv.ostream.mute()
        opt_drv.ostream.mute()

    opt_results = opt_drv.compute(molecule, basis)

    if mute_scf:
        grad_drv.ostream.unmute()
        opt_drv.ostream.unmute()

    optimized = Molecule.read_xyz_string(opt_results['final_geometry'])
    optimized.set_charge(molecule.get_charge())
    optimized.set_multiplicity(molecule.get_multiplicity())

    return optimized, opt_results


def compute_hessian(active_site,
                    atom_pairs=None,
                    mute_scf=True,
                    basis_set_label='def2-svp',
                    scf_drv=None,
                    xcfun='PBE0',
                    ostream=None,
                    comm=None):
    """
    Computes the nuclear Hessian of the active site.

    With atom_pairs the analytical Hessian is restricted to those blocks
    and everything else is left at zero, which is all the Seminario method
    needs when only the metal terms are being fitted. The diagonal blocks
    are added by the Hessian driver itself, so only the off-diagonal pairs
    need to be given.

    compute() passes the pairs of extract_pairs unless
    calculate_partial_hessian is off, in which case it asks for the whole
    Hessian instead.

    :param atom_pairs:
        The list of zero-based (i, j) tuples, typically from
        extract_pairs. None computes the full Hessian.

    :return:
        The Hessian as a (3N, 3N) numpy array in Hartree per Bohr squared.
    """

    molecule = active_site['molecule']

    scf_drv, basis = _get_scf_driver(molecule,
                                     basis_set_label=basis_set_label,
                                     scf_drv=scf_drv,
                                     xcfun=xcfun,
                                     comm=comm,
                                     ostream=ostream)

    assert_msg_critical(
        scf_drv.solvation_model is None, 'compute_hessian: ScfHessianDriver '
        'does not support a solvation model')

    # The Hessian driver reuses scf_drv.scf_results without checking
    # which geometry they belong to, so the SCF is run here for the
    # current one. This costs nothing: the driver would otherwise run the
    # same SCF itself.
    _print_muted_notice('the SCF for the Hessian',
                        mute_scf=mute_scf,
                        ostream=ostream)
    _run_scf(scf_drv, molecule, basis, mute_scf=mute_scf)

    _print_muted_notice('the Hessian', mute_scf=mute_scf, ostream=ostream)

    hessian_drv = ScfHessianDriver(scf_drv)
    # the numerical path ignores atom_pairs entirely and would silently
    # compute the full Hessian instead
    hessian_drv.numerical = False
    if atom_pairs is None:
        hessian_drv.atom_pairs = None
    else:
        hessian_drv.atom_pairs = [tuple(pair) for pair in atom_pairs]

    if mute_scf:
        hessian_drv.ostream.mute()

    hessian_drv.compute(molecule, basis)

    if mute_scf:
        hessian_drv.ostream.unmute()

    hessian = np.copy(hessian_drv.hessian)

    return hessian


def compute_resp_charges(active_site, mute_scf=True, comm=None, ostream=None):
    """
    Computes RESP charges for the active site.

    :return:
        The partial charges as an (N,) numpy array.
    """

    comm = MPI.COMM_WORLD if comm is None else comm
    ostream = _stream(ostream)

    molecule = active_site['molecule']

    _print_muted_notice('the RESP charge fit at Hartree-Fock/6-31G*',
                        mute_scf=mute_scf,
                        ostream=ostream)

    resp_drv = RespChargesDriver(comm, ostream)

    if mute_scf:
        resp_drv.ostream.mute()

    # Neither a basis nor SCF results are passed: the driver then defaults
    # to Hartree-Fock with 6-31G*, which is what RESP charges are meant to
    # be fitted to, and runs its own SCF. Handing it the active site's own
    # functional and basis would silently fit the charges at a level the
    # RESP parameters were never derived for.
    charges = resp_drv.compute(molecule)

    if mute_scf:
        resp_drv.ostream.unmute()

    charges = comm.bcast(charges, root=mpi_master())
    charges = np.array(charges)

    return charges


def d4_charges(active_site, ostream=None):
    """
    Returns D4 partial charges for the active site.

    The fallback wherever charges are wanted and none were fitted. They
    are cheap - a fraction of a millisecond - and they sum to the charge
    of the site exactly, but they are not RESP: on a binuclear zinc site
    they put about 0.3 e less on each metal. Good enough to relax a
    geometry on, and worth knowing about before they reach anything else.

    :param active_site:
        The active site.

    :return:
        The charges as an (N,) numpy array, capping hydrogens included.
    """

    ostream = _stream(ostream)

    molecule = active_site['molecule']
    charges = np.array(molecule.get_partial_charges(molecule.get_charge()))

    ostream.print_info(
        f'Using D4 partial charges: {charges.size} atoms summing to '
        f'{charges.sum():+.3f} e. No RESP charges were supplied.')
    ostream.flush()

    return charges


# ----------------------------------------------------------------------
# fitting
# ----------------------------------------------------------------------


def get_metal_keys(forcefield, active_site):
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


def _manual_bond_records_by_key(active_site, topology, binding_modes):
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
    caps = set(active_site['cap_indices'])
    site_of = {
        top_index: site_index
        for site_index, top_index in active_site['atom_map'].items()
        if site_index not in caps
    }

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


def manual_bond_keys(active_site, topology, binding_modes):
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

    return set(
        _manual_bond_records_by_key(active_site, topology, binding_modes))


def manual_bond_equilibria(active_site, topology, binding_modes):
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
        key: record['equilibrium'] for key, record in
        _manual_bond_records_by_key(active_site, topology,
                                    binding_modes).items()
        if record.get('equilibrium') is not None
    }


def build_forcefield(
        active_site,
        hessian=None,
        partial_charges=None,
        average_metal_terms=False,
        prune_weak_bridge_bonds=True,
        reparameterize_metal_angles=True,
        comm=None,
        ostream=None,
        default_metal_angle_force_constant=DEFAULT_METAL_ANGLE_FORCE_CONSTANT,
        default_metal_bond_force_constant=DEFAULT_METAL_BOND_FORCE_CONSTANT,
        metal_angle_equilibria=None,
        metal_bond_equilibria=None,
        protected_bonds=None,
        bond_equilibria=None,
        weak_bridge_tolerance=0.25):
    """
    Builds the active site force field and fits the metal terms.

    Without a Hessian the metal terms are seeded by _seed_metal_terms
    instead of fitted, which is the force field the crude pre-QM pass
    runs on: getting the equilibrium geometry roughly right matters far
    more than the stiffness at that stage.

    :param hessian:
        The Hessian as a (3N, 3N) numpy array. If None, default values for the metal terms are used instead of fitted.
    :param partial_charges:
        The partial charges fitted on the active site. The charge of the
        capping hydrogens is redistributed over the remaining atoms before
        they are applied, since the caps do not exist in the protein.
        If None, D4 charges are used.
    :param protected_bonds:
        The metal bond keys the weak bridge pruning must leave alone,
        which is what a bond added by hand needs: manual_bond_keys turns
        the records of the binding modes into them.
    :param bond_equilibria:
        Equilibrium distances in nanometers for individual metal bonds,
        from manual_bond_equilibria. Only the seeded pass reads them: with
        a Hessian the equilibria come from the geometry it was computed on.

    :return:
        The force field generator.
    """

    comm = MPI.COMM_WORLD if comm is None else comm
    ostream = _stream(ostream)

    assert_msg_critical(
        not isinstance(hessian, str), 'build_forcefield: the Hessian must be '
        'a numpy array.')

    molecule = active_site['molecule']
    n_atoms = molecule.number_of_atoms()

    forcefield = MMForceFieldGenerator(comm, ostream)
    forcefield.connectivity_matrix = np.asarray(
        active_site['connectivity_matrix'])
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
    if partial_charges is None:
        partial_charges = d4_charges(active_site, ostream=ostream)

    partial_charges = np.asarray(partial_charges)
    assert_msg_critical(
        partial_charges.shape == (n_atoms,), 'build_forcefield: expected '
        f'{n_atoms} partial charges, got {partial_charges.shape}')
    # The capping hydrogens stand in for alpha carbons and do not
    # exist anywhere the force field is used, so their charge is
    # folded into the rest of the active site before anything is written.
    partial_charges = redistribute_cap_charges(active_site, partial_charges)
    forcefield.partial_charges = partial_charges
    for index in range(n_atoms):
        forcefield.atoms[index]['charge'] = partial_charges[index]

    annotate_atoms(forcefield, active_site)

    bonds, angles = get_metal_keys(forcefield, active_site)

    # switching the angles off has to leave them at whatever the generator
    # guessed in both passes, so it is gated once, here
    if not reparameterize_metal_angles:
        angles = []

    if hessian is None:
        _seed_metal_terms(
            forcefield,
            active_site,
            bonds,
            angles,
            default_metal_angle_force_constant=(
                default_metal_angle_force_constant),
            default_metal_bond_force_constant=default_metal_bond_force_constant,
            metal_angle_equilibria=metal_angle_equilibria,
            metal_bond_equilibria=metal_bond_equilibria,
            bond_equilibria=bond_equilibria,
            ostream=ostream)
    else:
        hessian = np.asarray(hessian)
        assert_msg_critical(
            hessian.shape == (3 * n_atoms, 3 * n_atoms),
            'build_forcefield: Hessian shape '
            f'{hessian.shape} does not match {(3 * n_atoms, 3 * n_atoms)}')
        forcefield.ostream.mute()
        forcefield.reparameterize(hessian,
                                  reparameterize_keys=bonds + angles,
                                  average_metal_terms=average_metal_terms)
        forcefield.ostream.unmute()

        if prune_weak_bridge_bonds:
            bonds, angles = _prune_weak_bridges(
                forcefield,
                active_site,
                bonds,
                angles,
                weak_bridge_tolerance=weak_bridge_tolerance,
                protected=protected_bonds,
                ostream=ostream)

        _check_force_constants(forcefield,
                               active_site,
                               bonds,
                               angles,
                               hessian,
                               ostream=ostream)

    return forcefield


def annotate_atoms(forcefield, active_site):
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

            atom['comment'] = '; '.join(
                part for part in (comment, note) if part)


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


def _seed_metal_terms(
        forcefield,
        active_site,
        bonds,
        angles,
        default_metal_angle_force_constant=DEFAULT_METAL_ANGLE_FORCE_CONSTANT,
        default_metal_bond_force_constant=DEFAULT_METAL_BOND_FORCE_CONSTANT,
        metal_angle_equilibria=None,
        metal_bond_equilibria=None,
        bond_equilibria=None,
        ostream=None):
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
    :param bond_equilibria:
        Equilibrium distances in nanometers for individual bond keys. These
        beat both the measurement and the element table, since they name one
        bond rather than a kind of bond: manual_bond_equilibria turns what
        add_metal_bond was told into them.
    """

    ostream = _stream(ostream)

    labels = active_site['molecule'].get_labels()
    molecule = active_site['molecule']

    if bond_equilibria is None:
        bond_equilibria = {}

    for key in bonds:
        elements = tuple(labels[index] for index in key)
        equilibrium = bond_equilibria.get(tuple(sorted(key)))
        comment = 'asked for with the bond'

        if equilibrium is None:
            equilibrium = _lookup_equilibrium(metal_bond_equilibria, elements)
            comment = 'given equilibrium'

        if equilibrium is None:
            # the getters index atoms from one, and the force field keeps
            # its bond lengths in nanometers
            equilibrium = 0.1 * molecule.get_distance_in_angstroms(
                [index + 1 for index in key])
            comment = 'measured on the input geometry'
        forcefield.bonds[key]['equilibrium'] = equilibrium
        forcefield.bonds[key]['force_constant'] = (
            default_metal_bond_force_constant)
        forcefield.bonds[key]['comment'] = comment

    for key in angles:
        elements = tuple(labels[index] for index in key)
        equilibrium = _lookup_equilibrium(metal_angle_equilibria, elements)
        if equilibrium is None:
            equilibrium = molecule.get_angle_in_degrees(
                [index + 1 for index in key])
            comment = 'measured on the input geometry'
        else:
            comment = 'given equilibrium'
        forcefield.angles[key]['equilibrium'] = equilibrium
        forcefield.angles[key]['force_constant'] = (
            default_metal_angle_force_constant)
        forcefield.angles[key]['comment'] = comment

    ostream.flush()


def _prune_weak_bridges(forcefield,
                        active_site,
                        bonds,
                        angles,
                        weak_bridge_tolerance=0.25,
                        protected=None,
                        ostream=None):
    """
    Drops the long arm of a bridging residue that the fit gave no force
    constant.

    A residue that reaches two metals can hold one of them far more
    weakly than the other, whether it does the bridging through one atom
    that binds both (mu-1,1) or through two atoms of the same group that
    take one metal each (a mu-1,3 carboxylate). Both are the same
    situation seen from the residue, so the residue is what is looked at:
    the metal bonds of everything the metals leave connected together.

    Two independent things have to agree before an arm is dropped: the
    Hessian, by giving it no force constant at all, and the geometry, by
    holding it at least weak_bridge_tolerance further out than the
    shortest metal bond of that same residue. A zero on its own says
    nothing here - it can equally mean a geometry that is not stationary,
    or a Hessian that never covered the pair - which is why the distance
    has to agree.

    A bond with no stiffness is not the same thing as no bond: it leaves
    the pair at an equilibrium the dynamics never restores while its
    angles, torsions and exclusions all still act as though the two were
    bonded. So the bond goes, and with it every angle, torsion and
    improper that crosses it.

    Only the longer arms are ever dropped, so a bridging residue can
    never lose the contact it is held by.

    :param forcefield:
        The force field being built, whose terms are removed in place.
    :param active_site:
        The active site, for the metal indices and the geometry the fit
        was made on.
    :param bonds:
        The metal bond keys.
    :param angles:
        The metal angle keys.

    :return:
        The metal bond and angle keys that are left.
    """

    ostream = _stream(ostream)

    # a bond asked for by hand is a decision, and the two things this
    # reads - a zero force constant and a long distance - are exactly what
    # a hand-added bond looks like when the Hessian does not cover it, so
    # leaving it in reach of the heuristic would take it straight back out
    protected = {tuple(key) for key in (protected or ())}

    metals = set(active_site['metal_indices'])
    coordinates = active_site['molecule'].get_coordinates_in_angstrom()
    labels = active_site['molecule'].get_labels()
    fragments = _ligand_fragments(forcefield, metals)

    # the metal bonds of each residue, keyed by the fragment it is
    reached = {}
    for key in bonds:
        ligands = [index for index in key if index not in metals]
        if len(ligands) != 1:
            # a metal-metal bond belongs to no residue
            continue
        reached.setdefault(fragments[ligands[0]], []).append(key)

    removed = []
    for keys in reached.values():
        touched = {index for key in keys for index in key if index in metals}
        if len(touched) < 2:
            # one metal is a grip, however many atoms it is made with
            continue

        def length(key):
            return float(
                np.linalg.norm(coordinates[key[0]] - coordinates[key[1]]))

        lengths = {key: length(key) for key in keys}
        shortest = min(lengths.values())

        for key in keys:
            if key in protected:
                continue
            if forcefield.bonds[key]['force_constant'] != 0.0:
                continue
            if lengths[key] - shortest < weak_bridge_tolerance:
                continue
            removed.append((key, lengths[key], shortest))

    if not removed:
        return bonds, angles

    for key, length, shortest in removed:
        pair = set(key)
        names = '-'.join(labels[index] for index in key)

        crossing = {
            'angle': _terms_crossing(forcefield.angles, pair),
            'torsion': _terms_crossing(forcefield.dihedrals, pair),
            'improper': _terms_crossing(forcefield.impropers, pair, path=False),
        }

        del forcefield.bonds[key]
        for table, found in ((forcefield.angles, crossing['angle']),
                             (forcefield.dihedrals, crossing['torsion']),
                             (forcefield.impropers, crossing['improper'])):
            for crossed in found:
                del table[crossed]

        counts = ', '.join(
            f'{len(found)} {name}(s)' for name, found in crossing.items())
        ostream.print_info(
            f'Removed the bridging bond {key} {names}: it was fitted to '
            f'no force constant and is {length:.2f} A long against '
            f'{shortest:.2f} A for the shortest metal bond of the same '
            f'residue. {counts} crossing it were removed with it.')

    ostream.flush()

    # filtered rather than re-derived, so that an angle list the caller
    # emptied on purpose stays empty
    return ([key for key in bonds if key in forcefield.bonds],
            [key for key in angles if key in forcefield.angles])


def _ligand_fragments(forcefield, metals):
    """
    Numbers the residues a force field holds, as what the metals leave
    connected together.

    The bonds of the force field are walked rather than the connectivity
    matrix, so that what counts as one residue is what this force field
    is actually wired as. A metal is not part of any of them, which is
    the whole point: it is what would otherwise join two residues into
    one.

    :param forcefield:
        The force field.
    :param metals:
        The indices of the metal centers.

    :return:
        A dictionary from atom index to fragment number, holding every
        atom that is not a metal.
    """

    neighbors = {}
    for first, second in forcefield.bonds:
        if first in metals or second in metals:
            continue
        neighbors.setdefault(first, set()).add(second)
        neighbors.setdefault(second, set()).add(first)

    fragments = {}
    count = 0

    for atom in forcefield.atoms:
        if atom in metals or atom in fragments:
            continue

        stack = [atom]
        while stack:
            current = stack.pop()
            if current in fragments:
                continue
            fragments[current] = count
            stack.extend(neighbors.get(current, set()) - fragments.keys())

        count += 1

    return fragments


def _terms_crossing(table, pair, path=True):
    """
    Finds the terms of one table that act across a bond.

    A bond, an angle and a torsion are paths, so a term crosses the bond
    when the two atoms are neighbours in its key. An improper is not a
    path - its central atom is written first, with the substituents after
    it - so there it is enough that both atoms are in the key at all.

    :param table:
        The bonds, angles, dihedrals or impropers of a force field.
    :param pair:
        The two atoms of the bond, as a set.
    :param path:
        Whether the keys of the table are paths.

    :return:
        The keys that cross the bond.
    """

    if not path:
        return [key for key in table if pair <= set(key)]

    return [
        key for key in table if any({key[index], key[index + 1]} == pair
                                    for index in range(len(key) - 1))
    ]


def _check_force_constants(forcefield,
                           active_site,
                           bonds,
                           angles,
                           hessian=None,
                           ostream=None):
    """
    Warns about metal terms whose force constant was fitted to zero, and
    says which of the two reasons put it there.

    A term reads its own blocks of the Hessian and nothing else: a bond
    (i, j) reads the (i, j) block and an angle (i, j, k) reads the (i, j)
    and (j, k) blocks. So a zero comes from one of two places, and the fix
    is not the same:

    - **the Hessian does not cover the term.** compute_hessian is
      restricted to the pairs extract_pairs walks out of the connectivity,
      and everything else is left at zero. A bond that was not there when
      those pairs were taken has nothing to project, which is what a
      Hessian reused from a folder or supplied by hand runs into when the
      coordination has moved on since. Only recomputing it on this
      coordination fixes that.
    - **the projection was negative and Seminario clamped it**, which
      means the geometry is not stationary along that coordinate. On an
      unrelaxed structure this typically wipes out the long, strained
      metal-ligand bonds, which are exactly the ones that matter for a
      bridged binuclear site.

    Without a Hessian the two cannot be told apart and everything is
    reported as the second.

    :param bonds:
        The metal bond keys.
    :param angles:
        The metal angle keys.
    :param hessian:
        The Hessian the terms were fitted to, for telling an uncovered
        term from a clamped one.
    """

    ostream = _stream(ostream)

    labels = active_site['molecule'].get_labels()

    zero = [(key, 'bond')
            for key in bonds
            if forcefield.bonds[key]['force_constant'] == 0.0]
    zero += [(key, 'angle')
             for key in angles
             if forcefield.angles[key]['force_constant'] == 0.0]

    if not zero:
        return

    n_bonds = sum(1 for _, kind in zero if kind == 'bond')
    n_angles = len(zero) - n_bonds

    uncovered = [
        (key, kind) for key, kind in zero if not _hessian_covers(hessian, key)
    ]
    clamped = [pair for pair in zero if pair not in uncovered]

    def report(terms):
        for key, kind in terms:
            names = '-'.join(labels[index] for index in key)
            ostream.print_warning(
                f'  zero force constant: {kind} {key} {names}')

    ostream.print_warning(
        f'{n_bonds} of {len(bonds)} metal bond(s) and {n_angles} of '
        f'{len(angles)} metal angle(s) got a zero force constant.')

    if uncovered:
        ostream.print_warning(
            f'{len(uncovered)} of them are not covered by the Hessian at '
            'all: it was restricted to the atom pairs of a different '
            'coordination, so it holds no data for these terms. Recompute '
            'it on this one rather than reusing it.')
        report(uncovered)

    if clamped:
        ostream.print_warning(
            f'{len(clamped)} of them were clamped from a negative '
            'projection. The Hessian is most likely not evaluated at a '
            'stationary point; run optimize_active_site first.')
        report(clamped)

    ostream.flush()


def _hessian_covers(hessian, key):
    """
    Says whether a Hessian holds anything for one term.

    The blocks a term reads are the bonded pairs of its key, which is what
    extract_pairs walks out of the connectivity: (i, j) for a bond and
    (i, j) together with (j, k) for an angle. A block left at exactly zero
    was never computed, since compute_hessian fills only the pairs it is
    given.

    :param hessian:
        The Hessian, or None when there is none to look at.
    :param key:
        The bond or angle key.

    :return:
        True when every block the term reads holds something, and when
        there is no Hessian to say otherwise.
    """

    if hessian is None:
        return True

    hessian = np.asarray(hessian)
    pairs = [(key[index], key[index + 1]) for index in range(len(key) - 1)]

    for first, second in pairs:
        # the driver fills the pairs it is given, and a key is not stored
        # in any particular order, so both orientations are looked at
        block = hessian[3 * first:3 * first + 3, 3 * second:3 * second + 3]
        mirror = hessian[3 * second:3 * second + 3, 3 * first:3 * first + 3]
        if not np.any(block) and not np.any(mirror):
            return False

    return True


# ----------------------------------------------------------------------
# enzyme
# ----------------------------------------------------------------------


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
        len(rest) > 0, 'redistribute_cap_charges: the active site '
        'is nothing but capping hydrogens')

    cap_charge = sum(charges[index] for index in caps)
    charges[list(caps)] = 0.0
    charges[rest] += cap_charge / len(rest)

    return charges


def redistribute_charges(system,
                         topology,
                         active_site,
                         partial_charges,
                         ostream=None):
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
    shift = redistribute_backbone_charges(system,
                                          topology,
                                          active_site,
                                          charges,
                                          ostream=ostream)

    return cap_charge, shift


def redistribute_backbone_charges(system,
                                  topology,
                                  active_site,
                                  partial_charges,
                                  ostream=None):
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

    ostream = _stream(ostream)

    assert_msg_critical('openmm' in sys.modules,
                        'redistribute_backbone_charges: openmm '
                        'is required')

    charges = np.asarray(partial_charges)
    caps = set(active_site['cap_indices'])
    atom_map = active_site['atom_map']

    # the caps map to the alpha carbons they replaced, so they have to be
    # left out or CA would be counted as covered by the active site
    covered = {
        atom_map[index]: charges[index]
        for index in range(len(charges))
        if index not in caps
    }

    # The atom map indexes the topology the active site was extracted
    # from. Handing over a different one - most easily by running
    # prepare_protein after the extraction rather than before it - silently
    # writes the active site charges onto whatever atoms happen to hold those
    # indices, so the elements are checked before anything is modified.
    atoms = list(topology.atoms())
    labels = active_site['molecule'].get_labels()
    mismatched = [
        index for index in covered
        if index >= len(atoms) or atoms[index].element is None or
        atoms[index].element.symbol != labels[next(
            c for c, t in atom_map.items() if t == index and c not in caps)]
    ]

    assert_msg_critical(
        not mismatched, 'redistribute_backbone_charges: the '
        'atom map does not match this topology. Extract the active site '
        'from the same topology the system is built from; '
        'prepare_protein renumbers the atoms, so it must run first.')

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
        'redistribute_backbone_charges: the '
        f'coordination region is off by {difference:+.4f} e with no '
        'uncovered atom to absorb it')

    shift = difference / len(uncovered) if uncovered else 0.0

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

    ostream.print_blank()
    ostream.print_header('Charges written into the protein')
    ostream.print_header(31 * '-')
    ostream.print_header(
        _param('coordination region',
               f'{len(residue_indices)} residues, {len(region)} atoms'))
    ostream.print_header(
        _param('covered by the active site', f'{len(covered)} atoms'))
    ostream.print_header(
        _param('left to the protein', f'{len(uncovered)} atoms'))
    ostream.print_blank()
    ostream.print_header(
        _param('region charge, protein', f'{total_before:+.4f} e'))
    ostream.print_header(
        _param('region charge, fitted', f'{total_after:+.4f} e'))
    ostream.print_header(_param('difference to recover',
                                f'{difference:+.4f} e'))
    ostream.print_header(_param('shift per uncovered atom', f'{shift:+.4f} e'))
    ostream.print_blank()
    # a shift and not a scaling: after the caps are folded in the active site
    # carries its own formal charge exactly, so the amount to recover is
    # always minus the sum of the backbone charges and the scale factor
    # that achieves it is identically zero
    ostream.print_info(
        'The region is corrected as a whole, so the ligand-to-metal '
        'donation the fit captured is kept; balancing each residue to its '
        'formal charge separately would undo it.')
    ostream.print_blank()
    ostream.flush()

    return shift


def create_enzyme_system(topology,
                         active_site,
                         forcefield,
                         partial_charges=None,
                         forcefield_files=('amber14-all.xml',
                                           'amber14/tip3pfb.xml'),
                         ostream=None):
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

    ostream = _stream(ostream)

    assert_msg_critical('openmm' in sys.modules,
                        'create_enzyme_system: openmm is '
                        'required')

    openmm_ff = mmapp.ForceField(*forcefield_files)
    system = openmm_ff.createSystem(topology, nonbondedMethod=mmapp.NoCutoff)

    if partial_charges is None:
        partial_charges = d4_charges(active_site, ostream=ostream)

    redistribute_charges(system,
                         topology,
                         active_site,
                         partial_charges,
                         ostream=ostream)

    atom_map = active_site['atom_map']
    caps = set(active_site['cap_indices'])
    bonds, angles = get_metal_keys(forcefield, active_site)

    bond_force = None
    angle_force = None
    for force in system.getForces():
        if isinstance(force, mm.HarmonicBondForce):
            bond_force = force
        elif isinstance(force, mm.HarmonicAngleForce):
            angle_force = force

    assert_msg_critical(
        bond_force is not None and angle_force is not None,
        'create_enzyme_system: the protein '
        'system has no harmonic bond or angle force to extend')

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

    ostream.print_info(
        f'Added {sum(1 for term in added if term[0] == "bond")} metal '
        f'bond(s) and {sum(1 for term in added if term[0] == "angle")} '
        'metal angle(s) to the enzyme system.')
    ostream.flush()

    return system, added


# ----------------------------------------------------------------------
# persistence
# ----------------------------------------------------------------------


def _resolve_source(supplied, filename, folder=None, ostream=None):
    """
    Applies the precedence the three resolvers share: what the caller
    handed in beats what an earlier run left in the working folder, and
    neither means there is nothing to use.

    Stated once here rather than three times, so the rule cannot come to
    mean different things for a geometry, a Hessian and a set of charges.

    :param supplied:
        What the caller passed, or None.
    :param filename:
        The name the step writes its result under, one of the file name
        constants.
    :param folder:
        The working folder to fall back on.
    :param ostream:
        The output stream. The reuse is announced on it.

    :return:
        The tuple of the source and the flag that is True when it came
        from the folder rather than from the caller. The source is None
        when there is nothing to use.
    """

    if supplied is not None:
        return supplied, False

    source = _folder_file(filename, folder=folder)

    if source is None:
        return None, False

    ostream = _stream(ostream)
    ostream.print_info(f'Reusing {filename} from {folder}.')
    ostream.flush()

    return source, True


def _resolve_optimized_geometry(active_site,
                                folder=None,
                                optimized_geometry=None,
                                ostream=None):
    """
    Validates a geometry supplied through optimized_geometry, or left
    behind in the working folder by an earlier run. The element sequence is checked against the extracted active site
    :param active_site:
        The active site, to validate against.

    :return:
        The molecule, or None if there is nothing to use.
    """

    ostream = _stream(ostream)

    source, _ = _resolve_source(optimized_geometry,
                                GEOMETRY_FILE,
                                folder=folder,
                                ostream=ostream)

    if source is None:
        return None

    if isinstance(source, Molecule):
        molecule = Molecule(source)
    else:
        path = Path(source)
        assert_msg_critical(
            path.is_file(), '_resolve_optimized_geometry: the geometry file '
            f'{path} not found')
        molecule = Molecule.read_xyz_file(str(path))

    active_site = active_site['molecule']

    assert_msg_critical(
        molecule.number_of_atoms() == active_site.number_of_atoms(),
        '_resolve_optimized_geometry: the geometry has '
        f'{molecule.number_of_atoms()} atoms but the extracted active site '
        f'has {active_site.number_of_atoms()}')

    assert_msg_critical(
        list(molecule.get_labels()) == list(active_site.get_labels()),
        '_resolve_optimized_geometry: the elements of '
        'optimized_geometry do not match the extracted active site, so it '
        'describes a different structure')

    molecule.set_charge(active_site.get_charge())
    molecule.set_multiplicity(active_site.get_multiplicity())

    return molecule


def _resolve_hessian(active_site, folder=None, hessian=None, ostream=None):
    """
    Validates a Hessian supplied through the hessian setting, or left
    behind in the working folder by an earlier run.

    :param active_site:
        The active site, to validate the shape against.

    :return:
        The Hessian, or None if there is nothing to use.
    """

    ostream = _stream(ostream)

    source, reused = _resolve_source(hessian,
                                     HESSIAN_FILE,
                                     folder=folder,
                                     ostream=ostream)

    if source is None:
        return None

    assert_msg_critical(
        not isinstance(source, (str, Path)) or Path(source).is_file(),
        f'_resolve_hessian: hessian file {source} not found')

    if isinstance(source, (str, Path)):
        hessian = np.loadtxt(source)
    else:
        hessian = np.asarray(source)

    n_atoms = active_site['molecule'].number_of_atoms()
    expected = (3 * n_atoms, 3 * n_atoms)

    assert_msg_critical(
        hessian.shape == expected,
        f'_resolve_hessian: hessian has shape {hessian.shape} '
        f'but the extracted active site needs {expected}')

    # The shape and the elements match any site of the same composition, so
    # a file left behind by a run whose coordination differed passes both
    # and then fits zeros. A block the metal terms read that was never
    # filled is what says so, and is worth recomputing rather than warning
    # about: an explicitly supplied Hessian is an instruction, a file lying
    # in the folder is a guess.
    if reused and not _hessian_covers_site(hessian, active_site):
        ostream.print_warning(
            f'{HESSIAN_FILE} in {folder} holds nothing for some of the '
            'metal terms of this active site, so it was computed for a '
            'different coordination. Ignoring it.')
        ostream.flush()
        return None

    return hessian


def _hessian_covers_site(hessian, active_site):
    """
    Says whether a Hessian holds data for every metal term of a site.

    The blocks the metal terms read are the pairs extract_pairs walks out
    of the connectivity, which is exactly what compute_hessian fills.

    :param hessian:
        The Hessian.
    :param active_site:
        The active site whose metal terms have to be covered.

    :return:
        True when every pair the terms read holds something.
    """

    pairs, _ = extract_pairs(active_site['connectivity_matrix'],
                             active_site['metal_indices'],
                             bond_count=2)

    return all(_hessian_covers(hessian, pair) for pair in pairs)


def _resolve_partial_charges(active_site,
                             folder=None,
                             partial_charges=None,
                             ostream=None):
    """
    Validates charges supplied through the partial_charges setting, or
    left behind in the working folder by an earlier run.

    :param active_site:
        The active site, to validate the count against.

    :return:
        The partial charges, or None if there is nothing to use.
    """

    ostream = _stream(ostream)

    source, _ = _resolve_source(partial_charges,
                                CHARGES_FILE,
                                folder=folder,
                                ostream=ostream)

    if source is None:
        return None

    if isinstance(source, (str, Path)):
        assert_msg_critical(
            Path(source).is_file(),
            '_resolve_partial_charges: the charges file '
            f'{source} not found')
        charges = np.loadtxt(source)
    else:
        charges = np.asarray(source)

    n_atoms = active_site['molecule'].number_of_atoms()

    assert_msg_critical(
        charges.shape == (n_atoms,),
        f'_resolve_partial_charges: the charges have shape '
        f'{charges.shape} but the extracted active site has {n_atoms} atoms')

    total = float(np.sum(charges))
    expected = int(active_site['molecule'].get_charge())

    if abs(total - expected) > 1.0e-3:
        ostream.print_warning(
            f'The supplied partial charges sum to {total:+.3f}, but the '
            f'active site charge is {expected:+d}')
        ostream.flush()

    return charges


def save_forcefield(filename, forcefield):
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

    MMForceFieldGenerator.save_forcefield_as_json(forcefield, str(filename))


def load_forcefield(filename, active_site=None):
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
        f'load_forcefield: forcefield file {filename} not found')

    forcefield = MMForceFieldGenerator.load_forcefield_from_json_file(
        str(filename))

    if active_site is not None:
        _check_forcefield(forcefield, active_site)
        forcefield.molecule = active_site['molecule']

    forcefield.partial_charges = np.array(
        [atom['charge'] for atom in forcefield.atoms.values()])

    return forcefield


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


def _check_forcefield(forcefield, active_site):
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

    labels = active_site['molecule'].get_labels()
    elements = _forcefield_elements(forcefield)

    assert_msg_critical(
        len(elements) == len(labels),
        f'_check_forcefield: the force field has {len(elements)} '
        f'atoms but the extracted active site has {len(labels)}')

    assert_msg_critical(
        elements == labels,
        '_check_forcefield: the elements of the force field do '
        'not match the extracted active site, so it describes a different '
        'structure')


# ----------------------------------------------------------------------
# reporting
# ----------------------------------------------------------------------


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


def _print_param_list(label, items, value_width=20, ostream=None):
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

    ostream = _stream(ostream)

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
        ostream.print_header(_param(label if i == 0 else '', chunk))


def _print_binding_modes(binding_modes, ostream=None):
    """
    Prints the detected coordination sphere.
    """

    ostream = _stream(ostream)

    ostream.print_blank()
    ostream.print_header('Coordination sphere')
    ostream.print_header(19 * '-')

    for metal in binding_modes['metals']:
        ostream.print_header(
            _param(f'metal {metal["element"]} (index {metal["index"]})',
                   f'charge {metal["formal_charge"]:+d}'))

    ostream.print_blank()
    valstr = '{:>10} {:>9} | {:>18} | {:>16}'.format('residue', 'atoms',
                                                     'distances (A)', 'mode')
    ostream.print_header(valstr)
    ostream.print_header(60 * '-')

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
            distances = ', '.join(
                f'{d:.2f}' for ligand in row for d in ligand['distances'])
            modes = '/'.join(dict.fromkeys(ligand['mode'] for ligand in row))
            valstr = '{:>10} {:>9} | {:>18} | {:>16}'.format(
                row[0]['residue'], atoms, distances, modes)
            ostream.print_header(valstr)

    for note in binding_modes['notes']:
        ostream.print_warning(note)

    ostream.print_blank()
    ostream.flush()


def _print_binding_mode_update(changes,
                               largest_shift,
                               dropped_residues=None,
                               ostream=None):
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

    ostream = _stream(ostream)

    ostream.print_blank()
    ostream.print_header('Coordination update')
    ostream.print_header(19 * '-')
    ostream.print_header(_param('largest bond change',
                                f'{largest_shift:.2f} A'))

    if not changes:
        ostream.print_header(_param('coordination', 'unchanged'))
        ostream.print_blank()
        ostream.print_info(
            'The new geometry gives the same coordination sphere; the '
            'binding modes and the connectivity matrix are kept as they '
            'are.')
        ostream.print_blank()
        ostream.flush()
        return

    ostream.print_header(_param('contacts changed', len(changes)))
    ostream.print_blank()

    valstr = '{:>10} {:>12} | {:>46}'.format('change', 'atom', 'detail')
    ostream.print_header(valstr)
    ostream.print_header(72 * '-')

    for kind, atom, detail in changes:
        ostream.print_header('{:>10} {:>12} | {:>46}'.format(
            kind, atom, detail))

    ostream.print_blank()
    ostream.print_info(
        'The binding modes and the connectivity matrix were updated to '
        'the new geometry. Overwrite the ones you hold, or the fit will '
        'use a coordination the geometry no longer has.')

    for residue in dropped_residues or []:
        ostream.print_warning(
            f'{residue} no longer coordinates a metal, but it is still '
            'part of the truncated active site; extract the active site '
            'again to leave it out')

    ostream.print_blank()
    ostream.flush()


def _print_active_site(active_site, binding_modes, ostream=None):
    """
    Prints the composition of the truncated active site.
    """

    ostream = _stream(ostream)

    molecule = active_site['molecule']

    ostream.print_blank()
    ostream.print_header('Truncated active site')
    ostream.print_header(21 * '-')
    ostream.print_header(_param('atoms', molecule.number_of_atoms()))
    ostream.print_header(_param('charge', f'{int(molecule.get_charge()):+d}'))
    ostream.print_header(
        _param('multiplicity', int(molecule.get_multiplicity())))
    ostream.print_header(
        _param('capping hydrogens', len(active_site['cap_indices'])))
    ostream.print_header(
        _param('bonds', int(active_site['connectivity_matrix'].sum() // 2)))
    _print_param_list('residues', active_site['residues'], ostream=ostream)

    variants = sorted(binding_modes['variants'].values())
    _print_param_list('protonation', variants, ostream=ostream)

    ostream.print_blank()
    ostream.flush()


def _print_partial_charges(topology,
                           active_site,
                           partial_charges,
                           ostream=None):
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

    ostream = _stream(ostream)

    charges = np.asarray(partial_charges)
    caps = sorted(active_site['cap_indices'])
    metals = active_site['metal_indices']
    labels = active_site['molecule'].get_labels()
    n_atoms = len(charges)
    rest = [index for index in range(n_atoms) if index not in caps]
    cap_charge = float(sum(charges[index] for index in caps))

    ostream.print_blank()
    ostream.print_header('Partial charges')
    ostream.print_header(15 * '-')
    ostream.print_header(
        _param('active site charge',
               f'{int(active_site["molecule"].get_charge()):+d}'))
    ostream.print_header(_param('fitted total', f'{charges.sum():+.4f} e'))
    ostream.print_header(_param('on capping hydrogens', f'{cap_charge:+.4f} e'))
    ostream.print_header(
        _param(f'spread over {len(rest)} atoms',
               f'{cap_charge / len(rest):+.4f} e each'))
    ostream.print_blank()

    # group what is left by the residue each atom came from
    atoms = list(topology.atoms())
    by_residue = {}
    for index in rest:
        residue = atoms[active_site['atom_map'][index]].residue
        by_residue.setdefault(residue, []).append(index)

    corrected = redistribute_cap_charges(active_site, charges)

    valstr = '{:>16} {:>7} | {:>12}'.format('fragment', 'atoms', 'charge')
    ostream.print_header(valstr)
    ostream.print_header(45 * '-')

    for residue, indices in by_residue.items():
        total = sum(corrected[index] for index in indices)
        if len(indices) == 1 and indices[0] in metals:
            name = f'{labels[indices[0]]} (metal)'
        else:
            name = f'{residue.name}{residue.id}'
        valstr = '{:>16} {:>7} | {:>12.4f}'.format(name, len(indices), total)
        ostream.print_header(valstr)

    ostream.print_blank()
    ostream.flush()


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
        'given', 'measured', or 'mixed'.
    """

    given = {table[key].get('comment') == 'given equilibrium' for key in keys}

    if len(given) != 1:
        return 'mixed'

    return 'given' if given.pop() else 'measured'


def _print_mm_optimization(active_site,
                           forcefield,
                           relaxed,
                           frozen_indices,
                           max_iterations=0,
                           bond_change_warning=0.25,
                           ostream=None):
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
    :param max_iterations:
        The iteration limit the relaxation was given.
    :param bond_change_warning:
        How far a metal-ligand bond may move before it is reported.
    """

    ostream = _stream(ostream)

    labels = active_site['molecule'].get_labels()
    molecule = active_site['molecule']
    metals = sorted(active_site['metal_indices'])
    bonds, angles = get_metal_keys(forcefield, active_site)

    before = molecule.get_coordinates_in_angstrom()
    after = relaxed.get_coordinates_in_angstrom()
    shift = np.linalg.norm(after - before, axis=1)

    ostream.print_blank()
    ostream.print_header('Crude MM relaxation')
    ostream.print_header(19 * '-')

    ostream.print_header(_param('frozen atoms', len(frozen_indices)))
    ostream.print_header(
        _param('metal centers',
               'frozen' if set(metals) <= set(frozen_indices) else 'free'))
    ostream.print_header(
        _param('iteration limit',
               max_iterations if max_iterations > 0 else 'convergence'))
    ostream.print_header(
        _param(
            'metal bonds',
            f'{len(bonds)}, k = {_seeded_constant(forcefield.bonds, bonds)}'))
    ostream.print_header(
        _param(
            'metal angles', f'{len(angles)}, k = '
            f'{_seeded_constant(forcefield.angles, angles)}'
            if angles else 'left untouched'))
    ostream.print_header(
        _param('bond equilibria', _seeded_equilibria(forcefield.bonds, bonds)))
    if angles:
        ostream.print_header(
            _param('angle equilibria',
                   _seeded_equilibria(forcefield.angles, angles)))
    ostream.print_blank()

    valstr = '{:>12} {:>8} | {:>10} | {:>9} | {:>8}'.format(
        'atoms', 'elements', 'before (A)', 'after (A)', 'change')
    ostream.print_header(valstr)
    ostream.print_header(60 * '-')

    worst_bond = None
    for key in bonds:
        one_based = [index + 1 for index in key]
        was = molecule.get_distance_in_angstroms(one_based)
        now = relaxed.get_distance_in_angstroms(one_based)
        names = '-'.join(labels[index] for index in key)
        valstr = '{:>12} {:>8} | {:>10.2f} | {:>9.2f} | {:>+8.2f}'.format(
            str(key), names, was, now, now - was)
        ostream.print_header(valstr)
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
            ostream.print_header(valstr)

    ostream.print_blank()

    largest = int(np.argmax(shift))
    ostream.print_header(
        _param('largest shift', f'{shift[largest]:.2f} A on '
               f'{labels[largest]} {largest}'))
    ostream.print_header(_param('mean shift', f'{shift.mean():.2f} A'))

    # A ligand swinging around its metal moves a long way in Cartesian
    # terms while the coordination sphere itself is untouched, so what
    # the pass has to be held to is the bond lengths, not the shifts.
    if worst_bond is not None:
        ostream.print_header(
            _param('largest bond change',
                   f'{worst_bond[1]:+.2f} A on {worst_bond[0]}'))
        if abs(worst_bond[1]) > bond_change_warning:
            ostream.print_warning(
                f'The crude relaxation changed a {worst_bond[0]} bond by '
                f'{worst_bond[1]:+.2f} A. Check the metal terms it was '
                'given before trusting the geometry it produced.')

    ostream.print_blank()
    ostream.flush()


def _print_metal_parameters(active_site, forcefield, ostream=None):
    """
    Prints the fitted metal bonds and angles.
    """

    ostream = _stream(ostream)

    labels = active_site['molecule'].get_labels()
    metals = set(active_site['metal_indices'])
    coords = active_site['molecule'].get_coordinates_in_angstrom()
    bonds, angles = get_metal_keys(forcefield, active_site)

    ostream.print_blank()
    ostream.print_header('Metal bonds')
    ostream.print_header(11 * '-')
    valstr = '{:>12} {:>7} | {:>9} | {:>21}'.format('atoms', 'elements',
                                                    'r0 (A)',
                                                    'k (kcal/mol/A^2)')
    ostream.print_header(valstr)
    ostream.print_header(60 * '-')

    for key in bonds:
        params = forcefield.bonds[key]
        names = '-'.join(labels[index] for index in key)
        # kJ/mol/nm^2 to kcal/mol/A^2
        force_constant = params['force_constant'] / 100.0 / 4.184
        valstr = '{:>12} {:>7} | {:>9.3f} | {:>21.1f}'.format(
            str(key), names, params['equilibrium'] * 10.0, force_constant)
        ostream.print_header(valstr)

    ostream.print_blank()
    ostream.print_header('Metal angles')
    ostream.print_header(12 * '-')
    valstr = '{:>14} {:>9} | {:>12} | {:>19}'.format('atoms', 'elements',
                                                     'theta0 (deg)',
                                                     'k (kJ/mol/rad^2)')
    ostream.print_header(valstr)
    ostream.print_header(60 * '-')

    for key in angles:
        params = forcefield.angles[key]
        names = '-'.join(labels[index] for index in key)
        bridging = key[0] in metals and key[2] in metals
        # the marker gets a fixed-width field of its own, otherwise the
        # centering of print_header would shift the marked line
        valstr = '{:>14} {:>9} | {:>12.1f} | {:>19.1f} {:<9}'.format(
            str(key), names, params['equilibrium'], params['force_constant'],
            'bridging' if bridging else '')
        ostream.print_header(valstr)

    ostream.print_blank()

    for metal_a in sorted(metals):
        for metal_b in sorted(metals):
            if metal_a >= metal_b:
                continue
            distance = np.linalg.norm(coords[metal_a] - coords[metal_b])
            ostream.print_header(
                _param(f'{labels[metal_a]}-{labels[metal_b]} distance',
                       f'{distance:.3f} A'))

    ostream.print_blank()
    ostream.flush()

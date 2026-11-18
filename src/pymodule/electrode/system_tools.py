"""Everything on top of the graphite: solute placement, solvation (VeloxChem
SolvationBuilder), the OpenMM system, the uniform field, FEP and MD helpers.

The graphite itself comes from graphite_builder.GraphiteBuilder; this module only calls
its outputs (molecule/positions, forcefield_xml, apply_to_system, add_water_correction).
build_system.py is the working script that uses these steps.
"""
from __future__ import annotations

import datetime
import io
import itertools
import json
import os
import shutil
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import openmm as mm
from openmm import app, unit
from scipy.spatial import cKDTree

from ..environment import get_data_path  # noqa: E402

from .graphite import GraphiteBuilder, nonbonded_force  # noqa: E402
from .graphite.graphite import _min_image, graphite_forcefield_xml  # noqa: E402
from .graphite.params import (CARBON_SETS, FORCE_GROUPS,  # noqa: E402
                                     INTERLAYER, check_status)

NM = unit.nanometer
FARADAY_KJ_PER_MOL_PER_V = 96.485332    # e * 1 V in kJ/mol
DA_PER_NM3_TO_KG_PER_M3 = 1.66053906660  # 1 Da/nm^3 = 1.6605 kg/m^3
SOLUTE_RESNAME = 'MOL'
FIELD_FORCE_NAME = 'electrode uniform field'


# =========================================================== electrolyte registry

# Hydroxide models (two-site; the O-H bond is constrained via HBonds).
HYDROXIDE_MODELS = {
    'legacy_gaff': dict(
        qO=-1.186738, qH=0.186738,
        sigO=0.3242871334030835, epsO=0.389112,
        sigH=0.053792464601313685, epsH=0.0196648,
        r_OH=0.09902539941651334, status='legacy',
        source='GAFF oh/ho types + RESP charges, as in electrode_interface_electrolyte V2 '
               '(read from graphene_r000_E0p5/system.xml; there the bond was flexible)',
        note='Not fitted to any KOH solution property. Use for tests and continuity only.'),
    'bonthuis2016': dict(
        # H has no LJ in the paper (sigma 0); 0.1 nm with eps 0 is the same physics and
        # keeps sigma-based soft-core terms (FEP driver) finite.
        qO=-1.0, qH=0.0, sigO=0.381, epsO=0.05, sigH=0.1, epsH=0.0, r_OH=0.100,
        status='literature',
        source='Bonthuis, Mamatkulov & Netz, J. Chem. Phys. 144, 104503 (2016), '
               'doi 10.1063/1.4942771: Table IV (p. 104503-8) O sigma 0.381 nm, eps 0.05 kJ/mol, '
               'q -1.0; H q 0, no LJ; O-H 1.0 A constrained (p. 104503-3); SPC/E water',
        note='Read from the paper 2026-10-23. Fitted to the OH- solvation free energy and Na+OH- '
             'activity derivatives up to 1.5 mol/l. The paper combines sigma GEOMETRICALLY '
             '(Eq. 19); OpenMM NonbondedForce averages sigma (OH- O / SPC/E O: 0.3488 vs '
             '0.3473 nm, +0.4 %). Na+-OH- needed sigma_+- = 1.5 * sqrt(sigma_+ sigma_-) '
             '(Eq. 22) with Smith-Dang Na+ (sigma 0.261 nm, eps 0.4186 kJ/mol); K+-OH- was '
             'not fitted, so K+ pairing with this OH- is unvalidated (see '
             'cation_anion_sigma_scale in assemble.build).'),
}

# Hydronium models (rigid: O-H via HBonds, H-H added as constraints, so no angle term).
HYDRONIUM_MODELS = {
    'bonthuis2016': dict(
        qO=-1.4, qH=0.8, sigO=0.31, epsO=0.8, sigH=0.1, epsH=0.0, r_OH=0.098, r_HH=0.1619,
        status='literature',
        source='Bonthuis, Mamatkulov & Netz, J. Chem. Phys. 144, 104503 (2016), '
               'doi 10.1063/1.4942771: Table IV (p. 104503-8) O sigma 0.31 nm, eps 0.8 kJ/mol, '
               'q -1.4; H q +0.8, no LJ; O-H 0.98 A and H-H 0.1619 nm (H-O-H 111.4 deg), all '
               'distances constrained (p. 104503-3); SPC/E water',
        note='Read from the paper 2026-10-30. Fitted to the H3O+ solvation free energy (with '
             'Smith-Dang Cl- as reference ion) and H3O+Cl- activity derivatives, with geometric '
             'sigma combination (OpenMM averages sigma). H has no LJ in the paper; 0.1 nm with '
             'eps 0 is the same physics. Non-reactive: no proton hopping (Grotthuss), so H3O+ '
             'mobility is too low. Pairing with HSO4- was never fitted.'),
}

OPENMM_SPCE = 'amber14/spce.xml'   # SPC/E + Joung-Cheatham ions (frcmod.ionsjc_spce)
ELECTROLYTE_INPUTS = get_data_path() / 'electrode' / 'electrolyte'

# Solvent species for VeloxChem's SolvationBuilder. 'ffxml' is an OpenMM force-field
# file (or 'hydroxide:<model>' for a generated one); 'resname' and 'atom_names' must
# match its residue template. Water must be named HOH so OpenMM makes it rigid.
SPECIES = {
    'water_spce': dict(smiles='O', resname='HOH', atom_names=['O', 'H1', 'H2'],
                       ffxml=OPENMM_SPCE, status='literature',
                       source='SPC/E (Berendsen 1987), rigid via rigidWater; OpenMM amber14/spce.xml'),
    'K_jc_spce': dict(smiles='[K+]', resname='K', atom_names=['K'], ffxml=OPENMM_SPCE,
                      status='literature', source='Joung & Cheatham 2008, SPC/E set (amber14/spce.xml)'),
    'Na_jc_spce': dict(smiles='[Na+]', resname='NA', atom_names=['NA'], ffxml=OPENMM_SPCE,
                       status='literature', source='Joung & Cheatham 2008, SPC/E set (amber14/spce.xml)'),
    'Cl_jc_spce': dict(smiles='[Cl-]', resname='CL', atom_names=['CL'], ffxml=OPENMM_SPCE,
                       status='literature', source='Joung & Cheatham 2008, SPC/E set (amber14/spce.xml)'),
    'OH_legacy_gaff': dict(smiles='[OH-]', resname='OHM', atom_names=['O', 'H'],
                           ffxml='hydroxide:legacy_gaff', status='legacy',
                           source=HYDROXIDE_MODELS['legacy_gaff']['source'],
                           note=HYDROXIDE_MODELS['legacy_gaff']['note']),
    'OH_bonthuis2016': dict(smiles='[OH-]', resname='OHM', atom_names=['O', 'H'],
                            ffxml='hydroxide:bonthuis2016', status='literature',
                            source=HYDROXIDE_MODELS['bonthuis2016']['source'],
                            note=HYDROXIDE_MODELS['bonthuis2016']['note']),
    'H3O_bonthuis2016': dict(smiles='[OH3+]', resname='H3O', atom_names=['O', 'H1', 'H2', 'H3'],
                             ffxml='hydronium:bonthuis2016', status='literature',
                             rigid_HH=HYDRONIUM_MODELS['bonthuis2016']['r_HH'],
                             source=HYDRONIUM_MODELS['bonthuis2016']['source'],
                             note=HYDRONIUM_MODELS['bonthuis2016']['note']),
    # S first: the anion centre for slot counting (anion_center = first atom).
    'HSO4_gaff2_resp': dict(smiles='S(=O)(=O)([O-])O', resname='HSO',
                            atom_names=['S', 'O1', 'O2', 'O3', 'O4', 'H'],
                            ffxml=str(ELECTROLYTE_INPUTS / 'HSO4.xml'), status='computed',
                            source='GAFF2 + RESP HF/6-31G* at the B3LYP-D3/def2-SVP geometry '
                                   '(inputs/electrolyte/make_hso4.py, HSO4_provenance.json)',
                            note='Same recipe as the solutes; not fitted to any H2SO4 solution '
                                 'property (check the bulk density against experiment).'),
}

# Counterion name (as VeloxChem's SolvationBuilder names it) -> species key. The molecular
# ones (H3O, HSO4) are supplied to the packer by make_packer; the acid neutralises with its
# own ions, so a charged solute in H2SO4 gets HSO4- (and a rounding excess H3O+/HSO4-).
COUNTERION_SPECIES = {'K': 'K_jc_spce', 'Na': 'Na_jc_spce', 'Cl': 'Cl_jc_spce',
                      'H3O': 'H3O_bonthuis2016', 'HSO4': 'HSO4_gaff2_resp'}
MOLECULAR_COUNTERIONS = ('H3O', 'HSO4')


# =========================================================== solute placement

def plane_frame(xyz, masses, plane=None):
    """Rotate so the smallest-variance axis (molecular-plane normal) is z; COM at origin.
    plane: solute atom indices that define the plane instead (the ring of a molecule with a
    side chain): their mass-weighted centre goes to the origin, the axes are fitted to them
    alone, and +z points to the side the other heavy atoms are on (side chain up)."""
    if plane is None:
        com = (xyz * masses[:, None]).sum(0) / masses.sum()
        x = xyz - com
        _, vecs = np.linalg.eigh((x * masses[:, None]).T @ x)
        rot = vecs[:, ::-1]
        if np.linalg.det(rot) < 0:
            rot[:, 2] *= -1
        return x @ rot
    sel = np.asarray(plane, int)
    m = masses[sel]
    x = xyz - (xyz[sel] * m[:, None]).sum(0) / m.sum()
    _, vecs = np.linalg.eigh((x[sel] * m[:, None]).T @ x[sel])
    rot = vecs[:, ::-1]
    if np.linalg.det(rot) < 0:
        rot[:, 2] *= -1
    y = x @ rot
    rest = np.setdiff1d(np.flatnonzero(masses > 2.0), sel)
    if len(rest) and y[rest, 2].mean() < 0:
        y[:, 1:] *= -1                       # 180 deg about x: side chain on the +z side
    return y


def place_solute(xyz, masses, slab, height_nm=0.40, flip=False, rotation=None, xy=None,
                 min_contact_nm=0.25, plane=None):
    """Solute coordinates (nm): flat, COM height_nm above the slab's top layer, at xy
    (with plane, the ring atom indices: the ring centre at height_nm, side chain up)
    (default: the box centre). Raises if a solute atom comes within min_contact_nm of
    the top layer plane or of any atom of a modified surface (groups, step H)."""
    slab = GraphiteBuilder.from_model(slab).model
    x = plane_frame(np.asarray(xyz, float), np.asarray(masses, float), plane)
    if flip:
        x[:, 1:] *= -1
    if rotation is not None:
        x = x @ np.asarray(rotation).T
    cx, cy = (slab.box_xy[0] / 2, slab.box_xy[1] / 2) if xy is None else xy
    out = x + np.array([cx, cy, slab.top_z + height_nm])
    if out[:, 2].min() - slab.top_z < min_contact_nm:
        raise ValueError(f'solute atom only {out[:, 2].min() - slab.top_z:.3f} nm above the top '
                         'layer; raise height_nm')
    others = slab.xyz[[e != 'C' or not s for e, s in zip(slab.elements, slab.skeleton)]]
    if len(others):
        dv = _min_image((out[:, None] - others[None]).reshape(-1, 3), slab.periodic, slab.cell)
        d = np.linalg.norm(dv, axis=1).min()
        if d < min_contact_nm:
            raise ValueError(f'solute atom {d:.3f} nm from a surface group or edge H; move it '
                             '(xy) or raise height_nm')
    return out


def combined_molecule(solute_labels, solute_xyz_nm, slab, solute_charge):
    """VeloxChem Molecule: solute atoms first, then the carbon structure (angstrom), for packing."""
    import veloxchem as vlx
    slab = GraphiteBuilder.from_model(slab).model
    labels = list(solute_labels) + list(slab.elements)
    coords = 10.0 * np.vstack([np.asarray(solute_xyz_nm), slab.xyz])
    mol = vlx.Molecule(labels, coords, 'angstrom')
    mol.set_charge(float(solute_charge) + float(round(slab.charges().sum(), 6)))
    mol.set_multiplicity(1)
    return mol


# =========================================================== solvation (SolvationBuilder)

NEEDED_SOLVATE_API = ('target_density', 'neutralize', 'strict')


def solvation_builder_class(force_vendored=False):
    """VeloxChem's SolvationBuilder if its custom_solvate has the arguments used here,
    else the vendored copy (solvationbuilder_vendored.py, the same code as the local
    install). Returns (class, description)."""
    import inspect
    import os
    import veloxchem as vlx
    have = inspect.signature(vlx.SolvationBuilder.custom_solvate).parameters
    missing = [p for p in NEEDED_SOLVATE_API if p not in have]
    forced = force_vendored or bool(os.environ.get('ELECTRODE_VENDORED_SOLVATION'))
    if not missing and not forced:
        return vlx.SolvationBuilder, f'veloxchem {getattr(vlx, "__version__", "?")} (installed)'
    from .solvationbuilder_vendored import SolvationBuilder
    why = f"installed one lacks {', '.join(missing)}" if missing else 'forced'
    return SolvationBuilder, f'vendored solvationbuilder_vendored.py ({why})'


def make_packer(n_solute, clearance_heavy=2.5, clearance_h=1.5, redraws=200, surface=False,
                threshold=None):
    """SolvationBuilder that keeps solvent out of the slab.

    n_solute: number of solute atoms at the start of the packing 'solute'
    molecule; the graphite atoms follow them (see graphite.combined_molecule).
    Clearances are in angstrom from the outer carbon planes, or (surface=True) from
    every carbon-model atom, periodic in x and y. threshold: the packer's clash distance
    between any two atoms (angstrom; None keeps SolvationBuilder's 1.8).
    """
    base, source = solvation_builder_class()

    class SlabAwareSolvationBuilder(base):

        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            if threshold is not None:
                self.threshold = threshold
            self.stats = {'slab_rejections': 0, 'periodic_rejections': 0,
                          'solvation_builder': source}

        def _slab_range(self):
            if not hasattr(self, '_slab_zr'):
                z = [atom[-1][2] for atom in self.system if atom[0] == 0][n_solute:]
                self._slab_zr = (min(z), max(z)) if z else None
            return self._slab_zr

        def _surface_tree(self):
            if not hasattr(self, '_stree'):
                from scipy.spatial import cKDTree
                pts = np.array([atom[-1] for atom in self.system if atom[0] == 0][n_solute:], float)
                box = np.asarray(self.box, float)
                lo = pts[:, 2].min() - 10.0
                w = pts.copy()
                w[:, :2] = np.mod(w[:, :2], box[:2])
                w[w[:, 0] >= box[0], 0] = 0.0
                w[w[:, 1] >= box[1], 1] = 0.0
                w[:, 2] -= lo
                self._stree = (cKDTree(w, boxsize=[box[0], box[1], w[:, 2].max() + 20.0]), lo, box)
            return self._stree

        def _hits_slab(self, labels, coords):
            if surface:
                tree, lo, box = self._surface_tree()
                for label, c in zip(labels, coords):
                    q = np.array([np.mod(c[0], box[0]) % box[0], np.mod(c[1], box[1]) % box[1], c[2] - lo])
                    if tree.query_ball_point(q, clearance_h if label == 'H' else clearance_heavy):
                        return True
                return False
            zr = self._slab_range()
            if zr is None:
                return False
            for label, c in zip(labels, coords):
                clr = clearance_h if label == 'H' else clearance_heavy
                if zr[0] - clr < c[2] < zr[1] + clr:
                    return True
            return False

        def _periodic_clash(self, coords, tree):
            if tree is None:
                return False
            box = np.asarray(self.box, float)
            thr = self.threshold
            for c in coords:
                options = []
                for d in range(3):
                    opts = [0.0]
                    if c[d] < thr:
                        opts.append(box[d])
                    if c[d] > box[d] - thr:
                        opts.append(-box[d])
                    options.append(opts)
                for shift in itertools.product(*options):
                    if shift == (0.0, 0.0, 0.0):
                        continue
                    if tree.query_ball_point(c + np.array(shift), thr):
                        return True
            return False

        def _counterion_molecules(self):
            if self.ion_name in MOLECULAR_COUNTERIONS:      # H3O+ / HSO4- of an acid electrolyte
                return species_molecules([COUNTERION_SPECIES[self.ion_name]], allow_unverified=True)[0][0]
            return super()._counterion_molecules()

        def _insert_molecule(self, new_molecule, tree):
            labels = new_molecule.get_labels()
            for _ in range(redraws):
                result = super()._insert_molecule(new_molecule, tree)
                if not result:
                    return result
                coords = [r[-1] for r in result]
                if self._hits_slab(labels, coords):
                    self.stats['slab_rejections'] += 1
                    continue
                if self._periodic_clash(coords, tree):
                    self.stats['periodic_rejections'] += 1
                    continue
                return result
            return None

    return SlabAwareSolvationBuilder()


def species_molecules(species_keys, allow_unverified=False):
    """VeloxChem molecules (charge, multiplicity set) for SPECIES keys, plus warnings."""
    import veloxchem as vlx
    mols, warnings = [], []
    for key in species_keys:
        spec = SPECIES[key]
        w = check_status('species', key, spec, allow_unverified)
        if w:
            warnings.append(w)
        m = vlx.Molecule.read_smiles(spec['smiles'])
        m.set_multiplicity(1)
        if len(m.get_labels()) != len(spec['atom_names']):
            raise ValueError(f'{key}: molecule has {len(m.get_labels())} atoms, template '
                             f'{len(spec["atom_names"])}')
        mols.append(m)
    return mols, warnings


def pure_water(water='water_spce'):
    """Species and proportions for pure water, the default electrolyte of build()."""
    return [water], [1.0]


def koh_mixture(molality, water='water_spce', hydroxide='OH_bonthuis2016', cation='K_jc_spce'):
    """Species keys and molecular proportions for KOH(aq) at a given molality (0: pure water).

    Proportions are molecule ratios (water : OH- : K+ = 55.51 : m : m), which is
    how SolvationBuilder.custom_solvate defines a mixture.
    """
    if molality == 0:
        return pure_water(water)
    return [water, hydroxide, cation], [55.508, float(molality), float(molality)]


def h2so4_mixture(molality, water='water_spce', cation='H3O_bonthuis2016', anion='HSO4_gaff2_resp'):
    """Species keys and molecular proportions for H2SO4(aq) at a given molality (mol H2SO4 per
    kg of water, as weighed in), modelled as H3O+ + HSO4- (first dissociation complete,
    second neglected). Each H3O+ takes one water: water : H3O+ : HSO4- = 55.51 - m : m : m.
    """
    if molality == 0:
        return pure_water(water)
    m = float(molality)
    return [water, cation, anion], [55.508 - m, m, m]


def molarity_to_molality(molarity, solution_density_kg_l, molar_mass_kg_mol=0.098079):
    """mol/l of solution -> mol/kg of water (H2SO4 molar mass by default)."""
    return molarity / (solution_density_kg_l - molarity * molar_mass_kg_mol)


ACID_PACKING_THRESHOLD = 1.5   # angstrom


def packer_kwargs(species, packer_options=None):
    """make_packer options for an electrolyte. Random insertion of H3O+/HSO4- jams below
    liquid density at the default 1.8 A clash distance (it forbids H-bond contacts); 1.5 A
    packs 3 M H2SO4 at 1280 kg/m^3, and minimisation removes the close contacts."""
    kw = {}
    if set(species) & {COUNTERION_SPECIES[n] for n in MOLECULAR_COUNTERIONS}:
        kw['threshold'] = ACID_PACKING_THRESHOLD
    kw.update(packer_options or {})
    return kw


def set_counterions(packer, counterion, species):
    """Which ions custom_solvate(neutralize=True) may add. An acid electrolyte (H3O+ and HSO4-
    in species) neutralises with its own ions; otherwise counterion ('K', 'Na' or 'Cl')."""
    keys = set(species)
    if {COUNTERION_SPECIES['H3O'], COUNTERION_SPECIES['HSO4']} <= keys:
        packer.pcharge, packer.ncharge = 'H3O', 'HSO4'
    elif counterion in ('K', 'Na'):
        packer.pcharge = counterion
    elif counterion == 'Cl':
        packer.ncharge = counterion
    elif counterion is not None:
        raise ValueError(f"counterion {counterion!r}: use 'K', 'Na' or 'Cl'")


def total_target_density(liquid_density, slab, solute_mass_da, box_z_nm, wall_gap_nm=0.30):
    """custom_solvate's target_density (kg/m^3 of the whole box, graphite included).

    SolvationBuilder targets the density of everything in the box. To get a
    liquid of liquid_density (kg/m^3) in the effective liquid volume
    A * (Lz - t - 2*wall_gap), add the graphite and solute masses and divide by
    the full box volume. It is a first estimate; calibrate afterwards.
    """
    lx, ly = slab.box_xy
    v_box = lx * ly * box_z_nm
    removed = slab.info.get('removed_layer_fraction', 0.0) * INTERLAYER   # steps (surface.add_step)
    v_liquid = lx * ly * (box_z_nm - slab.thickness + removed - 2 * wall_gap_nm)
    fixed = (slab.mass + solute_mass_da) * DA_PER_NM3_TO_KG_PER_M3   # kg/m^3 * nm^3
    return float((fixed + liquid_density * v_liquid) / v_box), float(v_liquid)


# =========================================================== forces on top of the graphite

def add_ion_pair_sigma_scaling(system, cations, anions, scale, cutoff, switch=None):
    """Cation-anion LJ as in Bonthuis et al. 2016, Eq. 22: sigma_+- = scale*sqrt(sigma_+ sigma_-),
    eps_+- = sqrt(eps_+ eps_-). Added as a correction on top of NonbondedForce (which
    uses the arithmetic mean sigma), so the net cation-anion LJ is the scaled one.
    anions: the LJ-carrying anion sites (e.g. hydroxide O).
    """
    nb = nonbonded_force(system)
    f = mm.CustomNonbondedForce(
        '4*e*((sg/r)^12-(sg/r)^6) - 4*e*((sl/r)^12-(sl/r)^6);'
        'sg = lam*sqrt(sig1*sig2); sl = 0.5*(sig1+sig2); e = sqrt(eps1*eps2)')
    f.setName('cation-anion sigma scaling')
    f.addGlobalParameter('lam', float(scale))
    f.addPerParticleParameter('sig')
    f.addPerParticleParameter('eps')
    for i in range(system.getNumParticles()):
        _, sig, eps = nb.getParticleParameters(i)
        f.addParticle([sig.value_in_unit(mm.unit.nanometer), eps.value_in_unit(mm.unit.kilojoule_per_mole)])
    for k in range(nb.getNumExceptions()):
        i, j = nb.getExceptionParameters(k)[:2]
        f.addExclusion(i, j)
    f.addInteractionGroup(list(cations), list(anions))
    f.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
    f.setCutoffDistance(cutoff)
    if switch:
        f.setUseSwitchingFunction(True)
        f.setSwitchingDistance(switch)
    f.setUseLongRangeCorrection(False)
    f.setForceGroup(FORCE_GROUPS['ion_pair'])
    system.addForce(f)
    return f



def add_uniform_field(system, field_v_per_nm):
    """-q*Ez*z on every charged particle (solute included), in its own force group.

    Ez is a global parameter in kJ/mol/nm/e: change it in a running Context with
    context.setParameter('Ez', E * 96.485332). E > 0 points along +z.
    The FEP driver in fep.py scales the solute's entries with lambda.
    """
    nb = nonbonded_force(system)
    f = mm.CustomExternalForce('-q*Ez*z')
    f.setName(FIELD_FORCE_NAME)
    f.addGlobalParameter('Ez', field_v_per_nm * FARADAY_KJ_PER_MOL_PER_V)
    f.addPerParticleParameter('q')
    for i in range(system.getNumParticles()):
        q = nb.getParticleParameters(i)[0].value_in_unit(mm.unit.elementary_charge)
        if q != 0.0:
            f.addParticle(i, [q])
    f.setForceGroup(FORCE_GROUPS['field'])
    system.addForce(f)
    return f



def add_height_restraint(system, solute, top_layer, masses, z0, k, kind='harmonic', name=None):
    """Restrain the solute COM height above the top graphite layer (relative, no anchor).

    kind: 'harmonic' 0.5*k*(h-z0)^2; 'upper' penalises only h > z0; 'lower' only h < z0.
    """
    h = '(z1-z2)'
    expr = {'harmonic': f'0.5*k_h*({h}-z0_h)^2',
            'upper': f'0.5*k_h*max(0,{h}-z0_h)^2',
            'lower': f'0.5*k_h*min(0,{h}-z0_h)^2'}[kind]
    f = mm.CustomCentroidBondForce(2, expr)
    f.setName(name or f'height restraint ({kind})')
    f.addPerBondParameter('k_h')
    f.addPerBondParameter('z0_h')
    f.addGroup(list(solute), [float(m) for m in masses])
    f.addGroup(list(top_layer), [1.0] * len(top_layer))
    f.addBond([0, 1], [k, z0])
    f.setForceGroup(FORCE_GROUPS['restraint'])
    system.addForce(f)
    return f




# =========================================================== OpenMM system

@dataclass
class ElectrodeSystem:
    system: mm.System
    topology: app.Topology
    positions: np.ndarray        # nm
    indices: dict
    metadata: dict
    ff_sources: list             # what the ForceField was built from (paths or XML text)
    solute_xml: str | None
    solute_pdb: str | None


def hydroxide_xml(model):
    m = HYDROXIDE_MODELS[model]
    if any(m[k] is None for k in ('qO', 'qH', 'sigO', 'epsO', 'sigH', 'epsH', 'r_OH')):
        raise ValueError(f"hydroxide model '{model}' has no parameters yet: {m['note']} "
                         f"Source: {m['source']}.")
    return f"""<ForceField>
 <AtomTypes>
  <Type name="EB-OHM-O" class="EB-OHM-O" element="O" mass="15.9994"/>
  <Type name="EB-OHM-H" class="EB-OHM-H" element="H" mass="1.008"/>
 </AtomTypes>
 <Residues><Residue name="OHM">
  <Atom name="O" type="EB-OHM-O"/><Atom name="H" type="EB-OHM-H"/><Bond atomName1="O" atomName2="H"/>
 </Residue></Residues>
 <HarmonicBondForce><Bond type1="EB-OHM-O" type2="EB-OHM-H" length="{m['r_OH']!r}" k="400000"/></HarmonicBondForce>
 <NonbondedForce coulomb14scale="0.8333333333333334" lj14scale="0.5">
  <Atom type="EB-OHM-O" charge="{m['qO']!r}" sigma="{m['sigO']!r}" epsilon="{m['epsO']!r}"/>
  <Atom type="EB-OHM-H" charge="{m['qH']!r}" sigma="{m['sigH']!r}" epsilon="{m['epsH']!r}"/>
 </NonbondedForce>
</ForceField>"""


def hydronium_xml(model):
    m = HYDRONIUM_MODELS[model]
    atoms = ''.join(f'<Atom name="H{i}" type="EB-H3O-H"/>' for i in (1, 2, 3))
    bonds = ''.join(f'<Bond atomName1="O" atomName2="H{i}"/>' for i in (1, 2, 3))
    return f"""<ForceField>
 <AtomTypes>
  <Type name="EB-H3O-O" class="EB-H3O-O" element="O" mass="15.9994"/>
  <Type name="EB-H3O-H" class="EB-H3O-H" element="H" mass="1.008"/>
 </AtomTypes>
 <Residues><Residue name="H3O">
  <Atom name="O" type="EB-H3O-O"/>{atoms}{bonds}
 </Residue></Residues>
 <HarmonicBondForce><Bond type1="EB-H3O-O" type2="EB-H3O-H" length="{m['r_OH']!r}" k="400000"/></HarmonicBondForce>
 <NonbondedForce coulomb14scale="0.8333333333333334" lj14scale="0.5">
  <Atom type="EB-H3O-O" charge="{m['qO']!r}" sigma="{m['sigO']!r}" epsilon="{m['epsO']!r}"/>
  <Atom type="EB-H3O-H" charge="{m['qH']!r}" sigma="{m['sigH']!r}" epsilon="{m['epsH']!r}"/>
 </NonbondedForce>
</ForceField>"""


def _ff_source(spec):
    src = spec['ffxml']
    if src.startswith('hydroxide:'):
        return ('xml', hydroxide_xml(src.split(':', 1)[1]))
    if src.startswith('hydronium:'):
        return ('xml', hydronium_xml(src.split(':', 1)[1]))
    return ('file', src)


KB_KJ_MOL_K = 0.0083144626


def add_height_window(system, solute, top_layer, masses, lo, hi, k, name='solute height window'):
    """Flat-bottom restraint on h, the solute COM height above the top carbon layer (nm):
    zero for lo <= h <= hi, 0.5*k*(distance outside the window)^2 otherwise. It stays on in
    every FEP window, so a decoupled solute cannot leave the window. Everything it uses sits
    inside the periodic box (slab faces at z = t/2 and Lz - t/2), so z1 - z2 needs no PBC.
    """
    f = mm.CustomCentroidBondForce(
        2, '0.5*k_win*(max(0, h - hi_win)^2 + min(0, h - lo_win)^2); h = z1 - z2')
    f.setName(name)
    for p in ('k_win', 'lo_win', 'hi_win'):
        f.addPerBondParameter(p)
    f.addGroup(list(solute), [float(m) for m in masses])
    f.addGroup(list(top_layer), [1.0] * len(top_layer))
    f.addBond([0, 1], [float(k), float(lo), float(hi)])
    f.setForceGroup(FORCE_GROUPS['restraint'])
    system.addForce(f)
    return f


def window_length(lo, hi, k, temperature=298.15):
    """Accessible COM height range (nm) of a non-interacting solute in the window: the flat
    part plus the two harmonic walls, (hi - lo) + 2*sqrt(pi*kT/(2k))."""
    kt = KB_KJ_MOL_K * temperature
    return float((hi - lo) + 2.0 * np.sqrt(np.pi * kt / (2.0 * k)))


_MISSING = object()


def physical_fingerprint(meta, folder=None):
    """Canonical description of a leg's physical system for pairing legs (review 2026-11-18, v2).
    Metadata: field, liquid composition, nonbonded settings, cation-anion sigma scale, the
    graphite builder record (carbon set and parameters, layers, mobility, face charges,
    surface features/groups), ghost particles, solute size, restraint strength (not its
    window), temperature if recorded. Files: SHA-256 of every force-field XML kept in the leg
    folder (solute, graphite, ions), so the same path with different content is caught."""
    import hashlib
    g = meta.get('graphite') or {}
    fp = {
        'field V/nm': (meta.get('field') or {}).get('V_per_nm'),
        'liquid species counts': (meta.get('liquid') or {}).get('species_counts'),
        'electrolyte species': (meta.get('liquid') or {}).get('electrolyte_species'),
        'nonbonded settings': meta.get('nonbonded'),
        'cation-anion sigma scale': meta.get('cation_anion_sigma_scale'),
        'graphite carbon set': g.get('carbon_set'),
        'graphite carbon parameters': g.get('carbon_parameters'),
        'graphite layers': g.get('layers'),
        'graphite mobility': g.get('mobility'),
        'graphite builder record': g.get('builder'),
        'ghost counterions': meta.get('ghost_counterions'),
        'ghost solute': meta.get('ghost_solute'),
        'solute atoms': (meta.get('solute') or {}).get('n_atoms'),
        'restraint strength k': (meta.get('restraint') or {}).get('k_kj_mol_nm2'),
        'temperature K': (meta.get('pre_fep_md') or {}).get('temperature_K', meta.get('temperature_K')),
    }
    fp = {k: json.loads(json.dumps(v, sort_keys=True, default=str)) for k, v in fp.items()}
    if folder is not None:
        for f in sorted(Path(folder).glob('*.xml')):
            if f.name.startswith('system'):          # system.xml holds the leg's own restraint
                continue
            fp[f'file:{f.name}'] = hashlib.sha256(f.read_bytes()).hexdigest()
    return fp


def adsorption_free_energy(surface_dir, mid_dir, temperature=298.15, delta_nm=1.0, max_dlz_nm=0.25):
    """Adsorption from two FEPs in the same slab box (build_system legs 'surface' and 'mid').

    With F the driver's solvation free energy of each leg (both with ion_corrections='none'),
    and l_s the accessible length of the surface window, the adsorption constant is
        K_ads = l_s * exp(-(F_surface - F_mid) / kT)          (nm)
    i.e. the integral of exp(-W(h)/kT) over the surface window, with W measured from the
    mid-liquid. The mid window cancels if W is ~constant across it (keep it narrow).
    dG_ads(delta) = -kT ln(K_ads / delta); for delta = l_s it equals F_surface - F_mid.
    """
    legs = {}
    for name, d in (('surface', surface_dir), ('mid', mid_dir)):
        d = Path(d)
        meta = json.loads((d / 'metadata.json').read_text())
        res = json.loads((d / 'fep_results.json').read_text())
        legs[name] = (meta, res)
    (ms, rs), (mm_, rm) = legs['surface'], legs['mid']
    # Same box: in z-NPT legs (box_nm_initial present) the start box must match and x, y stay
    # fixed; the final heights differ by the barostat's fluctuations (the FEP ran at each
    # leg's own final box), so only the start box is compared in z.
    key = 'box_nm_initial' if ('box_nm_initial' in ms and 'box_nm_initial' in mm_) else 'box_nm'
    if not np.allclose(ms[key], mm_[key]) or not np.allclose(ms['box_nm'][:2], mm_['box_nm'][:2]):
        raise ValueError(f'legs differ in {key}: {ms[key]} vs {mm_[key]}')
    for path in (('field', 'V_per_nm'), ('graphite', 'carbon_set'), ('solute', 'source_xml')):
        a, b = ms[path[0]][path[1]], mm_[path[0]][path[1]]
        if a != b:
            raise ValueError(f'legs differ in {"/".join(path)}: {a} vs {b}')
    # Physical fingerprint (review 2026-11-18, v2): everything that defines the Hamiltonian
    # and composition must match; only the restraint window (lo/hi) and coordinates may differ.
    fp_s, fp_m = physical_fingerprint(ms, surface_dir), physical_fingerprint(mm_, mid_dir)
    unverified = []
    for key in sorted(set(fp_s) | set(fp_m)):
        a, b = fp_s.get(key, _MISSING), fp_m.get(key, _MISSING)
        if a is _MISSING and b is _MISSING:
            continue
        if a is _MISSING or b is _MISSING:
            if key.startswith('file:'):          # a force-field file only one leg kept
                unverified.append(key); continue
            raise ValueError(f'legs differ in {key}: only one leg records it')
        if a != b:
            raise ValueError(f'legs differ in {key}: {a} vs {b}')
    files = [k for k in fp_s if k.startswith('file:') and k in fp_m]
    if not files:
        unverified.append('force-field files (none kept in both folders)')
    # each z-NPT leg runs its FEP windows at its own last MD box (fep_windows 'fixed'); the
    # final heights differ by barostat noise (SD ~0.04 nm in the campaign; no measurable effect
    # on F, review 2026-11-18). Record the difference; refuse a pair beyond max_dlz_nm.
    dlz = float(ms['box_nm'][2] - mm_['box_nm'][2])
    if abs(dlz) > max_dlz_nm:
        raise ValueError(f'final box heights differ by {dlz:+.3f} nm (> {max_dlz_nm} nm)')
    for name, (meta, res) in legs.items():
        if meta.get('restraint', {}).get('leg') != name:
            raise ValueError(f'{name} folder has restraint leg {meta.get("restraint", {}).get("leg")!r}')
        if res.get('electrode_ion_corrections', 'none').split()[0] not in ('none', 'unavailable'):
            raise ValueError(f'{name} leg used ion corrections {res.get("electrode_ion_corrections")!r}')
    r = ms['restraint']
    kt = KB_KJ_MOL_K * temperature
    l_s = window_length(r['lo_nm'], r['hi_nm'], r['k_kj_mol_nm2'], temperature)
    d_f = float(rs['free_energy']) - float(rm['free_energy'])
    k_ads = l_s * float(np.exp(-d_f / kt))
    unc = lambda r: float(r['uncertainty']) if r.get('uncertainty') is not None else float('nan')   # null -> NaN
    err = float(np.hypot(unc(rs), unc(rm)))
    return dict(F_surface_kj_mol=float(rs['free_energy']), F_mid_kj_mol=float(rm['free_energy']),
                dF_surface_minus_mid_kj_mol=d_f, uncertainty_kj_mol=err,
                surface_window_nm=[r['lo_nm'], r['hi_nm']], surface_window_length_nm=l_s,
                K_ads_nm=k_ads, delta_nm=delta_nm,
                dG_ads_kj_mol=float(-kt * np.log(k_ads / delta_nm)),
                dG_ads_window_reference_kj_mol=d_f, temperature_K=temperature,
                box_z_final_surface_minus_mid_nm=dlz,
                fingerprint=dict(checked=sorted(k for k in fp_s if k in fp_m), unverified=unverified),
                note='dG_ads(delta) = -kT ln(K_ads/delta); the window-reference value uses '
                     'delta = surface window length.')


def _read_solute(solute_pdb, solute_xml):
    pdb = app.PDBFile(solute_pdb)
    xyz = np.array(pdb.positions.value_in_unit(NM))
    masses = np.array([a.element.mass.value_in_unit(unit.dalton) for a in pdb.topology.atoms()])
    solo = app.ForceField(solute_xml).createSystem(pdb.topology, nonbondedMethod=app.NoCutoff)
    nb = nonbonded_force(solo)
    q = sum(nb.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge)
            for i in range(solo.getNumParticles()))
    if abs(q - round(q)) > 1e-3:
        raise ValueError(f'solute charge {q:.4f} is not an integer')
    labels = [a.element.symbol for a in pdb.topology.atoms()]
    return pdb.topology, xyz, masses, int(round(q)), labels


# ---------------------------------------------------------------------- build

def build(graphite, solute_pdb, solute_xml, species=None, proportions=None, *, box_z_nm,
          liquid_density=997.0, wall_gap_nm=0.30, height_nm=0.40, flip=False, rotation=None,
          field_V_per_nm=0.0, counterion='K', cutoff_nm=1.0, switch_nm=0.9,
          dispersion_correction=False, ewald_tolerance=5e-4, seed=None,
          allow_unverified=False, packer_options=None, cation_anion_sigma_scale=None,
          solute_xy_nm=None, height_window=None, plane_atoms=None):
    """graphite (a GraphiteBuilder, or its CarbonModel) + solute + packed solvent
    -> ElectrodeSystem. The graphite is used as built; nothing here changes it.

    species: SPECIES keys, proportions: molecule ratios for SolvationBuilder.custom_solvate.
    Default: pure water (packing.pure_water); pass e.g. packing.koh_mixture(1.0) for KOH.
    liquid_density (kg/m^3) is the first-estimate density of the liquid region;
    calibrate it with md.calibrate.
    counterion: element SolvationBuilder adds only to neutralise a charged solute or
    charged surface groups ('K', 'Na' or 'Cl'); in pure water these are the only ions.
    solute_xy_nm: lateral position of the solute COM (default: box centre).
    height_window: dict(leg, lo, hi, k) adds the flat-bottom window on the solute height
    above the top layer (add_height_window) for the adsorption legs.
    """
    species, proportions = _electrolyte(species, proportions)
    solute_pdb, solute_xml = str(solute_pdb), str(solute_xml)   # OpenMM readers need str paths
    gb = GraphiteBuilder.from_model(graphite)
    slab = gb.model
    if gb.kind != 'slab':
        raise ValueError(f'build() needs a periodic slab, not a {gb.kind}; '
                         'use export.write_carbon for flakes and nanotubes')
    if box_z_nm - slab.thickness < 2.0:
        raise ValueError('box_z_nm leaves less than 2 nm of liquid between the slab faces')
    if seed is not None:
        np.random.seed(seed)       # SolvationBuilder draws from numpy's global RNG
    sol_top, sol_xyz, sol_mass, sol_q, sol_labels = _read_solute(solute_pdb, solute_xml)
    placed = place_solute(sol_xyz, sol_mass, slab, height_nm, flip, rotation, xy=solute_xy_nm,
                          plane=plane_atoms)
    packing_solute = combined_molecule(sol_labels, placed, slab, sol_q)
    mols, warnings = species_molecules(species, allow_unverified)
    warnings += [w for w in gb.warnings() if w not in warnings]
    surface = gb.has_surface_features
    total_density, v_liquid = total_target_density(liquid_density, slab, sol_mass.sum(),
                                                   box_z_nm, wall_gap_nm)
    packer = make_packer(len(sol_labels), surface=surface, **packer_kwargs(species, packer_options))
    set_counterions(packer, counterion, species)
    box_a = [10 * slab.box_xy[0], 10 * slab.box_xy[1], 10 * box_z_nm]
    packer.custom_solvate(solute=packing_solute, solvents=mols, proportion=list(proportions),
                          box_size=box_a, target_density=total_density, neutralize=True,
                          strict=True)
    es = assemble_from_packer(packer, gb, sol_top, solute_xml, species,
                              field_V_per_nm=field_V_per_nm, cutoff_nm=cutoff_nm,
                              switch_nm=switch_nm, dispersion_correction=dispersion_correction,
                              ewald_tolerance=ewald_tolerance, allow_unverified=allow_unverified,
                              cation_anion_sigma_scale=cation_anion_sigma_scale)
    es.solute_pdb = str(solute_pdb)
    es.metadata['warnings'] = warnings + es.metadata['warnings']
    es.metadata['liquid'].update(
        target_liquid_density_kg_m3=liquid_density, wall_gap_nm=wall_gap_nm,
        effective_volume_nm3=v_liquid, packer_target_total_density_kg_m3=total_density,
        proportions=list(proportions))
    es.metadata['solute'].update(height_nm=height_nm, flip=flip, xy_nm=solute_xy_nm, plane_atoms=plane_atoms,
                                 rotation=None if rotation is None else np.asarray(rotation).tolist(),
                                 source_pdb=str(solute_pdb), source_xml=str(solute_xml))
    es.metadata['packing'] = dict(seed=seed, threshold_A=packer.threshold,
                                  **getattr(packer, 'stats', {}))
    if height_window:
        w = height_window
        masses = [es.system.getParticleMass(i).value_in_unit(unit.dalton) for i in es.indices['solute']]
        add_height_window(es.system, es.indices['solute'], es.indices['top_layer'], masses,
                          w['lo'], w['hi'], w['k'])
        es.metadata['restraint'] = dict(
            leg=w.get('leg'), kind='flat-bottom window on solute COM height above the top layer',
            lo_nm=w['lo'], hi_nm=w['hi'], k_kj_mol_nm2=w['k'],
            accessible_length_nm_298K=window_length(w['lo'], w['hi'], w['k']),
            force_group=FORCE_GROUPS['restraint'])
    return es


def build_reference(solute_pdb, solute_xml, species=None, proportions=None, *, box_nm=5.0,
                    liquid_density=997.0, counterion='K', cutoff_nm=1.0, switch_nm=0.9,
                    dispersion_correction=False, ewald_tolerance=5e-4, seed=None,
                    allow_unverified=False):
    """The solute in a cubic box of the electrolyte: no graphite, no field. This is the
    bulk-water reference FEP, the system the driver's 'auto' ion corrections are made for."""
    import veloxchem as vlx
    species, proportions = _electrolyte(species, proportions)
    solute_pdb, solute_xml = str(solute_pdb), str(solute_xml)
    if seed is not None:
        np.random.seed(seed)
    sol_top, sol_xyz, sol_mass, sol_q, sol_labels = _read_solute(solute_pdb, solute_xml)
    mol = vlx.Molecule(sol_labels, 10.0 * plane_frame(sol_xyz, sol_mass), 'angstrom')
    mol.set_charge(float(sol_q))
    mol.set_multiplicity(1)
    mols, warnings = species_molecules(species, allow_unverified)
    packer = make_packer(len(sol_labels), **packer_kwargs(species))
    set_counterions(packer, counterion, species)
    packer.custom_solvate(solute=mol, solvents=mols, proportion=list(proportions),
                          box_size=[10 * box_nm] * 3, target_density=liquid_density,
                          neutralize=True, strict=True)
    es = assemble_from_packer(packer, None, sol_top, solute_xml, species, field_V_per_nm=0.0,
                              cutoff_nm=cutoff_nm, switch_nm=switch_nm,
                              dispersion_correction=dispersion_correction,
                              ewald_tolerance=ewald_tolerance, allow_unverified=allow_unverified)
    es.solute_pdb = solute_pdb
    es.metadata['warnings'] = warnings + es.metadata['warnings']
    es.metadata['liquid'].update(target_liquid_density_kg_m3=liquid_density,
                                 proportions=list(proportions))
    es.metadata['solute'].update(source_pdb=solute_pdb, source_xml=solute_xml)
    es.metadata['reference'] = dict(kind='bulk', box_nm=box_nm)
    es.metadata['packing'] = dict(seed=seed, threshold_A=packer.threshold,
                                  **getattr(packer, 'stats', {}))
    return es


def _electrolyte(species, proportions):
    if species is None:
        if proportions is not None:
            raise ValueError('proportions given without species')
        return pure_water()
    if proportions is None:
        raise ValueError('species given without proportions')
    return list(species), list(proportions)


def build_bulk(species=None, proportions=None, *, box_nm=4.0, liquid_density=997.0, counterion='K',
               cutoff_nm=1.0, switch_nm=0.9, dispersion_correction=False,
               ewald_tolerance=5e-4, seed=None, allow_unverified=False,
               cation_anion_sigma_scale=None):
    """Bulk box of the same solvent model (no graphite, no solute), e.g. for a density
    reference. Default: pure water."""
    import veloxchem as vlx
    species, proportions = _electrolyte(species, proportions)
    if seed is not None:
        np.random.seed(seed)
    mols, warnings = species_molecules(species, allow_unverified)
    first = vlx.Molecule(mols[0].get_labels(), mols[0].get_coordinates_in_angstrom(), 'angstrom')
    first.set_charge(mols[0].get_charge())
    first.set_multiplicity(1)
    packer = make_packer(first.number_of_atoms(), **packer_kwargs(species))
    set_counterions(packer, counterion if counterion in ('K', 'Na') else None, species)
    packer.custom_solvate(solute=first, solvents=mols, proportion=list(proportions),
                          box_size=[10 * box_nm] * 3, target_density=liquid_density,
                          neutralize=True, strict=True)
    es = assemble_from_packer(packer, None, None, None, species, cutoff_nm=cutoff_nm,
                              switch_nm=switch_nm, dispersion_correction=dispersion_correction,
                              ewald_tolerance=ewald_tolerance, allow_unverified=allow_unverified,
                              cation_anion_sigma_scale=cation_anion_sigma_scale)
    es.metadata['warnings'] = warnings + es.metadata['warnings']
    es.metadata['liquid'].update(target_liquid_density_kg_m3=liquid_density,
                                 proportions=list(proportions))
    return es


# ---------------------------------------------------------------- assembly

def assemble_from_packer(packer, graphite, solute_top, solute_xml, species, *, field_V_per_nm=0.0,
                         cutoff_nm=1.0, switch_nm=0.9, dispersion_correction=False,
                         ewald_tolerance=5e-4, allow_unverified=False,
                         cation_anion_sigma_scale=None):
    """OpenMM system from a packed SolvationBuilder.

    The packer's 'solute' molecule is [solute atoms, graphite atoms] (graphite given),
    or one molecule of species[0] (bulk: graphite None and solute_top None).
    """
    gb = None if graphite is None else GraphiteBuilder.from_model(graphite)
    slab = gb.model if gb is not None else None
    coords = np.asarray(packer.system_molecule.get_coordinates_in_angstrom()) / 10.0
    labels = list(packer.system_molecule.get_labels())
    box = np.asarray(packer.box, float) / 10.0
    n_sol = 0 if solute_top is None else solute_top.getNumAtoms()
    n_c = 0 if slab is None else slab.n_atoms
    specs = [SPECIES[k] for k in species]

    # ---- molecule blocks in packer order: [mol 0] [species blocks] [counterions]
    blocks = []                     # (species key, natoms, list of coordinate arrays)
    off = n_sol + n_c
    if solute_top is None and slab is None:
        n0 = len(specs[0]['atom_names'])
        first = [coords[:n0]]
        off = n0
    for k, (key, mol, count) in enumerate(zip(species, packer.solvents, packer.added_solvent_counts)):
        na = mol.number_of_atoms()
        mols = [coords[off + i * na: off + (i + 1) * na] for i in range(count)]
        if k == 0 and solute_top is None and slab is None:
            mols = first + mols
        blocks.append((key, na, mols, mol))
        off += count * na
    n_ion = int(getattr(packer, 'added_counterions', 0) or 0)
    ion_key = None
    if n_ion:
        ion_name = getattr(packer, 'ion_name', None)
        if ion_name not in COUNTERION_SPECIES:
            ion_name = packer.counterion.get_labels()[0]
        ion_key = COUNTERION_SPECIES[ion_name]
        na = packer.counterion.number_of_atoms()
        blocks.append((ion_key, na, [coords[off + i * na: off + (i + 1) * na] for i in range(n_ion)],
                       packer.counterion))
        off += n_ion * na
    if off != len(coords):
        raise RuntimeError(f'packed atom bookkeeping mismatch: {off} vs {len(coords)}')

    # ---- place the slab mid-plane on z = 0 and wrap molecules whole
    shift = np.zeros(3)
    if slab is not None:
        shift[2] = slab.xyz[:, 2].mean() - coords[n_sol:n_sol + n_c, 2].mean()   # back to model z
    sol_xyz = coords[:n_sol] + shift
    if n_sol:
        com = sol_xyz.mean(0)
        sol_xyz[:, :2] -= np.floor(com[:2] / box[:2]) * box[:2]
    # every carbon-model atom inside [0, L): the bottom half of the slab, and anything on the
    # bottom face, sits just below Lz next to the liquid it touches. (At z < 0 a GPU platform
    # folds charged surface atoms back into the box, which shifts -q*E*z by q*E*Lz.)
    c_xyz = coords[n_sol:n_sol + n_c] + shift
    if slab is not None:
        c_xyz[:, 2] = slab.xyz[:, 2]     # exact model z: a layer at z = 0 must not round to Lz - 1e-16
    c_xyz = _wrap(c_xyz, box)
    wrapped = []
    for key, na, mols, mol in blocks:
        wrapped.append([(m + shift) - np.floor((m[0] + shift) / box) * box for m in mols])

    # ---- topology
    top = app.Topology()
    pos = []
    idx = {k: [] for k in ('solute', 'carbon', 'top_layer', 'bottom_layer', 'water', 'water_O',
                           'water_H', 'cation', 'cation_center', 'anion', 'anion_center', 'carbon_model',
                           'carbon_h', 'surface_groups')}
    if n_sol:
        ch = top.addChain('S')
        amap = {}
        for res in solute_top.residues():
            r = top.addResidue(res.name, ch)
            for a in res.atoms():
                amap[a.index] = top.addAtom(a.name, a.element, r)
                idx['solute'].append(amap[a.index].index)
        for b in solute_top.bonds():
            top.addBond(amap[b[0].index], amap[b[1].index])
        pos += list(sol_xyz)
        if labels[:n_sol] != [a.element.symbol for a in solute_top.atoms()]:
            raise RuntimeError('solute atom order changed during packing')
    if n_c:
        ch = top.addChain('G')
        elem = {'C': app.element.carbon, 'H': app.element.hydrogen, 'O': app.element.oxygen}
        skel = slab.skeleton
        for k, (p, lay, e, rn) in enumerate(zip(c_xyz, slab.layer, slab.elements, gb.residue_names())):
            a = top.addAtom(e, elem[e], top.addResidue(rn, ch))
            pos.append(p)
            idx['carbon_model'].append(a.index)
            if not skel[k]:
                idx['surface_groups'].append(a.index)
            elif e != 'C':
                idx['carbon_h'].append(a.index)
            else:
                idx['carbon'].append(a.index)
                if lay == 0:
                    idx['top_layer'].append(a.index)
                if lay == slab.layers - 1:
                    idx['bottom_layer'].append(a.index)
    counts = {}
    for (key, na, mols, mol), xyzs in zip(blocks, wrapped):
        spec = SPECIES[key]
        ch = top.addChain('I' if key == ion_key and mol is getattr(packer, 'counterion', None) else 'W')
        conn = np.asarray(mol.get_connectivity_matrix())
        elems = [app.element.get_by_symbol(s) for s in mol.get_labels()]
        q = mol.get_charge()
        for m in xyzs:
            r = top.addResidue(spec['resname'], ch)
            atoms = [top.addAtom(n, e, r) for n, e in zip(spec['atom_names'], elems)]
            for i in range(na):
                for j in range(i + 1, na):
                    if conn[i, j]:
                        top.addBond(atoms[i], atoms[j])
            pos += list(m)
            ids = [a.index for a in atoms]
            if spec['resname'] == 'HOH':
                idx['water'] += ids
                idx['water_O'].append(ids[0])
                idx['water_H'] += ids[1:]
            elif q > 0:
                idx['cation'] += ids
                idx['cation_center'].append(ids[0])
            elif q < 0:
                idx['anion'] += ids
                idx['anion_center'].append(ids[0])
        counts[key] = counts.get(key, 0) + len(xyzs)
    top.setPeriodicBoxVectors([mm.Vec3(box[0], 0, 0), mm.Vec3(0, box[1], 0),
                               mm.Vec3(0, 0, box[2])] * NM)
    positions = np.array(pos)

    # ---- force field and system
    sources, seen, warnings = [], set(), []
    if solute_xml:
        sources.append(('file', str(solute_xml)))
    if slab is not None:
        sources.append(('xml', gb.forcefield_xml()))
    species_keys = list(dict.fromkeys(k for k, *_ in blocks))
    for key in species_keys:
        s = _ff_source(SPECIES[key])
        if s not in seen:
            seen.add(s)
            sources.append(s)
    ff = app.ForceField(*[io.StringIO(v) if t == 'xml' else v for t, v in sources])
    system = ff.createSystem(top, nonbondedMethod=app.PME, nonbondedCutoff=cutoff_nm * NM,
                             constraints=app.HBonds, rigidWater=True,
                             ewaldErrorTolerance=ewald_tolerance, removeCMMotion=False)
    for res in top.residues():        # rigid ions (H3O+): H-H constraints on top of HBonds
        spec = next((SPECIES[k] for k in species_keys if SPECIES[k]['resname'] == res.name), None)
        if spec and spec.get('rigid_HH'):
            hs = [a.index for a in res.atoms() if a.element.symbol == 'H']
            for i, j in itertools.combinations(hs, 2):
                system.addConstraint(i, j, spec['rigid_HH'] * NM)
    nb = nonbonded_force(system)
    nb.setUseDispersionCorrection(dispersion_correction)
    if switch_nm:
        nb.setUseSwitchingFunction(True)
        nb.setSwitchingDistance(switch_nm * NM)
    for f in system.getForces():   # (SWIG proxies are new objects per call: test the type, not identity)
        f.setForceGroup(FORCE_GROUPS['nonbonded'] if isinstance(f, mm.NonbondedForce)
                        else FORCE_GROUPS['bonded'])
    mobility = gb.apply_to_system(system, idx['carbon_model'], positions) if n_c else None
    if gb is not None:
        gb.add_water_correction(system, idx['carbon'], idx['water_O'], idx['water_H'],
                                cutoff_nm, switch_nm)
    add_uniform_field(system, field_V_per_nm)
    hydroxides = {SPECIES[k]['ffxml'] for k, *_ in blocks if SPECIES[k]['ffxml'].startswith('hydroxide:')}
    if cation_anion_sigma_scale is not None and idx['cation'] and idx['anion_center']:
        add_ion_pair_sigma_scaling(system, idx['cation'], idx['anion_center'],
                                   cation_anion_sigma_scale, cutoff_nm, switch_nm)
    elif 'hydroxide:bonthuis2016' in hydroxides and idx['cation']:
        warnings.append('Bonthuis 2016 OH- with cations and plain Lorentz-Berthelot cation-anion LJ: '
                        'the paper needed sigma_+- = 1.5*sqrt(sigma_+ sigma_-) for Smith-Dang Na+ and '
                        'never fitted K+; cation-OH- association is unvalidated '
                        '(set cation_anion_sigma_scale, ideally fitted to KOH activity data).')
    if 'hydroxide:bonthuis2016' in hydroxides:
        warnings.append('Bonthuis 2016 OH- was fitted with geometric sigma combination; OpenMM '
                        'averages sigma (OH- O / water O 0.3488 vs 0.3473 nm).')

    total_q = sum(nb.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge)
                  for i in range(system.getNumParticles()))
    if abs(total_q) > 1e-4:
        raise RuntimeError(f'system is not neutral: {total_q:.6f} e')
    n_close = _close_contacts(top, positions, box, 0.14)
    if n_close:
        warnings.append(f'{n_close} atom pairs from different molecules closer than 0.14 nm; '
                        'minimise before dynamics')
    min_gap = min_dist = None
    if slab is not None:
        heavy = [i for i in idx['water_O'] + idx['cation_center'] + idx['anion_center']]
        dz = np.abs(positions[heavy, 2] - np.round(positions[heavy, 2] / box[2]) * box[2])
        min_gap = float(dz.min() - slab.top_z) if heavy else None
        if heavy:
            tree = cKDTree(_wrap(positions[idx['carbon_model']], box), boxsize=box)
            min_dist = float(tree.query(_wrap(positions[heavy], box))[0].min())

    try:
        import veloxchem as vlx
        vlx_version = getattr(vlx, '__version__', '?')
    except Exception:   # noqa: BLE001
        vlx_version = None
    meta = dict(
        builder='electrode_system (graphite from graphite_builder)', version=3,
        created=datetime.datetime.now().isoformat(timespec='seconds'),
        openmm=mm.__version__, veloxchem=vlx_version,
        box_nm=box.tolist(),
        graphite=dict(layers=slab.layers if slab else 0, cells=list(slab.cells) if slab else None,
                      n_carbon=slab.n_carbon if slab else 0, thickness_nm=slab.thickness if slab else None,
                      top_face_z=slab.top_z if slab else None,
                      bottom_face_z=float(box[2] - slab.top_z) if slab else None,
                      carbon_set=slab.carbon_set if slab else None,
                      carbon_parameters=CARBON_SETS[slab.carbon_set] if slab else None,
                      stacking=slab.info.get('stacking') if slab else None,
                      orientation=slab.info.get('orientation') if slab else None,
                      mobility=mobility,
                      nearest_solvent_heavy_atom_from_face_nm=min_gap,
                      nearest_solvent_heavy_atom_to_carbon_model_atom_nm=min_dist,
                      n_atoms=n_c,
                      builder=gb.metadata() if gb is not None else None),
        liquid=dict(mid_z=float(box[2] / 2), species_counts=counts,
                    n_water=len(idx['water_O']), counterion=ion_key, n_counterions=n_ion,
                    n_slots=len(idx['water_O']) + len(idx['anion_center']) + len(idx['cation_center']),
                    electrolyte_species=list(species)),
        solute=dict(present=bool(n_sol), n_atoms=n_sol,
                    resname=next(solute_top.residues()).name if n_sol else None,
                    alchemical_region='residue index 0' if n_sol else None),
        species={k: SPECIES[k] for k in dict.fromkeys(k for k, *_ in blocks)},
        field=dict(V_per_nm=field_V_per_nm, force_name=FIELD_FORCE_NAME,
                   Ez_kJ_per_mol_nm_e=field_V_per_nm * FARADAY_KJ_PER_MOL_PER_V,
                   convention='E > 0 points along +z: anions are pushed to -z. The top face '
                              '(z = +t/2, liquid above it) is anode-like, the bottom face '
                              '(z = Lz - t/2) cathode-like.'),
        force_groups=FORCE_GROUPS,
        cation_anion_sigma_scale=cation_anion_sigma_scale,
        nonbonded=dict(cutoff_nm=cutoff_nm, switch_nm=switch_nm,
                       dispersion_correction=dispersion_correction,
                       ewald_tolerance=ewald_tolerance, rigid_water=True, constraints='HBonds'),
        n_particles=system.getNumParticles(), total_charge=total_q,
        warnings=warnings,
        rules=['Run at fixed volume whenever graphite or a field is present, unless md.pressure_coupling=\'z\' '
               '(rigid slab at z = 0, field-scaling check on the platform in use).',
               'Mobile (restrained/flexible) carbon: keep the time step at or below 2 fs.',
               'Calibrate liquid_density against a bulk reference (md.calibrate) before production.',
               'Change one thing per labelled run.'])
    return ElectrodeSystem(system, top, positions, idx, meta,
                           [(t, v if t == 'file' else '<generated xml>') for t, v in sources],
                           str(solute_xml) if solute_xml else None, None)


def _wrap(x, box):
    """Coordinates in [0, box): np.mod can return box itself after rounding."""
    y = np.mod(x, box)
    return np.where(y >= box, 0.0, y)


def _close_contacts(top, positions, box, cutoff):
    mol_of = np.empty(top.getNumAtoms(), int)
    for k, res in enumerate(top.residues()):
        for a in res.atoms():
            mol_of[a.index] = k
    tree = cKDTree(_wrap(positions, box), boxsize=box)
    pairs = tree.query_pairs(cutoff, output_type='ndarray')
    if len(pairs) == 0:
        return 0
    carbon = np.array([a.residue.chain.id == 'G' for a in top.atoms()])   # the carbon model
    keep = (mol_of[pairs[:, 0]] != mol_of[pairs[:, 1]]) & ~(carbon[pairs[:, 0]] & carbon[pairs[:, 1]])
    return int(keep.sum())


# ------------------------------------------------------------------- output

def write_fep_inputs(es, outdir, overwrite=False):
    """Files for SolvationFepDriver.compute_solvation_from_openmm_files(
    solute_pdb='solute.pdb', solute_xml='solute.xml', system_pdb='system.pdb',
    system_xml='system.xml', other_xml_files=[])."""
    out = Path(outdir)
    if out.exists() and any(out.iterdir()) and not overwrite:
        raise FileExistsError(f'{out} is not empty')
    out.mkdir(parents=True, exist_ok=True)
    (out / 'system.xml').write_text(mm.XmlSerializer.serialize(es.system))
    with open(out / 'system.pdb', 'w') as fh:
        app.PDBFile.writeFile(es.topology, es.positions * NM, fh)
    if es.solute_pdb:
        shutil.copyfile(es.solute_pdb, out / 'solute.pdb')
    if es.solute_xml:
        shutil.copyfile(es.solute_xml, out / 'solute.xml')
    if es.metadata['graphite']['carbon_set']:
        (out / 'graphite.xml').write_text(graphite_forcefield_xml(es.metadata['graphite']['carbon_set']))
    for key, spec in es.metadata['species'].items():
        if spec['ffxml'].startswith('hydroxide:'):
            (out / f'{key}.xml').write_text(hydroxide_xml(spec['ffxml'].split(':', 1)[1]))
        elif spec['ffxml'].startswith('hydronium:'):
            (out / f'{key}.xml').write_text(hydronium_xml(spec['ffxml'].split(':', 1)[1]))
        elif Path(spec['ffxml']).is_file():          # generated species (HSO4-), for the record
            shutil.copyfile(spec['ffxml'], out / f'{key}.xml')
    (out / 'metadata.json').write_text(json.dumps(es.metadata, indent=1, default=str))
    (out / 'indices.json').write_text(json.dumps(es.indices))
    return out


# =========================================================== FEP

def check_gsc_sigma(system, max_sigma_nm=0.6):
    """Refuse an FEP window whose Gaussian soft-core force carries a placeholder sigma (review
    2026-11-18). OpenMM's amber14 SPC/E water H has sigma = 1 nm with epsilon 0; harmless for
    LJ, but a VeloxChem build whose GSC term copies NonbondedForce sigma then gives every
    water H a huge soft-core repulsion (zero overlap). The old 3-stage driver replaced such
    sigmas by 0.1 nm. Returns the number of GSC forces checked."""
    n = 0
    for f in system.getForces():
        if not isinstance(f, mm.CustomNonbondedForce) or f.getNumPerParticleParameters() != 1:
            continue
        e = f.getEnergyFunction()
        if not all(t in e for t in ('lambda', 'alpha', 'beta', 'exp(')):
            continue
        n += 1
        big = max(float(f.getParticleParameters(k)[0]) for k in range(f.getNumParticles()))
        if big > max_sigma_nm:
            raise RuntimeError(f'GSC force {f.getName()!r} has sigma {big:.3f} nm (> {max_sigma_nm}): '
                               'placeholder sigma reached the soft-core term; use the validated '
                               'driver (environment.sh) or the GSC fix')
    return n


FEP_DRIVER_ENV = 'ELECTRODE_FEP_DRIVER'                # path of a solvationfepdriver.py to use instead
FEP_DRIVER_SHA_ENV = 'ELECTRODE_FEP_DRIVER_SHA256'     # required SHA-256 of the driver source (pinning)


def fep_driver_base():
    """The SolvationFepDriver class to subclass, its source file and SHA-256 (review 2026-11-18).

    Default: the installed veloxchem.SolvationFepDriver. With ELECTRODE_FEP_DRIVER set, that
    file is loaded as a module of the veloxchem package instead (tests use the pinned
    production driver, third_party/veloxchem/solvationfepdriver.py). With
    ELECTRODE_FEP_DRIVER_SHA256 set, any other driver source is refused."""
    import hashlib
    import importlib.util
    import inspect
    import veloxchem as vlx
    path = os.environ.get(FEP_DRIVER_ENV)
    if path:
        name = 'veloxchem._pinned_solvationfepdriver'
        mod = sys.modules.get(name)
        if mod is None or getattr(mod, '__file__', None) != str(Path(path).resolve()):
            spec = importlib.util.spec_from_file_location(name, str(Path(path).resolve()))
            mod = importlib.util.module_from_spec(spec)
            mod.__package__ = 'veloxchem'
            sys.modules[name] = mod
            spec.loader.exec_module(mod)
        cls, src = mod.SolvationFepDriver, str(Path(path).resolve())
    else:
        cls = vlx.SolvationFepDriver
        src = inspect.getsourcefile(cls)
    sha = hashlib.sha256(Path(src).read_bytes()).hexdigest()
    want = os.environ.get(FEP_DRIVER_SHA_ENV)
    if want and sha != want.strip().lower():
        raise RuntimeError(f'FEP driver {src} has SHA-256 {sha}, expected {want} ({FEP_DRIVER_SHA_ENV})')
    return cls, src, sha


def electrode_fep_driver(expected_solute=None, carbon=None, mobility='frozen', pressure_coupling=None):
    """Return an instance of the electrode-aware SolvationFepDriver subclass.

    expected_solute / carbon: atom indices (e.g. from indices.json) to verify.
    mobility: the carbon mode in metadata['graphite']['mobility']['mode'].
    pressure_coupling: None (fixed volume: any barostat in a window is an error) or 'z'
    (metadata['pressure_coupling']: the windows keep the system's z-only barostat).
    """
    import veloxchem as vlx
    base, driver_src, driver_sha = fep_driver_base()

    class ElectrodeFepDriver(base):

        def _create_solvated_systems(self, lambdas):
            systems, topology, positions = super()._create_solvated_systems(lambdas)
            solute_ids, _ = self._get_alchemical_region(topology)
            solute_set = set(solute_ids)
            if expected_solute is not None and sorted(solute_ids) != sorted(expected_solute):
                raise RuntimeError('FEP alchemical region does not match the solute atoms')
            for lam, system in zip(lambdas, systems):
                baro = [f for f in system.getForces() if isinstance(
                    f, (mm.MonteCarloBarostat, mm.MonteCarloAnisotropicBarostat, mm.MonteCarloMembraneBarostat))]
                if pressure_coupling == 'z':
                    if len(baro) != 1 or not isinstance(baro[0], mm.MonteCarloAnisotropicBarostat) or (
                            baro[0].getScaleX() or baro[0].getScaleY() or not baro[0].getScaleZ()):
                        raise RuntimeError('z-NPT FEP windows need exactly one z-only anisotropic barostat')
                elif baro:
                    raise RuntimeError('electrode FEP windows must run at fixed volume')
                if carbon is not None and mobility == 'frozen' and any(
                        system.getParticleMass(i)._value != 0 for i in carbon):
                    raise RuntimeError('graphite is no longer frozen in an FEP window')
                if carbon is not None and mobility != 'frozen' and any(
                        system.getParticleMass(i)._value == 0 for i in carbon):
                    raise RuntimeError(f'{mobility} carbon has zero mass in an FEP window')
                n_gsc = check_gsc_sigma(system)
                # record how many soft-core forces were actually inspected (review v2 follow-up)
                self.electrode_gsc_checked = getattr(self, 'electrode_gsc_checked', {})
                self.electrode_gsc_checked.setdefault(int(self.stage), []).append(n_gsc)
                factor = 1.0 - float(lam) if self.stage == 1 else 0.0
                for f in system.getForces():
                    if isinstance(f, mm.CustomExternalForce) and f.getName() == FIELD_FORCE_NAME:
                        for k in range(f.getNumParticles()):
                            p, prm = f.getParticleParameters(k)
                            if p in solute_set:
                                f.setParticleParameters(k, p, [float(prm[0]) * factor])
            return systems, topology, positions

        def _create_vacuum_systems(self, lambdas):
            if self.solvent_name == 'omm_system_xml':
                self.solvent_name = 'omm_files'
                try:
                    return super()._create_vacuum_systems(lambdas)
                finally:
                    self.solvent_name = 'omm_system_xml'
            return super()._create_vacuum_systems(lambdas)

    d = ElectrodeFepDriver()
    d.electrode_driver_info = dict(source_file=driver_src, sha256=driver_sha,
                                   veloxchem_version=getattr(vlx, '__version__', None),
                                   openmm_version=mm.__version__)
    return d


def fep_window_systems(driver, lambdas):
    """driver._create_solvated_systems(lambdas) outside a full run, for checks. Some VeloxChem
    builds validate against an endpoint signature that only a full run sets; try with it
    unset, and raise RuntimeError('driver internals differ ...') if that is not enough.
    The production 3-stage driver (third_party/veloxchem, a84e9f76...) sets _active_ion_mode
    in compute_solvation_from_openmm_files before building windows; outside a run it is unset.
    Slab legs always run with ion corrections 'none' (ion_correction_settings), so the check
    sets that same value (review 2026-11-18, v2)."""
    if getattr(driver, '_active_ion_mode', None) is None:
        driver._active_ion_mode = 'none'
    try:
        return driver._create_solvated_systems(lambdas)
    except AttributeError as e:
        if '_endpoint_signature' not in str(e):
            raise
    driver._endpoint_signature = None
    try:
        return driver._create_solvated_systems(lambdas)
    except Exception as e:   # noqa: BLE001
        raise RuntimeError(f'driver internals differ in this VeloxChem build ({type(e).__name__}: {e}); '
                           'the full FEP run is not affected') from e


def ion_correction_settings(driver, meta, settings):
    """Set the driver's ion corrections for this system; returns (extra keyword arguments for
    compute_solvation_from_openmm_files, description). Slab: 'none'. Bulk reference in pure
    SPC/E: 'auto' with water_model='spce'; in an electrolyte: 'none'. Explicit settings win;
    builds without the feature say so."""
    if 'ion_corrections' in settings:
        return {}, f"{settings['ion_corrections']} (set explicitly)"
    if not hasattr(driver, 'ion_corrections'):
        return {}, 'unavailable in this VeloxChem build'
    if (meta.get('graphite') or {}).get('layers', 0) > 0:
        driver.ion_corrections = 'none'
        return {}, 'none (slab system: bulk-water presets do not apply)'
    mixture = (meta.get('liquid') or {}).get('electrolyte_species') or []
    if any(k != 'water_spce' for k in mixture):
        driver.ion_corrections = 'none'
        return {}, ('none (bulk reference in an electrolyte: the presets are for one ion in pure '
                    'water; the ions screen the solute charge)')
    if 'water_spce' in (meta.get('species') or {}):
        return {'water_model': 'spce'}, 'auto spce (bulk reference)'
    return {}, 'auto (bulk reference, no water preset declared: uncorrected)'


def run_fep(folder, **settings):
    """Run the four-stage solvation FEP on a folder written by assemble.write_fep_inputs.

    settings are driver attributes (num_steps, num_snapshots, num_equil_steps,
    temperature, platform, output_folder, save_trajectory_xtc, ...).
    """
    import os
    folder = Path(folder)
    idx = json.loads((folder / 'indices.json').read_text())
    meta = json.loads((folder / 'metadata.json').read_text())
    mode = ((meta.get('graphite') or {}).get('mobility') or {}).get('mode', 'frozen')
    driver = electrode_fep_driver(idx['solute'], idx['carbon'] or None, mode,
                                  pressure_coupling=meta.get('pressure_coupling'))
    for k, v in settings.items():
        if not hasattr(driver, k):
            raise AttributeError(f'SolvationFepDriver has no attribute {k!r}')
        setattr(driver, k, v)
    # Ion corrections (VeloxChem builds that have them): their presets are for a solute in
    # homogeneous bulk water in a cubic box. Slab systems: none (adsorption comes from two
    # legs in the same box instead). Bulk reference in SPC/E: auto with water_model='spce'.
    kwargs, applied = ion_correction_settings(driver, meta, settings)
    cwd = os.getcwd()
    os.chdir(folder)
    try:
        results = driver.compute_solvation_from_openmm_files(
            solute_pdb='solute.pdb', solute_xml='solute.xml', system_pdb='system.pdb',
            system_xml='system.xml', other_xml_files=[], **kwargs)
    finally:
        os.chdir(cwd)
    results['electrode_ion_corrections'] = applied
    results['fep_driver'] = getattr(driver, 'electrode_driver_info', None)    # source file + SHA-256
    (folder / 'fep_results.json').write_text(json.dumps(results, indent=2, default=float))
    return results


# =========================================================== MD helpers

def get_platform(name=None):
    names = [mm.Platform.getPlatform(i).getName() for i in range(mm.Platform.getNumPlatforms())]
    for n in ([name] if name else ['CUDA', 'OpenCL', 'CPU']):
        if n not in names:
            continue
        props = {'Precision': 'mixed'} if n in ('CUDA', 'OpenCL') else {}
        plat = mm.Platform.getPlatformByName(n)
        try:
            probe = mm.System()
            probe.addParticle(1.0)
            mm.Context(probe, mm.VerletIntegrator(0.001), plat, props)
        except Exception:   # noqa: BLE001
            continue
        return plat, props
    raise RuntimeError('no usable OpenMM platform')


RIGID_SLAB_FORCE_NAME = 'carbon slab as one molecule (zero-strength bonds)'


def make_slab_rigid_for_barostat(es):
    """Prepare a frozen-carbon slab for a z barostat (in place).

    OpenMM's Monte Carlo barostats rescale the centre of every molecule. Unbonded carbon
    atoms are each their own molecule (zero mass does not exempt them), so the slab would be
    rescaled with the box (its layer spacing follows Lz), and a slab stored across the periodic
    boundary would be torn open. Fix: join all carbon atoms into one molecule with
    zero-strength bonds (no energy, no force) and store them contiguous around z = 0, where
    the builder centres the slab; its centre is then ~0 and the barostat leaves it in place
    (tested on OpenCL: thickness 1.0050 nm unchanged at 1000 bar)."""
    carbon = list(es.indices.get('carbon') or [])
    if not carbon:
        return es
    if any(f.getName() == RIGID_SLAB_FORCE_NAME for f in es.system.getForces()):
        return es
    lz = es.metadata['box_nm'][2]
    x = np.array(es.positions, float)
    c = np.array(carbon)
    x[c[x[c, 2] > lz / 2], 2] -= lz
    es.positions = x
    hb = mm.HarmonicBondForce()
    hb.setName(RIGID_SLAB_FORCE_NAME)
    for i, j in zip(c[:-1], c[1:]):
        hb.addBond(int(i), int(j), 0.0, 0.0)
    hb.setForceGroup(FORCE_GROUPS['bonded'])
    es.system.addForce(hb)
    return es


def z_barostat(temperature=298.15, pressure_bar=1.0, frequency=25, seed=1):
    """One normal-pressure barostat: x and y fixed, z free (as in the water_4V sweeps)."""
    b = mm.MonteCarloAnisotropicBarostat(mm.Vec3(pressure_bar, pressure_bar, pressure_bar) * unit.bar,
                                         temperature * unit.kelvin, False, False, True, frequency)
    b.setRandomNumberSeed(seed)
    return b


def check_field_scaling(sim, es, scale=1.001, tol_kj_mol=0.05):
    """The field energy -q E z must change by -E sum q dz when the box is rescaled in z on the
    platform in use. GPU platforms evaluate it on wrapped coordinates, so a charged atom
    across the z boundary gives a jump of ~q E Lz and corrupts every barostat move. With the
    slab across the boundary the liquid never wraps; this checks it on the actual platform."""
    fg = es.metadata.get('force_groups', {}).get('field')
    E = es.metadata['field'].get('Ez_kJ_per_mol_nm_e', 0.0)
    if fg is None or not E:
        return None
    ctx = sim.context
    st = ctx.getState(getPositions=True, getEnergy=True, groups={fg})
    x = st.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    a, b, c = st.getPeriodicBoxVectors()
    e0 = st.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    field = [f for f in sim.system.getForces() if f.getName() == FIELD_FORCE_NAME][0]
    q = np.zeros(len(x))
    for k in range(field.getNumParticles()):
        i, prm = field.getParticleParameters(k)
        q[i] = prm[0]
    x1 = x.copy()
    x1[:, 2] *= scale
    ctx.setPeriodicBoxVectors(a, b, c * scale)
    ctx.setPositions(x1 * unit.nanometer)
    e1 = ctx.getState(getEnergy=True, groups={fg}).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    ctx.setPeriodicBoxVectors(a, b, c)
    ctx.setPositions(x * unit.nanometer)
    expected = -E * float((q * (x1[:, 2] - x[:, 2])).sum())
    if abs((e1 - e0) - expected) > tol_kj_mol:
        raise RuntimeError(f'field energy under z scaling: {e1 - e0:+.3f} kJ/mol on '
                           f'{ctx.getPlatform().getName()}, expected {expected:+.3f}: a charged atom wraps '
                           f'in z; do not use a barostat with the field on this system')
    return dict(platform=ctx.getPlatform().getName(), scale=scale, delta_kj_mol=e1 - e0, expected_kj_mol=expected)


def make_simulation(es, temperature=298.15, friction_per_ps=1.0, timestep_fs=2.0, seed=1,
                    platform=None, barostat=False, minimize=True, barostat_z=None):
    meta = es.metadata
    if barostat and (meta['graphite']['layers'] > 0 or meta['field']['V_per_nm'] != 0.0):
        raise ValueError('a barostat is only allowed for bulk boxes at zero field; '
                         'run slab or field systems at fixed volume (or barostat_z=...)')
    if barostat_z is not None and meta['graphite']['layers'] > 0 and not any(
            f.getName() == RIGID_SLAB_FORCE_NAME for f in es.system.getForces()):
        raise ValueError('call make_slab_rigid_for_barostat(es) before a z barostat on a slab')
    system = mm.XmlSerializer.deserialize(mm.XmlSerializer.serialize(es.system))
    if barostat:
        system.addForce(mm.MonteCarloBarostat(1.0 * unit.bar, temperature * unit.kelvin, 25))
    if barostat_z is not None:          # frequency: 0 = off until set (fixed-volume preparation)
        system.addForce(z_barostat(temperature, 1.0, barostat_z, seed))
    integ = mm.LangevinMiddleIntegrator(temperature * unit.kelvin, friction_per_ps / unit.picosecond,
                                        timestep_fs * unit.femtosecond)
    integ.setRandomNumberSeed(seed)
    plat, props = get_platform(platform)
    sim = app.Simulation(es.topology, system, integ, plat, props)
    sim.context.setPositions(es.positions * unit.nanometer)
    if minimize:
        sim.minimizeEnergy(maxIterations=5000)
    sim.context.setVelocitiesToTemperature(temperature * unit.kelvin, seed)
    return sim


def slot_centres(es):
    i = es.indices
    return np.array(i['water_O'] + i.get('cation_center', i['cation']) + i['anion_center'], dtype=int)


def slot_profile(sim, es, n_frames=50, stride=500, bin_nm=0.02, skip=0):
    """Slot density along z (cell coordinates) and its mid-liquid mean.

    Slab: mid region = at least 1 nm from both faces. Bulk: N_slots / <V>.
    """
    centres = slot_centres(es)
    meta = es.metadata
    if skip:
        sim.step(skip)
    lz = meta['box_nm'][2]
    edges = np.arange(0.0, lz + bin_nm, bin_nm)
    hist = np.zeros(len(edges) - 1)
    volumes = []
    for _ in range(n_frames):
        sim.step(stride)
        st = sim.context.getState(getPositions=True)
        pos = st.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        a, b, c = st.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.nanometer)
        volumes.append(a[0] * b[1] * c[2])
        hist += np.histogram(np.mod(pos[centres, 2], c[2]), bins=edges)[0]
    area = meta['box_nm'][0] * meta['box_nm'][1]
    z = 0.5 * (edges[1:] + edges[:-1])
    dens = hist / (n_frames * area * bin_nm)
    if meta['graphite']['layers'] > 0:
        top, bottom = meta['graphite']['top_face_z'], meta['graphite']['bottom_face_z']
        mask = (z > top + 1.0) & (z < bottom - 1.0)
        if not mask.any():
            raise ValueError('liquid too thin for a mid region 1 nm from both faces')
        mid = float(dens[mask].mean())
    else:
        mid = float(len(centres) / np.mean(volumes))
    return dict(z=z, density=dens, mid=mid, mean_volume=float(np.mean(volumes)))


def bulk_reference(species=None, proportions=None, *, box_nm=4.0, equil_ps=200, frames=100, stride=250,
                   platform=None, temperature=298.15, seed=1, **build_kwargs):
    """Slot density of the solvent model in a bulk NPT box at 1 bar, zero field
    (default: pure water)."""
    es = build_bulk(species, proportions, box_nm=box_nm, seed=seed, **build_kwargs)
    sim = make_simulation(es, temperature=temperature, barostat=True, platform=platform, seed=seed)
    prof = slot_profile(sim, es, n_frames=frames, stride=stride, skip=int(equil_ps / 0.002))
    return dict(slot_density=prof['mid'], volume_nm3=prof['mean_volume'],
                n_slots=len(slot_centres(es)), system=es)


def calibrate(build_fn, target_slot_density, liquid_density=997.0, iterations=3, tol=0.005,
              equil_ps=100, frames=50, stride=500, platform=None, temperature=298.15, log=print):
    """Adjust liquid_density until the slab's mid-liquid slot density matches the target.

    build_fn(liquid_density) -> ElectrodeSystem, e.g.
        lambda rho: build(graphite, ..., liquid_density=rho, field_V_per_nm=0.0)
    Run at zero field. Returns (liquid_density, history).
    """
    history = []
    rho = liquid_density
    for it in range(iterations):
        es = build_fn(rho)
        sim = make_simulation(es, temperature=temperature, platform=platform)
        prof = slot_profile(sim, es, n_frames=frames, stride=stride, skip=int(equil_ps / 0.002))
        err = prof['mid'] / target_slot_density - 1
        history.append(dict(iteration=it, liquid_density=rho, n_slots=es.metadata['liquid']['n_slots'],
                            mid_slot_density=prof['mid'], rel_error=err))
        log(f'calibration {it}: liquid_density {rho:.1f} kg/m3, slots {es.metadata["liquid"]["n_slots"]}, '
            f'mid {prof["mid"]:.3f}/nm3 vs {target_slot_density:.3f} ({100 * err:+.2f} %)')
        if abs(err) < tol:
            break
        rho *= target_slot_density / prof['mid']
    return rho, history

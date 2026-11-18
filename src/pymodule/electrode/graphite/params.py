"""Parameter registry for the graphite builder (carbon only; the electrolyte registry
lives in veloxchem.electrode.build_system).

Every parameter set carries a ``status``. 'literature' builds silently,
'legacy' and 'computed' (fitted by this module) build with a recorded note, anything else ('unverified', 'placeholder')
is refused unless ``allow_unverified=True``.

Units: nm, kJ/mol, elementary charge.
"""
from __future__ import annotations

KJ_PER_KCAL = 4.184

CC_BOND = 0.142        # graphite C-C bond, nm
INTERLAYER = 0.335     # Bernal AB interlayer spacing, nm
CARBON_MASS = 12.011

GRAPHITE_RESNAME = 'GRA'

# GAFF2 aromatic parameters (VeloxChem database/gaff-2.11.xml) used for flexible carbon.
# Equilibrium bond/angle are set to the lattice values (0.142 nm, 120 deg) so the ideal
# lattice is the minimum in its periodic box; force constants are GAFF2's.
GAFF2_CA_BOND_K = 316812.48          # kJ/mol/nm^2 (ca-ca, r0 0.1398 nm in GAFF2)
GAFF2_CA_ANGLE_K = 575.7184          # kJ/mol/rad^2 (ca-ca-ca, 120.02 deg)
GAFF2_CA_TORSION_K = 15.167          # kJ/mol, X-ca-ca-X, periodicity 2, phase pi
GAFF2_CH_BOND = (0.1086, 331121.76)  # ca-ha r0 (nm), k
GAFF2_CCH_ANGLE_K = 407.5216         # ca-ca-ha, 119.88 deg
GAFF2_CH_IMPROPER_K = 4.6024        # improper ca-X-X-ha (centre ca), periodicity 2, phase pi
GAFF2_HA = (0.262547852235958, 0.0673624)   # ha sigma (nm), epsilon (kJ/mol)

NBFIX_FORCE_NAME = 'graphite-water NBFIX correction'

# Graphite carbon LJ, applied through NonbondedForce (Lorentz-Berthelot with
# every partner, so the FEP driver's solute scaling also decouples
# solute-graphite). 'pairs' are optional C-water overrides added as a
# separate correction force that never involves the solute.
CARBON_SETS = {
    'steele': dict(
        sigma=0.340, epsilon=0.2328, pairs={}, status='literature',
        source='Steele, Surf. Sci. 36, 317 (1973): sigma 3.40 A, eps/k 28 K',
        note='Fits PAH desorption from graphite (TPD) within -20..+16 %; '
             'binds a water monomer ~35 % too weakly with LB (graphite_ff_validation).'),
    'wu_aluru': dict(
        sigma=0.340, epsilon=0.2328,
        pairs={'CO': (0.3436, 0.0850 * KJ_PER_KCAL), 'CH': (0.2690, 0.0383 * KJ_PER_KCAL)},
        status='literature',
        source='Wu & Aluru, J. Phys. Chem. B 117, 8802 (2013), doi 10.1021/jp402051t: '
               'RPA-based set, Table 1 (p. 8803) and Table 5 (p. 8807): sigma_C-OW 3.436 A, '
               'eps_C-OW 0.0850 kcal/mol, sigma_C-HW 2.690 A, eps_C-HW 0.0383 kcal/mol',
        note='Read from the paper 2026-10-23. Of their four sets (MP2, DFT-SAPT, RPA, CCSD(T)) '
             'the RPA set gives a 42.4 deg nanodroplet / 36.4 deg macroscopic contact angle '
             '(their exp. reference 42 +- 7 deg; fresh-HOPG measurements are 64-69 deg) and the '
             'CNT RBM shift. Their MD: rigid SPC/E, LJ cutoff 1.2 nm (contact angle depends on '
             'it, Fig. 3b), two fixed carbon layers. The paper gives C-water pairs only; Steele '
             'carbon for every other partner is our choice, not theirs.'),
    'werder': dict(
        sigma=0.340, epsilon=0.2328, pairs={'CO': (0.319, 0.0937 * KJ_PER_KCAL)},
        status='unverified',
        source='Werder et al., J. Phys. Chem. B 107, 1345 (2003) (SPC/E, ~86 deg)',
        note='C-O 3.190 A / 0.0937 kcal/mol as tabulated by Wu & Aluru 2013, Table 1 '
             '(fit to an 86 deg contact angle); the Werder paper itself has not been read.'),
    'builder_fixed': dict(
        sigma=0.3851 / 2 ** (1 / 6), epsilon=0.4396, pairs={}, status='legacy',
        source='VeloxChem EvbSystemBuilder after the 2026-10-17 sigma fix (UFF C_R)',
        note='Overbinds PAHs 1.2-1.6x vs TPD; for comparison with old runs only.'),
    'builder_old': dict(
        sigma=0.3851, epsilon=0.4396, pairs={}, status='legacy',
        source='VeloxChem EvbSystemBuilder before 2026-10-17 (UFF r_min used as sigma)',
        note='Used in every graphite run before 2026-10-17; reproduces old runs only.'),
}

# Oxygen groups for surface.add_groups. All closed-shell (see surface.py); 'site' is
# where they go ('basal' or 'edge'), 'charge' their formal charge. Parameters come from
# surface_templates.generate_templates (GAFF2 + RESP HF/6-31G* on an H-capped graphene
# cluster, VeloxChem MMForceFieldGenerator) and are cached in data/surface_groups.json.
_GROUP_METHOD = ('GAFF2 types/bonded terms and RESP charges (HF/6-31G*) from VeloxChem '
                 'MMForceFieldGenerator on an H-capped graphene cluster of the same local '
                 'geometry; group atoms relaxed with GAFF2, skeleton fixed (surface_templates.py)')
SURFACE_GROUPS = {
    'epoxide': dict(site='basal', charge=0, anchors=2, status='computed', source=_GROUP_METHOD,
                    note='O bridging a C-C bond on the basal plane; anchors stay flat (no pucker).'),
    'hydroxyl_pair': dict(site='basal', charge=0, anchors=2, status='computed', source=_GROUP_METHOD,
                          note='Two C-OH at para positions of one ring, same face (a single basal '
                               'OH would leave an unpaired pi electron). Basal OH stays protonated.'),
    'phenol': dict(site='edge', charge=0, anchors=1, status='computed', source=_GROUP_METHOD,
                   note='Edge C-OH. pKa ~10, so mostly phenolate at pH 14.'),
    'phenolate': dict(site='edge', charge=-1, anchors=1, status='computed', source=_GROUP_METHOD,
                      note='Edge C-O(-), the pH 14 form of phenol.'),
    'carboxyl': dict(site='edge', charge=0, anchors=1, status='computed', source=_GROUP_METHOD,
                     note='Edge COOH (syn). pKa ~4-5, so carboxylate at pH 14.'),
    'carboxylate': dict(site='edge', charge=-1, anchors=1, status='computed', source=_GROUP_METHOD,
                        note='Edge COO(-), the pH 14 form of carboxyl.'),
    'o_quinone': dict(site='edge', charge=0, anchors=2, status='computed', source=_GROUP_METHOD,
                      note='Two C=O on an armchair-edge carbon pair (o-quinone, the closed-shell '
                           'edge carbonyl). Armchair edges only.'),
}

# Force-group convention shared with electrode_system (which adds the field, ion-pair
# and solute-restraint groups).
FORCE_GROUPS = {
    'bonded': 0,           # bonds, angles, torsions (solute; flexible solvents)
    'nonbonded': 1,        # NonbondedForce: PME Coulomb + all LJ including graphite
    'nbfix': 3,            # optional graphite C x water correction
    'ion_pair': 6,         # optional cation x anion sigma-scaling correction
    'field': 4,            # uniform field -q*E*z
    'restraint': 5,        # optional solute height restraints
    'carbon_restraint': 7, # carbon position restraints (restrained / anchored flexible carbon)
}


def check_status(kind, name, entry, allow_unverified):
    """Return None (literature), a note to record (legacy; computed = fitted by this
    module, not taken from a paper), or raise."""
    status = entry.get('status')
    if status == 'literature':
        return None
    msg = (f"{kind} '{name}' has status '{status}': {entry.get('note', '')} "
           f"Source: {entry.get('source', '?')}.")
    if status in ('legacy', 'computed'):
        return msg
    if not allow_unverified:
        raise ValueError(msg + ' Pass allow_unverified=True to use it anyway.')
    return msg

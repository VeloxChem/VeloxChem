"""Carbon force-field terms added to an OpenMM system (graphite only).

Graphite LJ itself lives in NonbondedForce (Lorentz-Berthelot), so VeloxChem's
SolvationFepDriver decouples solute-graphite together with solute-solvent.
The optional NBFIX correction acts only between graphite carbons and water,
never the solute, so the FEP driver does not need to know about it.
"""
from __future__ import annotations

import openmm as mm

from .params import CARBON_SETS, FORCE_GROUPS, NBFIX_FORCE_NAME


def nonbonded_force(system):
    return next(f for f in system.getForces() if isinstance(f, mm.NonbondedForce))


def add_nbfix_correction(system, carbon, water_O, water_H, carbon_set, cutoff, switch=None):
    """E = LJ(pair override) - LJ(Lorentz-Berthelot) for graphite C x water O/H.

    Added on top of NonbondedForce, so the net C-water interaction is the override.
    Returns None if the carbon set has no overrides.
    """
    pairs = CARBON_SETS[carbon_set]['pairs']
    if not pairs or not water_O:
        return None
    nb = nonbonded_force(system)
    expr = ('4*ep*((sp/r)^12-(sp/r)^6) - 4*el*((sl/r)^12-(sl/r)^6);'
            'sp = select(delta(k-4), sCO, select(delta(k-5), sCH, sl));'
            'ep = select(delta(k-4), eCO, select(delta(k-5), eCH, el));'
            'sl = 0.5*(sig1+sig2); el = sqrt(eps1*eps2); k = kind1+kind2')
    f = mm.CustomNonbondedForce(expr)
    f.setName(NBFIX_FORCE_NAME)
    sco, eco = pairs.get('CO', (1.0, 0.0))
    sch, ech = pairs.get('CH', (1.0, 0.0))
    for name, value in (('sCO', sco), ('eCO', eco), ('sCH', sch), ('eCH', ech)):
        f.addGlobalParameter(name, value)
    for p in ('sig', 'eps', 'kind'):
        f.addPerParticleParameter(p)
    kind = {}
    kind.update({i: 3 for i in carbon})
    if 'CO' in pairs:
        kind.update({i: 1 for i in water_O})
    if 'CH' in pairs:
        kind.update({i: 2 for i in water_H})
    for i in range(system.getNumParticles()):
        _, sig, eps = nb.getParticleParameters(i)
        f.addParticle([sig.value_in_unit(mm.unit.nanometer),
                       eps.value_in_unit(mm.unit.kilojoule_per_mole), kind.get(i, 0)])
    for k in range(nb.getNumExceptions()):
        i, j = nb.getExceptionParameters(k)[:2]
        f.addExclusion(i, j)
    f.addInteractionGroup(list(carbon), list(water_O) + list(water_H))
    f.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
    f.setCutoffDistance(cutoff)
    if switch:
        f.setUseSwitchingFunction(True)
        f.setSwitchingDistance(switch)
    f.setUseLongRangeCorrection(False)
    f.setForceGroup(FORCE_GROUPS['nbfix'])
    system.addForce(f)
    return f


def add_carbon_mobility(system, model, atoms, positions, frozen_exclusion_nm=0.36):
    """Set up the carbon model's atoms (system indices `atoms`, in model order).

    Per-atom charges and group LJ (surface groups, charged flake edges) are written
    into NonbondedForce first. Then the skeleton (lattice C and edge H) follows
    model.mobility:

    frozen      mass 0; LJ between skeleton atoms closer than frozen_exclusion_nm is
                excluded (a constant, ~1e4 kJ/mol per carbon, that would swamp the energy).
    restrained  real masses; harmonic tether 0.5*k*d^2 of every skeleton atom to its built
                position (periodicdistance, so GPU coordinate wrapping is harmless);
                same-layer skeleton pairs within frozen_exclusion_nm excluded.
    flexible    real masses; GAFF2 aromatic bonds/angles/torsions (C-H impropers for
                edges) with the lattice values as equilibrium, C-H constrained;
                tethered like 'restrained' if model.restraint_k is set.

    Surface-group atoms always move: they get their template bonded terms (X-H bonds
    constrained). Pairs on the bond graph are treated as in AMBER/GAFF (1-2 and 1-3
    excluded, 1-4 LJ x 0.5 and Coulomb x 1/1.2): all pairs in flexible mode, pairs with
    at least one group atom otherwise. positions: system positions (nm), the tether
    reference. Returns a dict for the metadata.
    """
    import numpy as np
    from .graphite import _wrap_tree
    from . import params as P
    nb = nonbonded_force(system)
    atoms = np.asarray(atoms, int)
    mode = model.mobility
    periodic = system.usesPeriodicBoundaryConditions()
    if periodic:
        nb.setExceptionsUsePeriodicBoundaryConditions(True)
    skel = model.skeleton
    info = dict(mode=mode, restraint_k_kj_mol_nm2=model.restraint_k,
                charged_atoms=_set_atom_parameters(nb, model, atoms))

    if mode in ('frozen', 'restrained'):
        sk = np.flatnonzero(skel)
        tree = _wrap_tree(model.xyz[sk], model.periodic, model.cell)
        pairs = tree.query_pairs(frozen_exclusion_nm, output_type='ndarray')
        pairs = sk[pairs] if len(pairs) else np.zeros((0, 2), int)
        if mode == 'restrained':
            pairs = pairs[model.layer[pairs[:, 0]] == model.layer[pairs[:, 1]]]
        for i, j in pairs:
            nb.addException(int(atoms[i]), int(atoms[j]), 0.0, 1.0, 0.0)
        info.update(excluded_pairs=len(pairs), exclusion_radius_nm=frozen_exclusion_nm,
                    exclusion_scope='all skeleton pairs' if mode == 'frozen' else 'same layer')
        if mode == 'frozen':
            for i in atoms[skel]:
                system.setParticleMass(int(i), 0.0)
    if mode == 'flexible':
        info.update(_add_carbon_bonded(system, model, atoms, periodic))
    if model.terms:
        info.update(_add_group_terms(system, model, atoms, periodic))
    info['bond_graph_exceptions'] = _graph_exceptions(nb, model, atoms, all_pairs=mode == 'flexible')

    if model.restraint_k:
        dist = 'periodicdistance(x,y,z,x0,y0,z0)' if periodic else 'sqrt((x-x0)^2+(y-y0)^2+(z-z0)^2)'
        f = mm.CustomExternalForce(f'0.5*k_carbon*{dist}^2')
        f.setName('carbon position restraint')
        f.addGlobalParameter('k_carbon', float(model.restraint_k))
        for p in ('x0', 'y0', 'z0'):
            f.addPerParticleParameter(p)
        pos = np.asarray(positions)
        for i in atoms[skel]:
            f.addParticle(int(i), [float(v) for v in pos[i]])
        f.setForceGroup(P.FORCE_GROUPS['carbon_restraint'])
        system.addForce(f)
        info['restraint_force'] = f.getName()
    return info


def _set_atom_parameters(nb, model, atoms):
    """Charges and LJ overrides of the model into NonbondedForce; returns the number of
    charged atoms."""
    import numpy as np
    q = model.charges()
    lj = model.lj
    n = 0
    for k, i in enumerate(atoms):
        has_lj = lj is not None and not np.isnan(lj[k, 0])
        if q[k] == 0.0 and not has_lj:
            continue
        _, sig, eps = nb.getParticleParameters(int(i))
        sig, eps = (lj[k, 0], lj[k, 1]) if has_lj else (sig, eps)
        nb.setParticleParameters(int(i), float(q[k]), sig, eps)
        n += q[k] != 0.0
    return int(n)


def _graph_exceptions(nb, model, atoms, all_pairs):
    """AMBER/GAFF exclusions on the model's bond graph: 1-2, 1-3 excluded; 1-4 Coulomb
    x 1/1.2 and LJ x 0.5 (Lorentz-Berthelot). all_pairs=False: only pairs that involve
    a surface-group atom (skeleton pairs are then handled by distance exclusions)."""
    import numpy as np
    n = model.n_atoms
    nbrs = [[] for _ in range(n)]
    for i, j in model.bonds():
        nbrs[int(i)].append(int(j))
        nbrs[int(j)].append(int(i))
    sources = range(n) if all_pairs else np.flatnonzero(~model.skeleton)
    done = {}
    for s0 in sources:
        dist = {int(s0): 0}
        frontier = [int(s0)]
        for d in (1, 2, 3):
            nxt = []
            for a in frontier:
                for b in nbrs[a]:
                    if b not in dist:
                        dist[b] = d
                        nxt.append(b)
            frontier = nxt
        for b, d in dist.items():
            if d and (min(s0, b), max(s0, b)) not in done:
                done[(min(s0, b), max(s0, b))] = d
    for (i, j), d in done.items():
        ai, aj = int(atoms[i]), int(atoms[j])
        if d < 3:
            nb.addException(ai, aj, 0.0, 1.0, 0.0)
        else:
            qi, si, ei = _plain(nb.getParticleParameters(ai))
            qj, sj, ej = _plain(nb.getParticleParameters(aj))
            nb.addException(ai, aj, qi * qj / 1.2, 0.5 * (si + sj), 0.5 * (ei * ej) ** 0.5)
    return len(done)


def _plain(prm):
    u = mm.unit
    q, sig, eps = prm
    return (q.value_in_unit(u.elementary_charge), sig.value_in_unit(u.nanometer),
            eps.value_in_unit(u.kilojoule_per_mole))


def _add_group_terms(system, model, atoms, periodic):
    """Bonded terms of the surface groups (model.terms, model indices); bonds to H
    become constraints, as with HBonds for the solute."""
    a = lambda i: int(atoms[int(i)])            # noqa: E731
    el = model.elements
    t = model.terms
    bf, af, tf = mm.HarmonicBondForce(), mm.HarmonicAngleForce(), mm.PeriodicTorsionForce()
    rb = mm.RBTorsionForce()
    n_con = 0
    for i, j, r0, k in t.get('bonds', []):
        if 'H' in (el[int(i)], el[int(j)]):
            system.addConstraint(a(i), a(j), float(r0))
            n_con += 1
        else:
            bf.addBond(a(i), a(j), float(r0), float(k))
    for i, j, k, th, kk in t.get('angles', []):
        af.addAngle(a(i), a(j), a(k), float(th), float(kk))
    for i, j, k, l, per, ph, kk in t.get('torsions', []):
        tf.addTorsion(a(i), a(j), a(k), a(l), int(per), float(ph), float(kk))
    for i, j, k, l, *c in t.get('rb', []):
        rb.addTorsion(a(i), a(j), a(k), a(l), *[float(x) for x in c])
    added = {}
    for f, name in ((bf, 'surface group bonds'), (af, 'surface group angles'),
                    (tf, 'surface group torsions'), (rb, 'surface group RB torsions')):
        n = (f.getNumBonds() if f is bf else f.getNumAngles() if f is af else f.getNumTorsions())
        if not n:
            continue
        f.setName(name)
        f.setUsesPeriodicBoundaryConditions(periodic)
        f.setForceGroup(FORCE_GROUPS['bonded'])
        system.addForce(f)
        added[name] = n
    return dict(group_terms=added, group_constraints=n_con,
                n_group_atoms=int((~model.skeleton).sum()))


def _add_carbon_bonded(system, model, atoms, periodic):
    """GAFF2 'ca' bonded terms on the carbon skeleton (see params.GAFF2_*)."""
    import math
    from . import params as P
    el = model.elements
    skel = model.skeleton
    bonds = [(int(i), int(j)) for i, j in model.bonds() if skel[i] and skel[j]]
    nbrs = [[] for _ in range(model.n_atoms)]
    for i, j in bonds:
        nbrs[i].append(j)
        nbrs[j].append(i)
    a = lambda i: int(atoms[i])            # noqa: E731  model -> system index
    bf, af, tf = mm.HarmonicBondForce(), mm.HarmonicAngleForce(), mm.PeriodicTorsionForce()
    n_con = 0
    for i, j in bonds:
        if el[i] == 'H' or el[j] == 'H':
            system.addConstraint(a(i), a(j), P.GAFF2_CH_BOND[0])
            n_con += 1
        else:
            bf.addBond(a(i), a(j), P.CC_BOND, P.GAFF2_CA_BOND_K)
    theta0 = 2 * math.pi / 3
    for j in range(model.n_atoms):
        for x in range(len(nbrs[j])):
            for y in range(x + 1, len(nbrs[j])):
                i, k = nbrs[j][x], nbrs[j][y]
                kk = P.GAFF2_CCH_ANGLE_K if 'H' in (el[i], el[k]) else P.GAFF2_CA_ANGLE_K
                af.addAngle(a(i), a(j), a(k), theta0, kk)
    for j, k in bonds:
        if el[j] == 'H' or el[k] == 'H':
            continue
        for i in nbrs[j]:
            for l in nbrs[k]:
                if i != k and l != j:
                    tf.addTorsion(a(i), a(j), a(k), a(l), 2, math.pi, P.GAFF2_CA_TORSION_K)
    n_imp = 0
    for c in range(model.n_atoms):
        hs = [h for h in nbrs[c] if el[h] == 'H']
        cs = [n for n in nbrs[c] if el[n] == 'C']
        if el[c] == 'C' and len(hs) == 1 and len(cs) == 2:
            # AMBER order: the central atom third (OpenMM's GAFF2 template lists it first)
            tf.addTorsion(a(cs[0]), a(cs[1]), a(c), a(hs[0]), 2, math.pi, P.GAFF2_CH_IMPROPER_K)
            n_imp += 1
    for f, name in ((bf, 'carbon bonds'), (af, 'carbon angles'), (tf, 'carbon torsions')):
        f.setName(name)
        f.setUsesPeriodicBoundaryConditions(periodic)
        f.setForceGroup(FORCE_GROUPS['bonded'])
        system.addForce(f)
    return dict(n_bonds=bf.getNumBonds(), n_ch_constraints=n_con, n_angles=af.getNumAngles(),
                n_torsions=tf.getNumTorsions() - n_imp, n_impropers=n_imp,
                bonded_source='GAFF2 ca/ha force constants (VeloxChem database/gaff-2.11.xml); '
                              'equilibrium C-C 0.142 nm and 120 deg (lattice), not GAFF2 0.1398 nm')

"""Parametrise the surface oxygen groups on H-capped graphene clusters (run once; QM).

HF/6-31G* RESP takes ~10-15 min per edge cluster and longer for the basal ones on
a laptop, so run production templates on the cluster, one job per template:

    python -m veloxchem.electrode.graphite.surface_templates --list           # the 11 keys
    python -m veloxchem.electrode.graphite.surface_templates --out tpl/ carboxylate/zigzag   # one job
    python -m veloxchem.electrode.graphite.surface_templates --merge tpl/*.json  # into data/surface_groups.json

Without --out the keys are generated in turn and merged into data/surface_groups.json
directly. --basis STO-3G --test writes cheap templates for the self-tests only, to
tests/data/surface_groups_test.json (never used by default).

For each group and site type (basal, armchair edge, zigzag edge):
1. cut a cluster from an ideal single graphene layer around a site of that type
   (basal: interior; armchair/zigzag: an H-terminated edge made by surface.add_step):
   every atom within SHELL bonds of the anchor(s) plus the atoms that close rings
   with them, pruned to >= 2 C neighbours, H-capped;
2. place the group with the ideal geometry of surface.ideal_group_atoms;
3. GAFF2 topology with zero charges (VeloxChem MMForceFieldGenerator, resp=False) and
   minimise the group atoms with the skeleton fixed;
4. RESP charges at HF/6-31G* (the MMForceFieldGenerator default) at that geometry,
   with GAFF2 types and bonded terms, through write_openmm_files;
5. store in the site frame (surface._site_frame): the group atoms (relaxed positions,
   RESP charges, GAFF2 LJ and types), the skeleton atoms within two bonds of the
   anchor(s) with their charges shifted equally so that the fragment carries the
   group's formal charge, every bonded term with a group atom, and the group bonds.
Results go to data/surface_groups.json, one entry per 'group/site'.
"""
from __future__ import annotations

import contextlib
import datetime
import json
import os
import sys
import tempfile
from pathlib import Path

import numpy as np

from .graphite import CarbonModel, make_slab
from .params import SURFACE_GROUPS
from .surface import (TEMPLATE_FILE, _hydrogenate, _neighbours, _prune, _shell, _site_frame,
                      _take, add_step, candidate_sites, ideal_group_atoms)

SHELL = 3           # bonds from the anchors kept in the cluster (plus ring closure)
                    # = the depth of the environment check in surface._apply_template
ZNUM = {'H': 1, 'C': 6, 'O': 8}


def site_keys():
    keys = []
    for g, spec in SURFACE_GROUPS.items():
        if spec['site'] == 'basal':
            keys.append(f'{g}/basal')
        elif g == 'o_quinone':
            keys.append(f'{g}/armchair')
        else:
            keys += [f'{g}/armchair', f'{g}/zigzag']
    return keys


def build_cluster(group, site):
    """(cluster CarbonModel with the group, anchors, origin, frame, group atom indices,
    reference atom indices (every non-cap skeleton atom), formal charge, atom names,
    cap H indices)."""
    sheet = make_slab('steele', 1, (3.0, 3.0))
    if site != 'basal':
        sheet = add_step(sheet, edge=site, keep_fraction=0.5)
    nbrs = _neighbours(sheet)
    centre = np.array([sheet.cell[0] / 2, sheet.cell[1] / 2, 0.0])
    best = None
    for c in candidate_sites(sheet, group, 'top', nbrs):
        if site == 'basal':
            anchors = c
        elif group == 'o_quinone':
            anchors = c[:2]
        else:
            if c[1] != site:
                continue
            anchors = c[:1]
        origin, _ = _site_frame(sheet, anchors, 'basal' if site == 'basal' else 'edge', 'top', nbrs)
        d = np.linalg.norm((origin - centre)[:2])
        if best is None or d < best[0]:
            best = (d, anchors, origin)
    _, anchors, origin = best
    # cut: atoms within SHELL bonds of the anchors plus the ring-closing atoms of the next
    # shell, so every reference atom (<= 2 bonds) keeps all its neighbours
    shell = _shell(nbrs, anchors, SHELL)
    ring = {j for i in shell for j in nbrs[i] if j not in shell
            and sum(1 for k in nbrs[j] if k in shell) >= 2}
    keep = np.array(sorted(set(shell) | ring))
    dv = sheet.xyz - origin
    for d in range(2):
        dv[:, d] -= np.round(dv[:, d] / sheet.cell[d]) * sheet.cell[d]
    clus = CarbonModel('flake', 'steele', [sheet.elements[i] for i in keep], origin + dv[keep],
                       np.zeros(len(keep), int), (False, False, False), (None, None, None),
                       info=dict(layers=1))
    anchor_pos = [origin + dv[a] for a in anchors]
    clus = _drop_orphan_h(_prune(clus, 0))
    n_real = clus.n_atoms
    clus = _hydrogenate(clus, 0)
    caps = set(range(n_real, clus.n_atoms))          # H that only cap the cut
    anchors = [int(np.argmin(np.linalg.norm(clus.xyz - p, axis=1))) for p in anchor_pos]
    nbrs = _neighbours(clus)
    kind = 'basal' if site == 'basal' else 'edge'
    origin, frame = _site_frame(clus, anchors, kind, 'top', nbrs)
    replaced = [j for a in anchors for j in nbrs[a] if clus.elements[j] == 'H'] if kind == 'edge' else []
    refs = [i for i in range(clus.n_atoms) if i not in caps and i not in set(replaced)]
    if group == 'o_quinone':
        anchor_xyz = [(clus.xyz[a], clus.xyz[h] - clus.xyz[a]) for a, h in zip(anchors, replaced)]
        anchor_xyz = [(a, o / np.linalg.norm(o)) for a, o in anchor_xyz]
    else:
        anchor_xyz = [clus.xyz[a] for a in anchors]
    new = ideal_group_atoms(group, origin, frame, anchor_xyz)
    keep = [i for i in range(clus.n_atoms) if i not in set(replaced)]
    remap = {o: n for n, o in enumerate(keep)}
    clus = _take(clus, keep)
    n0 = clus.n_atoms
    clus.elements = clus.elements + [e for e, _, _ in new]
    clus.xyz = np.vstack([clus.xyz, [x for _, _, x in new]])
    clus.layer = np.concatenate([clus.layer, np.zeros(len(new), int)])
    clus._bonds = None
    group_idx = list(range(n0, n0 + len(new)))
    names = [nm for _, nm, _ in new]
    return (clus, [remap[a] for a in anchors], origin, frame, group_idx, [remap[r] for r in refs],
            SURFACE_GROUPS[group]['charge'], names, sorted(remap[c] for c in caps))


def _drop_orphan_h(m):
    nbrs = _neighbours(m)
    orphan = [i for i in range(m.n_atoms) if m.elements[i] == 'H' and not nbrs[i]]
    return _take(m, [i for i in range(m.n_atoms) if i not in set(orphan)]) if orphan else m


def _vlx_molecule(elements, xyz_nm, charge):
    import veloxchem as vlx
    n_el = sum(ZNUM[e] for e in elements) - charge
    if n_el % 2:
        raise ValueError(f'cluster has {n_el} electrons: not closed-shell')
    mol = vlx.Molecule(list(elements), 10.0 * np.asarray(xyz_nm), 'angstrom')
    mol.set_charge(charge)
    mol.set_multiplicity(1)
    return mol


def _openmm_system(ff, workdir, stem):
    from openmm import app
    ff.write_openmm_files(str(Path(workdir) / stem), 'CLU')
    pdb = app.PDBFile(str(Path(workdir) / f'{stem}.pdb'))
    system = app.ForceField(str(Path(workdir) / f'{stem}.xml')).createSystem(
        pdb.topology, nonbondedMethod=app.NoCutoff, constraints=None, removeCMMotion=False)
    return system, pdb


@contextlib.contextmanager
def _cwd(path):
    old = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(old)


def generate_template(key, workdir=None, log=print, basis='6-31G*'):
    import openmm as mm
    import veloxchem as vlx
    from openmm import unit
    group, site = key.split('/')
    clus, anchors, origin, frame, gidx, refs, qf, names, caps = build_cluster(group, site)
    gset = set(gidx)
    workdir = Path(workdir or tempfile.mkdtemp(prefix=f'eb_tpl_{group}_{site}_'))
    formula = ''.join(f'{e}{clus.elements.count(e)}' for e in ('C', 'H', 'O') if e in clus.elements)
    log(f'{key}: cluster {formula} (charge {qf}), {len(refs)} reference atoms; {workdir}')

    # 1. GAFF2 with zero charges, relax the group atoms
    mol = _vlx_molecule(clus.elements, clus.xyz, qf)
    ff0 = vlx.MMForceFieldGenerator()
    ff0.ostream.mute()
    with _cwd(workdir):          # VeloxChem writes MOL*.h5 checkpoints into the current folder
        ff0.create_topology(mol, resp=False)
    sys0, pdb0 = _openmm_system(ff0, workdir, 'relax')
    for i in range(sys0.getNumParticles()):
        if i not in gset:
            sys0.setParticleMass(i, 0.0)
    ctx = mm.Context(sys0, mm.VerletIntegrator(0.001), mm.Platform.getPlatformByName('Reference'))
    ctx.setPositions(clus.xyz * unit.nanometer)
    mm.LocalEnergyMinimizer.minimize(ctx, 0.1, 10000)
    xyz = ctx.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    shift = float(np.abs(xyz - clus.xyz).max())
    log(f'{key}: GAFF2 relaxation of the group, max shift {shift:.4f} nm')

    # 2. RESP (HF/6-31G*) + GAFF2 at the relaxed geometry
    t0 = datetime.datetime.now()
    mol = _vlx_molecule(clus.elements, xyz, qf)
    ff = vlx.MMForceFieldGenerator()
    ff.ostream.mute()
    with _cwd(workdir):          # own folder per template: parallel tasks never share MOL_scf.h5
        ff.create_topology(mol, basis=vlx.MolecularBasis.read(mol, basis, ostream=None))
    system, pdb = _openmm_system(ff, workdir, 'cluster')
    log(f'{key}: RESP done in {(datetime.datetime.now() - t0).total_seconds():.0f} s')

    nb = next(f for f in system.getForces() if isinstance(f, mm.NonbondedForce))
    prm = [nb.getParticleParameters(i) for i in range(system.getNumParticles())]
    q = np.array([p[0].value_in_unit(unit.elementary_charge) for p in prm])
    sig = [p[1].value_in_unit(unit.nanometer) for p in prm]
    eps = [p[2].value_in_unit(unit.kilojoule_per_mole) for p in prm]
    types = [ff.atoms[i]['type'] for i in range(len(q))]

    # fragment charge: each cap H's charge goes to its carbon, so refs + group = formal
    qfix = q.copy()
    for c in caps:
        (k,) = [j for j in _neighbours(clus)[c]]
        qfix[k] += q[c]
    residual = qf - (qfix[refs].sum() + qfix[gidx].sum())
    if abs(residual) > 1e-4:
        raise RuntimeError(f'{key}: fragment charge off by {residual:.5f}')
    qfix[refs[0]] += residual

    u, v, n = frame
    R = np.array([u, v, n])                    # lab -> frame (rows)
    to_frame = lambda p: (R @ (np.asarray(p) - origin)).round(6).tolist()   # noqa: E731
    order = list(refs) + list(gidx)
    index = {a: k for k, a in enumerate(order)}

    terms = {'bonds': [], 'angles': [], 'torsions': [], 'rb': []}
    other = set()
    for f in system.getForces():
        if isinstance(f, mm.HarmonicBondForce):
            for k in range(f.getNumBonds()):
                i, j, r0, kk = f.getBondParameters(k)
                if {i, j} & gset:
                    terms['bonds'].append([index[i], index[j], r0.value_in_unit(unit.nanometer),
                                           kk.value_in_unit(unit.kilojoule_per_mole / unit.nanometer ** 2)])
        elif isinstance(f, mm.HarmonicAngleForce):
            for k in range(f.getNumAngles()):
                i, j, l, th, kk = f.getAngleParameters(k)
                if {i, j, l} & gset:
                    terms['angles'].append([index[i], index[j], index[l], th.value_in_unit(unit.radian),
                                            kk.value_in_unit(unit.kilojoule_per_mole / unit.radian ** 2)])
        elif isinstance(f, mm.PeriodicTorsionForce):
            for k in range(f.getNumTorsions()):
                i, j, l, m_, per, ph, kk = f.getTorsionParameters(k)
                if {i, j, l, m_} & gset:
                    terms['torsions'].append([index[i], index[j], index[l], index[m_], int(per),
                                              ph.value_in_unit(unit.radian),
                                              kk.value_in_unit(unit.kilojoule_per_mole)])
        elif isinstance(f, mm.RBTorsionForce):
            for k in range(f.getNumTorsions()):
                i, j, l, m_, *c = f.getTorsionParameters(k)
                if {i, j, l, m_} & gset:
                    terms['rb'].append([index[i], index[j], index[l], index[m_]] +
                                       [x.value_in_unit(unit.kilojoule_per_mole) for x in c])
        elif not isinstance(f, (mm.NonbondedForce, mm.CMMotionRemover)):
            other.add(type(f).__name__)
    if set(caps) & set(index):
        raise RuntimeError(f'{key}: a cap H is in the reference set')
    if other:
        raise RuntimeError(f'{key}: unhandled force types in the generated system: {other}')
    graph_bonds = [[index[b[0].index], index[b[1].index]] for b in pdb.topology.bonds()
                   if {b[0].index, b[1].index} & gset]

    import veloxchem
    cn = _neighbours(clus)
    rpos = {r: k for k, r in enumerate(refs)}
    capset = set(caps)
    tpl = dict(
        group=group, site=site, formal_charge=qf,
        skeleton_refs=[dict(element=clus.elements[r], pos=to_frame(xyz[r]), charge=float(qfix[r]),
                            gaff_type=types[r], resp_charge=float(q[r]),
                            neighbours=sorted(rpos[x] for x in cn[r] if x in rpos),
                            n_caps=sum(1 for x in cn[r] if x in capset)) for r in refs],
        group_atoms=[dict(element=clus.elements[g], name=names[k], pos=to_frame(xyz[g]),
                          charge=float(qfix[g]), sigma=float(sig[g]), epsilon=float(eps[g]),
                          gaff_type=types[g]) for k, g in enumerate(gidx)],
        terms=terms, graph_bonds=graph_bonds,
        provenance=dict(
            method=f'RESP HF/{basis} + GAFF2 (VeloxChem MMForceFieldGenerator.create_topology), '
                   'group atoms relaxed with GAFF2 at zero charge, skeleton fixed',
            basis=basis, test_only=basis != '6-31G*',
            cluster=formula, cluster_atoms=clus.n_atoms, cluster_charge=qf, shell_bonds=SHELL,
            n_reference_atoms=len(refs), n_cap_h=len(caps),
            cap_h_charge_moved_to_their_carbons=float(q[caps].sum()) if caps else 0.0,
            group_resp_sum=float(q[gidx].sum()),
            relax_max_shift_nm=shift, veloxchem=getattr(veloxchem, '__version__', '?'),
            created=datetime.datetime.now().isoformat(timespec='seconds')))
    (workdir / 'template.json').write_text(json.dumps(tpl, indent=1))
    return tpl


def generate_templates(keys=None, path=TEMPLATE_FILE, log=print, basis='6-31G*'):
    """Generate templates in turn and merge each into the JSON file as soon as it is done."""
    for key in keys or site_keys():
        merge_templates([generate_template(key, log=log, basis=basis)], path)
        log(f'{key}: saved to {path}')


def merge_templates(items, path=TEMPLATE_FILE):
    """Merge templates (dicts, or JSON files written with --out) into the template file."""
    path = Path(path)
    data = json.loads(path.read_text()) if path.exists() else {'version': 1, 'templates': {}}
    for it in items:
        tpl = it if isinstance(it, dict) else json.loads(Path(it).read_text())
        data['templates'][f"{tpl['group']}/{tpl['site']}"] = tpl
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=1))
    return data


TEST_TEMPLATE_FILE = TEMPLATE_FILE.parent / 'surface_groups_test.json'


def main(argv):
    import argparse
    ap = argparse.ArgumentParser(prog='python -m veloxchem.electrode.graphite.surface_templates')
    ap.add_argument('keys', nargs='*', help='group/site keys (default: all missing ones)')
    ap.add_argument('--list', action='store_true', help='print all keys and exit')
    ap.add_argument('--out', help='write one <group>__<site>.json per key into this folder')
    ap.add_argument('--merge', nargs='+', help='merge --out JSON files into the template file')
    ap.add_argument('--basis', default='6-31G*')
    ap.add_argument('--test', action='store_true', help=f'write to {TEST_TEMPLATE_FILE.name} (self-tests)')
    a = ap.parse_args(argv)
    path = TEST_TEMPLATE_FILE if a.test else TEMPLATE_FILE
    if a.basis != '6-31G*' and not a.test:
        ap.error('a basis other than 6-31G* is for --test templates only')
    if a.list:
        print('\n'.join(site_keys()))
        return
    if a.merge:
        merge_templates(a.merge, path)
        print(f'merged {len(a.merge)} templates into {path}')
        return
    done = set(json.loads(path.read_text())['templates']) if path.exists() else set()
    todo = a.keys or [k for k in site_keys() if k not in done]
    log = lambda s: print(s, flush=True)       # noqa: E731
    log(f'templates to generate ({a.basis}): {todo}')
    for key in todo:
        tpl = generate_template(key, log=log, basis=a.basis)
        if a.out:
            out = Path(a.out) / f"{key.replace('/', '__')}.json"
            out.parent.mkdir(parents=True, exist_ok=True)
            out.write_text(json.dumps(tpl, indent=1))
            log(f'{key}: written to {out}')
        else:
            merge_templates([tpl], path)
            log(f'{key}: saved to {path}')


if __name__ == '__main__':
    main(sys.argv[1:])

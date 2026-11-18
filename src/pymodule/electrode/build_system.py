"""Working script: build a graphite electrode system and, optionally, run the solvation FEP.

Same steps as prepare-system.ipynb, with the graphite built on its own first:

  1. solute force field   VeloxChem MMForceFieldGenerator from SMILES or xyz (RESP, GAFF2),
                          or an existing solute.pdb + solute.xml
  2. graphite             graphite_builder.GraphiteBuilder: the carbon electrode only
  3. solvent              VeloxChem SolvationBuilder, packed around graphite + solute
                          (pure water unless an electrolyte is given)
  4. OpenMM system        solute.xml + graphite terms + water/ion force field + uniform field
  5. minimise + short MD  starting structure for the FEP, at fixed volume (or z-NPT with
                          md.pressure_coupling='z', as water_4V); with md.report_ps,
                          initial.pdb + trajectory.dcd + observables.csv (MDOutput, as water_4V)
  6. files                system.xml, system.pdb, solute.pdb, solute.xml, metadata.json
  7. FEP                  SolvationFepDriver.compute_solvation_from_openmm_files(...)

Edit CONFIG below, or put any of its keys in a JSON file:

    python build_system.py                 # CONFIG as written
    python build_system.py my_run.json     # CONFIG updated with my_run.json

Adsorption of a (charged) solute = two FEPs in the same slab box, the solute held by a
flat-bottom window on its height above the top layer: restraint {'leg': 'surface'} and
{'leg': 'mid'}; then
    python build_system.py --adsorption runs/surface runs/mid
gives K_ads and dG_ads (system_tools.adsorption_free_energy). Slab runs use the driver's
ion_corrections='none'; the bulk reference (system='bulk_reference': solute in a cubic water
box) uses 'auto' with SPC/E.

Differences from the notebook: no barostat (fixed volume; a barostat with a field is
wrong on GPU), and no absolute-position restraint (the solute is placed at height_nm).
"""
from __future__ import annotations

import copy
import json
import sys
from pathlib import Path

import numpy as np
from openmm import unit

from . import system_tools as T  # noqa: E402
from .system_tools import GraphiteBuilder  # noqa: E402
from ..environment import get_data_path  # noqa: E402

CONFIG = dict(
    output='runs/example',
    system='electrode',            # or 'bulk_reference': the solute in a cubic box, no graphite
    bulk_box_nm=5.0,               # edge of the bulk-reference box
    # 1. solute: a folder with solute.pdb + solute.xml, or SMILES / xyz (+ charge for xyz)
    solute_dir=str(get_data_path() / 'electrode' / 'examples' / '14a'),
    solute_smiles=None,
    solute_xyz=None,
    solute_charge=0,
    # 2. graphite (GraphiteBuilder arguments) and optional surface features
    graphite=dict(carbon_set='wu_aluru', layers=3, lateral_nm=[5.0, 5.0], stacking='AB',
                  orientation='armchair_x', mobility='frozen'),
    surface=dict(step=None,        # e.g. {'edge': 'zigzag', 'keep_fraction': 0.5}
                 defects=None,     # e.g. {'vacancies': 1, 'stone_wales': 1}
                 groups=None,      # e.g. {'carboxylate': 2, 'epoxide': 1}
                 seed=1),
    # 3. solvent
    electrolyte=None,              # None: pure water; {'KOH': 1.0}: 1 mol/kg KOH;
                                   # {'H2SO4': 3.39}: 3.39 mol/kg H2SO4 as H3O+ + HSO4- (= 3 M)
    counterion='K',                # added only to neutralise a charged solute / surface
                                   # (an acid electrolyte neutralises with its own H3O+/HSO4-)
    box_z_nm=7.0,
    liquid_density=997.0,          # kg/m3; calibrate once per campaign (calibrate=True)
    calibrate=False,
    height_nm=0.40,                # solute COM above the top carbon layer
    solute_xy_nm=None,             # lateral position (default: box centre)
    flip=False,
    # starting orientation (slab systems): None = flat, ring plane parallel to the sheet; or
    #   {'flip': True, 'spin_deg': 90}            upside down, turned 90 deg in the plane
    #   {'tilt_deg': 90, 'bay': 'up'}             edge-on, N atoms pointing away from the sheet
    # (tilt is a turn about the molecule's long axis; raise height_nm for tilted starts)
    orientation=None,
    # solute atom indices (xyz order) of the ring that defines the molecular plane, for a
    # molecule with a side chain on one face: placement puts the ring centre at height_nm with
    # the chain up, and the tilt, ring height (observables, pre-FEP MD) use these atoms.
    # None: the plane of all heavy atoms, placed by the COM (planar molecules).
    plane_atoms=None,
    seed=1,
    # adsorption legs: None (no restraint), or a flat-bottom window on the solute COM height h
    # above the top layer, on in every FEP window:
    #   {'leg': 'surface', 'h_min': 0.0, 'h_max': 0.8, 'k': 2000}   (nm, kJ/mol/nm^2)
    #   {'leg': 'mid', 'half_width': 0.15, 'k': 2000}               (placed and held mid-liquid)
    restraint=None,
    liquid_density_from=None,      # calibration.json of a calibration run: overrides liquid_density
    calibration=dict(bulk_box_nm=4.0, bulk_equil_ps=200, bulk_frames=100, bulk_stride=250,
                     iterations=3, tol=0.005, equil_ps=100, frames=50, stride=500),
    cutoff_nm=1.0,                 # LJ/PME real-space cutoff (Wu-Aluru fitted with 1.2 nm)
    switch_nm=0.9,
    # 4. field along +z in V/nm (E > 0: top face anode-like, bottom face cathode-like)
    field_V_per_nm=0.0,
    # 5. minimisation + MD before the FEP
    md=dict(minimize=True, steps=10000, timestep_fs=2.0, platform=None,
            report_ps=10.0,            # trajectory.dcd + observables.csv every 10 ps (None: off)
            dcd_ps=None,               # trajectory.dcd interval if longer (a multiple of report_ps)
            # None: fixed volume throughout. 'z': as the water_4V sweeps, x/y fixed and z free
            # at 1 bar during production MD and in every FEP window, after nvt_ps at fixed
            # volume and an optional high-pressure/temperature stage (hpht); the slab is made
            # one rigid molecule at z = 0 and the field energy is checked under z scaling.
            pressure_coupling=None, pressure_bar=1.0, barostat_frequency=25, nvt_ps=1000.0,
            hpht=None,                 # e.g. {'ps': 4000, 'pressure_bar': 1000, 'temperature_k': 500}
            # FEP windows after z-NPT MD: 'fixed' = fixed volume at the last MD box (exact for
            # VeloxChem's MBAR re-evaluation, which restores positions but not the box), or
            # 'z' = keep the z barostat in every window, as the water_4V sweeps did (each
            # snapshot is then re-evaluated in the window's starting box)
            fep_windows='fixed'),
    # 7. FEP (SolvationFepDriver settings)
    fep=dict(run=False, num_steps=50000, num_snapshots=500, save_trajectory_xtc=True,
             save_system_xml=True, platform=None),
)


def solute_files(cfg, out):
    """Step 1: solute.pdb + solute.xml, from solute_dir or generated with MMForceFieldGenerator."""
    if cfg.get('solute_dir'):
        d = Path(cfg['solute_dir'])
        return d / 'solute.pdb', d / 'solute.xml'
    import veloxchem as vlx
    stem = Path(out) / 'solute_ff' / 'solute'
    if not stem.with_suffix('.xml').exists():
        if cfg.get('solute_smiles'):
            mol = vlx.Molecule.read_smiles(cfg['solute_smiles'])
        elif cfg.get('solute_xyz'):
            mol = vlx.Molecule.read_xyz_file(cfg['solute_xyz'])
            mol.set_charge(cfg.get('solute_charge', 0))
        else:
            raise ValueError('give solute_dir, solute_smiles or solute_xyz')
        mol.set_multiplicity(1)
        stem.parent.mkdir(parents=True, exist_ok=True)
        ff = vlx.MMForceFieldGenerator()
        ff.create_topology(mol)                  # GAFF2 + RESP (HF/6-31G*), as in the notebook
        ff.write_openmm_files(str(stem), T.SOLUTE_RESNAME)
    return stem.with_suffix('.pdb'), stem.with_suffix('.xml')


def make_graphite(cfg):
    """Step 2: the carbon electrode, with its optional surface features."""
    gb = GraphiteBuilder(**cfg['graphite'])
    s = cfg.get('surface') or {}
    if s.get('step'):
        gb.add_step(**s['step'])
    if s.get('defects'):
        gb.add_defects(seed=s.get('seed'), **s['defects'])
    if s.get('groups'):
        gb.add_groups(s['groups'], seed=s.get('seed'))
    if s.get('face_charge_uC_cm2') is not None:          # constant-charge electrode: [top, bottom] or top
        fc = s['face_charge_uC_cm2']
        fc = list(fc) if isinstance(fc, (list, tuple)) else [fc, None]
        gb.set_face_charge(*fc)
    return gb


def solute_orientation(cfg, solute_pdb, solute_xml):
    """Step 2b: the solute's starting orientation above the sheet -> (flip, rotation or None),
    as place_solute applies them to the plane frame (x long axis, z ring normal)."""
    o = cfg.get('orientation') or {}
    flip = bool(o.get('flip', cfg.get('flip', False)))
    tilt, spin = np.radians(o.get('tilt_deg', 0.0)), np.radians(o.get('spin_deg', 0.0))
    if o.get('bay') and tilt:
        _, xyz, mass, _, labels = T._read_solute(str(solute_pdb), str(solute_xml))
        x = T.plane_frame(xyz, mass, cfg.get('plane_atoms'))
        if flip:
            x[:, 1:] *= -1
        n_atoms = [i for i, e in enumerate(labels) if e == 'N']
        if not n_atoms:
            raise ValueError("orientation 'bay' needs N atoms in the solute")
        if (x[n_atoms, 1].mean() > 0) != (o['bay'] == 'up'):   # +tilt about x sends +y up
            tilt = -tilt
    if not tilt and not spin:
        return flip, None
    c, s = np.cos(tilt), np.sin(tilt)
    rx = np.array([[1, 0, 0], [0, c, -s], [0, s, c]])
    c, s = np.cos(spin), np.sin(spin)
    rz = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])
    return flip, rz @ rx


def electrolyte(cfg):
    """Step 3: species and proportions for SolvationBuilder.custom_solvate."""
    e = cfg.get('electrolyte')
    if not e:
        return T.pure_water()
    if set(e) == {'KOH'}:
        return T.koh_mixture(e['KOH'])
    if set(e) == {'H2SO4'}:
        return T.h2so4_mixture(e['H2SO4'])
    if set(e) == {'species', 'proportions'}:
        return list(e['species']), list(e['proportions'])
    raise ValueError("electrolyte: None, {'KOH': molality}, {'H2SO4': molality} or "
                     "{'species': [...], 'proportions': [...]}")


def leg_window(cfg, graphite):
    """Solute placement height and flat-bottom window for the restraint leg (or none)."""
    r = cfg.get('restraint')
    if not r:
        return cfg['height_nm'], None
    k = float(r.get('k', 2000.0))
    if r['leg'] == 'surface':
        return cfg['height_nm'], dict(leg='surface', lo=float(r.get('h_min', 0.0)),
                                      hi=float(r.get('h_max', 0.8)), k=k)
    if r['leg'] == 'mid':
        h_mid = 0.5 * (cfg['box_z_nm'] - graphite.thickness)
        w = float(r.get('half_width', 0.15))
        return h_mid, dict(leg='mid', lo=h_mid - w, hi=h_mid + w, k=k)
    raise ValueError("restraint leg must be 'surface' or 'mid'")


CONTACT_NM = 0.35                   # contact cutoff, as in the water_4V sweep


class MDOutput:
    """water_4V-style record of the pre-FEP MD, written into the leg folder:

      initial.pdb       the minimised start (frame 0)
      trajectory.dcd    all atoms, unwrapped, one frame per report interval (frames 1, 2, ...)
      observables.csv   one row per frame, starting with frame 0:
                        step, time; potential energy with the field term split out (its own
                        force group), kinetic energy, temperature; box and field x Lz;
                        solute COM z, height above the top layer, ring-plane tilt to z and
                        reorientation since frame 0; closest heavy atom to the top layer and
                        heavy atoms within 0.35 nm of it; water O, cation-centre and anion-O
                        atoms within 0.35 nm of solute heavy atoms; water O between the solute
                        and the sheet (within 0.5 nm laterally of the COM); the same distances
                        to the bottom face (the periodic image below the slab: for a free
                        solute, which may adsorb on either face)
    dcd_every: trajectory.dcd interval in steps (a multiple of every; default every).
    The final frame is system.pdb (the FEP start)."""

    def __init__(self, es, sim, out, every, timestep_fs, dcd_every=None, plane_atoms=None):
        from openmm import app
        self.es, self.sim, self.every, self.dt_fs = es, sim, int(every), float(timestep_fs)
        self.dcd_every = int(dcd_every or every)
        n = es.system.getNumParticles()
        mass = np.array([es.system.getParticleMass(i).value_in_unit(unit.dalton) for i in range(n)])
        idx = es.indices
        self.solute = np.array(idx['solute'], int)
        self.m_solute = mass[self.solute]
        self.heavy = self.solute[self.m_solute > 2.0]
        # ring atoms of a side-chain molecule: the plane (tilt) and an extra ring_height column
        self.ring = self.solute[np.asarray(plane_atoms, int)] if plane_atoms is not None else None
        self.m_ring = mass[self.ring] if self.ring is not None else None
        self.plane = self.ring if self.ring is not None else self.heavy
        self.top = np.array(idx.get('top_layer') or [], int)
        self.bottom = np.array(idx.get('bottom_layer') or [], int)
        anion = np.array(idx.get('anion') or [], int)
        anion_o = anion[(mass[anion] > 15.0) & (mass[anion] < 17.0)] if len(anion) else anion
        self.partners = dict(water_oxygen=np.array(idx.get('water_O') or [], int),
                             cation_center=np.array(idx.get('cation_center') or [], int),
                             anion_oxygen=anion_o if len(anion_o) else np.array(idx.get('anion_center') or [], int))
        self.field_group = es.metadata.get('force_groups', {}).get('field')
        self.E = es.metadata['field']['V_per_nm']
        import openmm as mm
        dof = 3 * int((mass > 0).sum()) - es.system.getNumConstraints()
        if any(isinstance(f, mm.CMMotionRemover) for f in es.system.getForces()):
            dof -= 3
        self.dof = dof
        out = Path(out)
        out.mkdir(parents=True, exist_ok=True)
        state = sim.context.getState(getPositions=True)
        x0 = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        with open(out / 'initial.pdb', 'w') as fh:
            app.PDBFile.writeFile(sim.topology, x0 * unit.nanometer, fh)
        self.normal0 = self._normal(x0)
        self.reporter = app.DCDReporter(str(out / 'trajectory.dcd'), self.dcd_every, enforcePeriodicBox=False)
        sim.reporters.append(self.reporter)
        self.fh = open(out / 'observables.csv', 'w')
        self.header = False
        self.box_z = []

    def _normal(self, x):
        h = x[self.plane]
        return np.linalg.svd(h - h.mean(0), full_matrices=False)[2][-1]

    @staticmethod
    def _contacts(a, b, box, cut):
        if not len(b):
            return 0
        d = b[None, :, :] - a[:, None, :]
        d -= box * np.round(d / box)
        return int(((d ** 2).sum(-1) <= cut * cut).any(0).sum())

    def record(self):
        sim = self.sim
        st = sim.context.getState(getPositions=True, getEnergy=True)
        x = st.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        box = np.diag(st.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.nanometer))
        pot = st.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        kin = st.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole)
        fe = (sim.context.getState(getEnergy=True, groups={self.field_group}).getPotentialEnergy()
              .value_in_unit(unit.kilojoule_per_mole) if self.field_group is not None else 0.0)
        if not (np.isfinite(x).all() and np.isfinite(pot + kin)):
            raise RuntimeError('non-finite coordinates or energy during the pre-FEP MD')
        self.box_z.append(float(box[2]))
        com = (x[self.solute] * self.m_solute[:, None]).sum(0) / self.m_solute.sum()
        n = self._normal(x)
        row = dict(step=sim.currentStep, time_ns=sim.currentStep * self.dt_fs / 1e6,
                   potential_energy_kJ_mol=pot, field_energy_kJ_mol=fe,
                   potential_energy_without_field_kJ_mol=pot - fe, kinetic_energy_kJ_mol=kin,
                   temperature_K=2 * kin / (self.dof * 0.00831446261815324),
                   box_x_nm=box[0], box_y_nm=box[1], box_z_nm=box[2], box_voltage_V=self.E * box[2],
                   solute_com_z_nm=float(com[2] % box[2]),
                   plane_tilt_to_z_deg=float(np.degrees(np.arccos(min(1.0, abs(n[2]))))),
                   plane_reorientation_deg=float(np.degrees(np.arccos(min(1.0, abs(n @ self.normal0))))))
        if len(self.top):
            ztop = x[self.top, 2].mean()
            dz = x[self.heavy, 2] - ztop
            w = self.partners['water_oxygen']
            lat = x[w, :2] - com[:2]
            lat -= box[:2] * np.round(lat / box[:2])
            below = (x[w, 2] > ztop) & (x[w, 2] < com[2]) & ((lat ** 2).sum(1) < 0.25)
            if self.ring is not None:
                row['ring_height_nm'] = float((x[self.ring, 2] * self.m_ring).sum() / self.m_ring.sum() - ztop)
            row.update(solute_height_nm=float(com[2] - ztop), solute_min_plane_distance_nm=float(dz.min()),
                       solute_heavy_plane_contacts=int((dz <= CONTACT_NM).sum()),
                       water_oxygen_between_solute_and_sheet=int(below.sum()))
        else:
            row.update(solute_height_nm='', solute_min_plane_distance_nm='', solute_heavy_plane_contacts='',
                       water_oxygen_between_solute_and_sheet='')
        if len(self.bottom):       # bottom face: below the liquid, at z = Lz - t/2 (mod Lz)
            zbot = float(np.mod(x[self.bottom, 2].mean(), box[2]))
            zc = float(com[2] % box[2])
            hz = np.mod(x[self.heavy, 2], box[2])
            w = self.partners['water_oxygen']
            wz = np.mod(x[w, 2], box[2])
            lat = x[w, :2] - com[:2]
            lat -= box[:2] * np.round(lat / box[:2])
            above = (wz < zbot) & (wz > zc) & ((lat ** 2).sum(1) < 0.25)
            if self.ring is not None:
                zr = float(np.mod((x[self.ring, 2] * self.m_ring).sum() / self.m_ring.sum(), box[2]))
                row['ring_height_bottom_nm'] = zbot - zr
            row.update(solute_height_bottom_nm=zbot - zc, solute_min_bottom_distance_nm=float((zbot - hz).min()),
                       solute_heavy_bottom_contacts=int(((zbot - hz) <= CONTACT_NM).sum()),
                       water_oxygen_between_solute_and_bottom=int(above.sum()))
        for k, ii in self.partners.items():
            row[f'{k}_contacts'] = self._contacts(x[self.heavy], x[ii], box, CONTACT_NM) if len(ii) else ''
        if not self.header:
            self.fh.write(','.join(row) + '\n')
            self.header = True
        self.fh.write(','.join(f'{v:.6g}' if isinstance(v, float) else str(v) for v in row.values()) + '\n')
        self.fh.flush()

    def close(self):
        if self.reporter in self.sim.reporters:
            self.sim.reporters.remove(self.reporter)
        self.fh.close()


def solute_heights(es, sim, steps, every=100, bin_nm=0.05, output=None, lz_range=1.0, plane_atoms=None):
    """Solute COM height above the top layer, sampled during MD (window check), and the
    cation/anion centre density along z (per nm^3, cell z) in the first and second half of
    the run: with an electrolyte and a field, the second half should match the first once
    the double layers have formed. Also the tilt of the solute's ring plane (angle between
    its normal and z: 0 flat, 90 edge-on). With plane_atoms (ring atoms of a side-chain molecule)
    the tilt uses the ring alone and the ring-centre heights are returned too.
    Returns (heights, profiles or None, tilts, ring heights or None)."""
    idx = es.indices
    m = np.array([es.system.getParticleMass(i)._value for i in idx['solute']])
    heavy = [i for i, mi in zip(idx['solute'], m) if mi > 2.0]
    ring = None if plane_atoms is None else np.asarray(idx['solute'])[np.asarray(plane_atoms, int)]
    m_ring = None if ring is None else m[np.asarray(plane_atoms, int)]
    plane = heavy if ring is None else ring
    tilts, ring_hs = [], []
    ions = {k: idx[k] for k in ('cation_center', 'anion_center') if idx.get(k)}
    lx, ly, lz = es.metadata['box_nm']
    edges = np.arange(0.0, lz * lz_range + bin_nm, bin_nm)    # lz_range > 1 when z may grow
    hist = {k: np.zeros((2, len(edges) - 1)) for k in ions}
    n = max(1, steps // every)
    hs = []
    for f in range(n):
        sim.step(every)
        if output and sim.currentStep % output.every == 0:
            output.record()
        x = sim.context.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        hs.append(float((x[idx['solute'], 2] * m).sum() / m.sum() - x[idx['top_layer'], 2].mean()))
        if ring is not None:
            ring_hs.append(float((x[ring, 2] * m_ring).sum() / m_ring.sum() - x[idx['top_layer'], 2].mean()))
        normal = np.linalg.svd(x[plane] - x[plane].mean(0))[2][-1]
        tilts.append(float(np.degrees(np.arccos(min(1.0, abs(normal[2]))))))
        for k, ii in ions.items():
            zi = x[ii, 2] if lz_range > 1 else np.mod(x[ii, 2], lz)    # the liquid never wraps in z
            hist[k][int(2 * f >= n)] += np.histogram(zi, bins=edges)[0]
    if not ions:
        return np.array(hs), None, np.array(tilts), (np.array(ring_hs) if ring is not None else None)
    frames = np.array([n - n // 2, max(n // 2, 1)])      # frames in the first and second half
    prof = {k: (h / (frames[:, None] * lx * ly * bin_nm)).round(4).tolist() for k, h in hist.items()}
    return np.array(hs), dict(z_nm=(0.5 * (edges[1:] + edges[:-1])).round(4).tolist(), bin_nm=bin_nm,
                              halves='[first half, second half] of the pre-FEP MD', **prof), np.array(tilts), \
        (np.array(ring_hs) if ring is not None else None)


def equilibrate_z_npt(sim, md, out, seed=1, temperature=298.15):
    """water_4V preparation before z-NPT production: nvt_ps at fixed volume (barostat
    frequency 0), then the z barostat on; optional hpht stage (raised pressure and temperature,
    then back to pressure_bar and 298.15 K with new velocities). Logs equilibration.csv and
    hpht.csv; production time starts at 0 afterwards."""
    import openmm as mm
    from openmm import app
    dt = md.get('timestep_fs', 2.0)
    every = max(1, round(md.get('report_ps', 10.0) * 1000 / dt))
    baro = [f for f in sim.system.getForces() if isinstance(f, mm.MonteCarloAnisotropicBarostat)]
    if len(baro) != 1:
        raise RuntimeError('z-NPT preparation needs exactly one z barostat')
    baro = baro[0]

    def run(name, steps):
        with open(Path(out) / f'{name}.csv', 'w') as fh:
            r = app.StateDataReporter(fh, every, step=True, time=True, potentialEnergy=True,
                                      temperature=True, volume=True, density=True, separator=',')
            sim.reporters.append(r)
            try:
                sim.step(steps)
            finally:
                sim.reporters.remove(r)

    log = {}
    nvt = round(md.get('nvt_ps', 1000.0) * 1000 / dt)
    if nvt:
        print(f'z-NPT preparation: {nvt:,} steps at fixed volume', flush=True)
        run('equilibration', nvt)
    log['nvt_ps'] = nvt * dt / 1000
    baro.setFrequency(int(md.get('barostat_frequency', 25)))
    sim.context.reinitialize(preserveState=True)
    p0 = float(md.get('pressure_bar', 1.0))
    pz, tk = mm.MonteCarloAnisotropicBarostat.PressureZ(), mm.MonteCarloAnisotropicBarostat.Temperature()
    sim.context.setParameter(pz, p0)
    h = md.get('hpht') or {}
    if h.get('ps'):
        steps = round(h['ps'] * 1000 / dt)
        ph, th = float(h.get('pressure_bar', 1000.0)), float(h.get('temperature_k', 500.0))
        print(f'high-pressure/temperature z-NPT: {steps:,} steps at {ph:g} bar, {th:g} K', flush=True)
        sim.context.setParameter(pz, ph)
        sim.context.setParameter(tk, th)
        sim.integrator.setTemperature(th * unit.kelvin)
        run('hpht', steps)
        sim.context.setParameter(pz, p0)
        sim.context.setParameter(tk, temperature)
        sim.integrator.setTemperature(temperature * unit.kelvin)
        sim.context.setVelocitiesToTemperature(temperature * unit.kelvin, seed)
        log.update(hpht_ps=steps * dt / 1000, hpht_pressure_bar=ph, hpht_temperature_k=th)
    sim.currentStep = 0
    sim.context.setTime(0.0)
    lz = sim.context.getState().getPeriodicBoxVectors()[2][2].value_in_unit(unit.nanometer)
    log.update(pressure_bar=p0, barostat_frequency=baro.getFrequency(), box_z_after_preparation_nm=lz)
    print(f'z-NPT preparation done: box z {lz:.3f} nm; production at {p0:g} bar, z free', flush=True)
    return log


def ghost_solute(es):
    """Pure-liquid reference (cfg ghost_solute: true; 2026-11-16, Wang 2025 phi0_H2O): the
    builder needs a solute, so a neutral one is built and then switched off: zero charges
    and LJ epsilon (NonbondedForce particles and its exceptions), removed from the field
    force. It keeps its mass and bonded terms, so it is an intact molecule that interacts
    with nothing (it still counts in the kinetic degrees of freedom)."""
    import openmm as mm
    nb = T.nonbonded_force(es.system)
    sol = set(es.indices.get('solute') or [])
    q0 = 0.0
    for i in sol:
        q, sig, eps = nb.getParticleParameters(i)
        q0 += q.value_in_unit(mm.unit.elementary_charge)
        nb.setParticleParameters(i, 0.0, sig, 0.0)
    if abs(q0) > 1e-6:
        raise ValueError(f'ghost_solute: use a neutral solute (charge {q0:+.3f} e)')
    n_ex = 0
    for k in range(nb.getNumExceptions()):
        i, j, qq, sig, eps = nb.getExceptionParameters(k)
        if i in sol or j in sol:
            nb.setExceptionParameters(k, i, j, 0.0, sig, 0.0); n_ex += 1
    for f in es.system.getForces():
        if isinstance(f, mm.CustomExternalForce) and f.getName() == T.FIELD_FORCE_NAME:
            for k in range(f.getNumParticles()):
                i, par = f.getParticleParameters(k)
                if i in sol:
                    f.setParticleParameters(k, i, [0.0])
    es.metadata['ghost_solute'] = dict(n_atoms=len(sol), exceptions_zeroed=n_ex)
    print(f'ghost solute: {len(sol)} atoms, {n_ex} exceptions switched off (pure-liquid reference)')
    return es


def ghost_counterions(es):
    """Switch off the neutralising counterion(s) (cfg ghost_counterion: true; 2026-11-16, user:
    charged solute in water without a mobile ion). The builder must add one to pass its
    neutrality check; here its charge and LJ epsilon are set to 0 in NonbondedForce, it is
    removed from the field force and given zero mass (it stays put and interacts with
    nothing). The box then carries the solute's net charge; OpenMM PME treats that with a
    uniform neutralising background (no force on any particle). Only for pure water, where
    the counterion is the only ion ('anion' or 'cation' indices)."""
    import openmm as mm
    nb = T.nonbonded_force(es.system)
    ions = list(es.indices.get('anion') or []) + list(es.indices.get('cation') or [])
    if not ions or len(ions) > 4:
        raise ValueError(f'ghost_counterion: expected 1-4 counterion atoms, found {len(ions)}')
    q0 = 0.0
    for i in ions:
        q, sig, eps = nb.getParticleParameters(i)
        q0 += q.value_in_unit(mm.unit.elementary_charge)
        nb.setParticleParameters(i, 0.0, sig, 0.0)
        es.system.setParticleMass(i, 0.0)
    for f in es.system.getForces():
        if isinstance(f, mm.CustomExternalForce) and f.getName() == T.FIELD_FORCE_NAME:
            for k in range(f.getNumParticles()):
                i, par = f.getParticleParameters(k)
                if i in ions:
                    f.setParticleParameters(k, i, [0.0])
    es.metadata['ghost_counterions'] = dict(indices=ions, charge_removed_e=q0,
                                            net_charge_e=-q0, note='PME uniform background')
    print(f'ghost counterion(s) {ions}: charge {q0:+.3f} e switched off; box net charge {-q0:+.3f} e')
    return es


def main(cfg):
    out = Path(cfg['output'])
    out.mkdir(parents=True, exist_ok=True)
    solute_pdb, solute_xml = solute_files(cfg, out)                      # 1
    species, proportions = electrolyte(cfg)                              # 3
    bulk = cfg.get('system', 'electrode') == 'bulk_reference'
    graphite = None if bulk else make_graphite(cfg)                      # 2
    height, window = (None, None) if bulk else leg_window(cfg, graphite)
    flip, rotation = (False, None) if bulk else solute_orientation(cfg, solute_pdb, solute_xml)
    if bulk:
        print(f"bulk reference: {cfg['bulk_box_nm']} nm cubic box; electrolyte: {species}")
    else:
        print(f"graphite: {graphite.kind}, {graphite.n_atoms} atoms, {graphite.carbon_set}, "
              f"{graphite.model.mobility}; electrolyte: {species}; restraint: {window}")

    def make(rho, field):                                                # 3 + 4
        if bulk:
            return T.build_reference(solute_pdb, solute_xml, species, proportions,
                                     box_nm=cfg['bulk_box_nm'], liquid_density=rho,
                                     counterion=cfg.get('counterion', 'K'), seed=cfg.get('seed'),
                                     cutoff_nm=cfg.get('cutoff_nm', 1.0),
                                     switch_nm=cfg.get('switch_nm', 0.9))
        return T.build(graphite, solute_pdb, solute_xml, species, proportions,
                       box_z_nm=cfg['box_z_nm'], liquid_density=rho, height_nm=height,
                       solute_xy_nm=cfg.get('solute_xy_nm'), flip=flip, rotation=rotation,
                       field_V_per_nm=field, counterion=cfg.get('counterion', 'K'),
                       seed=cfg.get('seed'), height_window=window, plane_atoms=cfg.get('plane_atoms'),
                       cutoff_nm=cfg.get('cutoff_nm', 1.0), switch_nm=cfg.get('switch_nm', 0.9))

    rho, calibration = cfg['liquid_density'], None
    md = cfg.get('md') or {}
    if cfg.get('liquid_density_from'):                    # a calibration run's calibration.json
        cal = json.loads(Path(cfg['liquid_density_from']).read_text())
        key = 'bulk_density_kg_m3' if bulk else 'liquid_density'
        rho = float(cal[key])
        print(f"liquid_density {rho:.1f} kg/m3 from {cfg['liquid_density_from']} ({key})")
    if cfg.get('calibrate') and not bulk:
        c = cfg.get('calibration') or {}
        ref = T.bulk_reference(species, proportions, platform=md.get('platform'), liquid_density=rho,
                               box_nm=c.get('bulk_box_nm', 4.0), equil_ps=c.get('bulk_equil_ps', 200),
                               frames=c.get('bulk_frames', 100), stride=c.get('bulk_stride', 250))
        bulk_system = ref['system'].system
        mass = sum(bulk_system.getParticleMass(i).value_in_unit(unit.dalton)
                   for i in range(bulk_system.getNumParticles()))
        bulk_density = mass / ref['volume_nm3'] * T.DA_PER_NM3_TO_KG_PER_M3
        rho, history = T.calibrate(lambda r: make(r, 0.0), ref['slot_density'], rho,
                                   iterations=c.get('iterations', 3), tol=c.get('tol', 0.005),
                                   equil_ps=c.get('equil_ps', 100), frames=c.get('frames', 50),
                                   stride=c.get('stride', 500), platform=md.get('platform'))
        calibration = dict(bulk_slot_density=ref['slot_density'], bulk_density_kg_m3=bulk_density,
                           bulk_volume_nm3=ref['volume_nm3'], liquid_density=rho, history=history,
                           settings=c)
        print(f'calibrated liquid_density = {rho:.1f} kg/m3 (reuse it for every field); '
              f'bulk NPT density {bulk_density:.1f} kg/m3')
    es = make(rho, 0.0 if bulk else cfg['field_V_per_nm'])
    es.metadata['solute']['orientation'] = cfg.get('orientation')
    if cfg.get('ghost_counterion'):
        ghost_counterions(es)
    if cfg.get('ghost_solute'):
        ghost_solute(es)

    npt = md.get('pressure_coupling') == 'z'
    if md.get('pressure_coupling') not in (None, 'z'):
        raise ValueError("md.pressure_coupling must be None or 'z'")
    if md.get('minimize', True) or md.get('steps'):                     # 5
        prep = None
        if npt:
            T.make_slab_rigid_for_barostat(es)                          # no-op without graphite
            if md.get('fep_windows', 'fixed') not in ('fixed', 'z'):
                raise ValueError("md.fep_windows must be 'fixed' or 'z'")
        sim = T.make_simulation(es, timestep_fs=md.get('timestep_fs', 2.0), seed=cfg.get('seed') or 1,
                                platform=md.get('platform'), minimize=md.get('minimize', True),
                                barostat_z=0 if npt else None)
        if npt:
            check = T.check_field_scaling(sim, es)
            if check:
                print(f"field energy under z scaling OK on {check['platform']}: "
                      f"{check['delta_kj_mol']:+.3f} vs {check['expected_kj_mol']:+.3f} kJ/mol")
            prep = equilibrate_z_npt(sim, md, out, seed=cfg.get('seed') or 1)
            prep['field_scaling_check'] = check
        heights = profiles = tilts = ring_heights = None
        output = None
        if md.get('steps') and md.get('report_ps'):              # water_4V-style trajectory + observables
            every = max(1, round(md['report_ps'] * 1000 / md.get('timestep_fs', 2.0)))
            every = max(100, 100 * round(every / 100))           # a multiple of the height sampling stride
            dcd_every = every * max(1, round((md.get('dcd_ps') or md['report_ps']) / md['report_ps']))
            output = MDOutput(es, sim, out, every, md.get('timestep_fs', 2.0), dcd_every=dcd_every,
                              plane_atoms=cfg.get('plane_atoms'))
            output.record()                                      # frame 0 = initial.pdb
        try:
            if md.get('steps') and window:
                heights, profiles, tilts, ring_heights = solute_heights(
                    es, sim, int(md['steps']), output=output, lz_range=1.15 if npt else 1.0,
                    plane_atoms=cfg.get('plane_atoms'))
            elif md.get('steps') and output:
                left = int(md['steps'])
                while left > 0:
                    k = min(output.every, left)
                    sim.step(k)
                    left -= k
                    if sim.currentStep % output.every == 0:
                        output.record()
            elif md.get('steps'):
                sim.step(int(md['steps']))
        finally:
            if output:
                output.close()
        state = sim.context.getState(getPositions=True)
        es.positions = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        es.metadata['pre_fep_md'] = dict(minimized=md.get('minimize', True), steps=md.get('steps', 0),
                                         timestep_fs=md.get('timestep_fs', 2.0),
                                         ensemble='z-NPT (x, y fixed)' if npt else 'NVT')
        if npt:                     # carry the final box into system.xml/pdb; FEP windows keep z-NPT
            box = state.getPeriodicBoxVectors()
            es.topology.setPeriodicBoxVectors(box)
            es.system.setDefaultPeriodicBoxVectors(*box)
            lz_end = box[2][2].value_in_unit(unit.nanometer)
            es.metadata['box_nm_initial'] = list(es.metadata['box_nm'])
            es.metadata['box_nm'] = [es.metadata['box_nm'][0], es.metadata['box_nm'][1], lz_end]
            lz_mean = float(np.mean(output.box_z)) if output and output.box_z else lz_end
            g = es.metadata.get('graphite') or {}
            if g.get('top_face_z') is not None:        # profiles average over the moving box
                g['bottom_face_z_initial'] = g['bottom_face_z']
                g['bottom_face_z'] = lz_mean - g['top_face_z']
            fep_windows = md.get('fep_windows', 'fixed')
            if fep_windows == 'z':
                es.system.addForce(T.z_barostat(298.15, md.get('pressure_bar', 1.0),
                                                md.get('barostat_frequency', 25), cfg.get('seed') or 1))
                es.metadata['pressure_coupling'] = 'z'
            es.metadata['pre_fep_md']['z_npt'] = dict(
                prep, box_z_final_nm=lz_end, box_z_mean_nm=lz_mean,
                fep_windows=('z-NPT, same barostat' if fep_windows == 'z'
                             else 'fixed volume at the last MD box'))
        if output:
            es.metadata['pre_fep_md']['outputs'] = dict(
                report_interval_steps=output.every, dcd_interval_steps=output.dcd_every,
                report_ps=output.every * md.get("timestep_fs", 2.0) / 1000,
                files=['initial.pdb', 'trajectory.dcd (frames 1..N, unwrapped)', 'observables.csv (frames 0..N)'],
                contact_cutoff_nm=CONTACT_NM, field_force_group=output.field_group)
        if heights is not None:
            outside = float(((heights < window['lo']) | (heights > window['hi'])).mean())
            es.metadata['pre_fep_md']['solute_height_nm'] = dict(
                min=float(heights.min()), mean=float(heights.mean()), max=float(heights.max()),
                fraction_outside_window=outside)
            if ring_heights is not None:       # ring centre of a side-chain molecule (contact check)
                es.metadata['pre_fep_md']['ring_height_nm'] = dict(
                    min=float(ring_heights.min()), mean=float(ring_heights.mean()), max=float(ring_heights.max()),
                    plane_atoms=list(cfg['plane_atoms']))
                print(f"ring centre height {ring_heights.min():.2f}..{ring_heights.max():.2f} nm "
                      f"(mean {ring_heights.mean():.2f})")
            if profiles:
                es.metadata['pre_fep_md']['ion_profiles'] = profiles
            tail = tilts[-max(1, len(tilts) // 10):]
            es.metadata['pre_fep_md']['solute_tilt_deg'] = dict(
                first=float(tilts[0]), mean=float(tilts.mean()), last_10pct_mean=float(tail.mean()),
                max=float(tilts.max()), fraction_below_30=float((tilts < 30).mean()))
            print(f"solute tilt {tilts[0]:.0f} deg at start, {tail.mean():.0f} deg over the last 10 % "
                  f"(0 flat, 90 edge-on)")
            print(f"solute height {heights.min():.2f}..{heights.max():.2f} nm (mean "
                  f"{heights.mean():.2f}); window {window['lo']:.2f}..{window['hi']:.2f}; "
                  f"outside {100 * outside:.0f} % of frames")

    folder = T.write_fep_inputs(es, out, overwrite=True)                 # 6
    (folder / 'config.json').write_text(json.dumps(cfg, indent=1, default=str))
    if calibration:
        (folder / 'calibration.json').write_text(json.dumps(calibration, indent=1, default=float))
    liq = es.metadata['liquid']
    print(f"wrote {folder}: {es.metadata['n_particles']} particles, {liq['species_counts']}")
    for w in es.metadata['warnings']:
        print('WARNING:', w[:200])

    fep = dict(cfg.get('fep') or {})                                     # 7
    if fep.pop('run', False):
        results = T.run_fep(folder, **{k: v for k, v in fep.items() if v is not None})
        print('FEP free energy (kJ/mol):', results['free_energy'])
    return es


def adsorption(surface_dir, mid_dir, temperature=298.15):
    res = T.adsorption_free_energy(surface_dir, mid_dir, temperature)
    (Path(surface_dir) / 'adsorption.json').write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))
    return res


def load_config(argv):
    cfg = copy.deepcopy(CONFIG)
    if argv:
        user = json.loads(Path(argv[0]).read_text())
        for k, v in user.items():
            if isinstance(v, dict) and isinstance(cfg.get(k), dict):
                cfg[k].update(v)
            else:
                cfg[k] = v
    return cfg


if __name__ == '__main__':
    if sys.argv[1:2] == ['--adsorption']:
        adsorption(*sys.argv[2:4], *[float(t) for t in sys.argv[4:5]])
    else:
        main(load_config(sys.argv[1:]))

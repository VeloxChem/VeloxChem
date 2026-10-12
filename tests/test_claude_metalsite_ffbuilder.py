from pathlib import Path
from copy import deepcopy
from types import SimpleNamespace
import sys
import numpy as np
import random
import pytest

from veloxchem.molecule import Molecule
from veloxchem.outputstream import OutputStream
from veloxchem.metalsite.metalsiteffbuilder import MetalSiteForceFieldBuilder, Stage
from veloxchem.metalsite import util
from veloxchem.metalsite import qm
from veloxchem.metalsite.builder import ActiveSiteBuilder
from veloxchem.metalsite.qm import QmParameterizer


class TestMetalSiteForceFieldBuilder:

    @pytest.fixture(autouse=True)
    def _placed_hydrogens(self):
        """
        Starts every test from the same random state.

        protonate places its hydrogens through Modeller.addHydrogens, which
        draws from Python's global random stream: seven random() calls before
        it move a hydrogen by 1.9 A. Left alone, that makes a result depend on
        how many draws every earlier test happened to make, so a test that
        passes on its own fails in the suite and back again. Seeding here does
        not make two calls inside one test identical -- the stream still
        advances -- it makes the whole file reproducible from run to run.
        """

        random.seed(20250821)

    @pytest.fixture(autouse=True)
    def _kept_stdout(self, monkeypatch):
        """
        Keeps an output stream from closing the stdout pytest handed it.

        OutputStream.close() closes whatever stream it holds unless it is the
        current sys.stdout, and pytest swaps sys.stdout per test, so a stream
        collected after that swap -- one held alive by the traceback of a
        failing test, say -- closes the capture of every test after it.
        """

        monkeypatch.setattr(OutputStream, 'close', lambda self: None)

    @staticmethod
    def coordination(binding_modes):
        """
        What a set of binding modes says, without the atom indices that the
        protonation renumbers.

        :param binding_modes:
            The binding modes.

        :return:
            A set of (residue, atom, mode) tuples.
        """

        return {(ligand['residue'], ligand['atom'], ligand['mode'])
                for ligand in binding_modes['ligands']}

    @staticmethod
    def printing():
        """
        A stream that reaches the capture of a test that reads output.

        :return:
            The output stream.
        """

        return OutputStream(sys.stdout)

    @staticmethod
    def shells(ostream=None, **settings):
        """
        The phase classes, silent unless a stream is given, on the builder's
        default settings unless a test names one.

        A phase class holds nothing but a communicator, a stream and its
        settings, so a test builds one where it needs it rather than sharing
        any -- the way the interface does, through its own factories, so
        the defaults are the interface's and stated once.

        :param ostream:
            The output stream, for a test that reads output.
        :param settings:
            Settings to assign on the builder the shells are made from. The
            pruning is off unless a test asks for it: the Hessians here are
            made up, so the zeros they produce say nothing about a bridge.

        :return:
            A namespace holding one instance of each phase class.
        """

        if ostream is None:
            ostream = OutputStream(None)

        builder = MetalSiteForceFieldBuilder(ostream=ostream)
        builder.prune_weak_bridge_bonds = False
        for key, value in settings.items():
            setattr(builder, key, value)

        return SimpleNamespace(sites=builder._sites(),
                               qm=builder._qm(),
                               enzyme=builder._enzyme())

    # Metal-ligand distances in Angstrom of the reference structure, which is
    # a designed binuclear zinc site truncated to its coordinating residues.
    expected_contacts = {
        ('GLU73', 'OE2'): (2.09, 'bridging_single'),
        ('CYS77', 'SG'): (2.25, 'monodentate'),
        ('ASP130', 'OD1'): (2.50, 'bidentate'),
        ('ASP130', 'OD2'): (2.57, 'bidentate'),
        ('ASP134', 'OD2'): (1.96, 'monodentate'),
        ('HIS187', 'NE2'): (2.07, 'monodentate'),
        ('HIS191', 'NE2'): (2.17, 'monodentate'),
    }

    @staticmethod
    def structure_file():

        here = Path(__file__).parent

        return str(here / 'data' / 'zn_site_ECDDHH.pdb')

    @staticmethod
    def random_hessian(n_atoms, seed=0):
        """
        A symmetric matrix of the right shape. The force constants it produces
        are meaningless; these tests check which terms are fitted and what
        equilibrium values they get, not the numbers themselves.
        """

        size = 3 * n_atoms
        rng = np.random.default_rng(seed)
        amat = rng.normal(0.0, 0.05, (size, size))

        return amat + amat.T + np.eye(size) * 0.5

    def seeded_forcefield(self,
                          active_site,
                          partial_charges=None,
                          ostream=None,
                          bond_equilibria=None,
                          **settings):
        """
        The seeded force field of a site: metal terms off the geometry, on
        D4 charges unless a test says otherwise.

        :param active_site:
            The active site.
        :param partial_charges:
            The charges, D4 when None.
        :param ostream:
            The output stream, for a test that reads output.
        :param bond_equilibria:
            Per-bond equilibria, as build_forcefield takes them.
        :param settings:
            Settings of the builder the seeded build reads.

        :return:
            The force field generator.
        """

        if partial_charges is None:
            partial_charges = util.d4_charges(active_site)

        return self.shells(ostream, **settings).sites.build_forcefield(
            active_site, partial_charges, bond_equilibria=bond_equilibria)

    def fitted_forcefield(self,
                          active_site,
                          hessian,
                          partial_charges=None,
                          ostream=None,
                          protected_bonds=None,
                          **settings):
        """
        A seeded force field fitted against a Hessian, with the weak bridge
        pruning off unless a test asks for it: the Hessians here are made
        up.

        :param active_site:
            The active site.
        :param hessian:
            The Hessian.
        :param partial_charges:
            The charges, D4 when None.
        :param ostream:
            The output stream, for a test that reads output.
        :param protected_bonds:
            The bonds the pruning leaves alone, as fit_forcefield takes them.
        :param settings:
            Settings of the builder either pass reads.

        :return:
            The force field generator.
        """

        shells = self.shells(ostream, **settings)

        forcefield = shells.sites.build_forcefield(
            active_site,
            util.d4_charges(active_site)
            if partial_charges is None else partial_charges)

        return shells.qm.fit_forcefield(active_site,
                                        forcefield,
                                        hessian,
                                        protected_bonds=protected_bonds)

    @staticmethod
    def builder(tmp_path=None, name='run', **settings):
        """
        A builder that writes where the test says and never runs QM.

        The output stream goes nowhere: an OutputStream around the stdout
        pytest hands out closes it when it is collected, which breaks the
        capture of every test after this one.

        :param tmp_path:
            The folder to write into. A folder of its own under it is used.
        :param name:
            The name of that folder.
        :param settings:
            Settings to assign on the builder.

        :return:
            The builder.
        """

        builder = MetalSiteForceFieldBuilder(ostream=OutputStream(None))

        # the Hessians in these tests are made up, so the pruning that reads
        # their zeros is off unless a test asks for it
        builder.prune_weak_bridge_bonds = False
        # no QM at all: the geometry is whatever the crude pass leaves, and
        # any attempt to run a driver would fail on the basis set
        builder.do_qm_optimization = False
        builder.basis_set_label = 'no-such-basis'
        # the basis above does not stop RESP, which is given only the
        # molecule so that it runs its own HF/6-31G* -- it ignores
        # basis_set_label entirely and was costing 98 s per fit, about
        # seven eighths of the whole suite. D4 is the other branch of
        # calculate_partial_charges and takes 0.2 ms; a test that is about
        # charges passes partial_charges in or turns this back on.
        builder.do_resp = False
        if tmp_path is not None:
            builder.folder = str(Path(tmp_path) / name)

        for key, value in settings.items():
            setattr(builder, key, value)

        return builder

    def derive(self,
               topology,
               positions,
               request=None,
               coordinating_residues=None,
               **settings):
        """
        The coordination of a structure at the builder's detection settings,
        under an empty request unless one is given.

        :return:
            The binding modes.
        """

        if request is None:
            request = util.empty_request()

        return self.shells(**settings).sites.derive_binding_modes(
            topology,
            positions,
            request,
            coordinating_residues=coordinating_residues)

    @staticmethod
    def request_of(binding_modes):
        """
        The request a derived set of binding modes was produced under.

        derive_binding_modes copies the decisions into what it returns, so
        this reads them back out to feed the next derivation.
        """

        return {
            'manual_bonds':
            deepcopy(binding_modes.get('manual_bonds', [])),
            'variants':
            dict(binding_modes.get('variants', {})),
            'coordinating_residues':
            list(binding_modes.get('coordinating_residues', [])),
            'extra_residues':
            list(binding_modes.get('extra_residues', [])),
            'excluded_residues':
            list(binding_modes.get('excluded_residues', [])),
        }

    def bonded(self, binding_modes, topology, positions, *args, **kwargs):
        """
        Adds a metal bond and derives the coordination again, which is what
        the builder does: the edit is recorded in the request and the
        coordination is a function of it.

        :return:
            The new binding modes.
        """

        request = self.shells().sites.add_metal_bond(
            self.request_of(binding_modes), binding_modes, topology, positions,
            *args, **kwargs)

        return self.shells().sites.derive_binding_modes(topology, positions,
                                                        request)

    def unbonded(self, binding_modes, topology, positions, *args, **kwargs):
        """
        Removes a metal bond and derives the coordination again.

        :return:
            The new binding modes.
        """

        request = self.shells().sites.remove_metal_bond(
            self.request_of(binding_modes), binding_modes, *args, **kwargs)

        return self.shells().sites.derive_binding_modes(topology, positions,
                                                        request)

    def protonated(self, topology, positions, binding_modes, **settings):
        """
        Protonates and works the coordination out again on the result.

        Nothing is remapped across the renumbering: what survives it is the
        request, which is keyed on residue indices and atom names, and the
        coordination is derived afresh on the protonated structure. This is
        the two-step the pipeline itself runs, kept in one place so the tests
        do not each spell it out.

        :return:
            The protonated topology, its positions and the coordination.
        """

        request = {
            'manual_bonds': binding_modes.get('manual_bonds', []),
            'variants': {},
            'coordinating_residues': binding_modes.get('coordinating_residues',
                                                       []),
            'extra_residues': binding_modes.get('extra_residues', []),
            'excluded_residues': binding_modes.get('excluded_residues', []),
        }

        topology, positions, variants, _ = self.shells().sites.protonate(
            topology, positions, binding_modes, None)
        request['variants'] = variants

        binding_modes = self.shells(**settings).sites.derive_binding_modes(
            topology, positions, request)

        return topology, positions, binding_modes

    def prepared(self, **settings):
        """
        Runs the structural part of the pipeline, which needs no QM.

        The core functions keep no state, so every product of the pipeline is
        returned for the tests to pass back in explicitly.

        :param settings:
            Settings of the builder the coordination detection reads.

        :return:
            A dictionary with everything the structural steps produced.
        """

        sites = self.shells(prepare_protein=False, **settings).sites
        topology, positions = sites.load_and_prepare_protein(
            self.structure_file())
        request = util.empty_request()
        binding_modes = sites.derive_binding_modes(topology, positions, request)
        # the coordination is derived again on the protonated structure
        # rather than remapped onto it
        topology, positions, variants, _ = sites.protonate(
            topology, positions, binding_modes, None)
        request['variants'] = variants
        binding_modes = sites.derive_binding_modes(topology, positions, request)
        active_site = sites.extract_active_site(topology, positions,
                                                binding_modes)

        return {
            'topology': topology,
            'positions': positions,
            'request': request,
            'binding_modes': binding_modes,
            'active_site': active_site,
        }

    @staticmethod
    def pushed_structure(target,
                         resid='130',
                         names=('OD1', 'OD2', 'CG'),
                         push=1.6):
        """
        Moves one sidechain away from its metal, far enough that no cutoff
        reaches it any more.

        :param target:
            Where to write the altered structure.
        :param resid:
            The residue to move.
        :param names:
            The atoms of it to move.
        :param push:
            How far to move them, in Angstrom.

        :return:
            The path written.
        """

        import openmm.app as mmapp
        import openmm.unit as mmunit

        source = TestMetalSiteForceFieldBuilder.structure_file()
        pdb = mmapp.PDBFile(source)
        positions = pdb.getPositions(asNumpy=True).value_in_unit(
            mmunit.angstrom)

        metals = [
            atom.index for atom in pdb.topology.atoms()
            if atom.element is not None and atom.element.symbol == 'Zn'
        ]

        for atom in pdb.topology.atoms():
            if atom.residue.id != resid or atom.name not in names:
                continue
            nearest = min(metals,
                          key=lambda metal: np.linalg.norm(positions[atom.index]
                                                           - positions[metal]))
            direction = positions[atom.index] - positions[nearest]
            positions[
                atom.index] += direction / np.linalg.norm(direction) * push

        target = Path(target)
        with target.open('w') as fh:
            mmapp.PDBFile.writeFile(pdb.topology,
                                    positions * mmunit.angstrom,
                                    fh,
                                    keepIds=True)

        return str(target)

    def test_a_residue_can_be_asked_for_however_far_away_it_is(self, tmp_path):
        """
        A design whose ligand has not been placed properly yet has that ligand
        outside every cutoff, and no amount of reviewing the distances will
        bring it in. coordinating_residues says it coordinates anyway, through its
        closest sidechain donor, and the rest of the pipeline treats it as the
        ligand it was declared to be.
        """

        pytest.importorskip('openmm')

        structure = self.pushed_structure(tmp_path / 'far.pdb')

        builder = MetalSiteForceFieldBuilder()
        builder.ostream.mute()
        topology, positions = self.shells(
            prepare_protein=False).sites.load_and_prepare_protein(structure)

        # too far for either cutoff to find it
        found = self.derive(topology, positions)
        assert 'ASP130' not in {
            ligand['residue']
            for ligand in found['ligands']
        }

        binding_modes = self.derive(topology, positions,
                                    coordinating_residues=['130'])

        forced = [
            ligand for ligand in binding_modes['ligands']
            if ligand['residue'] == 'ASP130'
        ]
        assert len(forced) == 1
        assert forced[0]['atom'] == 'OD1'
        assert len(forced[0]['metals']) == 1
        assert forced[0]['distances'][0] > builder._sites().report_cutoff
        assert forced[0]['mode'] == 'monodentate'

        # the request is recorded, and said out loud
        assert binding_modes['coordinating_residues'] == [2]
        assert any('coordinating_residues' in note
                   for note in binding_modes['notes'])

        # and it reaches the active site as a bonded ligand
        topology, positions, binding_modes = self.protonated(
            topology, positions, binding_modes)
        active_site = self.shells().sites.extract_active_site(
            topology, positions, binding_modes)
        assert 'ASP130' in active_site['residues']

        # protonate works on a copy and remaps the indices, so the ligand has
        # to be read back out of what it returned
        remapped = next(ligand for ligand in binding_modes['ligands']
                        if ligand['residue'] == 'ASP130')
        site_index = next(site for site, top in active_site['atom_map'].items()
                          if top == remapped['index'])
        assert any(active_site['connectivity_matrix'][site_index][metal]
                   for metal in active_site['metal_indices'])

    def test_a_residue_asked_for_is_not_dropped_by_a_relaxation(self, tmp_path):
        """
        A distance is not what made it a ligand, so a geometry that moves it
        further away does not unmake it. update_binding_modes reads the
        request back out of the dictionary it is handed.
        """

        pytest.importorskip('openmm')

        structure = self.pushed_structure(tmp_path / 'far.pdb')

        builder = MetalSiteForceFieldBuilder()
        builder.ostream.mute()
        topology, positions = self.shells(
            prepare_protein=False).sites.load_and_prepare_protein(structure)
        binding_modes = self.derive(topology, positions,
                                    coordinating_residues=[130])
        topology, positions, binding_modes = self.protonated(
            topology, positions, binding_modes)
        active_site = self.shells().sites.extract_active_site(
            topology, positions, binding_modes)
        updated, _, changed = self.shells().sites.update_binding_modes(topology,
                                                                       active_site['molecule'],
                                                                       active_site,
                                                                       binding_modes)

        assert not changed
        assert 'ASP130' in {ligand['residue'] for ligand in updated['ligands']}
        assert updated['coordinating_residues'] == binding_modes[
            'coordinating_residues']

    def test_a_residue_asked_for_has_to_exist(self, tmp_path):
        """
        A residue that is named but not there is a typo, and a typo that
        quietly does nothing is worse than one that stops the run.
        """

        pytest.importorskip('openmm')

        structure = self.pushed_structure(tmp_path / 'far.pdb')

        builder = MetalSiteForceFieldBuilder()
        builder.ostream.mute()
        topology, positions = self.shells(
            prepare_protein=False).sites.load_and_prepare_protein(structure)

        # the id, the number and the label all name the same residue
        for request in ('130', 130, 'ASP130'):
            binding_modes = self.derive(topology, positions,
                                        coordinating_residues=[request])
            assert binding_modes['coordinating_residues'] == [2]

        with pytest.raises(Exception, match='not a residue'):
            self.derive(topology, positions, coordinating_residues=['ASP999'])

    def detected(self, **settings):
        """
        The structure and the coordination as suggested, before protonation.

        The manual edits are made where a review is made: on what
        derive_binding_modes proposed, before protonate has renumbered
        anything.

        :param settings:
            Settings of the builder the coordination detection reads.

        :return:
            A dictionary with the topology, the positions and the binding
            modes.
        """

        sites = self.shells(prepare_protein=False, **settings).sites
        topology, positions = sites.load_and_prepare_protein(
            self.structure_file())
        binding_modes = sites.derive_binding_modes(topology, positions,
                                                   util.empty_request())

        return {
            'topology': topology,
            'positions': positions,
            'binding_modes': binding_modes,
        }

    @staticmethod
    def metal_of(binding_modes, residue):
        """
        The metal a residue of the fixture coordinates.

        The two zinc ions are told apart by which residues reach them, so that
        nothing here depends on the atom indices the file happens to give.

        :param binding_modes:
            The binding modes.
        :param residue:
            The residue label, as 'CYS77'.

        :return:
            The atom index of the metal.
        """

        for ligand in binding_modes['ligands']:
            if ligand['residue'] == residue:
                return ligand['metals'][0]

        raise AssertionError(f'{residue} coordinates nothing')

    @staticmethod
    def contacts_of(binding_modes, residue):
        """
        The contacts a residue makes, as (atom, metals, mode) tuples.
        """

        return [(ligand['atom'], tuple(ligand['metals']), ligand['mode'])
                for ligand in binding_modes['ligands']
                if ligand['residue'] == residue]

    def test_a_bond_added_by_hand_reclassifies_the_residue(self):
        """
        A second oxygen on the metal makes the carboxylate bidentate, which
        the caller should not have to work out and write down itself.
        """

        pytest.importorskip('openmm')

        detected = self.detected()
        binding_modes = detected['binding_modes']
        metal = self.metal_of(binding_modes, 'ASP134')

        assert self.contacts_of(binding_modes,
                                'ASP134') == [('OD2', (metal, ), 'monodentate')]

        edited = self.bonded(binding_modes,
                             detected['topology'],
                             detected['positions'],
                             134,
                             metal,
                             atom='OD1')

        assert self.contacts_of(edited, 'ASP134') == [
            ('OD2', (metal, ), 'bidentate'),
            ('OD1', (metal, ), 'bidentate'),
        ]

        added = next(
            ligand for ligand in edited['ligands']
            if ligand['atom'] == 'OD1' and ligand['residue'] == 'ASP134')

        # everything a detected contact carries, not only the two keys the
        # bond itself is made of
        assert added['res_name'] == 'ASP'
        assert added['chain'] == binding_modes['ligands'][0]['chain']
        assert len(added['distances']) == len(added['metals'])
        assert added['distances'][0] == pytest.approx(3.35, abs=0.01)

        metal_entry = next(entry for entry in edited['metals']
                           if entry['index'] == metal)

        # the bond is written down in terms nothing downstream renumbers
        assert edited['manual_bonds'] == [{
            'res_index':
            added['res_index'],
            'atom':
            'OD1',
            'metal_res_index':
            metal_entry['res_index'],
            'action':
            'add',
        }]

        # and the dictionary that was handed in is untouched
        assert self.contacts_of(binding_modes,
                                'ASP134') == [('OD2', (metal, ), 'monodentate')]
        assert binding_modes['manual_bonds'] == []

    def test_an_added_bond_outranks_the_asymmetry_rule(self):
        """
        The far oxygen of a lopsided carboxylate is normally read as pointing
        away, and that rule reading an unrelaxed structure must not undo an
        edit made on purpose. What the chemistry forbids still applies.
        """

        pytest.importorskip('openmm')

        detected = self.detected()
        binding_modes = detected['binding_modes']
        metal = self.metal_of(binding_modes, 'ASP134')

        # 1.96 A against 3.35 A, far past the 0.75 A the rule allows
        edited = self.bonded(binding_modes,
                             detected['topology'],
                             detected['positions'],
                             'ASP134',
                             metal,
                             atom='OD1')

        assert len(self.contacts_of(edited, 'ASP134')) == 2
        assert any('a manual edit asked for them' in note
                   for note in edited['notes'])

    def test_an_added_bond_can_make_a_bridge(self):
        """
        One donor atom reaching both metals is a bridge, and a carboxylate is
        allowed to make one.
        """

        pytest.importorskip('openmm')

        detected = self.detected()
        binding_modes = detected['binding_modes']
        near = self.metal_of(binding_modes, 'ASP134')
        far = self.metal_of(binding_modes, 'ASP130')

        edited = self.bonded(binding_modes,
                             detected['topology'],
                             detected['positions'],
                             134,
                             far,
                             atom='OD2')

        contacts = self.contacts_of(edited, 'ASP134')
        assert contacts == [('OD2', tuple(sorted([near,
                                                  far])), 'bridging_single')]

        contact = next(ligand for ligand in edited['ligands']
                       if ligand['residue'] == 'ASP134')

        # the distances stay alongside the metals they belong to
        assert contact['distances'][contact['metals'].index(near)] == (
            pytest.approx(1.96, abs=0.01))

    def test_the_donor_atom_is_resolved_when_the_choice_does_not_matter(self):
        """
        A cysteine has one donor, and the two oxygens of a carboxylate that
        binds nothing yet are interchangeable, so neither has to be named.
        """

        pytest.importorskip('openmm')

        detected = self.detected()
        binding_modes = detected['binding_modes']
        metal = self.metal_of(binding_modes, 'CYS77')

        # taken off and put back on: the sulfur is the only thing it could be
        without = self.unbonded(binding_modes, detected['topology'],
                                detected['positions'], 77)
        assert self.contacts_of(without, 'CYS77') == []

        again = self.bonded(without, detected['topology'],
                            detected['positions'], 77, metal)
        assert self.contacts_of(again,
                                'CYS77') == [('SG', (metal, ), 'monodentate')]

        # a carboxylate that coordinates nothing yet: either oxygen would do,
        # so the nearer one is taken without asking
        free = self.unbonded(binding_modes, detected['topology'],
                             detected['positions'], 130)
        metal_130 = self.metal_of(binding_modes, 'ASP130')
        restored = self.bonded(free, detected['topology'],
                               detected['positions'], 130, metal_130)

        assert self.contacts_of(restored, 'ASP130') == [('OD1', (metal_130, ),
                                                         'monodentate')]

    def test_an_ambiguous_donor_atom_is_asked_for(self):
        """
        Which ring nitrogen of a histidine binds is a chemical choice, and so
        is the second oxygen of a carboxylate that already binds through the
        first. Guessing either would be guessing at the force field.
        """

        pytest.importorskip('openmm')

        detected = self.detected()
        binding_modes = detected['binding_modes']
        metal = self.metal_of(binding_modes, 'HIS187')

        with pytest.raises(Exception, match='chemical choice') as excinfo:
            self.bonded(binding_modes, detected['topology'],
                        detected['positions'], 187, metal)

        # and it says what to pass instead
        assert 'ND1' in str(excinfo.value)
        assert 'NE2' in str(excinfo.value)

        with pytest.raises(Exception, match='chemical choice'):
            self.bonded(binding_modes, detected['topology'],
                        detected['positions'], 130,
                        self.metal_of(binding_modes, 'ASP130'))

    def test_a_bond_the_chemistry_forbids_is_refused(self):
        """
        An imidazole nitrogen has one lone pair in the ring plane, and the
        truncation is sidechain-only. Refusing beats accepting the edit and
        trimming it back out further down, where the caller would not see it.
        """

        pytest.importorskip('openmm')

        detected = self.detected()
        binding_modes = detected['binding_modes']
        other = self.metal_of(binding_modes, 'CYS77')

        with pytest.raises(Exception, match='can bridge'):
            self.bonded(binding_modes,
                        detected['topology'],
                        detected['positions'],
                        187,
                        other,
                        atom='ND1')

        with pytest.raises(Exception, match='sidechain-only'):
            self.bonded(binding_modes,
                        detected['topology'],
                        detected['positions'],
                        187,
                        other,
                        atom='N')

        with pytest.raises(Exception, match='lone pair'):
            self.bonded(binding_modes,
                        detected['topology'],
                        detected['positions'],
                        187,
                        other,
                        atom='CB')

    def test_the_atom_and_the_metal_have_to_be_what_they_are_called(self):
        """
        A name or an index that means something else is a typo, and a typo
        that silently bonds the wrong atom is the worst outcome of the three.
        """

        pytest.importorskip('openmm')

        detected = self.detected()
        binding_modes = detected['binding_modes']
        metal = self.metal_of(binding_modes, 'HIS187')

        with pytest.raises(Exception, match='no atom named'):
            self.bonded(binding_modes,
                        detected['topology'],
                        detected['positions'],
                        187,
                        metal,
                        atom='OD9')

        cysteine = next(ligand for ligand in binding_modes['ligands']
                        if ligand['residue'] == 'CYS77')

        with pytest.raises(Exception, match='is not part of HIS187'):
            self.bonded(binding_modes,
                        detected['topology'],
                        detected['positions'],
                        187,
                        metal,
                        atom=cysteine['index'])

        with pytest.raises(Exception, match='no metal center at atom index'):
            self.bonded(binding_modes, detected['topology'],
                        detected['positions'], 187, cysteine['index'])

        with pytest.raises(Exception, match='not a residue'):
            self.bonded(binding_modes, detected['topology'],
                        detected['positions'], 999, metal)

    def test_adding_a_bond_that_is_already_there_changes_nothing(self):
        """
        Asking twice is not an error, and the second call has nothing to do.
        """

        pytest.importorskip('openmm')

        detected = self.detected()
        binding_modes = detected['binding_modes']

        edited = self.bonded(binding_modes, detected['topology'],
                             detected['positions'], 77,
                             self.metal_of(binding_modes, 'CYS77'))

        # the request gains no record, so the derivation says exactly what it
        # said before. It is not the same object -- deriving always builds a
        # new one -- so the comparison is on what it holds
        assert edited['manual_bonds'] == []
        assert self.coordination(edited) == self.coordination(binding_modes)

    def test_removing_a_bond_needs_the_metal_only_when_it_bridges(self):
        """
        A residue holding one metal has one bond to remove; one holding two
        does not, and picking for the caller would pick wrong half the time.
        """

        pytest.importorskip('openmm')

        detected = self.detected()
        binding_modes = detected['binding_modes']

        with pytest.raises(Exception, match='bridges the metals'):
            self.unbonded(binding_modes, detected['topology'],
                          detected['positions'], 73)

        kept = self.metal_of(binding_modes, 'ASP130')
        dropped = self.metal_of(binding_modes, 'CYS77')

        edited = self.unbonded(binding_modes,
                               detected['topology'],
                               detected['positions'],
                               73,
                               metal=dropped)
        assert self.contacts_of(edited,
                                'GLU73') == [('OE2', (kept, ), 'monodentate')]

        # and a metal it does not reach is a mistake worth stopping on
        with pytest.raises(Exception, match='is not bonded to the metal'):
            self.unbonded(binding_modes,
                          detected['topology'],
                          detected['positions'],
                          130,
                          metal=dropped)

    def test_a_removal_can_drop_one_arm_or_the_whole_grip(self):
        """
        Naming an atom takes a bidentate back to a monodentate; leaving it out
        removes the residue's hold on that metal altogether.
        """

        pytest.importorskip('openmm')

        detected = self.detected()
        binding_modes = detected['binding_modes']
        metal = self.metal_of(binding_modes, 'ASP130')

        one_arm = self.unbonded(binding_modes,
                                detected['topology'],
                                detected['positions'],
                                130,
                                atom='OD1')
        assert self.contacts_of(one_arm,
                                'ASP130') == [('OD2', (metal, ), 'monodentate')]

        whole = self.unbonded(binding_modes, detected['topology'],
                              detected['positions'], 130)
        assert self.contacts_of(whole, 'ASP130') == []
        records = {(record['atom'], record['action'])
                   for record in whole['manual_bonds']}
        assert records == {('OD1', 'remove'), ('OD2', 'remove')}

        with pytest.raises(Exception, match='coordinates no metal'):
            self.unbonded(binding_modes, detected['topology'],
                          detected['positions'], 999)

    def test_a_manual_edit_survives_the_rest_of_the_pipeline(self):
        """
        protonate renumbers the atoms and update_binding_modes reads the
        coordination off the geometry again, so an edit that only rewrote the
        contact list would be undone before anything was fitted. The records
        name their atoms in terms that neither step touches.
        """

        pytest.importorskip('openmm')

        detected = self.detected()
        binding_modes = detected['binding_modes']
        metal = self.metal_of(binding_modes, 'ASP134')

        edited = self.bonded(binding_modes,
                             detected['topology'],
                             detected['positions'],
                             134,
                             metal,
                             atom='OD1')
        edited = self.unbonded(edited, detected['topology'],
                               detected['positions'], 191)

        topology, positions, edited = self.protonated(detected['topology'],
                                                      detected['positions'],
                                                      edited)
        active_site = self.shells().sites.extract_active_site(
            topology, positions, edited)

        # the residue that lost its only bond is not part of the site at all
        assert 'HIS191' not in active_site['residues']

        updated, _, changed = self.shells().sites.update_binding_modes(topology,
                                                                       active_site['molecule'],
                                                                       active_site, edited)

        assert not changed
        assert len(self.contacts_of(updated, 'ASP134')) == 2
        assert self.contacts_of(updated, 'HIS191') == []
        assert updated['manual_bonds'] == edited['manual_bonds']

    def test_a_manual_edit_survives_a_re_derivation(self):
        """
        An edit is a decision, so it lives in the request and every later
        derivation replays it. Deriving twice from the same request has to
        give the same coordination, or the edit would quietly come undone
        the next time anything asked.
        """

        pytest.importorskip('openmm')

        detected = self.detected()
        binding_modes = detected['binding_modes']

        edited = self.bonded(binding_modes,
                             detected['topology'],
                             detected['positions'],
                             134,
                             self.metal_of(binding_modes, 'ASP134'),
                             atom='OD1')

        again = self.derive(detected['topology'], detected['positions'],
                            self.request_of(edited))

        assert again['manual_bonds'] == edited['manual_bonds']
        assert self.contacts_of(again,
                                'ASP134') == self.contacts_of(edited, 'ASP134')

    def test_removing_a_bond_clears_the_request_that_would_restore_it(
            self, tmp_path):
        """
        coordinating_residues makes a residue a ligand whatever its distance, so a
        removal that left the request in place would be undone by the next
        re-detection with nothing said about it.
        """

        pytest.importorskip('openmm')

        structure = self.pushed_structure(tmp_path / 'far.pdb')

        builder = MetalSiteForceFieldBuilder()
        builder.ostream.mute()
        topology, positions = self.shells(
            prepare_protein=False).sites.load_and_prepare_protein(structure)
        binding_modes = self.derive(topology, positions,
                                    coordinating_residues=[130])

        assert binding_modes['coordinating_residues'] == [2]

        edited = self.unbonded(binding_modes, topology, positions, 130)
        assert edited['coordinating_residues'] == []
        assert self.contacts_of(edited, 'ASP130') == []

        topology, positions, edited = self.protonated(topology, positions,
                                                      edited)
        active_site = self.shells().sites.extract_active_site(
            topology, positions, edited)
        updated, _, _ = self.shells().sites.update_binding_modes(topology,
                                                                 active_site['molecule'],
                                                                 active_site, edited)

        assert self.contacts_of(updated, 'ASP130') == []

    def test_binding_modes_match_the_structure(self):

        pytest.importorskip('openmm')

        prepared = self.prepared()
        binding_modes = prepared['binding_modes']

        assert len(binding_modes['metals']) == 2
        assert all(metal['element'] == 'Zn'
                   for metal in binding_modes['metals'])

        found = {
            (ligand['residue'], ligand['atom']): ligand
            for ligand in binding_modes['ligands']
        }

        assert set(found) == set(self.expected_contacts)

        for key, (distance, mode) in self.expected_contacts.items():
            assert min(found[key]['distances']) == pytest.approx(distance,
                                                                 abs=0.01)
            assert found[key]['mode'] == mode

    def test_bridging_ligand_binds_both_metals(self):

        pytest.importorskip('openmm')

        prepared = self.prepared()

        bridging = [
            ligand for ligand in prepared['binding_modes']['ligands']
            if ligand['mode'] == 'bridging_single'
        ]

        assert len(bridging) == 1
        assert bridging[0]['residue'] == 'GLU73'
        # 2.79 A is past a naive 2.6 A cutoff, but it is a real bond
        assert max(bridging[0]['distances']) == pytest.approx(2.79, abs=0.01)

    @staticmethod
    def contact(res_name, atom, distances, res_index=0):
        """
        Builds one ligand contact as derive_binding_modes records it.

        _assign_binding_modes only reads the residue context and the metals a
        contact reached, so the coordination rules can be checked without a
        structure behind them.

        :param res_name:
            The residue name, which decides what the residue is allowed to do.
        :param atom:
            The name of the donor atom.
        :param distances:
            One distance per metal the atom is within the primary cutoff of.
        :param res_index:
            The residue index, which groups contacts of the same sidechain.

        :return:
            The ligand contact.
        """

        return {
            'residue': f'{res_name}1',
            'res_name': res_name,
            'res_index': res_index,
            'chain': 'A',
            'atom': atom,
            'index': 0,
            'metals': list(range(len(distances))),
            'distances': list(distances),
        }

    def test_a_carboxylate_oxygen_may_reach_two_metals(self):
        """
        A carboxylate is the reference case for a real bridge.

        Both contacts have to survive, since the second one becomes a metal
        bond of the force field.
        """

        ligands = [self.contact('GLU', 'OE2', [2.09, 2.79])]
        notes = []

        self.shells().sites._assign_binding_modes(ligands, notes)

        assert ligands[0]['mode'] == 'bridging_single'
        assert ligands[0]['metals'] == [0, 1]
        assert notes == []

    def test_a_histidine_nitrogen_is_trimmed_to_one_metal(self):
        """
        An imidazole nitrogen has one lone pair, in the ring plane, so it
        cannot point at two metals however short the second contact looks.

        The far contact is dropped rather than warned about, because anything
        left in 'metals' becomes a bond that the force field then has to hold.
        A note records what went, so the decision can be reversed by hand.
        """

        ligands = [self.contact('HIS', 'NE2', [2.13, 2.98])]
        notes = []

        self.shells().sites._assign_binding_modes(ligands, notes)

        assert ligands[0]['mode'] == 'monodentate'
        assert ligands[0]['metals'] == [0]
        assert ligands[0]['distances'] == [2.13]
        assert len(notes) == 1
        assert 'cannot bridge' in notes[0]

    def test_a_histidine_reaching_two_metals_keeps_one_contact(self):
        """
        The same rule applies to a histidine whose two ring nitrogens sit near
        different metals: that is a bridge as well, so the residue is cut back
        to the one contact it can actually make.

        Keeping both left the residue labelled monodentate while it still
        reported, and bonded, two metals.
        """

        ligands = [
            self.contact('HIS', 'ND1', [2.88]),
            self.contact('HIS', 'NE2', [1.95]),
        ]
        # the second nitrogen reached the other metal
        ligands[1]['metals'] = [1]
        notes = []

        self.shells().sites._assign_binding_modes(ligands, notes)

        assert len(ligands) == 1
        assert ligands[0]['atom'] == 'NE2'
        assert ligands[0]['mode'] == 'monodentate'
        assert ligands[0]['metals'] == [1]
        assert ligands[0]['distances'] == [1.95]
        assert any('cannot bridge' in note for note in notes)

    def test_a_bridging_oxygen_beside_a_second_contact_is_a_double_bridge(self):
        """
        One oxygen holding both metals while the other reaches one of them is
        a single ligand binding twice over. Left alone the residue came out
        as bridging_single and monodentate at the same time, which is not a
        description of anything.
        """

        ligands = [
            self.contact('GLU', 'OE2', [2.09, 2.35]),
            self.contact('GLU', 'OE1', [2.20]),
        ]
        notes = []

        self.shells().sites._assign_binding_modes(ligands, notes)

        assert [ligand['mode'] for ligand in ligands] == ['bridging_double'] * 2
        # both contacts survive: all three are real bonds
        assert ligands[0]['metals'] == [0, 1]
        assert ligands[1]['metals'] == [0]
        assert notes == []

    def test_a_lopsided_carboxylate_is_monodentate(self):
        """
        A carboxylate is only chelating while both oxygens are at comparable
        distances. Past BIDENTATE_ASYMMETRY the far oxygen is pointing away,
        and bonding it anyway would bend the group open.
        """

        ligands = [
            self.contact('ASP', 'OD1', [1.96]),
            self.contact('ASP', 'OD2', [2.90]),
        ]
        notes = []

        self.shells().sites._assign_binding_modes(ligands, notes)

        assert len(ligands) == 1
        assert ligands[0]['atom'] == 'OD1'
        assert ligands[0]['mode'] == 'monodentate'
        assert len(notes) == 1
        assert 'not a chelating bidentate' in notes[0]

    def test_an_even_carboxylate_stays_bidentate(self):

        ligands = [
            self.contact('ASP', 'OD1', [2.50]),
            self.contact('ASP', 'OD2', [2.57]),
        ]
        notes = []

        self.shells().sites._assign_binding_modes(ligands, notes)

        assert [ligand['mode'] for ligand in ligands] == ['bidentate'] * 2
        assert notes == []

    def test_the_bidentate_asymmetry_is_a_fixed_rule(self):
        """
        Two carboxylate oxygens are read as a chelating bidentate only while
        the far one sits within BIDENTATE_ASYMMETRY of the near one; the
        threshold is a rule of the detection, not a setting.
        """

        def classify(far):
            ligands = [
                self.contact('GLU', 'OE1', [2.10]),
                self.contact('GLU', 'OE2', [far]),
            ]
            self.shells().sites._assign_binding_modes(ligands, [])
            return [ligand['mode'] for ligand in ligands]

        inside = 2.10 + util.BIDENTATE_ASYMMETRY - 0.05
        outside = 2.10 + util.BIDENTATE_ASYMMETRY + 0.05
        assert classify(inside) == (['bidentate'] * 2)
        assert classify(outside) == ['monodentate']

    def test_the_asymmetry_rule_leaves_other_residues_alone(self):
        """
        The rule is about the two oxygens of one carboxylate, so a pair of
        contacts on any other residue is not touched by it.
        """

        ligands = [
            self.contact('HIS', 'ND1', [1.96]),
            self.contact('HIS', 'NE2', [2.90]),
        ]
        notes = []

        self.shells().sites._assign_binding_modes(ligands, notes)

        assert len(ligands) == 2
        assert [ligand['mode'] for ligand in ligands] == ['bidentate'] * 2
        assert any('more than one ring nitrogen' in note for note in notes)

    def test_two_carboxylate_oxygens_on_two_metals_bridge(self):
        """
        The counterpart of the histidine case: a carboxylate reaching two
        metals through both of its oxygens is the mu-1,3 bridge, which is
        common enough in binuclear sites that it has its own mode.
        """

        ligands = [
            self.contact('ASP', 'OD1', [2.77]),
            self.contact('ASP', 'OD2', [2.65]),
        ]
        ligands[1]['metals'] = [1]
        notes = []

        self.shells().sites._assign_binding_modes(ligands, notes)

        assert [ligand['mode'] for ligand in ligands] == ['bridging_mu13'] * 2
        assert notes == []

    def test_protonation_variants_and_cluster_charge(self):

        pytest.importorskip('openmm')

        prepared = self.prepared()

        variants = sorted(prepared['binding_modes']['variants'].values())

        assert variants == ['ASP', 'ASP', 'CYX', 'GLU', 'HID', 'HID']

        # two Zn(2+) with CYX, two ASP, GLU and two neutral HID
        assert prepared['active_site']['molecule'].get_charge() == 0
        assert prepared['active_site']['molecule'].number_of_atoms() == 55
        assert len(prepared['active_site']['cap_indices']) == 6

    def test_labels_name_what_the_edit_methods_take(self):
        """
        The labels exist so that Molecule.show tells the user what to pass to
        add_metal_bond and remove_metal_bond, so every one of them has to be
        accepted by those methods rather than merely look right.
        """

        pytest.importorskip('openmm')

        prepared = self.prepared()
        topology = prepared['topology']
        positions = prepared['positions']
        binding_modes = prepared['binding_modes']
        active_site = prepared['active_site']

        labels = active_site['labels']
        elements = active_site['molecule'].get_labels()

        assert len(labels) == len(elements)

        # hydrogens, capping ones included, say nothing
        for index, element in enumerate(elements):
            if element == 'H':
                assert labels[index] == ''

        # a metal is named by the atom index _resolve_metal wants
        for index in active_site['metal_indices']:
            metal = int(labels[index])
            assert self.shells().sites._resolve_metal(binding_modes, metal)['index'] == \
                active_site['atom_map'][index]

        # a beta carbon carries its residue's resid, and the rest of the
        # sidechain carries the atom names the topology uses
        atoms = list(topology.atoms())
        for index, label in enumerate(labels):
            if index in active_site['metal_indices'] or label == '':
                continue
            atom = atoms[active_site['atom_map'][index]]
            if index in active_site['beta_carbon_indices']:
                assert label == str(atom.residue.id)
                assert self.shells().sites._resolve_residue(topology, label, None) is atom.residue
            else:
                assert label == atom.name

        assert sorted(labels[index]
                      for index in active_site['beta_carbon_indices']) == \
            ['130', '134', '187', '191', '73', '77']

        # and the labels of one residue drive an edit end to end
        beta_carbon = active_site['beta_carbon_indices'][2]
        donor = labels.index('OD1', beta_carbon)
        metal = int(labels[active_site['metal_indices'][1]])

        removed = self.unbonded(binding_modes,
                                prepared['topology'],
                                prepared['positions'],
                                labels[beta_carbon],
                                metal=metal,
                                atom=labels[donor])
        restored = self.bonded(removed,
                               topology,
                               positions,
                               labels[beta_carbon],
                               metal,
                               atom=labels[donor])

        def bonds(modes):
            return sorted((ligand['index'], tuple(sorted(ligand['metals'])))
                          for ligand in modes['ligands'])

        assert bonds(removed) != bonds(binding_modes)
        assert bonds(restored) == bonds(binding_modes)

    def test_show_active_site_draws_the_site_own_bonds_and_labels(self):
        """
        Two things about a truncated metal site cannot be recovered from the
        geometry, so showing it has to hand them over rather than let
        Molecule.show work them out: the bonds, since a metal-ligand bond is
        a decision nothing perceives and OpenMM's hydrogens sit past the C-H
        threshold, and the labels, which name an atom the way the edit
        methods want to be told about it.

        Everything else the caller passes goes straight through.
        """

        pytest.importorskip('openmm')

        active_site = self.prepared()['active_site']
        bonds = util.connectivity_bonds(active_site['connectivity_matrix'])

        captured = {}
        original = Molecule.show
        Molecule.show = lambda self, **kwargs: captured.update(kwargs) or 'ok'
        try:
            assert self.shells().sites.show_active_site(active_site,
                                                        width=640,
                                                        atom_indices=True) == 'ok'
        finally:
            Molecule.show = original

        assert captured['atom_labels'] == active_site['labels']
        assert captured['bonds'] == bonds

        # the rest is forwarded untouched
        assert captured['width'] == 640
        assert captured['atom_indices'] is True

        # one label per atom, and they say what an edit takes: the atom index
        # on a metal, the resid on the beta carbon standing for its residue,
        # and nothing on a hydrogen
        labels = active_site['labels']
        elements = active_site['molecule'].get_labels()

        assert len(labels) == active_site['molecule'].number_of_atoms()

        # a metal is labelled with its index in the protonated topology,
        # which is the number add_metal_bond takes and the only one the user
        # is ever shown -- not its index within the cluster
        metal = active_site['metal_indices'][0]
        assert labels[metal] == str(active_site['atom_map'][metal])

        # and a beta carbon with the resid of the residue it stands for
        assert labels[active_site['beta_carbon_indices'][0]].isdigit()
        assert all(labels[i] == '' for i, element in enumerate(elements)
                   if element == 'H')

    def test_show_active_site_lets_the_caller_override(self):
        """
        The bonds and the labels are defaults, not a fixed choice: a caller
        drawing a forming bond or a labelling of its own should not have to
        go around the function.
        """

        pytest.importorskip('openmm')

        active_site = self.prepared()['active_site']

        captured = {}
        original = Molecule.show
        Molecule.show = lambda self, **kwargs: captured.update(kwargs)
        try:
            self.shells().sites.show_active_site(active_site,
                                                 bonds=[(0, 1)],
                                                 atom_labels=False)
        finally:
            Molecule.show = original

        assert captured['bonds'] == [(0, 1)]
        assert captured['atom_labels'] is False

    def test_the_bonds_to_draw_are_derived_from_the_matrix(self):
        """
        The bonds exist so that Molecule.show draws the site's own bonds
        rather than perceiving them by distance, which here would lose every
        metal-ligand bond and invent C-H ones. They are derived from the
        matrix where they are drawn rather than stored beside it, so there is
        no second copy that can fall behind.
        """

        pytest.importorskip('openmm')

        prepared = self.prepared()
        active_site = prepared['active_site']
        matrix = active_site['connectivity_matrix']
        bonds = util.connectivity_bonds(matrix)

        assert 'bond_labels' not in active_site
        assert len(bonds) == int(matrix.sum() // 2)

        # zero-indexed plain ints, in the order and the form show wants
        for first, second in bonds:
            assert isinstance(first, int) and isinstance(second, int)
            assert first < second
            assert matrix[first, second] == 1

        # the metal-ligand bonds are in there, and nothing perceives those:
        # there is one for every contact the binding modes decided on
        metals = set(active_site['metal_indices'])
        metal_bonds = [pair for pair in bonds if set(pair) & metals]
        contacts = sum(
            len(ligand['metals'])
            for ligand in prepared['binding_modes']['ligands'])

        assert len(metal_bonds) == contacts > 0

        # a relaxation that changes the coordination changes both together
        coordinates = self.moved_ligand_geometry(prepared, 'ASP134', 'OD2', 1.6)
        _, updated_site, changed = self.shells().sites.update_binding_modes(
            prepared['topology'], coordinates, active_site,
            prepared['binding_modes'])

        assert changed
        assert len(util.connectivity_bonds(
            updated_site['connectivity_matrix'])) == len(bonds) - 1

        # and the site that was passed in still holds its own
        assert util.connectivity_bonds(matrix) == bonds

    def test_connectivity_comes_from_the_topology(self):
        """
        The atom typing collapses if the hydrogens are perceived by distance,
        because OpenMM places them at about 1.19 Angstrom and the C-H
        threshold of Molecule.get_connectivity_matrix is 1.183 Angstrom.
        """

        pytest.importorskip('openmm')

        prepared = self.prepared()
        connectivity_matrix = prepared['active_site']['connectivity_matrix']
        labels = prepared['active_site']['molecule'].get_labels()
        molecule = prepared['active_site']['molecule']

        for index, label in enumerate(labels):
            if label == 'H':
                assert connectivity_matrix[index].sum() == 1

        # distance perception disagrees, which is the whole point
        perceived = molecule.get_connectivity_matrix()

        assert not np.array_equal(perceived, connectivity_matrix)

        forcefield = self.seeded_forcefield(prepared['active_site'])

        # a sidechain CB must be sp3, not the sp carbon that mis-perception
        # would give
        cb_index = labels.index('C')

        assert forcefield.atoms[cb_index]['type'] == 'c3'

    @staticmethod
    def moved_ligand_geometry(prepared, residue, atom_name, shift):
        """
        Returns the active site coordinates with one ligand atom moved along
        the line to the metal it binds.

        :param prepared:
            The prepared pipeline products.
        :param residue:
            The residue label of the atom to move.
        :param atom_name:
            The name of the atom to move.
        :param shift:
            The distance in Angstrom to move it by, away from the metal.

        :return:
            The coordinates in Angstrom.
        """

        active_site = prepared['active_site']
        coordinates = active_site['molecule'].get_coordinates_in_angstrom()
        site_of = {
            top_index: site_index
            for site_index, top_index in active_site['atom_map'].items()
        }

        ligand = [
            entry for entry in prepared['binding_modes']['ligands']
            if entry['residue'] == residue and entry['atom'] == atom_name
        ][0]

        moved = site_of[ligand['index']]
        metal = site_of[ligand['metals'][0]]
        direction = coordinates[moved] - coordinates[metal]
        direction /= np.linalg.norm(direction)
        coordinates[moved] += direction * shift

        return coordinates

    def test_an_unchanged_geometry_leaves_the_binding_modes_alone(self):

        pytest.importorskip('openmm')

        prepared = self.prepared()

        binding_modes, active_site, updated = (self.shells().sites.update_binding_modes(
            prepared['topology'], prepared['active_site']['molecule'],
            prepared['active_site'], prepared['binding_modes']))

        assert updated is False
        # the arguments themselves come back, so the caller can overwrite
        # unconditionally without losing what it holds
        assert binding_modes is prepared['binding_modes']
        assert active_site is prepared['active_site']

    def test_a_stretched_contact_is_dropped_from_the_binding_modes(self):

        pytest.importorskip('openmm')

        prepared = self.prepared()
        old_modes = prepared['binding_modes']
        old_matrix = prepared['active_site']['connectivity_matrix']

        # 2.57 A plus 1.0 A puts OD2 past the primary cutoff, which leaves
        # ASP130 monodentate through OD1
        coordinates = self.moved_ligand_geometry(prepared, 'ASP130', 'OD2', 1.0)

        binding_modes, active_site, updated = (self.shells().sites.update_binding_modes(
            prepared['topology'], coordinates, prepared['active_site'],
            old_modes))
        connectivity_matrix = active_site['connectivity_matrix']

        assert updated is True

        found = {
            (ligand['residue'], ligand['atom']): ligand
            for ligand in binding_modes['ligands']
        }

        assert ('ASP130', 'OD2') not in found
        assert found[('ASP130', 'OD1')]['mode'] == 'monodentate'
        # every other contact survives untouched
        assert set(found) == set(self.expected_contacts) - {('ASP130', 'OD2')}

        assert connectivity_matrix.sum() == old_matrix.sum() - 2

        # the arguments are not written to
        assert len(old_modes['ligands']) == len(self.expected_contacts)
        assert old_matrix.sum() != connectivity_matrix.sum()

    def test_a_new_contact_is_added_to_the_binding_modes(self):

        pytest.importorskip('openmm')

        prepared = self.prepared()
        active_site = prepared['active_site']

        # GLU73 OE1 sits at 3.17 A, between the two cutoffs; pulling it in by
        # 1.2 A makes it a second contact of the bridging carboxylate
        coordinates = self.moved_ligand_geometry(prepared, 'GLU73', 'OE2', 0.0)
        site_of = {
            top_index: site_index
            for site_index, top_index in active_site['atom_map'].items()
        }
        oe1 = [
            atom for atom in prepared['topology'].atoms()
            if atom.name == 'OE1' and atom.index in site_of
        ][0]
        metal = active_site['metal_indices'][0]
        direction = coordinates[site_of[oe1.index]] - coordinates[metal]
        direction /= np.linalg.norm(direction)
        coordinates[site_of[oe1.index]] = coordinates[metal] + direction * 2.0

        binding_modes, updated_site, updated = (self.shells().sites.update_binding_modes(
            prepared['topology'], coordinates, active_site,
            prepared['binding_modes']))
        connectivity_matrix = updated_site['connectivity_matrix']

        assert updated is True

        found = {
            (ligand['residue'], ligand['atom']): ligand
            for ligand in binding_modes['ligands']
        }

        assert ('GLU73', 'OE1') in found
        assert (connectivity_matrix.sum() ==
                prepared['active_site']['connectivity_matrix'].sum() + 2)
        # the covalent bonds are untouched: only the metal rows differ
        difference = (connectivity_matrix -
                      prepared['active_site']['connectivity_matrix'])
        changed = set(np.argwhere(difference).flatten().tolist())

        assert changed <= set(
            active_site['metal_indices']) | {site_of[oe1.index]}

    def test_the_binding_mode_update_is_printed(self, capsys):

        pytest.importorskip('openmm')

        # unmuted, since the point is what reaches stdout
        prepared = self.prepared()

        coordinates = self.moved_ligand_geometry(prepared, 'ASP130', 'OD2', 1.0)

        capsys.readouterr()
        self.shells(self.printing()).sites.update_binding_modes(prepared['topology'],
                                                                coordinates,
                                                                prepared['active_site'],
                                                                prepared['binding_modes'])
        captured = capsys.readouterr().out

        assert 'Coordination update' in captured
        assert 'lost' in captured
        assert 'ASP130 OD2' in captured
        assert 'bidentate -> monodentate' in captured

    def test_a_geometry_of_the_wrong_size_is_rejected(self):

        pytest.importorskip('openmm')

        prepared = self.prepared()
        coordinates = prepared['active_site']['molecule']
        coordinates = coordinates.get_coordinates_in_angstrom()[:-1]

        with pytest.raises(Exception, match='while the active site has'):
            self.shells().sites.update_binding_modes(prepared['topology'], coordinates,
                                                     prepared['active_site'],
                                                     prepared['binding_modes'])

    def test_extract_pairs_is_smaller_than_a_full_hessian(self):

        pytest.importorskip('openmm')

        prepared = self.prepared()
        connectivity_matrix = prepared['active_site']['connectivity_matrix']
        metals = prepared['active_site']['metal_indices']

        pairs, atoms = util.extract_pairs(connectivity_matrix,
                                          metals,
                                          bond_count=2)

        n_atoms = prepared['active_site']['molecule'].number_of_atoms()

        assert len(atoms) < n_atoms / 2
        assert len(pairs) < n_atoms * (n_atoms + 1) // 20

        # every metal-ligand bond must be present or its force constant is
        # lost
        metal_bonds = {(min(i, j), max(i, j))
                       for i in metals
                       for j in np.where(connectivity_matrix[i])[0]}

        assert metal_bonds <= set(pairs)

    def test_bond_count_one_gives_only_the_metal_bonds(self):

        pytest.importorskip('openmm')

        prepared = self.prepared()
        metals = set(prepared['active_site']['metal_indices'])

        pairs, atoms = util.extract_pairs(
            prepared['active_site']['connectivity_matrix'],
            metals,
            bond_count=1)

        assert all(metals & set(pair) for pair in pairs)

    def test_the_partial_hessian_covers_an_unbonded_donor(self):
        """
        A donor the geometry left just outside metal_bond_cutoff is the one
        add_metal_bond is reached for afterwards, and a block that was never
        computed costs the whole Hessian again to fill in. So the pairs the
        partial Hessian is computed over cover every donor within
        partial_hessian_cutoff, bonded or not.
        """

        pytest.importorskip('openmm')

        prepared = self.prepared()
        active_site = prepared['active_site']
        matrix = np.asarray(active_site['connectivity_matrix'])
        metals = active_site['metal_indices']

        molecule = active_site['molecule']
        coordinates = molecule.get_coordinates_in_angstrom()
        labels = molecule.get_labels()

        strict, _ = util.extract_pairs(matrix, metals, bond_count=2)
        forgiving, _ = self.shells(
            partial_hessian_cutoff=3.5).qm.hessian_pairs(active_site)

        # nothing the metal terms read may be lost by widening the walk
        assert set(strict) <= set(forgiving)

        gained = set(forgiving) - set(strict)
        assert gained

        # every pair gained is a near contact of a metal, or something the
        # walk reached through one
        near = {
            (min(metal, index), max(metal, index))
            for metal in metals
            for index, label in enumerate(labels)
            if label in util.DONOR_ELEMENTS and not matrix[metal, index]
            and np.linalg.norm(coordinates[index] - coordinates[metal]) <= 3.5
        }

        assert near
        assert near <= gained

    def test_the_partial_hessian_cutoff_can_be_switched_off(self):

        pytest.importorskip('openmm')

        prepared = self.prepared()
        active_site = prepared['active_site']

        strict, strict_atoms = util.extract_pairs(
            active_site['connectivity_matrix'],
            active_site['metal_indices'],
            bond_count=2)
        pairs, atoms = self.shells(
            partial_hessian_cutoff=None).qm.hessian_pairs(active_site)

        assert pairs == strict
        assert atoms == strict_atoms

    def test_the_forgiving_walk_leaves_the_active_site_alone(self):
        """
        What is computed is widened, never what is bonded: the site keeps the
        coordination it was perceived with, so the fit keeps fitting that.
        """

        pytest.importorskip('openmm')

        prepared = self.prepared()
        active_site = prepared['active_site']
        before = np.array(active_site['connectivity_matrix'], copy=True)

        self.shells(partial_hessian_cutoff=5.0).qm.hessian_pairs(active_site)

        assert np.array_equal(active_site['connectivity_matrix'], before)

    def test_fitted_equilibria_reproduce_the_input_geometry(self):
        """
        The equilibrium values must come from the geometry the Hessian was
        computed on, and topologically equivalent ligands sitting at different
        distances must stay distinct.
        """

        pytest.importorskip('openmm')

        prepared = self.prepared()
        molecule = prepared['active_site']['molecule']
        coords = molecule.get_coordinates_in_angstrom()

        forcefield = self.fitted_forcefield(prepared['active_site'],
                                            self.random_hessian(molecule.number_of_atoms()))

        bonds, angles = util.get_metal_keys(forcefield, prepared['active_site'])

        assert len(bonds) == 8
        assert len(angles) == 24

        for key in bonds:
            distance = np.linalg.norm(coords[key[0]] - coords[key[1]])
            # the equilibrium values are in nm
            assert forcefield.bonds[key]['equilibrium'] * 10 == pytest.approx(
                distance, abs=0.01)

        # the two histidines are topologically equivalent but 0.1 A apart
        labels = prepared['active_site']['molecule'].get_labels()
        his_bonds = [
            key for key in bonds
            if sorted(labels[index] for index in key) == ['N', 'Zn']
        ]

        assert len(his_bonds) == 2

        equilibria = [forcefield.bonds[key]['equilibrium'] for key in his_bonds]

        assert abs(equilibria[0] - equilibria[1]) > 0.005

    def test_averaging_metal_terms_merges_the_two_histidines(self):

        pytest.importorskip('openmm')

        prepared = self.prepared()
        molecule = prepared['active_site']['molecule']

        forcefield = self.fitted_forcefield(prepared['active_site'],
                                            self.random_hessian(molecule.number_of_atoms()),
                                            average_metal_terms=True)

        labels = prepared['active_site']['molecule'].get_labels()
        bonds, angles = util.get_metal_keys(forcefield, prepared['active_site'])
        his_bonds = [
            key for key in bonds
            if sorted(labels[index] for index in key) == ['N', 'Zn']
        ]

        equilibria = [forcefield.bonds[key]['equilibrium'] for key in his_bonds]

        assert equilibria[0] == pytest.approx(equilibria[1])

    def test_bridging_angle_term_exists(self):

        pytest.importorskip('openmm')

        prepared = self.prepared()
        molecule = prepared['active_site']['molecule']

        forcefield = self.fitted_forcefield(prepared['active_site'],
                                            self.random_hessian(molecule.number_of_atoms()))

        metals = set(prepared['active_site']['metal_indices'])
        bridging = [
            key for key in forcefield.angles
            if key[0] in metals and key[2] in metals
        ]

        assert len(bridging) == 1
        # the main geometric handle on the metal-metal distance
        assert forcefield.angles[bridging[0]]['equilibrium'] == pytest.approx(
            98.0, abs=1.0)

    def test_build_forcefield_rejects_the_xtb_string(self):
        """
        Passing 'xtb' would re-optimize the cluster without constraints.
        """

        pytest.importorskip('openmm')

        prepared = self.prepared()

        with pytest.raises(Exception):
            self.fitted_forcefield(prepared['active_site'], 'xtb')

    def test_build_forcefield_without_a_hessian_measures_the_equilibria(self):
        """
        The seeded pass takes its equilibrium values off the input geometry.

        Before any QM has run, the structure file is the only description of
        the site there is, and a designed binuclear site is strained enough
        that a table of textbook distances would pull the two metals onto
        values neither of them has.
        """

        pytest.importorskip('openmm')

        prepared = self.prepared()
        molecule = prepared['active_site']['molecule']
        forcefield = self.seeded_forcefield(prepared['active_site'])

        bonds, angles = util.get_metal_keys(forcefield, prepared['active_site'])

        assert bonds and angles

        for key in bonds:
            # the getters index atoms from one, the force field uses nanometers
            measured = 0.1 * molecule.get_distance_in_angstroms(
                [index + 1 for index in key])
            assert forcefield.bonds[key]['equilibrium'] == pytest.approx(
                measured)
            assert forcefield.bonds[key]['force_constant'] == (
                util.DEFAULT_METAL_BOND_FORCE_CONSTANT)
            assert 'measured' in forcefield.bonds[key]['comment']

        for key in angles:
            measured = molecule.get_angle_in_degrees(
                [index + 1 for index in key])
            assert forcefield.angles[key]['equilibrium'] == pytest.approx(
                measured)
            assert forcefield.angles[key]['force_constant'] == (
                util.DEFAULT_METAL_ANGLE_FORCE_CONSTANT)

        # and the strain is real: at least one bond is far enough from the
        # literature value that seeding from a table would have moved it
        literature = util.LITERATURE_METAL_BONDS
        labels = prepared['active_site']['molecule'].get_labels()
        deviations = []
        for key in bonds:
            elements = tuple(labels[index] for index in key)
            expected = literature.get(elements, literature.get(elements[::-1]))
            if expected is not None:
                deviations.append(
                    abs(forcefield.bonds[key]['equilibrium'] - expected))
        assert max(deviations) > 0.005

    @staticmethod
    def guessed_covalent_terms(forcefield, active_site):
        """
        The bonds and angles that found no parameters, metals left out.

        A term holding a metal is meant to be guessed -- that is what the
        seeding and the Seminario fit are for. A covalent one is a gap.

        :param forcefield:
            The force field.
        :param active_site:
            The active site, for the metal indices.

        :return:
            The keys, as (kind, atom types) pairs.
        """

        metals = set(active_site['metal_indices'])
        terms = []
        for kind in ('bonds', 'angles'):
            for key, term in getattr(forcefield, kind).items():
                if metals & set(key):
                    continue
                if term.get('comment') == util.UNPARAMETERIZED_COMMENT:
                    types = '-'.join(forcefield.atoms[index]['type'].strip()
                                     for index in key)
                    terms.append((kind, types))

        return terms

    def test_metal_blind_typing_parameterizes_the_residues(self):
        """
        The coordinating residues are typed as the amino acids they are.

        GAFF types an atom from what it is bonded to, so with the metal
        bonds in the perception a carboxylate oxygen is typed as an ether
        one, the carbon behind it loses the carbonyl type, and the covalent
        terms of the residue fall back to a flat guessed constant that the
        fit -- restricted to the metal keys -- never repairs. Looking past
        the metal bonds gives the residues their own parameters back.
        """

        pytest.importorskip('openmm')

        prepared = self.prepared()
        active_site = prepared['active_site']

        bonded = self.seeded_forcefield(active_site, metal_blind_typing=False)
        blind = self.seeded_forcefield(active_site, metal_blind_typing=True)

        # the site is a designed binuclear zinc center with a bridging
        # glutamate, a bidentate aspartate, a cysteine and two histidines,
        # and every one of them loses parameters to the coordination
        assert self.guessed_covalent_terms(bonded, active_site)
        assert self.guessed_covalent_terms(blind, active_site) == []

        # the bridging oxygen reaches two metals, which leaves the
        # perception with no branch at all and drops it to UFF
        assert any(atom['comment'] == util.UFF_TYPE_COMMENT
                   for atom in bonded.atoms.values())
        assert not any(atom['comment'] == util.UFF_TYPE_COMMENT
                       for atom in blind.atoms.values())

    def test_metal_blind_typing_leaves_the_topology_alone(self):
        """
        Only the typing looks past the metal bonds; the bonds are still made.

        The switch is about which parameters a term gets, not about which
        terms there are: a metal bond the coordination asked for still has
        to reach the force field, with the value the seeding gave it.
        """

        pytest.importorskip('openmm')

        prepared = self.prepared()
        active_site = prepared['active_site']
        metals = set(active_site['metal_indices'])

        bonded = self.seeded_forcefield(active_site, metal_blind_typing=False)
        blind = self.seeded_forcefield(active_site, metal_blind_typing=True)

        for kind in ('bonds', 'angles', 'dihedrals'):
            assert set(getattr(bonded, kind)) == set(getattr(blind, kind))

        metal_bonds = [key for key in blind.bonds if metals & set(key)]
        assert metal_bonds
        for key in metal_bonds:
            assert (blind.bonds[key]['force_constant'] == bonded.bonds[key]
                    ['force_constant'])
            assert (blind.bonds[key]['equilibrium'] == bonded.bonds[key]
                    ['equilibrium'])

    def test_the_builder_types_metal_blind_by_default(self, tmp_path):
        """
        A site is a set of amino acids, so the shell turns the switch on,
        and the ActiveSiteBuilder it makes reads it off the shell.
        """

        builder = self.builder(tmp_path)

        assert builder.metal_blind_typing is True
        assert builder._sites().metal_blind_typing is True

        builder.metal_blind_typing = False

        assert builder._sites().metal_blind_typing is False

    def test_metal_equilibria_can_be_overridden(self):
        """
        The override tables replace the measured values per element
        combination, in either key order.
        """

        pytest.importorskip('openmm')

        prepared = self.prepared()
        labels = prepared['active_site']['molecule'].get_labels()

        # deliberately mixed key order: a bond reads the same either way
        equilibria = {
            ('Zn', 'O'): 0.111,
            ('N', 'Zn'): 0.222,
            ('S', 'Zn'): 0.333,
        }

        forcefield = self.seeded_forcefield(
            prepared['active_site'], default_metal_bond_equilibria=equilibria)
        bonds, angles = util.get_metal_keys(forcefield, prepared['active_site'])

        for key in bonds:
            elements = set(labels[index] for index in key)
            expected = {
                frozenset(('Zn', 'O')): 0.111,
                frozenset(('Zn', 'N')): 0.222,
                frozenset(('Zn', 'S')): 0.333,
            }[frozenset(elements)]
            assert forcefield.bonds[key]['equilibrium'] == expected
            assert forcefield.bonds[key]['comment'] == 'given equilibrium'

        # the same for one angle, keyed backwards
        triple = tuple(labels[index] for index in angles[0])

        forcefield = self.seeded_forcefield(
            prepared['active_site'],
            default_metal_angle_equilibria={triple[::-1]: 77.0})
        _, angles = util.get_metal_keys(forcefield, prepared['active_site'])

        matched = 0
        for key in angles:
            elements = tuple(labels[index] for index in key)
            if elements in (triple, triple[::-1]):
                assert forcefield.angles[key]['equilibrium'] == 77.0
                matched += 1
            else:
                assert forcefield.angles[key]['equilibrium'] != 77.0
        assert matched > 0

    def test_the_seeded_force_constants_are_the_flat_defaults(self):
        """
        Nothing before a Hessian says anything about the stiffness of a
        metal term, so every seeded one carries the same flat constant.
        """

        pytest.importorskip('openmm')

        prepared = self.prepared()

        forcefield = self.seeded_forcefield(prepared['active_site'])
        bonds, angles = util.get_metal_keys(forcefield, prepared['active_site'])

        assert all(forcefield.bonds[key]['force_constant'] ==
                   util.DEFAULT_METAL_BOND_FORCE_CONSTANT for key in bonds)
        assert all(forcefield.angles[key]['force_constant'] ==
                   util.DEFAULT_METAL_ANGLE_FORCE_CONSTANT for key in angles)

    def test_metal_angles_can_be_switched_off(self):
        """
        With reparameterize_metal_angles off the metal angles keep whatever
        the generator guessed, in the seeded pass and in the Hessian fit
        alike, while the bonds are still parameterized in both.
        """

        pytest.importorskip('openmm')

        prepared = self.prepared()
        n_atoms = prepared['active_site']['molecule'].number_of_atoms()
        hessian = self.random_hessian(n_atoms)

        seeded = self.seeded_forcefield(prepared['active_site'],
                                        reparameterize_metal_angles=False)
        fitted = self.fitted_forcefield(prepared['active_site'],
                                        hessian,
                                        reparameterize_metal_angles=False)

        bonds, angles = util.get_metal_keys(seeded, prepared['active_site'])
        assert bonds and angles

        # neither pass touched the angles, so both agree with each other
        for key in angles:
            assert (seeded.angles[key]['force_constant'] == fitted.angles[key]
                    ['force_constant'])
            assert (seeded.angles[key]['equilibrium'] == fitted.angles[key]
                    ['equilibrium'])

        # the gate is specific to the angles: the bonds were still fitted
        assert any(seeded.bonds[key]['force_constant'] != fitted.bonds[key]
                   ['force_constant'] for key in bonds)

        # and with the flag on, the seeded pass does move them
        touched = self.seeded_forcefield(prepared['active_site'])

        assert any(touched.angles[key]['force_constant'] != seeded.angles[key]
                   ['force_constant'] for key in angles)

    def test_build_forcefield_applies_the_partial_charges(self):
        """
        Checks that supplied charges reach the force field, cap charge folded
        in.

        They once did not reach it at all. Injecting a custom connectivity
        matrix requires setting topology_update_flag, and create_topology
        responds to that flag by resetting partial_charges to None, after which
        resp=False fills them with zeros. Charges assigned before that call
        were discarded without a word.

        The redistribution happens here and not only on transfer, because the
        force field is written to disk and used to minimize the cluster. An
        exported force field that puts charge on capping hydrogens - atoms that
        exist nowhere it is applied - is wrong wherever it is read.
        """

        pytest.importorskip('openmm')

        prepared = self.prepared()
        active_site = prepared['active_site']
        n_atoms = active_site['molecule'].number_of_atoms()
        charges = np.linspace(-0.5, 0.5, n_atoms)

        forcefield = self.seeded_forcefield(active_site,
                                            partial_charges=charges)

        expected = util.redistribute_cap_charges(active_site, charges)

        assert np.allclose(forcefield.partial_charges, expected)

        # both the array and the per-atom entries, since the OpenMM export
        # reads the per-atom entries
        for index in range(n_atoms):
            assert forcefield.atoms[index]['charge'] == pytest.approx(
                expected[index])

        # the caps carry nothing and the total is unchanged
        for index in active_site['cap_indices']:
            assert forcefield.atoms[index]['charge'] == 0.0
        assert sum(forcefield.atoms[index]['charge']
                   for index in range(n_atoms)) == pytest.approx(charges.sum())

    def test_redistribute_cap_charges_is_idempotent(self):
        """
        build_forcefield and redistribute_charges both redistribute, and either can
        be handed charges that have already been through it. The second pass
        has to be a no-op, or the cap charge would be moved twice.
        """

        pytest.importorskip('openmm')

        prepared = self.prepared()
        active_site = prepared['active_site']
        n_atoms = active_site['molecule'].number_of_atoms()

        once = util.redistribute_cap_charges(active_site,
                                             np.linspace(-0.5, 0.5, n_atoms))
        twice = util.redistribute_cap_charges(active_site, once)

        assert np.allclose(once, twice)

    def test_build_forcefield_rejects_wrong_number_of_charges(self):

        pytest.importorskip('openmm')

        prepared = self.prepared()

        with pytest.raises(Exception):
            self.seeded_forcefield(prepared['active_site'],
                                   partial_charges=np.zeros(3))

    def test_minimize_active_site_freezes_the_beta_carbons(self):

        pytest.importorskip('openmm')

        prepared = self.prepared()
        molecule = prepared['active_site']['molecule']
        coords = molecule.get_coordinates_in_angstrom()

        forcefield = self.seeded_forcefield(prepared['active_site'])
        minimized = self.shells().sites.minimize_active_site(
            prepared['active_site'], forcefield,
            util.constrained_indices(prepared['active_site'], False))

        displacement = np.linalg.norm(minimized - coords, axis=1)
        frozen = util.constrained_indices(prepared['active_site'], False)
        free = [
            index for index in range(len(displacement)) if index not in frozen
        ]

        # the coordinates make a round trip through a PDB file, so they carry
        # a rounding of 0.001 Angstrom that is not motion
        assert displacement[frozen].max() < 2.0e-3
        # and the minimization has to have done something
        assert max(displacement[index] for index in free) > 1.0e-2

    def test_an_untemplated_improper_still_reaches_the_relaxation(
            self, monkeypatch):
        """
        The relaxation runs on a system OpenMM builds from the generator's
        own XML, and a residue template can express an improper only as a
        central atom and three atoms bonded to it. The one a bidentate
        carboxylate gets names the metal, which its central carbon is not
        bonded to, so that path drops it without a word. It is put on the
        system by index instead, where nothing has to be bonded to anything.
        """

        pytest.importorskip('openmm')

        import tempfile

        import openmm as mm
        import openmm.app as mmapp

        prepared = self.prepared()
        active_site = prepared['active_site']
        forcefield = self.seeded_forcefield(active_site)

        bonded = {}
        for first, second in forcefield.bonds:
            bonded.setdefault(first, set()).add(second)
            bonded.setdefault(second, set()).add(first)

        untemplated = [
            key for key in forcefield.impropers
            if not set(key[1:]) <= bonded.get(key[0], set())
        ]
        assert untemplated, (
            'the fixture holds a bidentate carboxylate, whose improper names '
            'a metal its central carbon is not bonded to')

        def torsions_of(system):
            force = next(one for one in system.getForces()
                         if isinstance(one, mm.PeriodicTorsionForce))
            return [
                tuple(force.getTorsionParameters(index)[:4])
                for index in range(force.getNumTorsions())
            ]

        # what the residue template alone makes of the force field
        with tempfile.TemporaryDirectory() as temp_dir:
            stem = str(Path(temp_dir) / 'active_site')
            forcefield.write_openmm_files(stem)
            pdb = mmapp.PDBFile(f'{stem}.pdb')
            templated = mmapp.ForceField(f'{stem}.xml').createSystem(
                pdb.topology, nonbondedMethod=mmapp.NoCutoff)

        covered = {frozenset(key) for key in torsions_of(templated)}
        for key in untemplated:
            assert frozenset(key) not in covered

        # and what the relaxation runs on
        seen = {}
        original = ActiveSiteBuilder._add_untemplated_impropers

        def recording(sites, system, ff):
            seen['system'] = system
            seen['added'] = original(sites, system, ff)
            return seen['added']

        monkeypatch.setattr(ActiveSiteBuilder, '_add_untemplated_impropers',
                            recording)
        self.shells().sites.minimize_active_site(
            active_site, forcefield, util.constrained_indices(active_site, False))

        assert seen['added'] == untemplated

        built = torsions_of(seen['system'])
        for key in untemplated:
            # the central atom goes third, where OpenMM puts it itself
            assert (key[1], key[2], key[0], key[3]) in built

        # the ones the template does express are not added a second time
        for key in forcefield.impropers:
            if key in untemplated:
                continue
            assert sum(frozenset(one) == frozenset(key) for one in built) == 1

    def test_constrained_indices_default_to_the_beta_carbons(self):

        pytest.importorskip('openmm')

        prepared = self.prepared()

        beta_carbons = prepared['active_site']['beta_carbon_indices']
        caps = prepared['active_site']['cap_indices']

        # one of each per truncated sidechain
        assert len(beta_carbons) == 6
        assert len(caps) == 6
        assert not set(beta_carbons) & set(caps)

        labels = prepared['active_site']['molecule'].get_labels()
        assert all(labels[index] == 'C' for index in beta_carbons)
        assert all(labels[index] == 'H' for index in caps)

        assert util.constrained_indices(
            prepared['active_site'], False) == sorted(beta_carbons)

        assert util.constrained_indices(
            prepared['active_site'], True) == sorted(beta_carbons + caps)

    def test_minimize_can_also_freeze_the_capping_hydrogens(self):

        pytest.importorskip('openmm')

        prepared = self.prepared()
        coords = prepared['active_site'][
            'molecule'].get_coordinates_in_angstrom()

        forcefield = self.seeded_forcefield(prepared['active_site'])
        frozen = util.constrained_indices(prepared['active_site'], True)
        minimized = self.shells().sites.minimize_active_site(
            prepared['active_site'], forcefield, frozen)

        displacement = np.linalg.norm(minimized - coords, axis=1)

        assert len(frozen) == 12
        assert displacement[frozen].max() < 2.0e-3

    def test_mm_optimize_relaxes_the_site_and_freezes_the_beta_carbons(self):
        """
        The crude pass has to move the free atoms, hold the beta carbons, and
        hand back a new molecule rather than edit the active site in place.
        """

        pytest.importorskip('openmm')

        prepared = self.prepared()
        active_site = prepared['active_site']
        molecule = active_site['molecule']
        before = molecule.get_coordinates_in_angstrom()

        relaxed = self.shells().sites.mm_optimize_active_site(
            active_site, self.seeded_forcefield(active_site))

        after = relaxed.get_coordinates_in_angstrom()
        displacement = np.linalg.norm(after - before, axis=1)
        frozen = util.constrained_indices(active_site, False)
        free = [
            index for index in range(len(displacement)) if index not in frozen
        ]

        # the coordinates make a round trip through a PDB file, so they carry
        # a rounding of 0.001 Angstrom that is not motion
        assert displacement[frozen].max() < 2.0e-3
        assert max(displacement[index] for index in free) > 1.0e-2

        assert relaxed.get_labels() == molecule.get_labels()
        assert relaxed.get_charge() == molecule.get_charge()
        assert relaxed.get_multiplicity() == molecule.get_multiplicity()

        # the active site is an argument, not a place to write to
        assert active_site['molecule'] is molecule
        assert np.allclose(molecule.get_coordinates_in_angstrom(), before)

    def test_mm_optimize_can_freeze_the_metals(self, tmp_path):
        """
        Holding the metals keeps the frame of the site exactly as it was and
        lets the ligands relax around it. Off by default, since the drift of
        a free metal is a reading on the seeded terms.
        """

        pytest.importorskip('openmm')

        prepared = self.prepared()
        active_site = prepared['active_site']
        metals = sorted(active_site['metal_indices'])
        before = active_site['molecule'].get_coordinates_in_angstrom()

        forcefield = self.seeded_forcefield(active_site)

        free = self.shells().sites.mm_optimize_active_site(active_site, forcefield)
        held = self.shells(mm_constrain_metals=True).sites.mm_optimize_active_site(
            active_site, forcefield)

        # the coordinates make a round trip through a PDB file, so they carry
        # a rounding of 0.001 Angstrom that is not motion
        held_shift = np.linalg.norm(held.get_coordinates_in_angstrom()[metals] -
                                    before[metals],
                                    axis=1)
        assert held_shift.max() < 2.0e-3

        free_shift = np.linalg.norm(free.get_coordinates_in_angstrom()[metals] -
                                    before[metals],
                                    axis=1)
        assert free_shift.max() > 2.0e-3

        # the metal-metal distance is what freezing them really preserves
        pair = [metals[0] + 1, metals[1] + 1]
        assert held.get_distance_in_angstroms(pair) == pytest.approx(
            active_site['molecule'].get_distance_in_angstroms(pair), abs=2.0e-3)

        # and the builder's setting reaches the same place. The site is put
        # there by hand rather than built, so the stage has to be said as
        # well: what an operation is allowed to do is a function of the stage
        # alone, and no longer of which field happens to be filled in.
        builder = self.builder(tmp_path, mm_constrain_metals=True)
        builder._active_site = active_site
        builder._stage = Stage.ACTIVE_SITE
        by_setting = builder.mm_optimize_active_site()
        setting_shift = np.linalg.norm(
            by_setting.get_coordinates_in_angstrom()[metals] - before[metals],
            axis=1)
        assert setting_shift.max() < 2.0e-3

    def test_mm_optimize_honors_the_equilibrium_overrides(self):
        """
        The overrides have to reach the geometry, not just the parameters: a
        metal bond forced to an absurd length has to come out of the
        relaxation stretched towards it.
        """

        pytest.importorskip('openmm')

        prepared = self.prepared()
        active_site = prepared['active_site']
        molecule = active_site['molecule']
        labels = molecule.get_labels()
        metals = set(active_site['metal_indices'])

        forcefield = self.seeded_forcefield(active_site)
        bonds, _ = util.get_metal_keys(forcefield, active_site)
        # a sulfur bond: one ligand only, so nothing else fights the stretch
        key = next(key for key in bonds
                   if 'S' in [labels[index] for index in key])
        one_based = [index + 1 for index in key]
        before = molecule.get_distance_in_angstroms(one_based)

        elements = tuple(labels[index] for index in key)

        stretched = self.seeded_forcefield(
            active_site,
            default_metal_bond_equilibria={elements: 0.1 * (before + 0.5)},
            # the angles hold the ligand in place, so let the bond win
            reparameterize_metal_angles=False)
        relaxed = self.shells().sites.mm_optimize_active_site(active_site, stretched)

        after = relaxed.get_distance_in_angstroms(one_based)

        assert after > before + 0.1
        assert metals & set(key)

    def test_build_active_site_runs_the_crude_pass_without_any_qm(
            self, tmp_path):
        """
        build_active_site relaxes the site on the crude force field before any
        QM is asked for, and writes it where it can be looked at.
        """

        pytest.importorskip('openmm')

        builder = self.builder(tmp_path)
        molecule = builder.build_active_site(self.structure_file())

        assert builder.scf_drv is None

        written = tmp_path / 'run' / util.MM_GEOMETRY_FILE
        assert written.is_file()

        relaxed = Molecule.read_xyz_file(str(written))
        assert relaxed.number_of_atoms() == molecule.number_of_atoms()
        assert np.allclose(relaxed.get_coordinates_in_angstrom(),
                           molecule.get_coordinates_in_angstrom(),
                           atol=1.0e-3)

    def test_the_crude_pass_can_be_switched_off(self, tmp_path):

        pytest.importorskip('openmm')

        # prepare_protein off, so that this runs the same steps prepared()
        # does and the two geometries are comparable
        builder = self.builder(tmp_path, prepare_protein=False)

        # the hydrogens are placed off the global random stream, so the two
        # extractions only agree when they start from the same state
        random.seed(20250821)
        extracted = self.prepared()['active_site']['molecule']
        random.seed(20250821)
        molecule = builder.build_active_site(self.structure_file(),
                                             mm_opt=False)

        assert not (tmp_path / 'run' / util.MM_GEOMETRY_FILE).exists()
        assert np.allclose(molecule.get_coordinates_in_angstrom(),
                           extracted.get_coordinates_in_angstrom())

    def test_a_geometry_given_to_build_forcefield_replaces_the_site(
            self, tmp_path):
        """
        A geometry handed to build_forcefield is the geometry the fit
        describes, so it has to reach the active site rather than sit beside
        it.
        """

        pytest.importorskip('openmm')

        molecule = self.prepared()['active_site']['molecule']
        n_atoms = molecule.number_of_atoms()

        builder = self.builder(tmp_path)
        builder.build_active_site(self.structure_file())
        builder.build_forcefield(opt_geometry=molecule,
                                 hessian=self.random_hessian(n_atoms),
                                 partial_charges=np.zeros(n_atoms))

        assert np.allclose(
            builder.active_site_molecule.get_coordinates_in_angstrom(),
            molecule.get_coordinates_in_angstrom())

    def test_the_reporting_distance_follows_the_bonding_one(self, tmp_path):
        """
        A raised bonding cutoff has to reach further, and used not to.

        The scan collects candidates out to the reporting distance and
        decides bonding inside the bonding one, so a bonding cutoff raised
        past a fixed reporting distance bonded nothing new: the contacts it
        would have made were never collected. The reporting distance is
        therefore a margin above the bonding one rather than a setting of its
        own.
        """

        pytest.importorskip('openmm')

        builder = self.builder(tmp_path)

        assert builder._sites().report_cutoff == pytest.approx(
            builder.metal_bond_cutoff + util.REPORT_CUTOFF_MARGIN)

        builder.metal_bond_cutoff = 5.0

        assert builder._sites().report_cutoff > 5.0

        # and it reaches: the fixture holds contacts between the default
        # cutoff and this one that were beyond the old fixed 3.5 A
        near = self.detected(metal_bond_cutoff=3.0)
        far = self.detected(metal_bond_cutoff=5.0)

        def bonded(found):
            return {(ligand['residue'], ligand['atom'])
                    for ligand in found['binding_modes']['ligands']
                    if ligand['metals']}

        assert bonded(far) > bonded(near)

    def test_the_generator_is_not_heard_from_unless_it_is_asked_for(
            self, capsys):
        """
        The force field generator names every parameter it looks up, which
        buries what this module has to say about the site -- and a
        shoehorning, which rebuilds after every edit, prints it all again
        each time. It is silent by default and can be turned back on.
        """

        pytest.importorskip('openmm')

        prepared = self.prepared()

        self.seeded_forcefield(prepared['active_site'],
                               ostream=self.printing())
        quiet = capsys.readouterr().out

        self.seeded_forcefield(prepared['active_site'],
                               ostream=self.printing(),
                               mute_forcefield_generator=False)
        loud = capsys.readouterr().out

        assert 'GAFF' not in quiet
        assert 'GAFF' in loud
        assert len(loud) > len(quiet)

    def test_the_cutoffs_reach_the_re_detection(self, tmp_path, monkeypatch):
        """
        The detection settings have to reach every detection, not only the
        first one.

        Adopting a new geometry re-detects the coordination on it, and that
        re-detection is what the Hessian and the fit read their connectivity
        from. It used to be given the asymmetry alone, so a builder with a
        tightened metal_bond_cutoff honored it while extracting the site and
        silently fell back to the core default of 3.0 the moment a geometry
        was adopted.
        """

        pytest.importorskip('openmm')

        seen = {}
        original = ActiveSiteBuilder.update_binding_modes

        def recording(self, topology, geometry, active_site, binding_modes):
            seen.update(metal_bond_cutoff=self.metal_bond_cutoff,
                        report_cutoff=self.report_cutoff)
            return original(self, topology, geometry, active_site,
                            binding_modes)

        monkeypatch.setattr(ActiveSiteBuilder, 'update_binding_modes',
                            recording)

        # the reporting distance follows the bonding one rather than being
        # set beside it
        builder = self.builder(tmp_path,
                               metal_bond_cutoff=2.6,
                               prepare_protein=False)
        builder.build_active_site(self.structure_file())
        builder._adopt_geometry(builder.active_site_molecule)

        assert seen['metal_bond_cutoff'] == 2.6
        assert seen['report_cutoff'] == 3.1

    def test_given_input_skips_every_qm_step(self, tmp_path):
        """
        With a geometry, a Hessian and charges all given, the pipeline has to
        run through without touching a QM driver.
        """

        pytest.importorskip('openmm')

        molecule = self.prepared()['active_site']['molecule']
        n_atoms = molecule.number_of_atoms()

        builder = self.builder(tmp_path)
        builder.build_active_site(self.structure_file())
        forcefield = builder.build_forcefield(
            opt_geometry=molecule,
            hessian=self.random_hessian(n_atoms),
            partial_charges=np.zeros(n_atoms))

        assert builder.hessian.shape == (3 * n_atoms, 3 * n_atoms)
        assert np.allclose(builder.partial_charges, 0.0)
        assert builder.scf_drv is None
        assert forcefield is builder.active_site_forcefield

        bonds, angles = util.get_metal_keys(forcefield, builder._active_site)
        assert len(bonds) == 8

    def test_the_hessian_is_restricted_unless_the_whole_one_is_asked_for(
            self, tmp_path, monkeypatch):
        """
        Seminario reads one block of the Hessian per metal bond and two per
        metal angle, so everything else the analytical Hessian would hold is
        thrown away by the fit. That is what calculate_partial_hessian is for,
        and turning it off has to reach the driver as no restriction at all
        rather than as a wider set of pairs.
        """

        pytest.importorskip('openmm')

        molecule = self.prepared()['active_site']['molecule']
        n_atoms = molecule.number_of_atoms()

        def run(partial, monkeypatch):
            asked = []

            builder = self.builder(tmp_path, name=f'run_{partial}')
            builder.calculate_partial_hessian = partial

            def hessian(qm_self, active_site, atom_pairs, *qm_settings):
                asked.append(atom_pairs)
                return self.random_hessian(
                    active_site['molecule'].number_of_atoms())

            monkeypatch.setattr(QmParameterizer, 'compute_hessian', hessian)
            builder.build_active_site(self.structure_file())
            builder.build_forcefield(opt_geometry=molecule,
                                     partial_charges=np.zeros(n_atoms))

            assert len(asked) == 1

            return asked[0]

        pairs = run(True, monkeypatch)

        assert pairs is not None
        assert all(len(pair) == 2 for pair in pairs)

        # the full Hessian is the one the driver computes when it is handed
        # nothing to restrict itself to
        assert run(False, monkeypatch) is None

    def test_supplied_geometry_must_match_the_cluster(self):

        pytest.importorskip('openmm')

        prepared = self.prepared()
        # a molecule with the right size but the wrong elements
        wrong = Molecule(
            ['H'] * prepared['active_site']['molecule'].number_of_atoms(),
            prepared['active_site']['molecule'].get_coordinates_in_angstrom())

        with pytest.raises(Exception, match='elements'):
            self.shells().qm._resolve_optimized_geometry(
                prepared['active_site']['molecule'], optimized_geometry=wrong, folder=None)

    def test_supplied_hessian_must_match_the_cluster(self):

        pytest.importorskip('openmm')

        prepared = self.prepared()

        with pytest.raises(Exception, match='shape'):
            self.shells().qm._resolve_hessian(prepared['active_site'],
                                              hessian=np.zeros((6, 6)), folder=None)

    def test_supplied_charges_must_match_the_cluster(self):

        pytest.importorskip('openmm')

        prepared = self.prepared()

        with pytest.raises(Exception, match='shape'):
            self.shells().qm._resolve_partial_charges(prepared['active_site'],
                                                      partial_charges=np.zeros(3), folder=None)

    def test_zero_force_constants_are_reported(self, capsys):
        """
        A Hessian that is not evaluated at a stationary point clamps metal
        force constants to zero, which is the failure most likely to go
        unnoticed. The warning has to survive the muting of the force field
        generator, which shares the output stream of the builder.
        """

        pytest.importorskip('openmm')

        # unmuted, since the point is what reaches stdout
        prepared = self.prepared()
        active_site = prepared['active_site']

        n_atoms = active_site['molecule'].number_of_atoms()

        # every block holds something, and every projection comes out
        # negative, which is what a geometry that is not stationary looks like
        capsys.readouterr()
        forcefield = self.fitted_forcefield(active_site,
                                            np.ones((3 * n_atoms, 3 * n_atoms)),
                                            ostream=self.printing())
        captured = capsys.readouterr().out

        bonds, angles = util.get_metal_keys(forcefield, active_site)

        assert all(forcefield.bonds[key]['force_constant'] == 0.0
                   for key in bonds)
        assert 'zero force constant' in captured
        assert 'stationary point' in captured
        assert 'not covered by the Hessian' not in captured

    def test_a_term_the_hessian_does_not_cover_is_named_as_such(self, capsys):
        """
        compute_hessian is restricted to the pairs extract_pairs walks out of
        the connectivity, so a Hessian reused from a folder whose
        coordination has moved on since holds nothing at all for the bond
        that closed in the meantime. That reads as a zero force constant
        exactly like a geometry that is not stationary does, and the fix is
        not the same one, so the two have to be told apart.
        """

        pytest.importorskip('openmm')

        prepared = self.prepared()
        active_site = prepared['active_site']

        n_atoms = active_site['molecule'].number_of_atoms()
        matrix = np.asarray(active_site['connectivity_matrix'])
        metal = int(active_site['metal_indices'][0])
        ligand = int(np.where(matrix[metal])[0][0])

        hessian = np.ones((3 * n_atoms, 3 * n_atoms))
        # blank one metal bond, exactly as a Hessian restricted to another
        # set of atom pairs would have left it
        for first, second in ((metal, ligand), (ligand, metal)):
            hessian[3 * first:3 * first + 3, 3 * second:3 * second + 3] = 0.0

        capsys.readouterr()
        forcefield = self.fitted_forcefield(active_site,
                                            hessian,
                                            ostream=self.printing())
        captured = capsys.readouterr().out

        bonds, _ = util.get_metal_keys(forcefield, active_site)
        blanked = [key for key in bonds if set(key) == {metal, ligand}]

        assert len(blanked) == 1
        assert not self.shells().qm._hessian_covers(hessian, blanked[0])
        assert all(
            self.shells().qm._hessian_covers(hessian, key) for key in bonds
            if key != blanked[0])

        assert 'not covered by the Hessian' in captured
        assert f'zero force constant: bond {blanked[0]}' in captured
        # the rest are zero for the other reason, and still say so
        assert 'stationary point' in captured

    @staticmethod
    def bridging_bonds(forcefield, active_site):
        """
        The metal bonds of every atom that reaches more than one metal, with
        the length of each, measured on the geometry the fit was made on.

        :return:
            A dictionary of {bridging atom: {bond key: length in Angstrom}}.
        """

        metals = set(active_site['metal_indices'])
        coordinates = active_site['molecule'].get_coordinates_in_angstrom()
        bonds, _ = util.get_metal_keys(forcefield, active_site)

        reached = {}
        for key in bonds:
            ligands = [index for index in key if index not in metals]
            if len(ligands) != 1:
                continue
            reached.setdefault(ligands[0], {})[key] = float(
                np.linalg.norm(coordinates[key[0]] - coordinates[key[1]]))

        return {
            ligand: lengths
            for ligand, lengths in reached.items() if len(lengths) > 1
        }

    def test_the_unfitted_long_half_of_a_bridge_is_dropped(self, capsys):
        """
        An atom bridging two metals reaches both through the same lone pair
        and can hold one of them far more weakly. Where the fit says so by
        giving that half no force constant and the geometry agrees by holding
        it further out, keeping it as a bond with no stiffness is worse than
        not having it: the pair sits at an equilibrium nothing restores while
        its angles, torsions and exclusions all still act. So it goes, and
        everything that crosses it goes with it.
        """

        pytest.importorskip('openmm')

        prepared = self.prepared()
        active_site = prepared['active_site']

        n_atoms = active_site['molecule'].number_of_atoms()
        # every block holds something and every projection comes out
        # negative, so the whole metal site is fitted to zero
        hessian = np.ones((3 * n_atoms, 3 * n_atoms))

        kept = self.fitted_forcefield(active_site, hessian)

        bridges = self.bridging_bonds(kept, active_site)

        assert len(bridges) == 1

        lengths = next(iter(bridges.values()))
        long_key = max(lengths, key=lengths.get)
        short_key = min(lengths, key=lengths.get)

        # the fixture separates the two halves by 0.7 A
        assert lengths[long_key] - lengths[short_key] > util.WEAK_BRIDGE_TOLERANCE
        assert long_key in kept.bonds

        capsys.readouterr()
        pruned = self.fitted_forcefield(active_site,
                                        hessian,
                                        prune_weak_bridge_bonds=True,
                                        ostream=self.printing())
        captured = capsys.readouterr().out

        assert long_key not in pruned.bonds
        assert short_key in pruned.bonds

        pair = set(long_key)

        assert not self.shells().qm._terms_crossing(pruned.angles, pair)
        assert not self.shells().qm._terms_crossing(pruned.dihedrals, pair)
        assert not self.shells().qm._terms_crossing(pruned.impropers, pair, path=False)

        # what was removed is not then reported as a zero force constant
        assert f'Removed the bridging bond {long_key}' in captured
        assert f'zero force constant: bond {long_key}' not in captured

    def test_a_pruned_generator_describes_one_topology(self, capsys):
        """
        The generator carries the connectivity it was built on and the 1-4
        pairs derived from it beside the terms themselves, so dropping a bond
        out of the tables alone leaves it describing two topologies at once:
        the next create_topology would put the bond straight back, and the
        pair list would go on scaling a 1-4 interaction between atoms that
        are no longer 1-4 through anything.
        """

        pytest.importorskip('openmm')

        prepared = self.prepared()
        active_site = prepared['active_site']

        n_atoms = active_site['molecule'].number_of_atoms()
        hessian = np.ones((3 * n_atoms, 3 * n_atoms))

        kept = self.fitted_forcefield(active_site, hessian)
        long_key = max(next(
            iter(self.bridging_bonds(kept, active_site).values())).items(),
                       key=lambda item: item[1])[0]

        capsys.readouterr()
        pruned = self.fitted_forcefield(active_site,
                                        hessian,
                                        prune_weak_bridge_bonds=True)

        assert long_key not in pruned.bonds
        assert not pruned.connectivity_matrix[long_key[0], long_key[1]]
        assert not pruned.connectivity_matrix[long_key[1], long_key[0]]

        # the matrix the generator ends up with is the one that would have
        # produced these bonds had the pruning never been needed
        assert (sorted(pruned.bonds) == util.connectivity_bonds(
            pruned.connectivity_matrix))

        # A pair reached by a second dihedral path is still 1-4, so the pair
        # list is re-derived rather than subtracted from. The invariant that
        # says so is that a pruned generator cannot be told apart from one
        # built on a site that never had the bond, which is what
        # connectivity_from_forcefield makes it possible to ask.
        without = self.shells().sites.connectivity_from_forcefield(active_site, pruned)
        fresh = self.fitted_forcefield(without, hessian)

        assert set(pruned.pairs) == set(fresh.pairs)
        assert set(pruned.pairs) < set(kept.pairs)

    @staticmethod
    def mu_1_3_forcefield():
        """
        A carboxylate bridging two metals through a different oxygen each,
        which is the mu-1,3 mode the fixture does not hold.

        Zn0 is held by O2 at 2.40 A with no force constant, Zn1 by O3 at
        1.95 A with one: the same weak arm as a mu-1,1 bridge, except that no
        single atom reaches both metals.

        :return:
            The tuple of the molecule, the force field and the metal indices.
        """

        from veloxchem.mmforcefieldgenerator import MMForceFieldGenerator

        # Zn0  O2 - C4 - O3  Zn1, with C4 carrying a methyl at C5
        labels = ['Zn', 'Zn', 'O', 'O', 'C', 'C']
        coordinates = np.array([
            [-2.40, 0.00, 0.00],
            [3.16, 0.00, 0.00],
            [0.00, 0.00, 0.00],
            [1.21, 0.00, 0.00],
            [0.61, 1.10, 0.00],
            [0.61, 2.60, 0.00],
        ])

        forcefield = MMForceFieldGenerator()
        forcefield.atoms = {
            index: {
                'name': f'{label}{index}',
                'mass': 65.4 if label == 'Zn' else 12.0,
            }
            for index, label in enumerate(labels)
        }
        forcefield.bonds = {
            (0, 2): {
                'force_constant': 0.0
            },
            (1, 3): {
                'force_constant': 5.0e4
            },
            (2, 4): {
                'force_constant': 3.0e5
            },
            (3, 4): {
                'force_constant': 3.0e5
            },
            (4, 5): {
                'force_constant': 2.0e5
            },
        }
        forcefield.angles = {
            (0, 2, 4): {
                'force_constant': 0.0
            },
            (1, 3, 4): {
                'force_constant': 4.0e2
            },
            (2, 4, 3): {
                'force_constant': 6.0e2
            },
            (2, 4, 5): {
                'force_constant': 6.0e2
            },
        }
        forcefield.dihedrals = {
            (0, 2, 4, 3): {
                'barrier': 1.0
            },
            (0, 2, 4, 5): {
                'barrier': 1.0
            },
            (1, 3, 4, 5): {
                'barrier': 1.0
            }
        }
        forcefield.impropers = {(4, 2, 3, 5): {'barrier': 1.0}}

        # a real generator carries the connectivity it was built on beside
        # the terms, and the pruning keeps the two agreeing, so a hand-made
        # one has to bring it too
        forcefield.connectivity_matrix = np.zeros((len(labels), len(labels)),
                                                  dtype=int)
        for first, second in forcefield.bonds:
            forcefield.connectivity_matrix[first, second] = 1
            forcefield.connectivity_matrix[second, first] = 1

        molecule = Molecule(labels, coordinates, 'angstrom')

        return molecule, forcefield, [0, 1]

    def test_a_mu_1_3_bridge_is_pruned_like_a_single_atom_one(self):
        """
        A carboxylate can bridge two metals through two oxygens rather than
        one, and the weak arm of that is the same situation: what makes a
        bridge is the residue reaching both metals, not one atom doing it.
        The comparison is then against the shortest metal bond of that
        residue, whichever of its atoms carries it.
        """

        molecule, forcefield, metals = self.mu_1_3_forcefield()

        builder = MetalSiteForceFieldBuilder()
        builder.ostream.mute()
        active_site = {'molecule': molecule, 'metal_indices': metals}

        bonds, angles = util.get_metal_keys(forcefield, active_site)

        assert set(bonds) == {(0, 2), (1, 3)}

        # the two oxygens are one residue, so the metals leave them connected
        fragments = self.shells().qm._ligand_fragments(forcefield, set(metals))

        assert fragments[2] == fragments[3]
        assert 0 not in fragments and 1 not in fragments

        bonds, angles = self.shells().qm._prune_weak_bridges(forcefield, active_site, bonds,
                                                             angles, None)

        assert bonds == [(1, 3)]
        assert (0, 2) not in forcefield.bonds
        assert (0, 2, 4) not in forcefield.angles
        assert (0, 2, 4, 3) not in forcefield.dihedrals
        assert (0, 2, 4, 5) not in forcefield.dihedrals

        # nothing that does not cross the dropped bond is touched
        assert (1, 3, 4) in forcefield.angles
        assert (1, 3, 4, 5) in forcefield.dihedrals
        assert (4, 2, 3, 5) in forcefield.impropers

    def test_a_bridge_on_one_metal_is_not_a_bridge(self):
        """
        A bidentate carboxylate holds one metal with two oxygens, which is a
        grip rather than a bridge: there is no second metal for the residue
        to be weakly attached to, so nothing is dropped however the two arms
        compare.
        """

        molecule, forcefield, metals = self.mu_1_3_forcefield()

        # both arms onto Zn0, which is what makes it bidentate
        forcefield.bonds[(0, 3)] = forcefield.bonds.pop((1, 3))
        forcefield.angles[(0, 3, 4)] = forcefield.angles.pop((1, 3, 4))
        forcefield.dihedrals[(0, 3, 4, 5)] = forcefield.dihedrals.pop(
            (1, 3, 4, 5))

        builder = MetalSiteForceFieldBuilder()
        builder.ostream.mute()
        active_site = {'molecule': molecule, 'metal_indices': metals}

        bonds, angles = util.get_metal_keys(forcefield, active_site)
        pruned, _ = self.shells().qm._prune_weak_bridges(forcefield, active_site, bonds,
                                                         angles, None)

        assert set(pruned) == {(0, 2), (0, 3)}
        assert (0, 2) in forcefield.bonds

    def test_a_bridge_the_geometry_does_not_separate_is_kept(
            self, monkeypatch):
        """
        A zero force constant on its own says nothing about a bridge: it can
        equally mean a geometry that is not stationary or a Hessian that
        never covered the pair. The geometry has to agree, so a tolerance the
        two halves do not clear keeps both of them.
        """

        pytest.importorskip('openmm')

        prepared = self.prepared()
        active_site = prepared['active_site']

        n_atoms = active_site['molecule'].number_of_atoms()
        hessian = np.ones((3 * n_atoms, 3 * n_atoms))

        # the tolerance is a constant the prune reads where it is defined
        monkeypatch.setattr(qm, 'WEAK_BRIDGE_TOLERANCE', 1.0)
        forcefield = self.fitted_forcefield(active_site,
                                            hessian,
                                            prune_weak_bridge_bonds=True)

        bridges = self.bridging_bonds(forcefield, active_site)
        lengths = next(iter(bridges.values()))

        assert len(lengths) == 2
        assert max(lengths.values()) - min(lengths.values()) < 1.0
        assert all(forcefield.bonds[key]['force_constant'] == 0.0
                   for key in lengths)

    def test_constants_live_in_the_core(self):

        assert util.SUPPORTED_METAL_ELEMENTS == ('Zn', )
        assert 'Zn' in util.METAL_ELEMENTS
        assert util.VARIANT_CHARGES['CYX'] == -1
        assert ('Zn', 'S') in util.LITERATURE_METAL_BONDS
        assert 'CB' not in util.BACKBONE_ATOM_NAMES

    def test_derive_binding_modes_rejects_a_non_zinc_metal(self, tmp_path):
        """
        The literature distances, formal charges and coordination rules have
        only been validated for zinc, so another metal has to be refused
        rather than silently treated as if it were zinc.
        """

        pytest.importorskip('openmm')

        lines = Path(self.structure_file()).read_text().splitlines(True)
        iron = [
            line.replace('ZN', 'FE').replace('Zn', 'Fe')
            if 'ZN' in line else line for line in lines
        ]
        iron_file = tmp_path / 'fe_site.pdb'
        iron_file.write_text(''.join(iron))

        builder = MetalSiteForceFieldBuilder()
        builder.ostream.mute()
        topology, positions = self.shells(
            prepare_protein=False).sites.load_and_prepare_protein(str(iron_file))

        with pytest.raises(Exception, match='only'):
            self.derive(topology, positions)

    def test_reported_functional_follows_the_scf_driver(self):
        """
        xcfun only applies when no SCF driver is given, so the header has to
        report what the assigned driver will actually use.
        """

        pytest.importorskip('openmm')

        from veloxchem.scfrestdriver import ScfRestrictedDriver

        builder = MetalSiteForceFieldBuilder()

        assert builder._functional_label() == 'PBE0'

        builder.xcfun = None
        assert builder._functional_label() == 'Hartree-Fock'

        # an assigned driver wins, whatever xcfun says
        builder.xcfun = 'PBE0'
        builder.scf_drv = ScfRestrictedDriver()
        assert builder._functional_label() == 'Hartree-Fock'

        builder.scf_drv.xcfun = 'B3LYP'
        assert builder._functional_label() == 'B3LYP'

    def test_histidine_tautomer_is_chosen_from_the_ring_nitrogens(self):
        """
        The tautomer must be decided between ND1 and NE2 only.

        A ring carbon can end up closer to the metal than either nitrogen in a
        distorted or unrelaxed structure. A carbon cannot donate, so that says
        the geometry is off, not that the coordination is different, and it
        must not be allowed to influence the choice. It is reported as a note
        instead.
        """

        pytest.importorskip('openmm')

        prepared = self.prepared()
        topology = prepared['topology']
        positions = prepared['positions']

        histidine = next(residue for residue in topology.residues()
                         if residue.name.startswith('HI'))
        atoms = {atom.name: atom.index for atom in histidine.atoms()}

        # a metal placed almost on top of CD2, nearer than either nitrogen
        metal_position = positions[atoms['CD2']] + np.array([0.1, 0.0, 0.0])

        def distance(name):
            return np.linalg.norm(positions[atoms[name]] - metal_position)

        assert distance('CD2') < distance('ND1')
        assert distance('CD2') < distance('NE2')

        variant, note = self.shells().sites._histidine_variant(histidine, positions,
                                                               [metal_position])

        # whichever nitrogen is closer coordinates, so it carries no hydrogen
        closer = 'ND1' if distance('ND1') < distance('NE2') else 'NE2'
        assert variant == ('HIE' if closer == 'ND1' else 'HID')

        assert note is not None
        assert 'CD2' in note

    def test_histidine_note_is_absent_for_a_clean_geometry(self):
        """
        The distortion note must not fire on a well-formed site.

        Both histidines of the reference structure coordinate through NE2 with
        no carbon in the way, so they type as HID and produce no notes.
        """

        pytest.importorskip('openmm')

        prepared = self.prepared()
        binding_modes = prepared['binding_modes']

        assert sorted(binding_modes['variants'].values()) == [
            'ASP', 'ASP', 'CYX', 'GLU', 'HID', 'HID'
        ]

        for note in binding_modes['notes']:
            assert 'distorted' not in note

    def test_folder_is_generated_and_can_be_overridden(self):
        """
        Each builder gets its own output folder name.

        The name carries the time the builder was created so that two runs
        never quietly overwrite one another's results, and it is a plain
        attribute so a caller can put the files wherever it likes.
        """

        builder = MetalSiteForceFieldBuilder()

        assert builder.folder.startswith('metal_site_')

        builder.folder = 'somewhere_else'
        assert builder.folder == 'somewhere_else'

    def test_a_run_writes_every_artifact(self, tmp_path):
        """
        Checks that a run can be picked up again from what it wrote.

        Every step writes its own result as soon as it has it, so what the
        folder holds at the end is the whole run. The Hessian and the charges
        are indexed to the truncated cluster, so each file is read back and
        checked against what went in rather than merely existing. The force
        field is loaded into OpenMM, since an XML that cannot be parsed is no
        use whatever its size.
        """

        pytest.importorskip('openmm')

        import openmm.app as mmapp

        n_atoms = self.prepared()['active_site']['molecule'].number_of_atoms()
        hessian = self.random_hessian(n_atoms)
        charges = np.linspace(-0.4, 0.4, n_atoms)

        builder = self.builder(tmp_path)
        builder.build_active_site(self.structure_file())
        forcefield = builder.build_forcefield(hessian=hessian,
                                              partial_charges=charges)

        folder = Path(builder.folder)
        expected = {
            'protonated.pdb',
            util.MM_GEOMETRY_FILE,
            util.GEOMETRY_FILE,
            util.HESSIAN_FILE,
            util.CHARGES_FILE,
            util.FORCEFIELD_FILE,
        }
        assert {path.name for path in folder.iterdir()} == expected

        # everything has to survive the round trip, not just exist
        assert np.allclose(np.loadtxt(folder / util.HESSIAN_FILE), hessian)
        # the charges file holds what the force field carries, which is the
        # fit with the capping hydrogens' share folded back into the site --
        # the two used to disagree, the file keeping the raw fit
        written = np.loadtxt(folder / util.CHARGES_FILE)
        assert np.allclose(written, forcefield.partial_charges)
        assert np.allclose(
            written, util.redistribute_cap_charges(builder._active_site,
                                                   charges))
        for index in builder._active_site['cap_indices']:
            assert written[index] == pytest.approx(0.0, abs=1e-12)

        cluster = Molecule.read_xyz_file(str(folder / util.GEOMETRY_FILE))
        assert cluster.number_of_atoms() == n_atoms

        # the JSON is what a template reads, and it is the only force field
        # the folder carries: an OpenMM export is written by whoever wants
        # to simulate the cluster on its own, from the force field the run
        # produced, rather than by every run whether or not anyone asks
        reloaded = util.load_forcefield(folder / util.FORCEFIELD_FILE)
        reloaded.molecule = builder._active_site['molecule']
        reloaded.write_openmm_files(str(tmp_path / 'exported'))
        openmm_ff = mmapp.ForceField(str(tmp_path / 'exported.xml'))
        cluster_pdb = mmapp.PDBFile(str(tmp_path / 'exported.pdb'))
        system = openmm_ff.createSystem(cluster_pdb.topology,
                                        nonbondedMethod=mmapp.NoCutoff)
        assert system.getNumParticles() == n_atoms

        # the XML is what OpenMM reads; the JSON is what a template reads
        reloaded_ff = util.load_forcefield(folder / util.FORCEFIELD_FILE)
        assert reloaded_ff.bonds == forcefield.bonds
        assert reloaded_ff.angles == forcefield.angles

    def test_the_structural_half_writes_what_it_has(self, tmp_path):
        """
        A run that stops before any QM is paid for still leaves behind what
        it built: the protonated structure and the relaxed cluster.

        The coordination is not among them. It is derived from those two and
        the request, so writing it out would be a copy of an answer that can
        be asked for again, and the printed table is where it is reviewed.
        """

        pytest.importorskip('openmm')

        builder = self.builder(tmp_path)
        builder.build_active_site(self.structure_file())

        written = {path.name for path in Path(builder.folder).iterdir()}
        assert written == {
            'protonated.pdb',
            util.MM_GEOMETRY_FILE,
        }

    @staticmethod
    def plain_system(topology, charges):
        """
        Builds a bare OpenMM system carrying the given charges.

        redistribute_charges only needs a system with a NonbondedForce, so the
        tests do not have to build a real protein force field. The reference
        structure is a set of isolated residues, which amber would refuse for
        want of terminal templates.

        :param topology:
            The OpenMM topology.
        :param charges:
            One charge per atom.

        :return:
            The system.
        """

        import openmm as mm

        system = mm.System()
        force = mm.NonbondedForce()

        for atom, charge in zip(topology.atoms(), charges):
            system.addParticle(atom.element.mass)
            force.addParticle(charge, 0.3, 0.5)

        system.addForce(force)

        return system

    def test_d4_charges_stand_in_when_nothing_is_fitted(self, tmp_path, capsys):
        """
        A run without RESP still has charges: the force field, the enzyme
        system and partial_charges.txt all get the D4 ones, rather than the
        force field quietly using them while the rest of the run carries none.

        The basis set is unusable on purpose, so any QM step would fail rather
        than pass unnoticed.
        """

        pytest.importorskip('openmm')
        pytest.importorskip('pdbfixer')

        molecule = self.prepared()['active_site']['molecule']
        n_atoms = molecule.number_of_atoms()

        builder = self.builder(tmp_path, do_resp=False)
        builder.ostream = self.printing()
        builder.build_active_site(self.structure_file())
        builder.build_forcefield(opt_geometry=molecule,
                                 hessian=self.random_hessian(n_atoms))
        system, topology = builder.create_enzyme_system()

        charges = builder.partial_charges
        assert charges is not None
        assert len(charges) == n_atoms
        # D4 charges sum to the charge of the site
        assert charges.sum() == pytest.approx(float(molecule.get_charge()),
                                              abs=1.0e-6)

        # they reach the file and the enzyme system, not only the force field
        written = tmp_path / 'run' / util.CHARGES_FILE
        assert written.is_file()
        assert np.allclose(np.loadtxt(written), charges)
        assert system.getNumParticles() == topology.getNumAtoms()

        # and the run says where they came from, without listing them per
        # atom. The atom is picked out of a residue with more than one, since
        # a metal is a residue of its own and its per-residue total in the
        # summary is its own charge whether or not anything lists atoms.
        site = builder._active_site
        singles = set(site['metal_indices']) | set(site['cap_indices'])
        inside = next(i for i in range(n_atoms) if i not in singles)

        printed = capsys.readouterr().out
        assert 'D4 partial charges' in printed
        assert f'{charges[inside]:.4f}' not in printed

    def test_the_enzyme_system_carries_the_charges_it_is_given(self):
        """
        create_enzyme_system puts the charges it is handed -- the ones the
        force field beside it was built on -- onto the coordination sphere,
        rather than leaving the protein force field's own there.
        """

        pytest.importorskip('openmm')
        pytest.importorskip('pdbfixer')

        import openmm as mm

        builder = MetalSiteForceFieldBuilder()
        builder.ostream.mute()

        topology, positions = self.shells().sites.load_and_prepare_protein(
            self.structure_file())
        binding_modes = self.derive(topology, positions)
        topology, positions, binding_modes = self.protonated(
            topology, positions, binding_modes)
        active_site = self.shells().sites.extract_active_site(
            topology, positions, binding_modes)
        forcefield = self.seeded_forcefield(active_site)

        def charges_of(system):
            force = next(one for one in system.getForces()
                         if isinstance(one, mm.NonbondedForce))
            return [
                force.getParticleParameters(index)[0]._value
                for index in range(system.getNumParticles())
            ]

        system, _ = self.shells().enzyme.create_enzyme_system(
            topology, active_site, forcefield, util.d4_charges(active_site))

        # they are not simply amber's own charges left alone
        import openmm.app as mmapp

        bare = mmapp.ForceField(
            *self.shells().enzyme.protein_forcefield_files).createSystem(
            topology, nonbondedMethod=mmapp.NoCutoff)
        assert not np.allclose(charges_of(system), charges_of(bare))

    def test_metal_impropers_reach_the_enzyme_system(self):
        """
        The metal terms of the force field are bonds, angles and impropers,
        and all three have to be put onto the protein system: nothing in the
        protein force field describes a metal, so an improper left behind is
        simply gone.

        They are added by index rather than through a residue template,
        which is what a metal improper needs -- OpenMM builds one out of a
        template only for a central atom and three atoms bonded to it, and
        the improper a bidentate carboxylate gets names the metal, which its
        central carbon is not bonded to.
        """

        pytest.importorskip('openmm')
        pytest.importorskip('pdbfixer')

        import openmm as mm
        import openmm.app as mmapp

        topology, positions = self.shells().sites.load_and_prepare_protein(
            self.structure_file())
        binding_modes = self.derive(topology, positions)
        topology, positions, binding_modes = self.protonated(
            topology, positions, binding_modes)
        active_site = self.shells().sites.extract_active_site(
            topology, positions, binding_modes)
        forcefield = self.seeded_forcefield(active_site)

        impropers = util.get_metal_impropers(forcefield, active_site)
        assert impropers, 'the fixture should give the metals impropers'

        atom_map = active_site['atom_map']
        metals = {atom_map[index] for index in active_site['metal_indices']}

        def torsions_of(system):
            force = next(one for one in system.getForces()
                         if isinstance(one, mm.PeriodicTorsionForce))
            return {
                tuple(force.getTorsionParameters(index)[:4])
                for index in range(force.getNumTorsions())
            }

        bare = mmapp.ForceField('amber14-all.xml',
                                'amber14/tip3pfb.xml').createSystem(
                                    topology, nonbondedMethod=mmapp.NoCutoff)

        # nothing amber writes reaches a metal, so every one of these is
        # either transferred here or lost
        assert not any(metals & set(key) for key in torsions_of(bare))

        system, added = self.shells().enzyme.create_enzyme_system(
            topology, active_site, forcefield, util.d4_charges(active_site))
        built = torsions_of(system)

        for key in impropers:
            # the central atom of the key goes third, the way OpenMM orders
            # the same improper when it builds one from a residue template
            first, second, central, third = key[1], key[2], key[0], key[3]
            assert (atom_map[first], atom_map[second], atom_map[central],
                    atom_map[third]) in built

        assert [key for kind, key in added
                if kind == 'improper'] == list(impropers)

    def test_cap_charges_are_redistributed(self):
        """
        The capping hydrogens do not exist in the protein, so their charge has
        to go somewhere rather than simply vanish.

        On this site they carry about half an electron between them, which is
        far too much to discard quietly.
        """

        pytest.importorskip('openmm')

        prepared = self.prepared()
        active_site = prepared['active_site']
        caps = active_site['cap_indices']

        n_atoms = active_site['molecule'].number_of_atoms()
        charges = np.linspace(-0.5, 0.5, n_atoms)
        cap_charge = sum(charges[index] for index in caps)

        corrected = util.redistribute_cap_charges(active_site, charges)

        # the caps are emptied and the total is untouched
        assert all(corrected[index] == 0.0 for index in caps)
        assert corrected.sum() == pytest.approx(charges.sum())

        # and the charge is spread evenly over everything else
        rest = [i for i in range(n_atoms) if i not in caps]
        expected = cap_charge / len(rest)
        for index in rest:
            assert corrected[index] - charges[index] == pytest.approx(expected)

        # the input is left alone
        assert charges[caps[0]] != 0.0

    def test_redistribute_charges_preserves_the_region_charge(self):
        """
        Replacing the sidechain charges must not change the net charge.

        The fit only covers each coordinating residue from CB outwards, so the
        residues it cuts no longer carry what the protein force field intended.
        The difference is absorbed by the atoms the cluster does not cover, and
        the coordination region as a whole has to come out where it started -
        otherwise the protein's net charge moves, which is not allowed to
        happen.
        """

        pytest.importorskip('openmm')

        import openmm as mm

        prepared = self.prepared()
        topology = prepared['topology']
        active_site = prepared['active_site']

        n_atoms = active_site['molecule'].number_of_atoms()
        rng = np.random.default_rng(4)
        protein_charges = rng.normal(0.0, 0.3, topology.getNumAtoms())
        system = self.plain_system(topology, protein_charges)

        cluster_charges = rng.normal(0.0, 0.3, n_atoms)

        caps = set(active_site['cap_indices'])
        covered = {
            active_site['atom_map'][index]
            for index in range(n_atoms) if index not in caps
        }
        residues = {
            atom.residue.index
            for atom in topology.atoms() if atom.index in covered
        }
        region = [
            atom.index for residue in topology.residues()
            if residue.index in residues for atom in residue.atoms()
        ]

        before = sum(protein_charges[index] for index in region)

        cap_charge, shift = self.shells().enzyme.redistribute_charges(system, topology,
                                                                      active_site,
                                                                      cluster_charges)

        force = next(f for f in system.getForces()
                     if isinstance(f, mm.NonbondedForce))

        def charge_of(index):
            return force.getParticleParameters(index)[0].value_in_unit(
                mm.unit.elementary_charge)

        after = sum(charge_of(index) for index in region)

        assert after == pytest.approx(before, abs=1.0e-9)

        # the covered atoms carry the fitted charges, caps folded in
        corrected = util.redistribute_cap_charges(active_site, cluster_charges)
        for cluster_index in range(n_atoms):
            if cluster_index in caps:
                continue
            topology_index = active_site['atom_map'][cluster_index]
            assert charge_of(topology_index) == pytest.approx(
                corrected[cluster_index])

        # and the uncovered atoms all moved by the same amount
        uncovered = [index for index in region if index not in covered]
        assert uncovered

        # the alpha carbons are among them: the caps only stood in for them,
        # so the protein force field still owns the real atoms
        alpha_carbons = {
            atom.index
            for residue in topology.residues() if residue.index in residues
            for atom in residue.atoms() if atom.name == 'CA'
        }
        assert alpha_carbons <= set(uncovered)

        for index in uncovered:
            assert charge_of(index) - protein_charges[index] == pytest.approx(
                shift)

        assert cap_charge == pytest.approx(
            sum(cluster_charges[index] for index in caps))

    def test_redistribute_charges_rejects_a_mismatched_topology(self):
        """
        The atom map belongs to one topology and only that one.

        prepare_protein renumbers atoms, so running it after the extraction
        rather than before leaves the map pointing at the wrong atoms. That
        writes the cluster charges onto whatever happens to hold those indices
        and still balances the books, so it produces a plausible-looking system
        that is quietly wrong. It has to be refused.
        """

        pytest.importorskip('openmm')

        prepared = self.prepared()
        topology = prepared['topology']
        active_site = prepared['active_site']

        n_atoms = active_site['molecule'].number_of_atoms()
        charges = np.zeros(topology.getNumAtoms())
        system = self.plain_system(topology, charges)

        # shift every index by one, as a renumbering would
        shifted = dict(active_site)
        shifted['atom_map'] = {
            cluster: topology_index + 1
            for cluster, topology_index in active_site['atom_map'].items()
        }

        with pytest.raises(Exception, match='atom map'):
            self.shells().enzyme.redistribute_charges(system, topology, shifted,
                                                      np.zeros(n_atoms))

    def test_backbone_correction_is_separable_from_the_cap_correction(self):
        """
        The two corrections are separate methods because they apply in
        different places.

        The cap correction is expressible on the cluster alone, so
        build_forcefield applies it and the exported force field carries it.
        The backbone correction adjusts atoms the cluster deliberately leaves
        out, so it needs a protein system and can only run once one exists.
        redistribute_charges is the wrapper that does both, and has to agree
        with applying them in turn.
        """

        pytest.importorskip('openmm')

        import openmm as mm

        prepared = self.prepared()
        topology = prepared['topology']
        active_site = prepared['active_site']
        n_atoms = active_site['molecule'].number_of_atoms()

        rng = np.random.default_rng(7)
        protein_charges = rng.normal(0.0, 0.3, topology.getNumAtoms())
        cluster_charges = rng.normal(0.0, 0.3, n_atoms)

        def charges_of(system):
            force = next(f for f in system.getForces()
                         if isinstance(f, mm.NonbondedForce))
            return [
                force.getParticleParameters(index)[0].value_in_unit(
                    mm.unit.elementary_charge)
                for index in range(system.getNumParticles())
            ]

        # the wrapper
        combined = self.plain_system(topology, protein_charges)
        cap_charge, shift = self.shells().enzyme.redistribute_charges(combined, topology,
                                                                      active_site,
                                                                      cluster_charges)

        # the two steps by hand
        stepwise = self.plain_system(topology, protein_charges)
        folded = util.redistribute_cap_charges(active_site, cluster_charges)
        step_shift = self.shells().enzyme.redistribute_backbone_charges(stepwise, topology,
                                                                        active_site, folded)

        assert step_shift == pytest.approx(shift)
        assert cap_charge == pytest.approx(
            sum(cluster_charges[index] for index in active_site['cap_indices']))
        assert np.allclose(charges_of(combined), charges_of(stepwise))

    def test_backbone_correction_is_idempotent(self):
        """
        Running it twice must not move the charge twice.

        The first pass leaves the region where the protein force field had it,
        so the second finds nothing to recover and applies a zero shift.
        """

        pytest.importorskip('openmm')

        prepared = self.prepared()
        topology = prepared['topology']
        active_site = prepared['active_site']
        n_atoms = active_site['molecule'].number_of_atoms()

        rng = np.random.default_rng(11)
        system = self.plain_system(topology,
                                   rng.normal(0.0, 0.3, topology.getNumAtoms()))
        charges = util.redistribute_cap_charges(active_site,
                                                rng.normal(0.0, 0.3, n_atoms))

        first = self.shells().enzyme.redistribute_backbone_charges(system, topology,
                                                                   active_site, charges)
        second = self.shells().enzyme.redistribute_backbone_charges(system, topology,
                                                                    active_site, charges)

        assert first != 0.0
        assert second == pytest.approx(0.0, abs=1.0e-12)

    def test_the_enzyme_system_waits_for_the_force_field(self, tmp_path):
        """
        create_enzyme_system triggers nothing expensive: it is the last step
        of the pipeline, and asking for it before the fit is a mistake worth
        naming rather than a reason to start a QM run.
        """

        pytest.importorskip('openmm')

        builder = self.builder(tmp_path)

        with pytest.raises(Exception, match='no active site'):
            builder.create_enzyme_system()

        builder.build_active_site(self.structure_file())

        with pytest.raises(Exception, match='no force field'):
            builder.create_enzyme_system()

    def test_the_active_site_property_exposes_the_full_dictionary(
            self, tmp_path):
        """
        active_site is the same dictionary active_site_molecule reads its
        molecule off, so MetalForceFieldManager can describe and match it
        without reaching into a private attribute.
        """

        pytest.importorskip('openmm')

        builder = self.editable(tmp_path)

        assert builder.active_site is None

        builder.build_active_site(self.structure_file(), mm_opt=False)

        assert builder.active_site is not None
        assert builder.active_site['molecule'] is builder.active_site_molecule
        assert 'connectivity_matrix' in builder.active_site
        assert 'metal_indices' in builder.active_site

    def test_print_active_site_needs_an_active_site(self, tmp_path):
        """
        print_active_site is a read of the same stage-gated state
        show_active_site draws, and is refused the same way before one
        exists.
        """

        pytest.importorskip('openmm')

        builder = self.editable(tmp_path)

        with pytest.raises(Exception, match='no active site'):
            builder.print_active_site()

        builder.build_active_site(self.structure_file(), mm_opt=False)
        builder.print_active_site()

    def test_adopt_forcefield_requires_an_active_site(self, tmp_path):
        """
        Adopting a transferred force field makes no sense before there is an
        active site for it to describe.
        """

        builder = self.editable(tmp_path)

        with pytest.raises(Exception, match='no active site'):
            builder.adopt_forcefield(object())

    def test_adopt_forcefield_moves_to_fitted_without_a_fit(self, tmp_path):
        """
        adopt_forcefield is how MetalForceFieldManager.build_ff_from_template
        hands a force field transferred from a template back to the builder:
        it takes the object as given, with no fit of its own, and moves the
        builder from ACTIVE_SITE straight to FITTED so create_enzyme_system
        becomes legal.
        """

        pytest.importorskip('openmm')

        builder = self.editable(tmp_path)
        builder.build_active_site(self.structure_file(), mm_opt=False)
        assert builder.stage == 'active_site'

        n_atoms = builder.active_site_molecule.number_of_atoms()
        forcefield = SimpleNamespace(
            partial_charges=np.linspace(-0.4, 0.4, n_atoms))

        adopted = builder.adopt_forcefield(forcefield)

        assert adopted is forcefield
        assert builder.stage == 'fitted'
        assert builder.active_site_forcefield is forcefield
        assert np.allclose(builder.partial_charges, forcefield.partial_charges)

    def test_adopt_forcefield_lets_create_enzyme_system_run(self, tmp_path):
        """
        The end-to-end promise: a force field fitted for one builder and
        adopted onto another, structurally identical one is enough for
        create_enzyme_system to run on the second, without that builder ever
        calling its own build_forcefield.
        """

        pytest.importorskip('openmm')
        pytest.importorskip('pdbfixer')

        source = self.builder(tmp_path, name='source', do_resp=False)
        source.build_active_site(self.structure_file())
        source.build_forcefield(
            hessian=self.covering_hessian(source._active_site))

        target = self.builder(tmp_path, name='target', do_resp=False)
        target.build_active_site(self.structure_file())
        assert target.stage == 'active_site'

        target.adopt_forcefield(source.active_site_forcefield)
        assert target.stage == 'fitted'

        system, topology = target.create_enzyme_system()

        assert system.getNumParticles() == topology.getNumAtoms()

    def test_the_relaxation_runs_on_the_fitted_terms_once_there_are_any(
            self, tmp_path, monkeypatch):
        """
        The crude pass seeds its metal equilibria by measuring them on the
        geometry in front of it, so a contact the cutoffs did not close is
        already at its own minimum and the relaxation cannot move it. Once a
        fit exists its metal terms say more than that -- a per-bond
        equilibrium and a per-bond force constant -- so the relaxation runs
        on those instead of building the seeded ones again.
        """

        pytest.importorskip('openmm')
        pytest.importorskip('pdbfixer')

        builder = self.builder(tmp_path, do_resp=False)
        builder.build_active_site(self.structure_file())
        builder.build_forcefield(
            hessian=self.covering_hessian(builder._active_site))
        fitted = builder.active_site_forcefield

        seeded = []
        monkeypatch.setattr(
            ActiveSiteBuilder, 'build_forcefield',
            lambda *args, **kwargs: seeded.append(kwargs) or fitted)

        relaxed_on = []
        original = ActiveSiteBuilder.mm_optimize_active_site
        monkeypatch.setattr(
            ActiveSiteBuilder, 'mm_optimize_active_site',
            lambda sites, site, forcefield, **kwargs: relaxed_on.append(
                forcefield) or original(sites, site, forcefield, **kwargs))

        builder.mm_optimize_active_site()

        assert relaxed_on == [fitted]
        assert seeded == []

    def test_the_relaxation_seeds_a_force_field_before_there_is_a_fit(
            self, tmp_path, monkeypatch):
        """
        The other side of it: with nothing fitted yet there is nothing to
        take, so the seeded force field is built here as it always was, and
        the manual equilibria reach it.
        """

        pytest.importorskip('openmm')

        builder = self.editable(tmp_path)
        builder.build_active_site(self.structure_file(), mm_opt=False)

        relaxed_on = []
        original = ActiveSiteBuilder.mm_optimize_active_site
        monkeypatch.setattr(
            ActiveSiteBuilder, 'mm_optimize_active_site',
            lambda sites, site, forcefield, **kwargs: relaxed_on.append(
                forcefield) or original(sites, site, forcefield, **kwargs))

        builder.mm_optimize_active_site()

        assert len(relaxed_on) == 1
        assert relaxed_on[0] is not None
        assert builder.stage == 'active_site'

    def test_a_relaxation_keeps_an_adopted_fit_and_drops_one_made_here(
            self, tmp_path):
        """
        Whether a new geometry invalidates the fit above it depends on where
        that fit came from. One made here was fitted to a Hessian computed on
        the geometry the relaxation just replaced, so it goes. One adopted
        from a template describes the template rather than this geometry, and
        relaxing toward its parameters brings the two closer together, so it
        stays and create_enzyme_system remains legal.
        """

        pytest.importorskip('openmm')
        pytest.importorskip('pdbfixer')

        source = self.builder(tmp_path, name='source', do_resp=False)
        source.build_active_site(self.structure_file())
        source.build_forcefield(
            hessian=self.covering_hessian(source._active_site))

        assert source.stage == 'fitted'
        # taken before the relaxation, which is what drops it from the
        # builder that made it
        fitted = source.active_site_forcefield

        source.mm_optimize_active_site()
        assert source.stage == 'active_site'
        assert source.active_site_forcefield is None

        target = self.builder(tmp_path, name='target', do_resp=False)
        target.build_active_site(self.structure_file())
        adopted = target.adopt_forcefield(fitted)

        before = target.active_site_molecule.get_coordinates_in_angstrom()
        target.mm_optimize_active_site()
        after = target.active_site_molecule.get_coordinates_in_angstrom()

        assert target.stage == 'fitted'
        assert target.active_site_forcefield is adopted
        assert not np.allclose(before, after)

    def test_intermediates_are_reused_from_the_folder(self):
        """
        A folder left behind by an earlier run is picked up automatically.

        This is what makes a failure part way through cheap: the geometry
        optimization and the Hessian each cost hours on a real site, and a
        crash in a later step must not throw away the earlier ones. Each step
        writes its result as soon as it has it, and finds it again next time.
        """

        pytest.importorskip('openmm')

        import tempfile

        prepared = self.prepared()
        active_site = prepared['active_site']
        n_atoms = active_site['molecule'].number_of_atoms()

        hessian = self.random_hessian(n_atoms)
        charges = np.linspace(-0.4, 0.4, n_atoms)
        geometry = active_site['molecule']

        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)

            # nothing there yet
            assert self.shells().qm._resolve_hessian(active_site, folder=path, hessian=None) is None
            assert self.shells().qm._resolve_partial_charges(active_site,
                                                             folder=path, partial_charges=None) is None
            assert self.shells().qm._resolve_optimized_geometry(
                active_site['molecule'], folder=path, optimized_geometry=None) is None

            np.savetxt(path / util.HESSIAN_FILE, hessian)
            np.savetxt(path / util.CHARGES_FILE, charges)
            geometry.write_xyz_file(str(path / util.GEOMETRY_FILE))

            # found in the folder, without being handed the values
            assert np.allclose(
                self.shells().qm._resolve_hessian(active_site, folder=path, hessian=None),
                hessian)
            assert np.allclose(
                self.shells().qm._resolve_partial_charges(active_site, folder=path, partial_charges=None),
                charges)
            reused = self.shells().qm._resolve_optimized_geometry(
                active_site['molecule'], folder=path, optimized_geometry=None)
            assert reused.number_of_atoms() == n_atoms

    def test_something_handed_in_beats_the_folder(self):
        """
        What a call is given has to win over whatever the folder happens to
        hold, or a stale run could not be overridden.
        """

        pytest.importorskip('openmm')

        import tempfile

        prepared = self.prepared()
        active_site = prepared['active_site']
        n_atoms = active_site['molecule'].number_of_atoms()

        from_folder = self.random_hessian(n_atoms, seed=1)
        explicit = self.random_hessian(n_atoms, seed=2)

        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)
            np.savetxt(path / util.HESSIAN_FILE, from_folder)

            assert np.allclose(
                self.shells().qm._resolve_hessian(active_site,
                                                  folder=path,
                                                  hessian=explicit), explicit)

    def test_the_forcefield_survives_a_json_round_trip(self, tmp_path):
        """
        A fitted force field can be written and read back on its own.

        The QM behind the fitted metal terms is what a run pays for, so the
        parameters have to be storable without the structure they came from.
        The file carries no geometry, so the caller puts the molecule of the
        active site back on after load_forcefield: without it the generator
        cannot be written out for OpenMM.
        """

        pytest.importorskip('openmm')

        import openmm.app as mmapp

        prepared = self.prepared()
        active_site = prepared['active_site']
        n_atoms = active_site['molecule'].number_of_atoms()

        charges = np.linspace(-0.4, 0.4, n_atoms)
        forcefield = self.fitted_forcefield(active_site,
                                            self.random_hessian(n_atoms),
                                            partial_charges=charges)

        path = tmp_path / 'metal_forcefield.json'
        util.save_forcefield(path, forcefield)
        reloaded = util.load_forcefield(path)
        reloaded.molecule = active_site['molecule']

        assert reloaded.bonds == forcefield.bonds
        assert reloaded.angles == forcefield.angles
        assert reloaded.dihedrals == forcefield.dihedrals
        assert reloaded.impropers == forcefield.impropers
        assert reloaded.atoms == forcefield.atoms

        # the fitted metal terms are the point of the file
        bonds, angles = util.get_metal_keys(reloaded, active_site)
        fitted_bonds, fitted_angles = util.get_metal_keys(
            forcefield, active_site)
        assert bonds == fitted_bonds
        assert angles == fitted_angles
        assert bonds and angles

        # the charges the file carries are the redistributed ones, so the
        # caps have already given theirs away
        expected = util.redistribute_cap_charges(active_site, charges)
        assert np.allclose(reloaded.partial_charges, expected)

        reloaded.write_openmm_files(str(tmp_path / 'reloaded'))
        openmm_ff = mmapp.ForceField(str(tmp_path / 'reloaded.xml'))
        cluster_pdb = mmapp.PDBFile(str(tmp_path / 'reloaded.pdb'))
        system = openmm_ff.createSystem(cluster_pdb.topology,
                                        nonbondedMethod=mmapp.NoCutoff)
        assert system.getNumParticles() == n_atoms

    def test_a_forcefield_from_another_structure_is_rejected(self, tmp_path):
        """
        The bonds and angles are keyed by plain atom indices, which fit any
        cluster of the same size, so a force field is validated against the
        elements of the extracted active site rather than trusted.
        """

        pytest.importorskip('openmm')

        prepared = self.prepared()
        active_site = prepared['active_site']
        n_atoms = active_site['molecule'].number_of_atoms()

        forcefield = self.seeded_forcefield(active_site)

        # a force field for a cluster with one atom too few
        short = tmp_path / 'short.json'
        util.save_forcefield(short, forcefield)
        del forcefield.atoms[n_atoms - 1]
        util.save_forcefield(short, forcefield)

        with pytest.raises(Exception, match='atoms'):
            util._check_forcefield(util.load_forcefield(short), active_site, 'the folder')

        # and one of the right size describing something else
        swapped = tmp_path / 'swapped.json'
        forcefield = self.seeded_forcefield(active_site)
        metal = active_site['metal_indices'][0]
        forcefield.atoms[metal]['name'] = 'C99'
        util.save_forcefield(swapped, forcefield)

        with pytest.raises(Exception, match='elements'):
            util._check_forcefield(util.load_forcefield(swapped), active_site, 'the folder')

    def test_a_folder_from_another_structure_is_rejected(self):
        """
        Reuse is keyed on the folder alone, so a folder belonging to a
        different site would otherwise be picked up silently. The stored
        intermediates are validated against the extracted cluster.
        """

        pytest.importorskip('openmm')

        import tempfile

        prepared = self.prepared()
        active_site = prepared['active_site']

        with tempfile.TemporaryDirectory() as folder:
            # a Hessian for some other, smaller cluster
            np.savetxt(Path(folder) / util.HESSIAN_FILE, np.zeros((9, 9)))

            with pytest.raises(Exception, match='shape'):
                self.shells().qm._resolve_hessian(active_site, folder=Path(folder), hessian=None)

    def test_a_hessian_for_another_coordination_is_ignored(self, capsys):
        """
        A file left in the folder by a run whose coordination differed passes
        the shape and element checks and then fits zeros, so the blocks the
        metal terms actually read are looked at as well.

        An explicitly supplied Hessian is an instruction and is left alone;
        one lying in a folder is a guess.
        """

        pytest.importorskip('openmm')

        import tempfile

        prepared = self.prepared()
        active_site = prepared['active_site']
        n_atoms = active_site['molecule'].number_of_atoms()

        hessian = self.random_hessian(n_atoms)
        metal = sorted(active_site['metal_indices'])[0]
        ligand = next(index for index in range(n_atoms)
                      if active_site['connectivity_matrix'][metal, index])

        # exactly what a Hessian computed for a coordination without this
        # bond leaves behind
        for first, second in ((metal, ligand), (ligand, metal)):
            hessian[3 * first:3 * first + 3, 3 * second:3 * second + 3] = 0.0

        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)
            np.savetxt(path / util.HESSIAN_FILE, hessian)

            capsys.readouterr()
            reused = self.shells(self.printing()).qm._resolve_hessian(active_site,
                                                                      folder=path, hessian=None)
            captured = capsys.readouterr().out

            assert reused is None
            assert 'different coordination' in captured

            # handed in rather than found, it is used as given
            given = self.shells().qm._resolve_hessian(active_site,
                                                      folder=path,
                                                      hessian=hessian)
            assert np.allclose(given, hessian)

    def test_an_edit_survives_a_rebuild(self, tmp_path):
        """
        The point of rebuilding without a path: the coordination the cutoffs
        proposed is corrected, and the site is extracted again from the
        correction rather than from the structure file.
        """

        pytest.importorskip('openmm')

        builder = self.builder(tmp_path)
        builder.build_active_site(self.structure_file(), mm_opt=False)

        before = self.coordination(builder.binding_modes)
        assert ('ASP130', 'OD1', 'bidentate') in before

        builder.remove_metal_bond('ASP130', atom='OD1')
        rebuilt = builder.build_active_site(mm_opt=False)

        after = self.coordination(builder.binding_modes)
        assert ('ASP130', 'OD1', 'bidentate') not in after
        assert ('ASP130', 'OD2', 'monodentate') in after

        # the site itself was extracted again, one metal bond lighter
        metals = set(builder._active_site['metal_indices'])
        matrix = builder._active_site['connectivity_matrix']
        assert rebuilt.number_of_atoms(
        ) == builder._active_site['molecule'].number_of_atoms()
        assert sum(matrix[metal].sum() for metal in metals) == sum(
            len(ligand['metals'])
            for ligand in builder.binding_modes['ligands'])

    def test_the_results_are_read_off_the_properties(self, tmp_path):
        """
        What a step produced is read back off the builder rather than carried
        around by the caller, and nothing writes to it from outside.
        """

        pytest.importorskip('openmm')

        n_atoms = self.prepared()['active_site']['molecule'].number_of_atoms()

        builder = self.builder(tmp_path)

        assert builder.active_site_molecule is None
        assert builder.active_site_forcefield is None
        assert builder.enzyme_system is None

        molecule = builder.build_active_site(self.structure_file())
        forcefield = builder.build_forcefield(
            hessian=self.random_hessian(n_atoms),
            partial_charges=np.zeros(n_atoms))
        system, topology = builder.create_enzyme_system()

        assert builder.active_site_molecule is molecule
        assert builder.active_site_forcefield is forcefield
        assert builder.enzyme_system is system
        assert builder.enzyme_topology is topology
        assert builder.enzyme_positions.shape == (topology.getNumAtoms(), 3)
        assert builder.hessian.shape == (3 * n_atoms, 3 * n_atoms)
        assert np.allclose(builder.partial_charges, 0.0)

        for name in ('active_site_molecule', 'active_site_forcefield',
                     'enzyme_system', 'hessian', 'partial_charges'):
            with pytest.raises(AttributeError):
                setattr(builder, name, None)

    def test_the_optimization_constraint_can_be_read_first(self, tmp_path):
        """
        What is frozen is what the fitted parameters end up describing, so it
        has to be readable without paying for the optimization that uses it,
        and it has to be the same string that optimization runs under.
        """

        pytest.importorskip('openmm')

        builder = self.builder(tmp_path)

        assert builder.optimization_constraints is None

        builder.build_active_site(self.structure_file(), mm_opt=False)
        active_site = builder._active_site

        beta_carbons = sorted(active_site['beta_carbon_indices'])
        expected = ','.join(str(index + 1) for index in beta_carbons)

        assert builder.optimization_constraints == f'freeze xyz {expected}'
        assert builder.optimization_constraints == util.freeze_constraints(
            util.constrained_indices(active_site, False))

        # the setting reaches it, and so does an empty selection
        builder.constrain_capping_hydrogens = True
        with_caps = builder.optimization_constraints
        assert len(with_caps.split(',')) == len(beta_carbons) + len(
            active_site['cap_indices'])
        assert util.freeze_constraints([]) is None

        with pytest.raises(AttributeError):
            builder.optimization_constraints = 'freeze xyz 1'

    def test_the_hessian_gets_an_scf_for_its_own_geometry(self, monkeypatch):
        """
        ScfHessianDriver runs an SCF only when scf_drv.scf_results is empty,
        and never checks which geometry the results it finds belong to. A
        driver that has been used before -- one assigned to scf_drv, or one
        that has just been through the optimization -- therefore hands the
        Hessian the densities of another structure, silently. compute_hessian
        runs the SCF itself to stop that.

        This once regressed by losing a single line while the comment above it
        survived, so it is pinned here rather than left to a QM run to notice.
        """

        pytest.importorskip('openmm')

        prepared = self.prepared()
        active_site = prepared['active_site']

        order = []

        class FakeScfDriver:

            def __init__(self):
                self.ostream = OutputStream(None)
                # results from some other geometry, as a used driver holds
                self.scf_results = {'stale': True}
                self.solvation_model = None

            def compute(self, molecule, basis):
                order.append(('scf', molecule.number_of_atoms()))
                self.scf_results = {'stale': False}

        class FakeHessianDriver:

            def __init__(self, scf_driver):
                self.scf_driver = scf_driver
                self.ostream = OutputStream(None)
                self.numerical = True
                self.atom_pairs = None
                self.hessian = None

            def compute(self, molecule, basis):
                order.append(('hessian', self.scf_driver.scf_results['stale']))
                n_atoms = molecule.number_of_atoms()
                self.hessian = np.zeros((3 * n_atoms, 3 * n_atoms))

        driver = FakeScfDriver()
        monkeypatch.setattr(QmParameterizer, '_get_scf_driver',
                            lambda qm_self, molecule, *args: (driver, None))
        monkeypatch.setattr(qm, 'ScfHessianDriver', FakeHessianDriver)

        hessian = self.shells().qm.compute_hessian(active_site, None)

        n_atoms = active_site['molecule'].number_of_atoms()
        assert hessian.shape == (3 * n_atoms, 3 * n_atoms)
        # the SCF ran first, on this geometry, and the Hessian saw its results
        assert order == [('scf', n_atoms), ('hessian', False)]

    # ------------------------------------------------------------------
    # metal coordination planarity impropers
    # ------------------------------------------------------------------

    def test_a_coordinating_histidine_gets_a_planarity_improper(self):
        """
        An imidazole nitrogen's coordinating lone pair sits in the ring
        plane, so a weak improper nudges the metal it grips into that
        plane too, keyed on the ring nitrogen against its two ring
        neighbours.
        """

        labels = ['Zn', 'N', 'C', 'C']
        coordinates = np.array([
            [-2.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.6, 1.1, 0.0],
            [0.6, -1.1, 0.0],
        ])
        molecule = Molecule(labels, coordinates, 'angstrom')
        active_site = {'molecule': molecule, 'metal_indices': [0]}

        matrix = np.zeros((4, 4), dtype=int)
        for first, second in [(0, 1), (1, 2), (1, 3)]:
            matrix[first, second] = 1
            matrix[second, first] = 1

        forcefield = SimpleNamespace(connectivity_matrix=matrix, impropers={})

        added = self.shells().sites._add_metal_planarity_impropers(forcefield, active_site)

        assert added == 1
        assert (1, 2, 3, 0) in forcefield.impropers
        term = forcefield.impropers[(1, 2, 3, 0)]
        assert term['type'] == 'Fourier'
        assert term['barrier'] == pytest.approx(
            util.DEFAULT_METAL_PLANARITY_FORCE_CONSTANT)
        assert term['phase'] == 180.0
        assert term['periodicity'] == 2

    def test_a_bidentate_carboxylate_gets_a_planarity_improper(self):
        """
        A carboxylate gripping one metal with both oxygens has both donor
        lone pairs in its own plane, so the same nudge applies, keyed on
        the carboxylate carbon against the two oxygens.
        """

        molecule, forcefield, metals = self.mu_1_3_forcefield()

        # both arms onto Zn0, which is what makes it bidentate rather than
        # the mu-1,3 bridge the fixture otherwise holds
        forcefield.bonds[(0, 3)] = forcefield.bonds.pop((1, 3))
        forcefield.connectivity_matrix[1, 3] = 0
        forcefield.connectivity_matrix[3, 1] = 0
        forcefield.connectivity_matrix[0, 3] = 1
        forcefield.connectivity_matrix[3, 0] = 1

        active_site = {'molecule': molecule, 'metal_indices': metals}

        added = self.shells().sites._add_metal_planarity_impropers(forcefield, active_site)

        assert added == 1
        assert (4, 2, 3, 0) in forcefield.impropers
        term = forcefield.impropers[(4, 2, 3, 0)]
        assert term['barrier'] == pytest.approx(
            util.DEFAULT_METAL_PLANARITY_FORCE_CONSTANT)

    def test_a_mu_1_3_carboxylate_gets_no_planarity_improper(self):
        """
        A carboxylate reaching two metals through two different oxygens has
        no single metal in its plane to nudge towards, so nothing is added
        -- the structural detection off connectivity_matrix alone is enough
        to tell this apart from a genuine bidentate grip.
        """

        molecule, forcefield, metals = self.mu_1_3_forcefield()
        active_site = {'molecule': molecule, 'metal_indices': metals}

        before = dict(forcefield.impropers)
        added = self.shells().sites._add_metal_planarity_impropers(forcefield, active_site)

        assert added == 0
        assert forcefield.impropers == before

    def test_build_forcefield_adds_planarity_impropers_to_both_passes(self):
        """
        build_forcefield is the one function both the seeded crude pass and
        the fitted pass call, so the same nudge has to show up whichever
        way it was reached.
        """

        pytest.importorskip('openmm')

        prepared = self.prepared()
        active_site = prepared['active_site']

        his_n = self.donor_index(active_site, 187, 'NE2')
        asp_c = self.donor_index(active_site, 130, 'CG')

        seeded = self.seeded_forcefield(active_site)
        fitted = self.fitted_forcefield(active_site,
                                        self.covering_hessian(active_site))

        def planarity_terms(forcefield, center):
            return [
                term for key, term in forcefield.impropers.items()
                if key[0] == center
                and term['comment'] == 'metal coordination planarity restraint'
            ]

        for forcefield in (seeded, fitted):
            # ASP130's carboxylate carbon also carries GAFF's own ordinary
            # ring-planarity improper against its beta carbon; only the one
            # keyed on the metal is ours to check for
            assert len(planarity_terms(forcefield, his_n)) == 1
            assert len(planarity_terms(forcefield, asp_c)) == 1

    def test_metal_planarity_impropers_can_be_switched_off(self):
        """
        add_metal_planarity_impropers is the builder-facing knob, and it
        has to actually reach build_forcefield.
        """

        pytest.importorskip('openmm')

        prepared = self.prepared()
        active_site = prepared['active_site']

        forcefield = self.seeded_forcefield(active_site,
                                            add_metal_planarity_impropers=False)

        assert not any(
            term['comment'] == 'metal coordination planarity restraint'
            for term in forcefield.impropers.values())

    # ------------------------------------------------------------------
    # editing the active site
    # ------------------------------------------------------------------

    def editable(self, tmp_path, name='edit', **settings):
        """
        A builder for the edit methods, which skips the pdbfixer repair.

        None of these tests builds an enzyme system, which is the only thing
        that needs the structure to match a protein force field's templates,
        and the repair costs seconds on every call.

        :return:
            The builder.
        """

        settings.setdefault('prepare_protein', False)

        return self.builder(tmp_path, name=name, **settings)

    @staticmethod
    def covering_hessian(active_site, seed=0):
        """
        A Hessian holding data for exactly the pairs this coordination
        reads, and nothing else -- what a reused file fails to be once the
        coordination has moved on. Deliberately the strict extract_pairs
        walk rather than the forgiving hessian_pairs one a real run makes,
        so that a bond added afterwards is genuinely uncovered.

        :param active_site:
            The active site whose pairs are to be covered.
        :param seed:
            The seed of the random numbers.

        :return:
            The Hessian.
        """

        n_atoms = active_site['molecule'].number_of_atoms()
        full = TestMetalSiteForceFieldBuilder.random_hessian(n_atoms, seed)
        hessian = np.zeros_like(full)

        pairs, _ = util.extract_pairs(active_site['connectivity_matrix'],
                                      active_site['metal_indices'],
                                      bond_count=2)

        for first, second in pairs:
            for i, j in ((first, second), (second, first)):
                hessian[3 * i:3 * i + 3,
                        3 * j:3 * j + 3] = (full[3 * i:3 * i + 3,
                                                 3 * j:3 * j + 3])

        for index in range(n_atoms):
            hessian[3 * index:3 * index + 3,
                    3 * index:3 * index + 3] = (full[3 * index:3 * index + 3,
                                                     3 * index:3 * index + 3])

        return hessian

    @staticmethod
    def metal_label(builder, position=0):
        """
        The metal atom index as the active site labels it, which is what the
        edit methods take.

        :param builder:
            The builder holding the site.
        :param position:
            Which metal center is wanted.

        :return:
            The atom index, as an int.
        """

        active_site = builder._active_site
        index = active_site['metal_indices'][position]

        return int(active_site['labels'][index])

    @staticmethod
    def donor_index(active_site, resid, name):
        """
        The active site index of one named atom of one residue, found the way
        the labels say it should be: the beta carbon carries the resid and the
        rest of its sidechain follows it.

        :param active_site:
            The active site.
        :param resid:
            The residue id the beta carbon is labelled with.
        :param name:
            The atom name.

        :return:
            The index into the active site.
        """

        labels = active_site['labels']
        beta = next(index for index in active_site['beta_carbon_indices']
                    if labels[index] == str(resid))

        return labels.index(name, beta)

    def fitted(self, tmp_path, name='fit', mm_opt=False, **settings):
        """
        A builder that has been through build_forcefield on a Hessian
        covering its own coordination, so that the terms of the fit mean
        something and an edit after it can be told apart from them.

        :return:
            The builder.
        """

        builder = self.editable(tmp_path, name=name, **settings)
        builder.build_active_site(self.structure_file(), mm_opt=mm_opt)
        builder.build_forcefield(
            hessian=self.covering_hessian(builder._active_site))

        return builder

    def test_the_protonation_of_a_residue_can_be_changed(self, tmp_path):
        """
        The variants are chosen from the coordination, and this is how one of
        them is overruled: a doubly protonated histidine carries one hydrogen
        more than the tautomer and one unit more of charge.
        """

        pytest.importorskip('openmm')

        builder = self.editable(tmp_path)
        builder.build_active_site(self.structure_file(), mm_opt=False)

        before = builder.active_site_molecule
        # the rules pick a neutral tautomer, whichever way round it is
        variants = builder.binding_modes['variants']
        assert variants[4] in ('HID', 'HIE')

        builder.update_protonation_state('187', 'HIP')

        after = builder.active_site_molecule
        assert builder.binding_modes['variants'][4] == 'HIP'
        assert after.get_charge() == before.get_charge() + 1
        assert after.number_of_atoms() == before.number_of_atoms() + 1

    def test_a_protonation_change_is_a_setting_and_survives_a_rebuild(
            self, tmp_path):
        """
        The request goes into protonation_overrides, which is the same thing
        the setting is, so nothing distinguishes a variant asked for by hand
        from one that was there before the run started.
        """

        pytest.importorskip('openmm')

        builder = self.editable(tmp_path)
        builder.build_active_site(self.structure_file(), mm_opt=False)
        builder.update_protonation_state('187', 'HIP')

        assert builder._protonation_overrides == {'HIS187': 'HIP'}

        builder.build_active_site(mm_opt=False)

        assert builder.binding_modes['variants'][4] == 'HIP'

    def test_a_variant_the_residue_does_not_have_is_refused(self, tmp_path):
        """
        Modeller refuses an illegal variant from deep inside addHydrogens,
        with a message naming neither the residue nor what it would accept.
        The check happens before it is called, so the error says both.
        """

        pytest.importorskip('openmm')

        builder = self.editable(tmp_path)
        builder.build_active_site(self.structure_file(), mm_opt=False)

        # a variant of another residue, a variant of nothing at all, and the
        # Amber name for the thiolate, which Modeller calls CYX
        for resid, variant in (('187', 'ASH'), ('187', 'HIZ'), ('77', 'CYM')):
            with pytest.raises(Exception) as error:
                builder.update_protonation_state(resid, variant)
            assert 'is not a protonation variant of' in str(error.value)

    def test_a_variant_with_no_charge_of_its_own_is_refused(
            self, tmp_path, monkeypatch):
        """
        A variant OpenMM can build but the charge table does not know would
        fail much later, on the charge of the cluster, so both halves are
        checked in the same place.
        """

        pytest.importorskip('openmm')

        builder = self.editable(tmp_path)
        builder.build_active_site(self.structure_file(), mm_opt=False)

        charges = dict(util.VARIANT_CHARGES)
        charges.pop('HIP')
        # check_variant reads the table where it is defined
        monkeypatch.setattr(util, 'VARIANT_CHARGES', charges)

        with pytest.raises(Exception) as error:
            builder.update_protonation_state('187', 'HIP')

        assert 'check_variant: no charge known for variant HIP' in str(
            error.value)

    def test_a_residue_can_join_the_site_without_coordinating_anything(
            self, tmp_path):
        """
        The cluster is otherwise exactly the residues that reach a metal.
        A second-shell residue is put in by saying so, and it stays a
        spectator: it is truncated and capped like the rest and makes no
        metal bond at all.

        The contrast with coordinating_residues, which puts the same residue
        in as a ligand of its nearest metal, is the whole point of having
        both.
        """

        pytest.importorskip('openmm')

        structure = self.pushed_structure(tmp_path / 'far.pdb')

        builder = self.editable(tmp_path)
        builder.build_active_site(structure, mm_opt=False)

        # too far out for any cutoff to reach it
        assert 'ASP130' not in builder._active_site['residues']
        before = builder.active_site_molecule.number_of_atoms()

        builder.include_residue('130')

        assert 'ASP130' in builder._active_site['residues']
        assert builder.active_site_molecule.number_of_atoms() > before

        # in the cluster, and bonded to no metal in it
        assert 'ASP130' not in {
            ligand['residue']
            for ligand in builder.binding_modes['ligands']
        }

        metals = set(builder._active_site['metal_indices'])
        matrix = builder._active_site['connectivity_matrix']
        assert sum(matrix[metal].sum() for metal in metals) == sum(
            len(ligand['metals'])
            for ligand in builder.binding_modes['ligands'])

    def test_an_included_residue_is_protonated_for_the_ph(self, tmp_path):
        """
        The variant rules are all about a sidechain that coordinates a metal,
        so a residue that coordinates nothing is left at the default for the
        pH -- and what that came out as is recorded rather than guessed at
        later from the residue name.
        """

        pytest.importorskip('openmm')

        structure = self.pushed_structure(tmp_path / 'far.pdb')

        builder = self.editable(tmp_path)
        builder.build_active_site(structure, mm_opt=False)
        builder.include_residue('130')

        variants = builder.binding_modes['variants']
        assert variants[2] == 'ASP'

        charged = builder.active_site_molecule.get_charge()

        builder.update_protonation_state('130', 'ASH')

        assert builder.binding_modes['variants'][2] == 'ASH'
        assert builder.active_site_molecule.get_charge() == charged + 1

    def test_a_residue_the_truncation_cannot_cut_is_refused(self):
        """
        The scheme cuts every sidechain at CA-CB, which a glycine has no CB
        for and which leaves a proline's CD hanging off nothing. A residue
        that got into the site by coordinating a metal is never one of these,
        but include_residue takes any residue at all.
        """

        class Residue:

            def __init__(self, name, id):
                self.name = name
                self.id = id

        for name in ('GLY', 'PRO'):
            with pytest.raises(Exception) as error:
                self.shells().sites.check_truncatable(Residue(name, '52'))
            assert f'{name}52' in str(error.value)

        # and the residues the site is normally made of pass
        for name in ('ASP', 'GLU', 'CYS', 'HIS', 'TYR'):
            self.shells().sites.check_truncatable(Residue(name, '52'))

    def test_a_residue_can_be_taken_out_of_the_site(self, tmp_path):
        """
        Taking a residue out takes its metal bonds with it, and records them
        as removed, so that re-detecting the coordination on a relaxed
        geometry does not put back what the distances still look like.
        """

        pytest.importorskip('openmm')

        builder = self.editable(tmp_path)
        builder.build_active_site(self.structure_file(), mm_opt=False)

        before = builder.active_site_molecule.number_of_atoms()
        assert 'ASP130' in builder._active_site['residues']

        builder.remove_residue('130')

        assert 'ASP130' not in builder._active_site['residues']
        assert builder.active_site_molecule.number_of_atoms() < before
        assert 'ASP130' not in {
            ligand['residue']
            for ligand in builder.binding_modes['ligands']
        }

        removed = {(record['res_index'], record['action'])
                   for record in builder._request['manual_bonds']}
        assert (2, 'remove') in removed

    def test_a_residue_taken_out_can_be_put_back(self, tmp_path):
        """
        The two edits are requests rather than a membership list, so they
        have to undo one another cleanly. What comes back is the residue, not
        its coordination: the bonds were removed and stay removed.
        """

        pytest.importorskip('openmm')

        builder = self.editable(tmp_path)
        builder.build_active_site(self.structure_file(), mm_opt=False)

        builder.remove_residue('130')
        builder.include_residue('130')

        assert 'ASP130' in builder._active_site['residues']
        assert 'ASP130' not in {
            ligand['residue']
            for ligand in builder.binding_modes['ligands']
        }

    def test_a_metal_cannot_be_left_without_any_ligands(self, tmp_path):
        """
        A metal with nothing bonded to it is not a coordination sphere, and
        the fit downstream has no term to make of it, so the removal that
        would produce one is refused rather than left to fail later.

        Both zincs of the fixture are held by several residues, so the site
        has to be taken down to the last ligand of one of them first -- which
        exercises the removals stacking up as well.
        """

        pytest.importorskip('openmm')

        builder = self.editable(tmp_path)
        builder.build_active_site(self.structure_file(), mm_opt=False)

        # everything on the second zinc except GLU73, which bridges both
        for resid in ('130', '187', '191'):
            builder.remove_residue(resid)

        assert builder._active_site['residues'] == ['GLU73', 'CYS77', 'ASP134']

        with pytest.raises(Exception) as error:
            builder.remove_residue('73')
        assert 'with no ligand at all' in str(error.value)

        # and the refusal changed nothing
        assert builder._active_site['residues'] == ['GLU73', 'CYS77', 'ASP134']

    def test_a_residue_that_is_not_there_cannot_be_removed(self, tmp_path):
        """
        A residue of another structure and one that has already been taken
        out are different mistakes, and each says which it is.
        """

        pytest.importorskip('openmm')

        builder = self.editable(tmp_path)
        builder.build_active_site(self.structure_file(), mm_opt=False)

        with pytest.raises(Exception) as error:
            builder.remove_residue('999')
        assert 'not a residue of this structure' in str(error.value)

        builder.remove_residue('130')

        with pytest.raises(Exception) as error:
            builder.remove_residue('130')
        assert 'is not part of the active site' in str(error.value)

    def test_the_edited_membership_survives_a_relaxation(self, tmp_path):
        """
        update_binding_modes reads the distances again, and the two requests
        are exactly the decisions the distances are not allowed to overrule,
        so they have to reach the dictionary it returns.
        """

        pytest.importorskip('openmm')

        structure = self.pushed_structure(tmp_path / 'far.pdb')

        builder = self.editable(tmp_path)
        builder.build_active_site(structure, mm_opt=False)
        builder.include_residue('130')
        builder.remove_residue('134')

        modes = builder.binding_modes
        assert modes['extra_residues'] == [2]
        assert modes['excluded_residues'] == [3]

        active_site = builder._active_site
        geometry = active_site['molecule'].get_coordinates_in_angstrom()

        updated, _, _ = self.shells().sites.update_binding_modes(builder.enzyme_topology,
                                                                 geometry, active_site, modes)

        assert updated['extra_residues'] == [2]
        assert updated['excluded_residues'] == [3]
        assert util.active_site_residues(updated) == \
            util.active_site_residues(modes)

    def test_an_edit_reaches_the_force_field_without_a_rebuild(self, tmp_path):
        """
        An edit is not a request to be remembered for later. The site on the
        builder describes the correction as soon as it is made, so that going
        straight on to build_forcefield cannot fit the coordination that was
        just corrected away.
        """

        pytest.importorskip('openmm')

        builder = self.editable(tmp_path)
        builder.build_active_site(self.structure_file(), mm_opt=False)

        before = self.coordination(builder.binding_modes)
        assert ('ASP130', 'OD1', 'bidentate') in before

        builder.remove_metal_bond('ASP130', atom='OD1')

        # no build_active_site() call in between
        after = self.coordination(builder.binding_modes)
        assert ('ASP130', 'OD1', 'bidentate') not in after

        metals = set(builder._active_site['metal_indices'])
        matrix = builder._active_site['connectivity_matrix']
        assert sum(matrix[metal].sum() for metal in metals) == sum(
            len(ligand['metals'])
            for ligand in builder.binding_modes['ligands'])

    def test_the_metal_is_named_the_way_the_active_site_labels_it(
            self, tmp_path):
        """
        The labels are the only place the user is shown a metal's number, and
        they index the protonated topology.

        There is one numbering now. The coordination is derived on whatever
        topology it is asked about rather than carried across the
        renumbering, so the number the labels show is the number an edit
        takes, with nothing in between to translate it.
        """

        pytest.importorskip('openmm')

        builder = self.editable(tmp_path)
        builder.build_active_site(self.structure_file(), mm_opt=False)

        labelled = self.metal_label(builder, 0)

        # what the builder reports is what the active site labels, because
        # both are derived on the protonated topology
        assert builder.binding_modes['metals'][0]['index'] == labelled

        builder.add_metal_bond('130', labelled, atom='OD1')

        with pytest.raises(Exception) as error:
            builder.add_metal_bond('130', 9999, atom='OD1')
        assert 'no metal center at atom index 9999' in str(error.value)
        assert 'which is what its labels show' in str(error.value)

    def test_the_stage_advances_with_the_workflow(self, tmp_path):
        """
        The stage is the one thing that says what the builder will do, so it
        has to say where the builder actually is.
        """

        pytest.importorskip('openmm')

        builder = self.editable(tmp_path)
        assert builder.stage == 'empty'

        builder.build_active_site(self.structure_file(), mm_opt=False)
        assert builder.stage == 'active_site'

        n_atoms = builder.active_site_molecule.number_of_atoms()
        builder.build_forcefield(hessian=self.random_hessian(n_atoms),
                                 partial_charges=np.zeros(n_atoms))
        assert builder.stage == 'fitted'

    def test_a_new_input_drops_the_fit_that_no_longer_matches(self, tmp_path):
        """
        A Hessian, a set of charges or a geometry computed after the fit
        describes something the fit does not.

        This used to pass silently: the four steps guarded only on there
        being an active site, so calculate_hessian() after build_forcefield()
        left active_site_forcefield answering with a force field fitted from
        a different Hessian than the one hessian now reports.
        """

        pytest.importorskip('openmm')

        builder = self.fitted(tmp_path)
        assert builder.stage == 'fitted'
        assert builder.active_site_forcefield is not None

        n_atoms = builder.active_site_molecule.number_of_atoms()
        builder.calculate_partial_charges()

        assert builder.stage == 'active_site'
        assert builder.active_site_forcefield is None

        # and the input that was just computed is the one still there, while
        # the Hessian beside it -- the other input -- is untouched
        assert builder.partial_charges is not None
        assert builder.hessian is not None
        assert builder.hessian.shape == (3 * n_atoms, 3 * n_atoms)

    def test_a_refused_operation_names_the_step_to_take(self, tmp_path):
        """
        The message is read by someone who called the wrong thing, so it has
        to name the next step from where they are rather than the one they
        asked for -- which would fail for the same reason.
        """

        pytest.importorskip('openmm')

        builder = self.editable(tmp_path)

        with pytest.raises(Exception) as error:
            builder.create_enzyme_system()
        assert 'no active site' in str(error.value)

        builder.build_active_site(self.structure_file(), mm_opt=False)

        with pytest.raises(Exception) as error:
            builder.create_enzyme_system()
        assert 'no force field yet' in str(error.value)

    def test_a_metal_can_be_named_by_the_label_drawn_on_it(self, tmp_path):
        """
        The number a caller has for a metal is the one show_active_site drew
        on it, and the labels are strings.

        Passing that string straight back has to name the same atom as the
        int does, or the picture and the edit methods disagree about what
        '85' means -- which produced the self-contradicting complaint that
        there was no metal at index 85 on a site holding Zn 85.
        """

        pytest.importorskip('openmm')

        builder = self.editable(tmp_path)
        builder.build_active_site(self.structure_file(), mm_opt=False)

        site = builder._active_site
        drawn = site['labels'][site['metal_indices'][1]]
        assert isinstance(drawn, str)

        # the string off the picture and the int mean the same metal.
        # One derivation is used for both, since binding_modes is worked out
        # afresh on every access and two calls would hand back two objects.
        modes = builder.binding_modes
        assert (self.shells().sites._resolve_metal(modes, drawn)
                is self.shells().sites._resolve_metal(modes, int(drawn)))

        # and it works through the public edit, which is where it is used
        builder.remove_metal_bond('187', drawn)
        assert self.contacts_of(builder.binding_modes, 'HIS187') == []

        # something that is not an index at all is still refused, by name
        with pytest.raises(Exception, match='no metal center at atom index'):
            builder.remove_metal_bond('73', 'NE2')

    def test_the_protonation_and_the_residues_are_fixed_once_fitted(
            self, tmp_path):
        """
        The charges and the parameters of the force field describe the
        protonation and the residues it was fitted to, so changing either
        afterwards would leave it describing a cluster that no longer exists.
        Coordination is the exception, and the error says so.
        """

        pytest.importorskip('openmm')

        builder = self.fitted(tmp_path)

        for call, args in (('update_protonation_state', ('187', 'HID')),
                           ('include_residue', ('130', )), ('remove_residue',
                                                            ('134', ))):
            with pytest.raises(Exception) as error:
                getattr(builder, call)(*args)
            message = str(error.value)
            assert 'force field has already been built' in message
            assert 'add_metal_bond' in message

        # and it is reopened by extracting the site again
        builder.build_active_site(mm_opt=False)
        builder.update_protonation_state('187', 'HID')

        assert builder.active_site_forcefield is None

    def test_a_bond_added_after_the_fit_is_fitted_again(self, tmp_path):
        """
        Rebuilding the site after the fit would throw away the geometry, the
        Hessian and the charges it was made from. The bond edit instead
        rewrites the coordination of the site that is there and fits the
        metal terms again on them, which costs nothing and leaves the
        angles crossing the edited bond right, since the generator derives
        every one of them from the connectivity.
        """

        pytest.importorskip('openmm')

        builder = self.fitted(tmp_path)
        forcefield = builder.active_site_forcefield
        bonds, angles = util.get_metal_keys(forcefield, builder._active_site)

        metal = self.metal_label(builder, 0)
        builder.add_metal_bond('134', metal, atom='OD1')

        refitted = builder.active_site_forcefield
        new_bonds, new_angles = util.get_metal_keys(refitted,
                                                    builder._active_site)

        assert len(new_bonds) == len(bonds) + 1
        assert len(new_angles) > len(angles)
        assert refitted is not forcefield

        # nothing was recomputed
        assert builder.scf_drv is None

        # and it goes back out again
        builder.remove_metal_bond('134', metal=metal, atom='OD1')
        final_bonds, final_angles = util.get_metal_keys(
            builder.active_site_forcefield, builder._active_site)

        assert len(final_bonds) == len(bonds)
        assert len(final_angles) == len(angles)

    def test_a_bond_the_hessian_does_not_cover_is_added_with_a_warning(
            self, tmp_path, capsys):
        """
        A Hessian is restricted to the pairs of the coordination it was
        computed for, so a bond added afterwards has nothing to project and
        is fitted to zero. That is worth saying loudly -- the cure is to
        recompute the Hessian -- but it is not a reason to refuse the edit.
        """

        pytest.importorskip('openmm')

        builder = self.fitted(tmp_path)
        builder.ostream = self.printing()
        capsys.readouterr()

        metal = self.metal_label(builder, 0)
        builder.add_metal_bond('134', metal, atom='OD1')

        metal_site = builder._active_site['metal_indices'][0]
        donor = self.donor_index(builder._active_site, '134', 'OD1')

        printed = capsys.readouterr().out
        assert 'not covered by the Hessian at all' in printed

        # the bond is there, and it is the one the warning is about
        forcefield = builder.active_site_forcefield
        added = tuple(sorted((metal_site, donor)))

        assert added in forcefield.bonds
        assert forcefield.bonds[added]['force_constant'] == 0.0
        assert f'bond {added}' in printed

    def test_the_enzyme_system_is_dropped_when_the_terms_are_refitted(
            self, tmp_path):
        """
        The enzyme system carries the metal terms that were just refitted, so
        holding on to it would hand back a system built from parameters that
        no longer exist.
        """

        pytest.importorskip('openmm')

        builder = self.fitted(tmp_path)
        builder._enzyme_system = 'stale'

        metal = self.metal_label(builder, 0)
        builder.add_metal_bond('134', metal, atom='OD1')

        assert builder.enzyme_system is None

    def test_an_added_bond_can_be_given_an_equilibrium_to_pull_to(
            self, tmp_path):
        """
        A bond is often added by hand because the sidechain is turned the
        wrong way round. Measuring its equilibrium on that geometry would
        only hold it there, so a distance can be named instead and the crude
        pass swings the residue into place before any QM is paid for.
        """

        pytest.importorskip('openmm')

        def relaxed(equilibrium):
            builder = self.editable(tmp_path, name=f'eq{equilibrium}')
            builder.build_active_site(self.structure_file(), mm_opt=False)

            metal = self.metal_label(builder, 0)
            builder.add_metal_bond('134',
                                   metal,
                                   atom='OD1',
                                   equilibrium=equilibrium)

            key = sorted(
                self.shells().sites.manual_bond_keys(builder._active_site,
                                                     builder.enzyme_topology,
                                                     builder.binding_modes))[0]
            pair = [key[0] + 1, key[1] + 1]

            before = builder.active_site_molecule.get_distance_in_angstroms(
                pair)
            builder.mm_optimize_active_site()
            after = builder.active_site_molecule.get_distance_in_angstroms(pair)

            return builder, before, after

        # left alone, the equilibrium is measured on the geometry, so the
        # relaxation keeps the contact roughly where the structure had it
        loose, before, unchanged = relaxed(None)
        assert self.shells().sites.manual_bond_equilibria(loose._active_site,
                                                          loose.enzyme_topology,
                                                          loose.binding_modes) == {}
        assert unchanged == pytest.approx(before, abs=0.2)

        # asked for, it is pulled in instead
        tight, _, pulled = relaxed(2.1)
        assert pulled == pytest.approx(2.1, abs=0.2)
        assert pulled < unchanged - 0.3

        # and it is kept in nanometers, which is what the force field wants
        equilibria = self.shells().sites.manual_bond_equilibria(tight._active_site,
                                                                tight.enzyme_topology,
                                                                tight.binding_modes)
        assert list(equilibria.values()) == pytest.approx([0.21])

    def test_an_equilibrium_survives_a_re_derivation(self, tmp_path):
        """
        The distance is part of the decision the bond records, so asking for
        the coordination again has to give it back -- otherwise the pull
        would quietly vanish the next time anything rebuilt the site.
        """

        pytest.importorskip('openmm')

        builder = self.editable(tmp_path)
        builder.build_active_site(self.structure_file(), mm_opt=False)
        builder.add_metal_bond('134',
                               self.metal_label(builder, 0),
                               atom='OD1',
                               equilibrium=2.1)

        record = next(entry for entry in builder._request['manual_bonds']
                      if entry['action'] == 'add')
        assert record['equilibrium'] == pytest.approx(0.21)

        # deriving again from the request reproduces it
        again = self.derive(builder.enzyme_topology, builder.enzyme_positions,
                            builder._request)
        assert again['manual_bonds'] == builder._request['manual_bonds']

    def test_an_equilibrium_has_to_be_a_positive_distance(self, tmp_path):
        """
        Zero or a negative distance is a typo, and one that would otherwise
        reach the force field as a bond length.
        """

        pytest.importorskip('openmm')

        builder = self.editable(tmp_path)
        builder.build_active_site(self.structure_file(), mm_opt=False)

        with pytest.raises(Exception, match='must be positive'):
            builder.add_metal_bond('134',
                                   self.metal_label(builder, 0),
                                   atom='OD1',
                                   equilibrium=0)

    def test_a_hand_added_bond_is_not_pruned_away_again(self, tmp_path):
        """
        The pruning drops a bridging arm the fit gave no force constant and
        the geometry holds further out than the rest of its residue. A bond
        added by hand looks exactly like that whenever the Hessian does not
        cover it, so without protecting it the fit would take straight back
        out what the edit just put in.

        Both halves are checked on the same site, since a test that only
        showed the bond surviving would pass just as well if the pruning had
        never looked at it.
        """

        pytest.importorskip('openmm')

        builder = self.editable(tmp_path)
        builder.build_active_site(self.structure_file(), mm_opt=False)

        # ASP134 grips one metal through OD2 at 1.95 A; reaching the other
        # through OD1, at 3.96 A, makes it the mu-1,3 bridge with a long weak
        # arm that the pruning exists to drop. Reaching the same metal twice
        # would only make it bidentate, which is a grip and is never touched.
        metal = self.metal_label(builder, 1)
        builder.add_metal_bond('134', metal, atom='OD1')

        active_site = builder._active_site
        protected = self.shells().sites.manual_bond_keys(active_site, builder.enzyme_topology,
                                                         builder.binding_modes)

        assert len(protected) == 1
        added = next(iter(protected))

        # every metal term fits to zero, so only the protection can save it
        n_atoms = active_site['molecule'].number_of_atoms()
        hessian = np.eye(3 * n_atoms) * 0.5

        unprotected = self.fitted_forcefield(active_site,
                                             hessian,
                                             prune_weak_bridge_bonds=True)
        kept = self.fitted_forcefield(active_site,
                                      hessian,
                                      prune_weak_bridge_bonds=True,
                                      protected_bonds=protected)

        assert added not in unprotected.bonds
        assert added in kept.bonds

        # protecting one arm is not protecting the residue: what the pruning
        # took out of the unprotected fit besides this bond it still takes
        assert set(kept.bonds) - set(unprotected.bonds) == {added}

        # and the whole thing again through the builder, which protects on
        # its own
        forcefield = builder.build_forcefield(hessian=hessian)

        assert added in forcefield.bonds

    def test_the_site_stops_claiming_a_bond_the_fit_pruned(self, tmp_path):
        """
        The pruning is the one step that can decide a metal contact is not a
        bond after all, and it decides it on the force field. The active site
        is what show_active_site draws, what the summary counts, what a later
        Hessian walks its pairs out of and what the manager describes a query
        site by, so a site left claiming the pruned bond makes the builder
        answer two different ways about what its own cluster is bonded like.
        """

        pytest.importorskip('openmm')

        builder = self.editable(tmp_path, prune_weak_bridge_bonds=True)
        builder.build_active_site(self.structure_file(), mm_opt=False)

        before = util.connectivity_bonds(
            builder.active_site['connectivity_matrix'])

        # every metal term fits to zero, which is what leaves the fixture's
        # own GLU73 bridge to be decided on its geometry alone
        n_atoms = builder.active_site_molecule.number_of_atoms()
        forcefield = builder.build_forcefield(hessian=np.ones((3 * n_atoms,
                                                               3 * n_atoms)))

        after = util.connectivity_bonds(
            builder.active_site['connectivity_matrix'])
        pruned = set(before) - set(after)

        assert len(pruned) == 1
        assert pruned.isdisjoint(forcefield.bonds)
        assert sorted(forcefield.bonds) == after

        # the covalent bonds are none of the pruning's business, and neither
        # is anything the fit did not remove
        metals = set(builder.active_site['metal_indices'])
        assert all(metals & set(key) for key in pruned)

    # ------------------------------------------------------------------
    # the force field xml
    # ------------------------------------------------------------------

    def enzyme_pair(self, tmp_path):
        """
        Builds one site both ways: the terms injected onto a system, and the
        system the emitted force field builds for itself.

        :return:
            A dictionary with both systems, what create_enzyme_forcefield
            returned, and the pieces they were built from.
        """

        from veloxchem.metalsite import openmmxml
        import openmm.app as mmapp

        topology, positions = self.shells().sites.load_and_prepare_protein(
            self.structure_file())
        binding_modes = self.derive(topology, positions)
        topology, positions, binding_modes = self.protonated(
            topology, positions, binding_modes)
        active_site = self.shells().sites.extract_active_site(
            topology, positions, binding_modes)
        forcefield = self.seeded_forcefield(active_site)
        charges = util.d4_charges(active_site)

        reference, _ = self.shells().enzyme.create_enzyme_system(
            topology, active_site, forcefield, charges)
        result = self.shells().enzyme.create_enzyme_forcefield(
            topology, positions, active_site, forcefield, charges)

        path = tmp_path / openmmxml.SITE_XML_FILE
        path.write_text(result['xml'])
        built = mmapp.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml',
                                 str(path)).createSystem(
                                     result['topology'],
                                     nonbondedMethod=mmapp.NoCutoff)

        return {
            'reference': reference,
            'built': built,
            'result': result,
            'path': path,
            'topology': topology,
            'active_site': active_site,
            'forcefield': forcefield,
        }

    def test_the_xml_builds_the_system_the_injection_builds(self, tmp_path):
        """
        The force field is the same parameters written down rather than a
        second set of them, so the system it builds has to be the system
        create_enzyme_system injects the terms into -- every bond, angle,
        torsion and particle identical.

        The exceptions are the one place the two differ, and deliberately:
        the metal-ligand contacts are real bonds here, so OpenMM excludes
        the 1-2 and 1-3 pairs they make and scales the 1-4 ones, which the
        injected system never did.
        """

        pytest.importorskip('openmm')
        pytest.importorskip('pdbfixer')

        import openmm as mm
        import openmm.unit as mmunit
        from veloxchem.metalsite import openmmxml

        pair = self.enzyme_pair(tmp_path)
        remap = pair['result']['restructured']['atom_index_map']

        def force(system, kind):
            return next(one for one in system.getForces()
                        if isinstance(one, kind))

        def bonds(system, index_of):
            handle = force(system, mm.HarmonicBondForce)
            return sorted((frozenset((index_of(p[0]), index_of(p[1]))),
                           round(p[2].value_in_unit(mmunit.nanometer), 9),
                           round(
                               p[3].value_in_unit(mmunit.kilojoule_per_mole /
                                                  mmunit.nanometer**2), 5))
                          for p in (handle.getBondParameters(i)
                                    for i in range(handle.getNumBonds())))

        def angles(system, index_of):
            handle = force(system, mm.HarmonicAngleForce)
            out = []
            for i in range(handle.getNumAngles()):
                p = handle.getAngleParameters(i)
                triple = (index_of(p[0]), index_of(p[1]), index_of(p[2]))
                out.append((min(triple, tuple(reversed(triple))),
                            round(p[3].value_in_unit(mmunit.radian), 9),
                            round(
                                p[4].value_in_unit(mmunit.kilojoule_per_mole /
                                                   mmunit.radian**2), 5)))
            return sorted(out)

        def torsions(system, index_of):
            handle = force(system, mm.PeriodicTorsionForce)
            out = []
            for i in range(handle.getNumTorsions()):
                p = handle.getTorsionParameters(i)
                barrier = round(p[6].value_in_unit(mmunit.kilojoule_per_mole),
                                6)
                # a zeroed torsion is not a term, and the wildcard ones
                # written across a metal bond are zeroed rather than removed
                if barrier == 0.0:
                    continue
                quad = tuple(index_of(index) for index in p[:4])
                out.append((min(quad, tuple(reversed(quad))), p[4],
                            round(p[5].value_in_unit(mmunit.radian),
                                  9), barrier))
            return sorted(out)

        def particles(system, index_of):
            handle = force(system, mm.NonbondedForce)
            return {
                index_of(i):
                (round(p[0].value_in_unit(mmunit.elementary_charge),
                       8), round(p[1].value_in_unit(mmunit.nanometer), 9),
                 round(p[2].value_in_unit(mmunit.kilojoule_per_mole), 9))
                for i, p in ((i, handle.getParticleParameters(i))
                             for i in range(handle.getNumParticles()))
            }

        def exceptions(system, index_of):
            handle = force(system, mm.NonbondedForce)
            out = {}
            for i in range(handle.getNumExceptions()):
                p = handle.getExceptionParameters(i)
                out[frozenset(
                    (index_of(p[0]),
                     index_of(p[1])))] = (round(
                         p[2].value_in_unit(mmunit.elementary_charge**2),
                         9), round(p[3].value_in_unit(mmunit.nanometer), 9),
                                          round(
                                              p[4].value_in_unit(
                                                  mmunit.kilojoule_per_mole),
                                              9))
            return out

        reference, built = pair['reference'], pair['built']
        moved = remap.__getitem__
        kept = int

        assert bonds(reference, moved) == bonds(built, kept)
        assert angles(reference, moved) == angles(built, kept)
        assert torsions(reference, moved) == torsions(built, kept)
        assert particles(reference, moved) == particles(built, kept)

        before = exceptions(reference, moved)
        after = exceptions(built, kept)

        assert not [pair_ for pair_ in before if pair_ not in after]
        assert not [
            pair_ for pair_ in before
            if pair_ in after and before[pair_] != after[pair_]
        ]

        # everything the metal bonds put within three bonds of each other,
        # and nothing else
        site_index = pair['result']['restructured']['site_atom_index']
        metal_bonds = [(site_index[key[0]], site_index[key[1]])
                       for key in openmmxml.metal_term_keys(
                           pair['forcefield'], pair['active_site'])[0]]
        bonded = {}
        for first, second in pair['result']['topology'].bonds():
            bonded.setdefault(first.index, set()).add(second.index)
            bonded.setdefault(second.index, set()).add(first.index)
        without = {
            index: set(neighbours)
            for index, neighbours in bonded.items()
        }
        for first, second in metal_bonds:
            without[first].discard(second)
            without[second].discard(first)

        new = [pair_ for pair_ in after if pair_ not in before]
        assert new
        for pair_ in new:
            first, second = tuple(pair_)
            assert util._bond_separation(bonded, first, second, 3) <= 3
            assert util._bond_separation(without, first, second, 3) > 3

    def test_the_metal_contacts_become_bonds_of_the_topology(self, tmp_path):
        """
        The whole point of restructuring: a metal-ligand contact is a real
        bond of the topology, which is what makes OpenMM write a term and an
        exclusion across it. The active site is moved into one residue and
        the backbone stays where it was, so the chain keeps its residues and
        nothing downstream sees a gap.
        """

        pytest.importorskip('openmm')
        pytest.importorskip('pdbfixer')

        from veloxchem.metalsite import openmmxml

        pair = self.enzyme_pair(tmp_path)
        original, restructured = pair['topology'], pair['result']['topology']
        site = pair['result']['restructured']['site_residue']

        metal_bonds, _, _ = openmmxml.metal_term_keys(pair['forcefield'],
                                                      pair['active_site'])

        assert restructured.getNumAtoms() == original.getNumAtoms()
        assert (len(list(restructured.bonds())) == len(list(original.bonds())) +
                len(metal_bonds))

        caps = set(pair['active_site']['cap_indices'])
        moved = len(pair['active_site']['atom_map']) - len(caps)
        assert len(list(site.atoms())) == moved

        # the atoms of the merged residue are the site's, under names no
        # two of which are the same and none of which a PDB cannot hold
        names = [atom.name for atom in site.atoms()]
        assert len(set(names)) == len(names)
        assert all(len(name) <= 4 for name in names)

        # every coordinating residue still exists, with its backbone
        stubs = pair['result']['restructured']['stub_residues']
        assert {f'{stub.name}{stub.id}'
                for stub in stubs.values()
                } == set(pair['active_site']['residues'])
        for stub in stubs.values():
            assert {atom.name for atom in stub.atoms()} >= {'N', 'CA', 'C', 'O'}

    def test_the_stub_templates_are_told_apart_by_the_matcher(self, tmp_path):
        """
        OpenMM matches a residue to a template on the bond graph alone, and
        the backbone stubs left behind are all isomorphic, so two of them
        carrying different charges is a hard error. The XML answers that
        itself, with a template matcher registered from an
        InitializationScript, rather than making its caller pass a
        residueTemplates mapping into every createSystem and every Modeller
        call.
        """

        pytest.importorskip('openmm')
        pytest.importorskip('pdbfixer')

        import xml.etree.ElementTree as ET
        import openmm.app as mmapp

        pair = self.enzyme_pair(tmp_path)

        tree = ET.parse(pair['path'])
        root = tree.getroot()
        root.remove(root.find('InitializationScript'))
        without = tmp_path / 'no_matcher.xml'
        ET.ElementTree(root).write(without)

        with pytest.raises(Exception, match='Multiple non-identical'):
            mmapp.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml',
                             str(without)).createSystem(
                                 pair['result']['topology'],
                                 nonbondedMethod=mmapp.NoCutoff)

    def test_the_restructured_topology_survives_a_pdb_round_trip(
            self, tmp_path):
        """
        The file is only half of what a simulation needs; the topology whose
        chain and residue ids the templates are keyed on is the other half.
        Writing the site as a residue of its own is what makes that topology
        keep its metal bonds through a PDB: the reader refuses a CONECT
        record between a metal and a residue it thinks it already knows the
        bonding of, and takes one naming a residue it does not.
        """

        pytest.importorskip('openmm')
        pytest.importorskip('pdbfixer')

        import openmm.app as mmapp
        import openmm.unit as mmunit

        pair = self.enzyme_pair(tmp_path)
        restructured = pair['result']['topology']

        path = tmp_path / 'metal_site.pdb'
        with path.open('w') as handle:
            mmapp.PDBFile.writeFile(restructured,
                                    pair['result']['positions'] *
                                    mmunit.angstrom,
                                    handle,
                                    keepIds=True)

        read_back = mmapp.PDBFile(str(path)).topology

        def bonds_of(topology):
            atoms = list(topology.atoms())
            return {
                frozenset((atoms[first.index].name, atoms[second.index].name))
                for first, second in topology.bonds()
                if first.residue.name == 'MSI' or second.residue.name == 'MSI'
            }

        assert read_back.getNumAtoms() == restructured.getNumAtoms()
        assert bonds_of(read_back) == bonds_of(restructured)

        system = mmapp.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml',
                                  str(pair['path'])).createSystem(
                                      read_back, nonbondedMethod=mmapp.NoCutoff)
        assert system.getNumParticles() == restructured.getNumAtoms()

    def test_the_wildcard_torsions_across_a_metal_bond_are_zeroed(
            self, tmp_path):
        """
        Once a metal-ligand contact is a real bond, the protein force field's
        wildcard torsions reach across it -- and two of the three it writes
        on this site hold the metal in the plane of an imidazole, which the
        fitted improper on the coordinating nitrogen already does, at four
        times the barrier and with nothing having fitted it. Nothing about
        the metal should come from a wildcard.
        """

        pytest.importorskip('openmm')
        pytest.importorskip('pdbfixer')

        import openmm as mm
        import openmm.unit as mmunit
        from veloxchem.metalsite import openmmxml

        pair = self.enzyme_pair(tmp_path)
        site_index = pair['result']['restructured']['site_atom_index']
        metal_bonds = {
            frozenset((site_index[key[0]], site_index[key[1]]))
            for key in openmmxml.metal_term_keys(pair['forcefield'],
                                                 pair['active_site'])[0]
        }

        torsions = next(one for one in pair['built'].getForces()
                        if isinstance(one, mm.PeriodicTorsionForce))

        crossing = []
        for index in range(torsions.getNumTorsions()):
            parameters = torsions.getTorsionParameters(index)
            quad = parameters[:4]
            if any(
                    frozenset(step) in metal_bonds
                    for step in zip(quad, quad[1:])):
                crossing.append(parameters)

        assert crossing
        kept = [
            parameters for parameters in crossing
            if parameters[6].value_in_unit(mmunit.kilojoule_per_mole) != 0.0
        ]

        # what is left with a barrier is the fit's own and nothing else: an
        # improper whose written order steps along a metal bond
        _, _, impropers = openmmxml.metal_term_keys(pair['forcefield'],
                                                    pair['active_site'])
        ours = {(site_index[key[1]], site_index[key[2]], site_index[key[0]],
                 site_index[key[3]])
                for key in impropers}
        assert {tuple(parameters[:4]) for parameters in kept} <= ours
        assert len(kept) < len(crossing)

    def test_the_enzyme_forcefield_waits_for_the_force_field(self, tmp_path):
        """
        Writing the force field triggers nothing expensive, and there is
        nothing to write until the metal terms have been fitted.
        """

        pytest.importorskip('openmm')

        from veloxchem.metalsite import openmmxml

        n_atoms = self.prepared()['active_site']['molecule'].number_of_atoms()
        builder = self.builder(tmp_path)

        with pytest.raises(Exception, match='build_active_site'):
            builder.create_enzyme_forcefield()

        builder.build_active_site(self.structure_file())

        with pytest.raises(Exception, match='build_forcefield'):
            builder.create_enzyme_forcefield()

        builder.build_forcefield(hessian=self.random_hessian(n_atoms),
                                 partial_charges=np.zeros(n_atoms))

        result = builder.create_enzyme_forcefield()

        assert builder.stage == 'enzyme'
        assert builder.enzyme_forcefield is result
        assert result['topology'].getNumAtoms() == len(result['positions'])

        folder = Path(builder.folder)
        assert (folder / openmmxml.SITE_XML_FILE).is_file()
        assert (folder / openmmxml.SITE_TOPOLOGY_FILE).is_file()

        # a bond edit invalidates the fit, and with it what was written from
        # it, the same way it invalidates an enzyme system
        builder.add_metal_bond('134', self.metal_label(builder, 0), atom='OD1')
        assert builder.enzyme_forcefield is None

    def test_the_torsions_across_a_metal_bond_are_zeroed_by_the_script(
            self, tmp_path):
        """
        The wildcard torsions the protein force field writes across a metal
        bond are zeroed by the file's own script, unconditionally: nothing
        about the metal should come from a wildcard.
        """

        pytest.importorskip('openmm')

        n_atoms = self.prepared()['active_site']['molecule'].number_of_atoms()
        builder = self.builder(tmp_path)
        builder.build_active_site(self.structure_file())
        builder.build_forcefield(hessian=self.random_hessian(n_atoms),
                                 partial_charges=np.zeros(n_atoms))

        xml = builder.create_enzyme_forcefield()['xml']

        assert 'PeriodicTorsionForce' in xml
        assert 'setTorsionParameters' in xml


class TestShellDecorators:
    """
    The two decorators every method of a phase class carries. One rank is
    enough for what is checked here: the guard against a collective inside a
    master-only body, and that a failure on the master is what the caller
    sees.
    """

    @staticmethod
    def shell():

        from veloxchem.metalsite.util import Shell, on_master, collective

        class Phase(Shell):

            @on_master
            def cheap(self, value):
                return value * 2

            @on_master
            def failing(self):
                raise ValueError('the master failed')

            @on_master
            def nested(self, value):
                # a master-only method calling another runs it inline, in
                # the section the outer one opened
                assert self.master_section.depth == 1
                return self.cheap(value) + 1

            @collective
            def expensive(self):
                return 'entered'

            @on_master
            def wrong(self):
                return self.expensive()

        return Phase(ostream=OutputStream(None))

    def test_a_master_only_step_hands_back_what_the_body_returned(self):

        shell = self.shell()

        assert shell.cheap(21) == 42
        assert shell.nested(20) == 41
        # the section is one object for every shell of the process
        assert shell.master_section is self.shell().master_section
        assert not shell.master_section.active

    def test_a_failure_on_the_master_is_raised_to_the_caller(self):

        shell = self.shell()

        with pytest.raises(ValueError, match='the master failed'):
            shell.failing()

        # and the section is closed again, so the next call is not inline
        assert shell.master_section.depth == 0
        assert shell.cheap(1) == 2

    def test_a_collective_may_not_run_inside_a_master_only_step(self):

        from veloxchem.errorhandler import VeloxChemError

        shell = self.shell()

        assert shell.expensive() == 'entered'

        with pytest.raises(VeloxChemError,
                           match='expensive is a collective step'):
            shell.wrong()

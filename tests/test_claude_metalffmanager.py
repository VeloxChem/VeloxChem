from pathlib import Path
from types import SimpleNamespace
import numpy as np
import random
import tempfile
import pytest

from veloxchem.molecule import Molecule
from veloxchem.outputstream import OutputStream
from veloxchem.metalsite.metalffmanager import MetalForceFieldManager
from veloxchem.metalsite.metalsiteffbuilder import MetalSiteForceFieldBuilder
from veloxchem.metalsite import util
from veloxchem.metalsite.matching import SiteMatcher
from veloxchem.metalsite.templates import TemplateLoader
from veloxchem.metalsite.shoehorn import Shoehorner


class TestMetalForceFieldManager:

    @staticmethod
    def shells(ostream=None, **settings):
        """
        The phase classes, silent unless a stream is given, on the builder's
        default settings unless a test names one. The matcher, the loader
        and the shoehorner carry no setting; the sites and qm shells are
        made the way the interface makes them, through its factories, so
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
                               matcher=SiteMatcher(ostream=ostream),
                               loader=TemplateLoader(ostream=ostream),
                               shoehorner=Shoehorner(ostream=ostream))

    def seeded_forcefield(self,
                          active_site,
                          partial_charges=None,
                          metal_bond_equilibria=None,
                          **settings):
        """
        The seeded force field of a site, on D4 charges unless given.

        :param active_site:
            The active site.
        :param partial_charges:
            The charges, D4 when None.
        :param metal_bond_equilibria:
            The per-call element table, as build_forcefield takes it.
        :param settings:
            Settings of the builder the seeded build reads.

        :return:
            The force field generator.
        """

        if partial_charges is None:
            partial_charges = util.d4_charges(active_site)

        return self.shells(**settings).sites.build_forcefield(
            active_site,
            partial_charges,
            metal_bond_equilibria=metal_bond_equilibria)

    def fitted_forcefield(self,
                          active_site,
                          hessian,
                          partial_charges=None,
                          **settings):
        """
        A seeded force field fitted against a Hessian, with the weak bridge
        pruning off unless a test asks for it.

        :param active_site:
            The active site.
        :param hessian:
            The Hessian.
        :param partial_charges:
            The charges, D4 when None.
        :param settings:
            Settings of the builder either pass reads.

        :return:
            The force field generator.
        """

        shells = self.shells(**settings)

        forcefield = shells.sites.build_forcefield(
            active_site,
            util.d4_charges(active_site)
            if partial_charges is None else partial_charges)

        return shells.qm.fit_forcefield(active_site, forcefield, hessian)

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

    @staticmethod
    def structure_file():

        here = Path(__file__).parent

        return str(here / 'data' / 'zn_site_ECDDHH.pdb')

    @staticmethod
    def random_hessian(n_atoms, seed=0):
        """
        A symmetric matrix of the right shape. The force constants it produces
        are meaningless; these tests check that parameters end up on the right
        atoms, not what the numbers are.
        """

        size = 3 * n_atoms
        rng = np.random.default_rng(seed)
        amat = rng.normal(0.0, 0.05, (size, size))

        return amat + amat.T + np.eye(size) * 0.5

    def prepared(self,
                 structure=None,
                 metal_bond_cutoff=None,
                 extra_residues=None):
        """
        Runs the structural half of the builder pipeline, which needs no QM.

        The protein repair is skipped, exactly as a query builder skips it in
        these tests, so a template and a query see the same truncation.

        :param metal_bond_cutoff:
            The distance a donor atom becomes a ligand within, or None for the
            builder's own. It is what decides denticity, so it is what tells a
            template gripping a metal with one oxygen apart from one gripping
            it with two.
        :param extra_residues:
            Residue indices kept in the cluster whatever their distance, the
            way include_residue keeps one. A residue in here that coordinates
            nothing is what a template loaded from the site is expected to
            discard again.
        """

        if structure is None:
            structure = TestMetalForceFieldManager.structure_file()

        settings = {'prepare_protein': False}
        if metal_bond_cutoff is not None:
            settings['metal_bond_cutoff'] = metal_bond_cutoff
        sites = self.shells(**settings).sites

        topology, positions = sites.load_and_prepare_protein(structure)
        request = util.empty_request()
        request['extra_residues'] = list(extra_residues or [])
        binding_modes = sites.derive_binding_modes(topology, positions, request)
        # derived again on the protonated structure rather than remapped
        topology, positions, variants, _ = sites.protonate(
            topology, positions, binding_modes, None)
        request['variants'] = variants
        binding_modes = sites.derive_binding_modes(topology, positions, request)
        active_site = sites.extract_active_site(topology, positions,
                                                binding_modes)

        return {
            'topology': topology,
            'positions': positions,
            'binding_modes': binding_modes,
            'active_site': active_site,
        }

    def template_folder(self,
                        folder,
                        structure=None,
                        seed=0,
                        mm_copy=False,
                        relax=False,
                        charges=True,
                        metal_bond_cutoff=None,
                        extra_residues=None):
        """
        Writes the two files a template is built from.

        :param folder:
            Where to write them.
        :param structure:
            The structure to build the template from.
        :param seed:
            The seed of the Hessian, so that two templates can differ.
        :param mm_copy:
            Whether to also write the MM geometry, as a copy of the one the
            template carries. That is what a run that never optimized leaves
            behind.
        :param relax:
            Whether to store the geometry relaxed on the crude force field,
            with the literature metal bonds the manager uses for its fallback,
            instead of the one the structure file holds.
        :param charges:
            Whether the force field carries charges at all. A run without RESP
            leaves every atom at zero.
        :param metal_bond_cutoff:
            The cutoff the coordination is detected with, which is what makes
            a template monodentate or bidentate.
        :param extra_residues:
            Residue indices kept in the cluster whatever their distance.

        :return:
            The tuple of the folder and the force field written into it.
        """

        prepared = self.prepared(structure, metal_bond_cutoff, extra_residues)
        active_site = prepared['active_site']
        n_atoms = active_site['molecule'].number_of_atoms()

        if relax:
            # the literature bonds are what the manager seeds its own
            # comparison relaxation with
            seeded = self.seeded_forcefield(
                active_site,
                metal_bond_equilibria=util.LITERATURE_METAL_BONDS)
            geometry = self.shells().sites.mm_optimize_active_site(
                active_site, seeded)
        else:
            geometry = active_site['molecule']

        partial_charges = np.linspace(-0.4, 0.4, n_atoms) if charges else None
        forcefield = self.fitted_forcefield(active_site,
                                            self.random_hessian(n_atoms, seed),
                                            partial_charges=partial_charges)

        folder = Path(folder)
        folder.mkdir(parents=True, exist_ok=True)
        util.save_forcefield(folder / util.FORCEFIELD_FILE, forcefield)
        geometry.write_xyz_file(str(folder / util.GEOMETRY_FILE))

        if mm_copy:
            geometry.write_xyz_file(str(folder / util.MM_GEOMETRY_FILE))

        return folder, forcefield

    @staticmethod
    def metal_bonds(forcefield, metal_indices):
        """
        The metal bonds of a force field, as unordered pairs.
        """

        metals = set(metal_indices)

        return {frozenset(key) for key in forcefield.bonds if metals & set(key)}

    @staticmethod
    def mapped_metal_bonds(template, mapping):
        """
        The metal bonds of a template, on the atoms of the site it matched.
        """

        metals = set(template['metal_indices'])

        return {
            frozenset((mapping[first], mapping[second]))
            for first, second in template['forcefield'].bonds
            if metals & {first, second}
        }

    @staticmethod
    def manager():

        manager = MetalForceFieldManager()
        manager.ostream.mute()

        return manager

    @staticmethod
    def active_site_builder(structure=None, metal_bond_cutoff=None,
                            **settings):
        """
        A MetalSiteForceFieldBuilder with an active site already built on it,
        for feeding to compare_active_site.

        The pdbfixer repair is skipped: the fixtures here are isolated
        residues with no protein to repair, exactly as the manager used to
        skip it under do_prepare_protein=False. The extracted site itself is
        left unrelaxed (mm_opt=False on build_active_site), matching what
        prepared() has always done for a query; compare_active_site's own
        mm_opt argument, which relaxes the site for the comparison rather
        than for the builder, is a separate decision made by the caller.

        :param structure:
            The structure to build the active site from. Defaults to
            structure_file().
        :param metal_bond_cutoff:
            The cutoff the coordination is detected with, which is what
            decides denticity of the query.
        :param settings:
            Further settings to assign on the builder before building.

        :return:
            The builder, with build_active_site already called on it.
        """

        if structure is None:
            structure = TestMetalForceFieldManager.structure_file()

        builder = MetalSiteForceFieldBuilder(ostream=OutputStream(None))
        builder.prepare_protein = False
        builder.do_resp = False
        # every step writes its result the moment it has it, wherever
        # folder points; left at the default that would be a
        # metal_site_<timestamp> folder in the repository root, so this
        # points it at the system temp directory instead
        builder.folder = tempfile.mkdtemp(prefix='metal_site_')
        if metal_bond_cutoff is not None:
            builder.metal_bond_cutoff = metal_bond_cutoff

        for key, value in settings.items():
            setattr(builder, key, value)

        builder.build_active_site(structure, mm_opt=False)

        return builder

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
    def permuted_structure(source, target, order):
        """
        Rewrites a structure with its coordinating residues in another order.

        The active site indices of a force field follow the order the residues
        appear in the file, so this is what tells a real mapping apart from the
        identity.

        :param source:
            The structure to permute.
        :param target:
            Where to write the permuted structure.
        :param order:
            The new order of the residues that are not metal centers.

        :return:
            The path written.
        """

        import openmm.app as mmapp
        import openmm.unit as mmunit

        pdb = mmapp.PDBFile(str(source))
        positions = pdb.getPositions(asNumpy=True).value_in_unit(
            mmunit.angstrom)
        residues = list(pdb.topology.residues())

        metals = [
            residue for residue in residues
            if any(atom.element is not None and atom.element.symbol == 'Zn'
                   for atom in residue.atoms())
        ]
        sidechains = [residue for residue in residues if residue not in metals]

        topology = mmapp.Topology()
        coordinates = []
        chains = {}

        for residue in [sidechains[index] for index in order] + metals:
            chain = chains.get(residue.chain.id)
            if chain is None:
                chain = topology.addChain(residue.chain.id)
                chains[residue.chain.id] = chain
            new_residue = topology.addResidue(residue.name, chain, residue.id)
            for atom in residue.atoms():
                topology.addAtom(atom.name, atom.element, new_residue)
                coordinates.append(positions[atom.index])

        target = Path(target)
        with target.open('w') as fh:
            mmapp.PDBFile.writeFile(topology,
                                    np.array(coordinates) * mmunit.angstrom,
                                    fh,
                                    keepIds=True)

        return str(target)

    def test_a_template_is_two_files(self, tmp_path):
        """
        A template is the force field and the geometry it was fitted on.

        Everything else a match needs is in the force field already: its bonds
        are the bonding topology of the active site, metal-ligand bonds
        included, and its zero charges are the capping hydrogens.
        """

        pytest.importorskip('openmm')

        folder, forcefield = self.template_folder(tmp_path / 'site_a')
        manager = self.manager()

        manager.load_template_from_folder(folder)
        template = manager.templates['site_a']

        assert list(manager.templates) == ['site_a']
        assert template['geometry_kind'] == 'qm_opt'

        assert template['fine_topology'].number_of_nodes() == len(
            forcefield.atoms)
        assert template['fine_topology'].number_of_edges() == len(
            forcefield.bonds)

        # one cap per truncated residue, and both zincs found
        assert len(template['cap_indices']) == 6
        assert len(template['metal_indices']) == 2

        bonds, angles = self.shells().matcher.metal_keys(template)
        assert bonds and angles

    def test_the_name_defaults_to_the_folder(self, tmp_path):
        """
        A folder is what a run leaves behind, so its name is what a template
        is called unless something better is given.
        """

        pytest.importorskip('openmm')

        folder, _ = self.template_folder(tmp_path / 'zn_ECDDHH')
        manager = self.manager()

        manager.load_template_from_folder(folder)
        manager.load_template_from_folder(folder, name='second')

        assert sorted(manager.templates) == ['second', 'zn_ECDDHH']

    def test_an_unoptimized_geometry_needs_the_fallback(self, tmp_path):
        """
        save_results writes the active site under the name of the optimized
        geometry whatever the run did, so a file identical to the MM geometry
        is the MM geometry and no QM optimization ever ran. Taking that for a
        fitted geometry would put a crude structure in a template silently.
        """

        pytest.importorskip('openmm')

        folder, _ = self.template_folder(tmp_path / 'crude', mm_copy=True)
        manager = self.manager()

        with pytest.raises(Exception, match='fallback'):
            manager.load_template_from_folder(folder)

        manager.load_template_from_folder(folder, fallback='mm_opt')

        assert manager.templates['crude']['geometry_kind'] == 'mm_opt'

    def test_a_geometry_that_moved_is_an_optimized_one(self, tmp_path):
        """
        The two files differ once an optimization has actually run, which is
        what makes the kind more than the file name.
        """

        pytest.importorskip('openmm')

        folder, _ = self.template_folder(tmp_path / 'optimized', mm_copy=True)

        # stand in for an optimization having moved the geometry
        path = folder / util.MM_GEOMETRY_FILE
        molecule = Molecule.read_xyz_file(str(path))
        coords = molecule.get_coordinates_in_angstrom()
        coords[0] += 0.2
        Molecule(molecule.get_labels(), coords,
                 'angstrom').write_xyz_file(str(path))

        manager = self.manager()
        manager.load_template_from_folder(folder)

        assert manager.templates['optimized']['geometry_kind'] == 'qm_opt'

    def test_a_missing_file_is_reported(self, tmp_path):
        """
        A folder that holds neither of the two files says which one it is
        missing, rather than failing somewhere later.
        """

        pytest.importorskip('openmm')

        empty = tmp_path / 'empty'
        empty.mkdir()
        manager = self.manager()

        with pytest.raises(Exception, match=util.FORCEFIELD_FILE):
            manager.load_template_from_folder(empty)

        with pytest.raises(Exception, match='not found'):
            manager.load_template_from_folder(tmp_path / 'nowhere')

        folder, _ = self.template_folder(tmp_path / 'no_geometry')
        (folder / util.GEOMETRY_FILE).unlink()

        with pytest.raises(Exception, match='no geometry'):
            manager.load_template_from_folder(folder)

    def test_load_templates_from_folders_loads_every_folder(self, tmp_path):
        """
        Loading several folders at once loads each of them under the name
        load_template_from_folder would give it on its own, so a set of
        templates can be assembled in one call rather than a loop written by
        hand.
        """

        pytest.importorskip('openmm')

        first, _ = self.template_folder(tmp_path / 'first', seed=0)
        second, _ = self.template_folder(tmp_path / 'second', seed=1)
        third, _ = self.template_folder(tmp_path / 'third', seed=2)

        manager = self.manager()
        manager.load_templates_from_folders([first, second, third])

        assert sorted(manager.templates) == ['first', 'second', 'third']
        for name, folder in (('first', first), ('second', second),
                             ('third', third)):
            assert manager.templates[name]['folder'] == str(folder)

    def test_a_residue_coordinating_no_metal_is_discarded(self, tmp_path):
        """
        A template exists for the coordination it was fitted for, and a
        residue reaching none of the metals is no part of one: it is a
        neighbour the truncation happened to keep. Nothing is ever
        transferred from it, while every structure compared against the
        template would have to hold one like it or map onto nothing, so
        loading drops it and renumbers what is left.
        """

        pytest.importorskip('openmm')

        # ASP130 reaches its zinc at 2.50 A, so a tighter cutoff leaves it
        # coordinating nothing at all; naming it an extra residue is what
        # keeps it in the cluster regardless, the way include_residue does
        folder, forcefield = self.template_folder(tmp_path / 'stray',
                                                  metal_bond_cutoff=2.4,
                                                  extra_residues=[2])
        manager = self.manager()

        manager.load_template_from_folder(folder)
        template = manager.templates['stray']

        # the aspartate sidechain cut at CB, its capping hydrogen included
        n_atoms = len(forcefield.atoms) - 7
        assert template['molecule'].number_of_atoms() == n_atoms
        assert len(template['cap_indices']) == 5
        assert len(template['metal_indices']) == 2

        coarse = template['coarse_topology']
        residues = self.shells().matcher.residue_nodes(coarse)
        assert len(residues) == 5
        assert all(coarse.degree(node) > 0 for node in residues)

        # every atom of the template is indexed the one way: the geometry,
        # the force field and the marked atoms are renumbered together
        assert list(template['molecule'].get_labels()) == (
            util._forcefield_elements(template['forcefield']))
        assert len(template['charges']) == n_atoms
        assert all(index < n_atoms
                   for key in template['forcefield'].bonds
                   for index in key)
        assert all(index < n_atoms for index in template['cap_indices'] +
                   template['beta_carbon_indices'])

    def test_discarding_a_residue_keeps_the_charge_whole(self, tmp_path):
        """
        The atoms that go take their charge with them, which leaves what
        stays off a whole number of electrons. The difference is put back on
        the remaining atoms, and the force field carries what the template
        reports: one array means the charges here as everywhere else.
        """

        pytest.importorskip('openmm')

        folder, _ = self.template_folder(tmp_path / 'stray',
                                         metal_bond_cutoff=2.4,
                                         extra_residues=[2])
        manager = self.manager()

        manager.load_template_from_folder(folder)
        template = manager.templates['stray']
        charges = template['charges']

        total = float(np.sum(charges))
        assert total == pytest.approx(round(total), abs=1.0e-9)

        # the caps are emptied by redistribute_cap_charges and stay empty
        assert all(charges[index] == 0.0 for index in template['cap_indices'])

        forcefield = template['forcefield']
        assert np.allclose(
            [atom['charge'] for atom in forcefield.atoms.values()], charges)
        assert np.allclose(forcefield.partial_charges, charges)

    def test_the_compensating_charge_misses_the_capping_hydrogens(self):
        """
        The charge that went missing is added in equal parts to every atom
        that stays, except the capping hydrogens: a cap stands for a bond to
        a protein the site does not hold, and charge put on one is taken off
        again by the next redistribution.
        """

        # -0.4 e in total, which is 0.4 e from the nearest whole number and
        # 0.2 e for each of the two atoms that may carry it
        charges, shift = self.shells().loader.compensate_charges([0.3, 0.0, -0.7, 0.0],
                                                                 [1, 3])

        assert shift == pytest.approx(0.2)
        assert charges == pytest.approx([0.5, 0.0, -0.5, 0.0])
        assert float(np.sum(charges)) == pytest.approx(0.0)

    def test_comparing_needs_a_template(self):
        """
        Comparing against nothing is a mistake, not an empty answer. This is
        checked before anything about the active site, so no builder needs to
        be built to see it.
        """

        manager = self.manager()

        with pytest.raises(Exception, match='no templates'):
            manager.compare_active_site()

    def test_the_residue_order_does_not_matter(self, tmp_path):
        """
        The atom indices of a force field mean nothing outside the run that
        produced it: they follow the order the residues happen to appear in
        the structure file. Matching solves the mapping, so the same site
        written in another order is still the same site.

        This is the test that says the transfer is right rather than merely
        consistent.
        """

        pytest.importorskip('openmm')

        folder, _ = self.template_folder(tmp_path / 'site')
        permuted = self.permuted_structure(self.structure_file(),
                                           tmp_path / 'permuted.pdb',
                                           [3, 5, 0, 4, 1, 2])

        manager = self.manager()
        manager.load_template_from_folder(folder)
        template = manager.templates['site']

        manager.compare_active_site(self.active_site_builder(permuted),
                                    mm_opt=False)
        decision = manager._select_template(None)

        assert decision['score'] == pytest.approx(0.0, abs=1.0e-3)

        mapping = decision['entry']['mapping']
        assert any(key != value for key, value in mapping.items())

        forcefield = manager.build_ff_from_template()

        # every mapped atom keeps its element and its charge
        labels = manager._comparison['active_site']['molecule'].get_labels()
        template_labels = template['molecule'].get_labels()
        charges = forcefield.partial_charges

        for template_index, site_index in mapping.items():
            assert labels[site_index] == template_labels[template_index]
            assert charges[site_index] == pytest.approx(
                template['charges'][template_index])

        # and the metal bonds land on metal atoms of the permuted site
        metals = set(manager._comparison['active_site']['metal_indices'])
        bonds, _ = self.shells().matcher.metal_keys(template)

        for key in bonds:
            mapped = tuple(mapping[index] for index in key)
            assert metals & set(mapped)

    def test_the_hydrogens_are_left_out_of_both_measures(self, tmp_path):
        """
        The hydrogen placement noise reaches the internal coordinates as much
        as it reaches the cartesian RMSD, so one flag governs both. Over the
        heavy atoms a structure is exactly its own template; with the
        hydrogens in, it is not.
        """

        pytest.importorskip('openmm')

        folder, _ = self.template_folder(tmp_path / 'site')
        manager = self.manager()
        manager.load_template_from_folder(folder)
        template = manager.templates['site']
        site = 'active_site'

        builder = self.active_site_builder()
        manager.compare_active_site(builder, mm_opt=False)
        heavy = manager._comparison['templates']['site']['regions'][site][
            'ic_rmsd']

        assert heavy['bonds']['rms'] == pytest.approx(0.0, abs=1.0e-9)
        assert heavy['angles']['rms'] == pytest.approx(0.0, abs=1.0e-9)
        atoms = template['molecule'].number_of_atoms()
        assert len(self.shells().matcher.rmsd_indices(template, site,
                                                      True)) < atoms

        manager.compare_active_site(mm_opt=False, include_hydrogens=True)
        with_hydrogens = manager._comparison['templates']['site']['regions'][
            site]['ic_rmsd']

        assert with_hydrogens['angles']['rms'] > 0.1

        assert len(self.shells().matcher.rmsd_indices(template, site,
                                                      False)) == atoms

        # the setting is what include_hydrogens follows when it is not given
        manager.rmsd_heavy_atoms_only = False
        manager.compare_active_site(mm_opt=False)
        assert manager._comparison['include_hydrogens'] is True
        manager.rmsd_heavy_atoms_only = True
        manager.compare_active_site(mm_opt=False)
        assert manager._comparison['include_hydrogens'] is False

    def test_the_metal_shell_region_is_the_coordination_sphere(self, tmp_path):
        """
        Restricting an RMSD to the metals and what is within two bonds of them
        measures the part a transferred force field actually describes.

        On this site that region is every heavy atom of the coordination
        sphere and no hydrogen at all, which also puts it out of reach of the
        hydrogen placement noise: the same structure comes out at exactly zero
        against its own template, where the whole site does not.
        """

        pytest.importorskip('openmm')

        folder, _ = self.template_folder(tmp_path / 'site')
        manager = self.manager()
        manager.load_template_from_folder(folder)
        template = manager.templates['site']
        manager.rmsd_heavy_atoms_only = False

        whole = self.shells().matcher.region_indices(template, 'active_site')
        labels = template['molecule'].get_labels()
        assert whole == list(range(len(labels)))

        shell = self.shells().matcher.region_indices(template, 'metal_shell')

        assert set(template['metal_indices']) <= set(shell)
        assert len(shell) < len(whole)
        assert all(labels[index] != 'H' for index in shell)
        # the donor atoms are one bond out
        assert {'N', 'O', 'S'} <= {labels[index] for index in shell}

        # every region at once, so the two are measured on one mapping
        manager.compare_active_site(self.active_site_builder(),
                                    mm_opt=False,
                                    include_hydrogens=True)
        regions = manager._comparison['templates']['site']['regions']

        shell_found = regions['metal_shell']
        whole_found = regions['active_site']

        assert shell_found['rmsd'] == pytest.approx(0.0, abs=1.0e-9)
        assert shell_found['ic_rmsd']['angles']['rms'] == pytest.approx(
            0.0, abs=1.0e-9)

        # the same comparison over the whole site carries the hydrogen noise
        assert whole_found['rmsd'] > shell_found['rmsd']

    def test_the_metal_beta_carbon_region_is_the_frame_of_the_site(
            self, tmp_path):
        """
        The metals with the beta carbons are the frame the residues hang on:
        where each sidechain is anchored relative to the metals, with the
        conformation of the sidechain itself left out of it.

        Neither a force field nor a geometry says which atom is a CB. The
        truncation knows, and writes it into the comment of the atom, which is
        what survives into forcefield.json for a template to read back.
        """

        pytest.importorskip('openmm')

        folder, _ = self.template_folder(tmp_path / 'site')
        manager = self.manager()
        manager.load_template_from_folder(folder)
        template = manager.templates['site']
        region = self.shells().matcher.region_indices(template, 'metal_beta_carbons')

        # both zincs and one beta carbon per truncated residue
        assert len(region) == 2 + 6
        assert set(template['metal_indices']) <= set(region)
        labels = template['molecule'].get_labels()
        assert all(labels[index] == 'C'
                   for index in region
                   if index not in template['metal_indices'])

        # they are the atoms the builder marked, and no capping hydrogen is
        # in the region itself
        assert set(template['beta_carbon_indices']) < set(region)
        assert not set(region) & set(template['cap_indices'])

        marker = util.BETA_CARBON_COMMENT
        for index in template['beta_carbon_indices']:
            assert marker in template['forcefield'].atoms[index]['comment']

        # every atom of it is heavy, so leaving the hydrogens out changes
        # nothing
        assert self.shells().matcher.rmsd_indices(
            template, 'metal_beta_carbons', True) == region

        manager.compare_active_site(self.active_site_builder(), mm_opt=False)
        found = manager._comparison['templates']['site']['regions'][
            'metal_beta_carbons']

        assert found['rmsd'] == pytest.approx(0.0, abs=1.0e-9)
        assert found['atoms'] == len(region)

    @staticmethod
    def stretched_structure(source, target, push=0.8):
        """
        Pushes one oxygen of the bidentate ASP130 away from its metal, far
        enough that it stops counting as a ligand.

        The residue is the same residue in the same place either way; only how
        many of its atoms are within the cutoff changes, which is what a
        comparison is meant not to turn on.

        :param source:
            The structure to alter.
        :param target:
            Where to write the altered structure.
        :param push:
            How far to move the oxygen, in Angstrom.

        :return:
            The path written.
        """

        import openmm.app as mmapp
        import openmm.unit as mmunit

        pdb = mmapp.PDBFile(str(source))
        positions = pdb.getPositions(asNumpy=True).value_in_unit(
            mmunit.angstrom)

        metals = [
            atom.index
            for atom in pdb.topology.atoms()
            if atom.element is not None and atom.element.symbol == 'Zn'
        ]
        oxygens = [
            atom.index
            for atom in pdb.topology.atoms()
            if atom.name == 'OD2' and atom.residue.name == 'ASP' and
            atom.residue.id == '130'
        ]
        assert len(oxygens) == 1

        index = oxygens[0]
        nearest = min(metals,
                      key=lambda metal: np.linalg.norm(positions[index] -
                                                       positions[metal]))
        direction = positions[index] - positions[nearest]
        positions[index] += direction / np.linalg.norm(direction) * push

        target = Path(target)
        with target.open('w') as fh:
            mmapp.PDBFile.writeFile(pdb.topology,
                                    positions * mmunit.angstrom,
                                    fh,
                                    keepIds=True)

        return str(target)

    def test_the_residues_are_what_the_metals_leave_behind(self, tmp_path):
        """
        A template is a force field and a geometry, and neither says where one
        residue ends and the next begins. Taking the metals out of the graph
        does: what is left is one connected fragment per truncated sidechain,
        and that decomposition does not depend on how a residue grips a metal.
        """

        pytest.importorskip('openmm')

        folder, _ = self.template_folder(tmp_path / 'site')
        manager = self.manager()
        manager.load_template_from_folder(folder)
        template = manager.templates['site']

        coarse = template['coarse_topology']
        residues = [
            coarse.nodes[node] for node in coarse if node[0] == 'residue'
        ]
        metals = set(template['metal_indices'])
        heavy_atoms = [
            index
            for index, label in enumerate(template['molecule'].get_labels())
            if label != 'H'
        ]

        assert len(residues) == 6
        assert all(not metals & set(residue['heavy']) for residue in residues)
        assert sum(residue['heavy'].number_of_nodes()
                   for residue in residues) + len(metals) == len(heavy_atoms)

        keys = [residue['key'] for residue in residues]
        formulas = sorted(residue['formula'] for residue in residues)

        # the two aspartates and the two histidines are each one residue
        # twice, and everything else is its own
        assert len(set(keys)) == 4
        assert formulas == [
            'C2H3O2', 'C2H3O2', 'C3H5O2', 'C4H6N2', 'C4H6N2', 'CH3S'
        ]

        # the bridging glutamate is the one that reaches both metals
        bridging = [
            node for node in self.shells().matcher.residue_nodes(coarse)
            if coarse.degree(node) > 1
        ]
        assert len(bridging) == 1
        assert coarse.nodes[bridging[0]]['formula'] == 'C3H5O2'

    def test_denticity_does_not_decide_a_match(self, tmp_path):
        """
        A carboxylate that grips a metal with one oxygen in the structure and
        with two in the template is the same residue on the same metal. The
        coarse level says so, since it carries no denticity at all, and the
        contact that only one of the two makes is counted and reported rather
        than allowed to fail the comparison.
        """

        pytest.importorskip('openmm')

        folder, _ = self.template_folder(tmp_path / 'bidentate')
        monodentate = self.stretched_structure(self.structure_file(),
                                               tmp_path / 'mono.pdb')

        manager = self.manager()
        manager.load_template_from_folder(folder)
        template = manager.templates['bidentate']

        manager.compare_active_site(self.active_site_builder(monodentate),
                                    mm_opt=False)
        entry = manager._comparison['templates']['bidentate']

        # the site spec is the same site, whatever the metal bonds say
        assert entry['status'] == 'measured'
        assert self.shells().matcher.site_spec(
            manager._comparison['active_site']) == (self.shells().matcher.site_spec(
                template))

        bonds = entry['metal_bonds']
        assert bonds['template_only'] == 1
        assert bonds['query_only'] == 0
        assert bonds['shared'] == len(self.shells().matcher.metal_keys(template)[0]) - 1
        # the bonds they do share are the ones they had before
        assert bonds['deviation'] == pytest.approx(0.0, abs=1.0e-6)

        # and the geometry still decides, on a site that barely moved
        assert entry['regions']['metal_shell']['rmsd_heavy'] < 0.5

    def test_two_residues_on_one_metal_give_several_mappings(self, tmp_path):
        """
        Two histidines on the same metal can be lined up either way round, so
        the coarse level hands back both and the geometry chooses. The atom
        mappings are built residue by residue out of those, which is far fewer
        than an isomorphism of the whole site would enumerate.
        """

        pytest.importorskip('openmm')

        folder, _ = self.template_folder(tmp_path / 'site')
        manager = self.manager()
        manager.load_template_from_folder(folder)
        template = manager.templates['site']

        manager.compare_active_site(self.active_site_builder(), mm_opt=False)
        entry = manager._comparison['templates']['site']

        assert entry['n_coarse_mappings'] == 2
        assert 1 < entry['n_mappings'] < 100

        # the mapping that came out of it covers every atom, element for
        # element, and the metals map onto metals
        mapping = entry['mapping']
        atoms = template['molecule'].number_of_atoms()
        assert len(mapping) == atoms
        assert sorted(mapping.values()) == list(range(atoms))

        labels = manager._comparison['active_site']['molecule'].get_labels()
        template_labels = template['molecule'].get_labels()
        for index, image in mapping.items():
            assert labels[image] == template_labels[index]

        assert set(
            mapping[index] for index in template['metal_indices']) == (set(
                manager._comparison['active_site']['metal_indices']))

    def test_comparing_measures_every_template_over_every_region(
            self, tmp_path):
        """
        The comparison decides nothing: it measures every template over every
        region and hands the numbers back, so a decision the criteria cannot
        make can be made on them afterwards.
        """

        pytest.importorskip('openmm')

        folder, _ = self.template_folder(tmp_path / 'site')
        manager = self.manager()
        manager.load_template_from_folder(folder)

        manager.compare_active_site(self.active_site_builder(), mm_opt=False)
        results = manager._comparison

        assert results['geometry'] == 'input'
        assert results['include_hydrogens'] is False
        assert set(results['templates']) == {'site'}

        entry = results['templates']['site']
        assert entry['status'] == 'measured'
        assert entry['n_mappings'] > 1
        assert entry['metal_bonds']['deviation'] == pytest.approx(0.0,
                                                                  abs=1.0e-9)
        assert set(entry['regions']) == set(SiteMatcher.RMSD_REGIONS)

        for region, found in entry['regions'].items():
            assert found['atoms'] > 0
            assert found['rmsd_heavy'] == pytest.approx(0.0, abs=1.0e-9)
            assert found['ic_rmsd']['bonds']['rms'] == pytest.approx(
                0.0, abs=1.0e-9), region

    def test_comparing_reports_what_it_could_not_measure(self, tmp_path):
        """
        A template of another composition and one wired differently are both
        answers worth having, and each says which of the two it is rather
        than dropping out of the report.
        """

        pytest.importorskip('openmm')

        folder, _ = self.template_folder(tmp_path / 'site')

        # a site one coordinating residue short holds different atoms
        lines = Path(self.structure_file()).read_text().splitlines()
        kept = [
            line for line in lines
            if not (line.startswith('ATOM') and
                    line[17:26].strip().endswith('HIS A 191'))
        ]
        short = tmp_path / 'short.pdb'
        short.write_text('\n'.join(kept) + '\n')

        manager = self.manager()
        manager.load_template_from_folder(folder)

        manager.compare_active_site(self.active_site_builder(str(short)),
                                    mm_opt=False)
        entry = manager._comparison['templates']['site']

        assert entry['status'] == 'composition'
        assert entry['mapping'] is None
        assert entry['regions'] == {}

        # the same atoms coordinating differently: a residue taken off the
        # metal it binds, which the coarse level does notice
        manager = self.manager()
        manager.load_template_from_folder(folder)
        template = manager.templates['site']
        coarse = template['coarse_topology']
        fine = template['fine_topology']
        node = next(one for one in self.shells().matcher.residue_nodes(coarse)
                    if coarse.degree(one) == 1)
        metal = next(iter(coarse.neighbors(node)))[1]

        for atom in coarse.nodes[node]['heavy']:
            if fine.has_edge(metal, atom):
                fine.remove_edge(metal, atom)

        # a template is an active site, so it is re-described in place
        template.update(self.shells().matcher.describe(template, fine.edges()))

        manager.compare_active_site(self.active_site_builder(), mm_opt=False)

        assert manager._comparison['templates']['site']['status'] == 'spec'
        assert manager._comparison['templates']['site']['mapping'] is None

    def test_comparing_takes_the_hydrogens_on_request(self, tmp_path):
        """
        The hydrogens are out of every number by default, and the parameter
        puts them back for the call without touching the setting.
        """

        pytest.importorskip('openmm')

        folder, _ = self.template_folder(tmp_path / 'site')
        manager = self.manager()
        manager.load_template_from_folder(folder)

        builder = self.active_site_builder()
        manager.compare_active_site(builder, mm_opt=False)
        without = manager._comparison
        manager.compare_active_site(mm_opt=False, include_hydrogens=True)
        with_them = manager._comparison

        assert manager.rmsd_heavy_atoms_only is True
        assert with_them['include_hydrogens'] is True

        site = 'active_site'
        assert without['templates']['site']['regions'][site]['atoms'] < (
            with_them['templates']['site']['regions'][site]['atoms'])
        assert with_them['templates']['site']['regions'][site]['ic_rmsd'][
            'angles']['rms'] > 0.0

    def test_comparing_relaxes_the_structure_on_request(self, tmp_path):
        """
        One geometry is measured per call, and which one it was is recorded
        with the numbers. The relaxation is a comparison rather than a run, so
        it leaves no folder behind.
        """

        pytest.importorskip('openmm')

        folder, _ = self.template_folder(tmp_path / 'site')
        manager = self.manager()
        manager.load_template_from_folder(folder)
        manager.builder.folder = str(tmp_path / 'unwanted')

        manager.compare_active_site(self.active_site_builder())
        relaxed = manager._comparison

        assert relaxed['geometry'] == 'mm_relaxed'
        assert not (tmp_path / 'unwanted').exists()

        region = relaxed['templates']['site']['regions']['metal_shell']
        assert region['rmsd_heavy'] > 0.0

    def test_the_comparison_tables_are_printed(self, tmp_path, capsys):
        """
        Every number that comes back is also printed, since reading a table is
        how the numbers get looked at in the first place.
        """

        pytest.importorskip('openmm')

        folder, _ = self.template_folder(tmp_path / 'site')

        manager = MetalForceFieldManager()
        manager.load_template_from_folder(folder)
        capsys.readouterr()

        manager.compare_active_site(self.active_site_builder(), mm_opt=False)

        printed = capsys.readouterr().out

        assert 'Comparison against every template' in printed
        for region in SiteMatcher.RMSD_REGIONS:
            assert region in printed

        region, ic_type, measure = SiteMatcher.SELECTION_RANKED_ON
        assert f'Ranked on the {region} {ic_type} {measure}' in printed

    def test_a_template_that_predates_the_marker_says_so(self, tmp_path):
        """
        A force field written before the beta carbons were recorded cannot say
        where the sidechains were cut, and there is nothing to work it out
        from that is not a guess. The region says so rather than measuring
        over the wrong atoms.
        """

        pytest.importorskip('openmm')

        folder, forcefield = self.template_folder(tmp_path / 'old')

        marker = util.BETA_CARBON_COMMENT
        for atom in forcefield.atoms.values():
            atom['comment'] = atom['comment'].replace(marker, '').strip('; ')

        builder = MetalSiteForceFieldBuilder()
        builder.ostream.mute()
        util.save_forcefield(folder / util.FORCEFIELD_FILE, forcefield)

        manager = self.manager()
        manager.load_template_from_folder(folder)
        template = manager.templates['old']

        assert template['beta_carbon_indices'] == []

        with pytest.raises(Exception, match='beta carbons'):
            self.shells().matcher.region_indices(template, 'metal_beta_carbons')

    def test_the_capping_hydrogens_come_from_the_comments(self, tmp_path):
        """
        A template reads its capping hydrogens out of what the truncation
        wrote, not out of their charges.

        redistribute_cap_charges does leave them at exactly zero, and that is
        what the caps used to be found by. It is the correction doing its job
        rather than a marker: it says nothing in a force field whose charges
        are all zero, and it makes every hydrogen look like a cap there.
        """

        pytest.importorskip('openmm')

        folder, forcefield = self.template_folder(tmp_path / 'site')

        marker = util.CAP_COMMENT
        marked = [
            index for index, atom in forcefield.atoms.items()
            if marker in atom['comment']
        ]

        manager = self.manager()
        manager.load_template_from_folder(folder)
        template = manager.templates['site']

        assert template['cap_indices'] == marked
        # one per truncated residue, and all of them hydrogens
        assert len(template['cap_indices']) == 6
        labels = template['molecule'].get_labels()
        assert all(labels[index] == 'H' for index in template['cap_indices'])

        # empty every charge: the caps are still the same atoms
        for atom in forcefield.atoms.values():
            atom['charge'] = 0.0

        builder = MetalSiteForceFieldBuilder()
        builder.ostream.mute()
        util.save_forcefield(folder / util.FORCEFIELD_FILE, forcefield)

        manager.load_template_from_folder(folder)
        assert manager.templates['site']['cap_indices'] == marked

    def test_a_template_that_predates_the_cap_marker_says_so(self, tmp_path):
        """
        A force field written before the capping hydrogens were recorded
        cannot say which atoms they are, and guessing from the charges is
        what this replaced.
        """

        pytest.importorskip('openmm')

        folder, forcefield = self.template_folder(tmp_path / 'old')

        marker = util.CAP_COMMENT
        for atom in forcefield.atoms.values():
            atom['comment'] = atom['comment'].replace(marker, '').strip('; ')

        builder = MetalSiteForceFieldBuilder()
        builder.ostream.mute()
        util.save_forcefield(folder / util.FORCEFIELD_FILE, forcefield)

        manager = self.manager()

        with pytest.raises(Exception, match='capping hydrogens'):
            manager.load_template_from_folder(folder)

    def test_the_metal_shell_reaches_as_far_as_it_is_told(
            self, tmp_path, monkeypatch):
        """
        One bond out is the donor atoms alone, and every bond further pulls in
        another shell of the residues.
        """

        pytest.importorskip('openmm')

        folder, _ = self.template_folder(tmp_path / 'site')
        manager = self.manager()
        manager.load_template_from_folder(folder)
        template = manager.templates['site']

        sizes = []
        for bonds in (1, 2, 3):
            # the reach is a constant of the matcher, read where it is defined
            monkeypatch.setattr(SiteMatcher, 'METAL_SHELL_BONDS', bonds)
            sizes.append(
                len(
                    self.shells().matcher.region_indices(template,
                                                         'metal_shell')))

        assert sizes[0] < sizes[1] < sizes[2]
        # the metals and their donors
        assert sizes[0] == 2 + 7

    def test_hydrogen_positions_are_not_reproducible(self):
        """
        protonate places the hydrogens with Modeller.addHydrogens, which does
        not put them in the same place twice. Half the atoms of a truncated
        site are hydrogens, so an all-atom RMSD carries a noise floor of a few
        tenths of an Angstrom while the heavy atoms come out identical.

        That is why the RMSD is taken over the heavy atoms by default. This
        test is here so that flipping it back is a deliberate act.
        """

        pytest.importorskip('openmm')

        from veloxchem.superimpose import svd_superimpose

        first = self.prepared()['active_site']
        second = self.prepared()['active_site']

        assert first['molecule'].get_labels() == (
            second['molecule'].get_labels())

        one = first['molecule'].get_coordinates_in_angstrom()
        two = second['molecule'].get_coordinates_in_angstrom()
        heavy = [
            index for index, label in enumerate(first['molecule'].get_labels())
            if label != 'H'
        ]

        assert svd_superimpose(one[heavy],
                               two[heavy])[0] == pytest.approx(0.0, abs=1.0e-9)
        assert svd_superimpose(one, two)[0] > 0.1

        assert MetalForceFieldManager().rmsd_heavy_atoms_only is True

    def test_the_selection_criteria_are_what_the_class_says(self):
        """
        The two criteria sets are the contract of the selection, so the
        numbers themselves are written down here: a typo in them would be
        invisible everywhere else, since every site would still be measured
        and ranked exactly as before.
        """

        criteria = MetalForceFieldManager.SELECTION_CRITERIA

        assert criteria['tight'] == {
            'active_site': {
                'bonds': {
                    'rms': 0.5,
                    'max': 1.0
                }
            },
            'metal_shell': {
                'bonds': {
                    'rms': 0.5,
                    'max': 1.0
                }
            },
            'metal_beta_carbons': {
                'bonds': {
                    'rms': 0.5,
                    'max': 1.0
                }
            },
        }

        assert criteria['loose'] == {
            'active_site': {
                'bonds': {
                    'rms': 1.0,
                    'max': 2.0
                }
            },
            'metal_shell': {
                'bonds': {
                    'rms': 1.0,
                    'max': 1.5
                }
            },
            'metal_beta_carbons': {
                'bonds': {
                    'rms': 1.0,
                    'max': 1.5
                }
            },
        }

        # the tight set is the default, and the loose one is never tighter
        assert MetalForceFieldManager().selection_criteria == criteria['tight']

        for region, thresholds in criteria['tight'].items():
            for measure, limit in thresholds['bonds'].items():
                assert criteria['loose'][region]['bonds'][measure] >= limit

    def test_building_from_a_comparison_transfers_the_template(self, tmp_path):
        """
        The comparison measures and the selection decides: given the
        findings, the template that fits carries its metal terms and its
        charges onto the site, exactly as a match does, and the built force
        field is adopted onto the active site so create_enzyme_system can be
        called on it directly.
        """

        pytest.importorskip('openmm')

        folder, forcefield = self.template_folder(tmp_path / 'site')
        manager = self.manager()
        manager.load_template_from_folder(folder)
        template = manager.templates['site']

        builder = self.active_site_builder()
        manager.compare_active_site(builder, mm_opt=False)
        decision = manager._select_template(None)

        assert decision['name'] == 'site'
        assert decision['criteria_name'] == 'tight'
        assert decision['forced'] is False
        assert decision['candidates'] == ['site']
        assert manager._comparison['geometry'] == 'input'
        assert template['geometry_kind'] == 'qm_opt'

        # the site is its own template, so it is built exactly like it
        assert decision['score'] == pytest.approx(0.0, abs=1.0e-6)

        built_ff = manager.build_ff_from_template()

        # the force field is handed back to the active site's own builder
        assert manager.active_site is builder
        assert builder.stage == 'fitted'
        assert builder.active_site_forcefield is built_ff

        bonds, angles = self.shells().matcher.metal_keys(template)
        mapping = decision['entry']['mapping']

        for key in bonds:
            mapped = tuple(mapping[index] for index in key)
            found = built_ff.bonds.get(mapped, built_ff.bonds.get(mapped[::-1]))
            assert found['force_constant'] == pytest.approx(
                forcefield.bonds[key]['force_constant'])
            assert 'template site' in found['comment']

        for key in angles:
            mapped = tuple(mapping[index] for index in key)
            found = built_ff.angles.get(mapped,
                                        built_ff.angles.get(mapped[::-1]))
            assert found['force_constant'] == pytest.approx(
                forcefield.angles[key]['force_constant'])

        charges = built_ff.partial_charges
        for template_index, site_index in mapping.items():
            assert charges[site_index] == pytest.approx(
                template['charges'][template_index])

        # the caps carry nothing, here as in the template
        for index in manager._comparison['active_site']['cap_indices']:
            assert charges[index] == 0.0

    def test_active_site_tracks_the_live_builder_state(self, tmp_path):
        """
        compare_active_site tracks one builder as the active site and reads
        it back through the active_site property, so it can be edited in
        between comparisons. A bare call with no argument re-describes that
        builder from scratch, which is what lets an edit made after the first
        comparison be picked up by the next one without handing the builder
        in again.
        """

        pytest.importorskip('openmm')

        folder, _ = self.template_folder(tmp_path / 'site')
        manager = self.manager()
        manager.load_template_from_folder(folder)

        builder = self.active_site_builder()
        manager.compare_active_site(builder, mm_opt=False)

        assert manager.active_site is builder
        assert manager._comparison['templates']['site']['status'] == 'measured'

        # ASP134 grips the first metal through OD2 alone; adding a bond to
        # the second metal through OD1 makes it bridge both metals, which the
        # template's own coordination does not. The coarse level is exactly
        # which residue reaches which metal, so the two sites no longer line
        # up there and the match cannot be made any more.
        metal = self.metal_label(builder, 1)
        builder.add_metal_bond('134', metal, atom='OD1')

        manager.compare_active_site()

        assert manager.active_site is builder
        assert manager._comparison['templates']['site']['status'] != 'measured'

    def test_active_site_is_none_until_a_builder_is_given(self):
        """
        There is nothing to read back before compare_active_site has ever
        been handed a builder.
        """

        manager = self.manager()

        assert manager.active_site is None

    def test_building_needs_a_comparison_first(self):
        """
        There is no side-step from a bare structure to a built force field:
        build_ff_from_template always reads the last comparison, and there is
        a clear error rather than a stale one when none exists yet.
        """

        manager = self.manager()

        with pytest.raises(Exception, match='compare_active_site'):
            manager.build_ff_from_template()

    def test_a_template_outside_the_criteria_is_not_taken(self, tmp_path):
        """
        The criteria are a floor, not a preference: a template that is over
        them is refused however much the nearest thing on file it is,
        whether it was chosen automatically or named by hand. Which region
        failed, and by how much, is what the error names.
        """

        pytest.importorskip('openmm')

        folder, _ = self.template_folder(tmp_path / 'relaxed', relax=True)
        manager = self.manager()
        manager.load_template_from_folder(folder)

        # the structure has not been relaxed and the template has, so the
        # residues are anchored a few tenths of an Angstrom away from where
        # the template holds them
        manager.compare_active_site(self.active_site_builder(), mm_opt=False)
        beta = manager._comparison['templates']['relaxed']['regions'][
            'metal_beta_carbons']['ic_rmsd']['bonds']

        assert 0.1 < beta['rms'] < 0.5

        strict = {'metal_beta_carbons': {'bonds': {'rms': 0.1}}}

        manager.selection_criteria = strict
        decision = manager._select_template(None)

        assert decision['name'] is None
        assert decision['candidates'] == []
        assert 'metal_beta_carbons bonds rms' in decision['verdicts']['relaxed']

        with pytest.raises(Exception, match='criteria'):
            manager.build_ff_from_template()

        with pytest.raises(Exception, match='does not match'):
            manager.build_ff_from_template(template='relaxed')

        # the same site under the criteria it does meet
        manager.selection_criteria = 'loose'
        assert manager._select_template(None)['name'] == 'relaxed'
        assert manager.build_ff_from_template() is not None

    def test_a_criterion_that_was_not_measured_is_not_passed(self, tmp_path):
        """
        A region whose internal coordinates could not be measured is not a
        region that passed. Reading it as one would transfer a force field on
        the strength of a number that was never taken.
        """

        pytest.importorskip('openmm')

        folder, _ = self.template_folder(tmp_path / 'site')
        manager = self.manager()
        manager.load_template_from_folder(folder)

        manager.compare_active_site(self.active_site_builder(), mm_opt=False)

        assert manager._select_template(None)['name'] == 'site'

        entry = manager._comparison['templates']['site']
        entry['regions']['metal_shell']['ic_rmsd'].pop('bonds')

        decision = manager._select_template(None)

        assert decision['name'] is None
        assert decision['verdicts']['site'] == 'metal_shell no bonds measured'

        entry['regions']['metal_shell']['ic_rmsd'] = None

        assert manager._select_template(None)['verdicts']['site'] == (
            'metal_shell no internal '
            'coords')

    def test_the_closest_template_is_taken_when_several_pass(self, tmp_path):
        """
        Several templates describing one site is a question about the
        templates rather than about the structure, and the geometry is all
        there is to answer it with. The site that is built the most like the
        structure wins, which is the bonds of the whole active site rather
        than where it sits in space.
        """

        pytest.importorskip('openmm')

        near, _ = self.template_folder(tmp_path / 'near')
        far, _ = self.template_folder(tmp_path / 'far', relax=True)

        manager = self.manager()
        manager.load_template_from_folder(near)
        manager.load_template_from_folder(far)

        manager.compare_active_site(self.active_site_builder(), mm_opt=False)
        decision = manager._select_template(None)

        assert decision['candidates'] == ['near', 'far']
        assert decision['name'] == 'near'
        assert decision['scores']['near'] < decision['scores']['far']

        manager.build_ff_from_template()

        # the same comparison is still there, and it still picks near
        assert manager._select_template(None)['name'] == 'near'

    def test_a_named_template_bypasses_the_ranking(self, tmp_path):
        """
        Naming a template is a statement that it is the right one, so it does
        not have to win a ranking among several candidates -- but mapping
        onto every atom still decides whether it can be used at all, since a
        hole in the mapping cannot be forced shut.
        """

        pytest.importorskip('openmm')

        folder, _ = self.template_folder(tmp_path / 'site')
        manager = self.manager()
        manager.load_template_from_folder(folder)

        manager.compare_active_site(self.active_site_builder(), mm_opt=False)
        manager.selection_criteria = {
            'active_site': {
                'bonds': {
                    'rms': 1.0e-6
                }
            }
        }

        decision = manager._select_template('site')

        assert decision['name'] == 'site'
        assert decision['forced'] is True
        assert decision['criteria_name'] == 'custom'

        built_ff = manager.build_ff_from_template(template='site')
        assert built_ff is not None

        with pytest.raises(Exception, match='no template named'):
            manager._select_template("other")

    def test_a_template_must_map_onto_every_atom(self, tmp_path):
        """
        A template that holds other atoms describes another site, and no part
        of it is usable: it is not taken, and naming it is an error rather
        than a warning.
        """

        pytest.importorskip('openmm')

        folder, _ = self.template_folder(tmp_path / 'site')

        # a site one coordinating residue short is a different site
        lines = Path(self.structure_file()).read_text().splitlines()
        kept = [
            line for line in lines
            if not (line.startswith('ATOM') and
                    line[17:26].strip().endswith('HIS A 191'))
        ]
        short = tmp_path / 'short.pdb'
        short.write_text('\n'.join(kept) + '\n')

        manager = self.manager()
        manager.load_template_from_folder(folder)

        manager.compare_active_site(self.active_site_builder(str(short)),
                                    mm_opt=False)

        assert manager._comparison['templates']['site']['status'] == (
            'composition')

        decision = manager._select_template(None)

        assert decision['name'] is None
        assert decision['verdicts']['site'] == 'different atoms'

        with pytest.raises(Exception, match='criteria'):
            manager.build_ff_from_template()

        with pytest.raises(Exception, match='every atom'):
            manager.build_ff_from_template(template='site')

    def test_a_bidentate_template_builds_a_bidentate_site(self, tmp_path):
        """
        Denticity is a distance cutoff on an unrelaxed structure rather than
        chemistry, which is why the matching refuses to have an opinion about
        it. What matches on those terms has to be built on the template's
        terms too: a template bond the site does not make has no atoms to
        land on, and before the coordination was forced that was a fatal
        error at the transfer.

        The fixture holds ASP130 at 2.50 and 2.57 A from Zn1, with the next
        contact 0.2 A further out, so a cutoff either side of 2.55 makes the
        same site monodentate or bidentate and changes nothing else.
        """

        pytest.importorskip('openmm')

        folder, _ = self.template_folder(tmp_path / 'site',
                                         metal_bond_cutoff=2.6)

        manager = self.manager()
        manager.load_template_from_folder(folder)
        template = manager.templates['site']

        manager.compare_active_site(
            self.active_site_builder(metal_bond_cutoff=2.55), mm_opt=False)
        entry = manager._comparison['templates']['site']

        # the site is measured as it was detected, and it is one bond short
        assert entry['status'] == 'measured'
        assert entry['metal_bonds']['template_only'] == 1
        assert entry['metal_bonds']['query_only'] == 0

        active_site = manager._comparison['active_site']
        before = np.array(active_site['connectivity_matrix'])

        forcefield, changes, built_site = manager._loader.build_forcefield_from_template(
            template, entry['mapping'], active_site,
            manager.builder._sites())

        assert len(changes['added']) == 1
        assert changes['removed'] == []

        wanted = self.mapped_metal_bonds(template, entry['mapping'])
        built = self.metal_bonds(forcefield, built_site['metal_indices'])

        assert built == wanted

        # the added bond is the second carboxylate oxygen of ASP130
        first, second = changes['added'][0]
        labels = built_site['molecule'].get_labels()

        assert sorted((labels[first], labels[second])) == ['O', 'Zn']

        # what compare_active_site measured is left as it was measured
        assert np.array_equal(before,
                              np.array(active_site['connectivity_matrix']))
        assert not np.array_equal(
            before, np.array(built_site['connectivity_matrix']))

    def test_a_monodentate_template_builds_a_monodentate_site(self, tmp_path):
        """
        The other way round is the quiet half of the same problem: a contact
        the site makes alone has no fitted parameters to land on it, so
        before the coordination was forced it kept the seeded guess and said
        nothing about it. Every metal term of the built force field has to
        come from the template.
        """

        pytest.importorskip('openmm')

        folder, _ = self.template_folder(tmp_path / 'site',
                                         metal_bond_cutoff=2.55)

        manager = self.manager()
        manager.load_template_from_folder(folder)
        template = manager.templates['site']

        manager.compare_active_site(
            self.active_site_builder(metal_bond_cutoff=2.6), mm_opt=False)
        entry = manager._comparison['templates']['site']

        assert entry['metal_bonds']['query_only'] == 1

        active_site = manager._comparison['active_site']
        forcefield, changes, built_site = manager._loader.build_forcefield_from_template(
            template, entry['mapping'], active_site,
            manager.builder._sites())

        assert changes['added'] == []
        assert len(changes['removed']) == 1

        metals = set(built_site['metal_indices'])
        wanted = self.mapped_metal_bonds(template, entry['mapping'])

        assert self.metal_bonds(forcefield, metals) == wanted

        seeded = [
            key for table in (forcefield.bonds, forcefield.angles)
            for key, params in table.items()
            if (metals & set(key)) and 'template site' not in params['comment']
        ]

        assert seeded == []

    def test_a_site_wired_like_the_template_is_left_alone(self, tmp_path):
        """
        Forcing the coordination is for the sites whose coordination differs.
        One detected the same way is handed on as it came, rather than
        rebuilt into an equal copy of itself.
        """

        pytest.importorskip('openmm')

        folder, _ = self.template_folder(tmp_path / 'site',
                                         metal_bond_cutoff=2.6)

        manager = self.manager()
        manager.load_template_from_folder(folder)
        template = manager.templates['site']

        manager.compare_active_site(
            self.active_site_builder(metal_bond_cutoff=2.6), mm_opt=False)
        entry = manager._comparison['templates']['site']
        active_site = manager._comparison['active_site']

        forcefield, changes, built_site = manager._loader.build_forcefield_from_template(
            template, entry['mapping'], active_site,
            manager.builder._sites())

        assert changes == {'added': [], 'removed': []}
        assert built_site is active_site

    def test_the_selection_criteria_are_validated(self, tmp_path):
        """
        A criteria set that names a region, a coordinate type or a measure
        that does not exist would otherwise be a criterion that is silently
        never checked.
        """

        pytest.importorskip('openmm')

        folder, _ = self.template_folder(tmp_path / 'site')
        manager = self.manager()
        manager.load_template_from_folder(folder)

        manager.compare_active_site(self.active_site_builder(), mm_opt=False)

        bad = [
            ('tight', 'strict'),
            ('not a region', {
                'everything': {
                    'bonds': {}
                }
            }),
            ('internal', {
                'active_site': {
                    'lengths': {}
                }
            }),
            ('measure', {
                'active_site': {
                    'bonds': {
                        'mean': 0.5
                    }
                }
            }),
        ]

        for message, criteria in bad:
            manager.selection_criteria = criteria
            with pytest.raises(Exception, match=message):
                manager._select_template(None)

        manager.selection_criteria = 'tight'

        fresh = self.manager()
        with pytest.raises(Exception, match='compare_active_site'):
            fresh._select_template(None)

    def test_the_selection_table_is_printed(self, tmp_path, capsys):
        """
        The whole field is printed rather than the winner alone: whether the
        others were near misses or a long way off is what says how much the
        chosen template is worth.
        """

        pytest.importorskip('openmm')

        near, _ = self.template_folder(tmp_path / 'near')
        far, _ = self.template_folder(tmp_path / 'far', relax=True)

        manager = MetalForceFieldManager()
        manager.load_template_from_folder(near)
        manager.load_template_from_folder(far)

        manager.compare_active_site(self.active_site_builder(), mm_opt=False)
        capsys.readouterr()

        # deciding and printing the decision are separate, so that a
        # comparison that decides only to answer a boolean prints nothing
        decision = manager._select_template(None)
        manager._matcher.print_selection(manager._comparison, decision)

        printed = capsys.readouterr().out

        assert 'Choosing a template on the tight criteria' in printed
        for region in MetalForceFieldManager.SELECTION_CRITERIA['tight']:
            assert region[:13] in printed
        assert 'within the criteria' in printed
        assert 'Taking near' in printed

    def test_the_selection_table_is_printed_once(self, tmp_path, capsys):
        """
        compare_active_site decides only to answer whether anything matched,
        and build_ff_from_template decides again to act on it. Printing the
        whole field for both only invited the reader to look for what had
        changed between the two tables, and nothing had.
        """

        pytest.importorskip('openmm')

        folder, _ = self.template_folder(tmp_path / 'site')

        manager = MetalForceFieldManager()
        manager.load_template_from_folder(folder)
        capsys.readouterr()

        assert manager.compare_active_site(self.active_site_builder(),
                                           mm_opt=False) is True
        printed = capsys.readouterr().out
        assert 'Choosing a template on the' not in printed

        manager.build_ff_from_template()
        printed = capsys.readouterr().out
        assert printed.count('Choosing a template on the') == 1

    def test_show_template_draws_one_template(self, tmp_path):
        """
        show_template does not need ActiveSiteBuilder.show_active_site: a template has no
        'labels' or 'connectivity_matrix' key to hand it, since it never had
        a topology behind it, so the bonds it draws come from the force field
        instead -- the same source every other template computation reads
        its bonds from.
        """

        pytest.importorskip('openmm')
        pytest.importorskip('py3Dmol')

        folder, _ = self.template_folder(tmp_path / 'site')
        manager = self.manager()
        manager.load_template_from_folder(folder)

        manager.show_template('site')

        with pytest.raises(Exception, match='no template'):
            manager.show_template('missing')

    def test_show_all_templates_draws_every_template(self, tmp_path):
        """
        show_all_templates draws every loaded template in one grid, each with
        its own explicit metal-ligand bonds.
        """

        pytest.importorskip('openmm')
        pytest.importorskip('py3Dmol')

        near, _ = self.template_folder(tmp_path / 'near')
        far, _ = self.template_folder(tmp_path / 'far', relax=True)

        manager = self.manager()
        manager.load_template_from_folder(near)
        manager.load_template_from_folder(far)

        manager.show_all_templates()

        with pytest.raises(Exception, match='no templates'):
            self.manager().show_all_templates()

    def test_shoehorn_needs_an_active_site(self, tmp_path):
        """
        shoehorn edits the tracked active site, so it needs one to edit.
        """

        pytest.importorskip('openmm')

        folder, _ = self.template_folder(tmp_path / 'site')
        manager = self.manager()
        manager.load_template_from_folder(folder)

        with pytest.raises(Exception, match='no active site'):
            manager.shoehorn('site')

    def test_shoehorn_needs_a_known_template(self, tmp_path):
        """
        Naming a template that was never loaded is an error naming the
        mistake rather than a silent no-op.
        """

        pytest.importorskip('openmm')

        folder, _ = self.template_folder(tmp_path / 'site')
        manager = self.manager()
        manager.load_template_from_folder(folder)
        manager.compare_active_site(self.active_site_builder(), mm_opt=False)

        with pytest.raises(Exception, match='no template named'):
            manager.shoehorn('missing')

    @staticmethod
    def structure_without(residue_name, path):
        """
        The fixture with every residue of one name taken out of it, for a
        site that cannot be given what a template holds.

        :param residue_name:
            The residue name to drop.
        :param path:
            Where to write the structure.

        :return:
            The path, as a string.
        """

        source = Path(TestMetalForceFieldManager.structure_file())
        kept = [
            line for line in source.read_text().splitlines(True)
            if not (line.startswith(('ATOM', 'HETATM'))
                    and line[17:20].strip() == residue_name)
        ]
        path.write_text(''.join(kept))

        return str(path)

    @staticmethod
    def coordination(builder):
        """
        Which residue binds which metal, by residue label and atom name, for
        comparing a site before and after an edit.
        """

        return sorted(
            (ligand['residue'], ligand['atom'], len(ligand['metals']))
            for ligand in builder.binding_modes['ligands'])

    def test_shoehorn_leaves_a_matching_site_alone(self, tmp_path):
        """
        A site already built the way the template is built is one there is
        nothing to do to.
        """

        pytest.importorskip('openmm')

        folder, _ = self.template_folder(tmp_path / 'site')
        manager = self.manager()
        manager.load_template_from_folder(folder)

        builder = self.active_site_builder()
        manager.compare_active_site(builder, mm_opt=False)

        before = self.coordination(builder)

        assert manager.shoehorn('site') is True
        assert manager.active_site is builder
        assert self.coordination(builder) == before

    def test_shoehorn_picks_up_a_residue_the_site_is_missing(self, tmp_path):
        """
        A residue the template is made of and the site is not is looked for
        outward from the metal centers and included, and the coordination
        stage then bonds it to the center that wants it.
        """

        pytest.importorskip('openmm')

        folder, _ = self.template_folder(tmp_path / 'site')
        manager = self.manager()
        manager.load_template_from_folder(folder)

        builder = self.active_site_builder()
        builder.remove_residue('191')

        assert 'HIS191' not in builder.active_site['residues']

        manager.compare_active_site(builder, mm_opt=False)

        assert manager.shoehorn('site') is True
        assert 'HIS191' in builder.active_site['residues']

        # and the site it left behind is one the template describes
        assert manager.compare_active_site(mm_opt=False) is True

    def test_shoehorn_measures_the_site_it_leaves_behind(self, tmp_path):
        """
        The edits are exactly what makes the last comparison stale, so the
        site is compared again before shoehorn returns rather than by a
        caller who has to remember to. build_ff_from_template reads the last
        comparison, and would otherwise be answered out of numbers measured
        on the site as it stood before the edits.
        """

        pytest.importorskip('openmm')

        folder, _ = self.template_folder(tmp_path / 'site')
        manager = self.manager()
        manager.load_template_from_folder(folder)

        builder = self.active_site_builder()
        builder.remove_residue('191')
        manager.compare_active_site(builder, mm_opt=False)

        stale = manager._comparison
        assert 'HIS191' not in stale['active_site']['residues']

        assert manager.shoehorn('site') is True

        assert manager._comparison is not stale
        assert 'HIS191' in manager._comparison['active_site']['residues']

    def test_build_ff_from_template_falls_back_to_the_shoehorned_one(
            self, tmp_path):
        """
        Walking a site onto a template is the same statement that naming one
        makes, so an unnamed call after a shoehorning is not really unnamed.
        It is also the case the criteria are least able to judge: they
        measure a geometry whose coordination sphere is still open, and what
        would close it is the force field being asked for.
        """

        pytest.importorskip('openmm')

        folder, _ = self.template_folder(tmp_path / 'site')
        manager = self.manager()
        manager.load_template_from_folder(folder)

        builder = self.active_site_builder()
        builder.remove_residue('191')
        manager.compare_active_site(builder, mm_opt=False)

        assert manager.shoehorn('site') is True

        # no measurement of a relaxed site passes this
        manager.selection_criteria = {'active_site': {'bonds': {'rms': 0.0}}}
        assert manager._select_template(None)['name'] is None

        forcefield = manager.build_ff_from_template()

        assert forcefield is not None
        assert manager.active_site.stage == 'fitted'

    def test_the_shoehorned_template_beats_the_one_ranked_first(
            self, tmp_path):
        """
        A shoehorning decides which template is built from even when the
        criteria have something of their own to say.

        A shoehorned site does not pass against its template alone: it
        passes against every template of that family, and the ranking then
        separates near-identical siblings on thousandths of an Angstrom. So
        the site would be walked onto one template and built from another,
        leaving a geometry from one and parameters from the other.
        """

        pytest.importorskip('openmm')

        exact, _ = self.template_folder(tmp_path / 'exact')
        relaxed, _ = self.template_folder(tmp_path / 'relaxed', relax=True)

        manager = self.manager()
        manager.load_template_from_folder(exact)
        manager.load_template_from_folder(relaxed)
        manager.selection_criteria = 'loose'

        manager.compare_active_site(self.active_site_builder(), mm_opt=False)

        # the site is the one 'exact' was built from, so it ranks first,
        # while 'relaxed' is within the criteria alongside it
        decision = manager._select_template(None)
        assert decision['name'] == 'exact'
        assert decision['verdicts']['relaxed'] is None

        assert manager.shoehorn('relaxed') is True

        manager.build_ff_from_template()

        assert manager.built_from == 'relaxed'

    def test_the_shoehorn_fallback_belongs_to_the_site_it_was_made_on(
            self, tmp_path):
        """
        A shoehorning is a statement about one site, not about the manager,
        so handing in another site takes it back and an unnamed call is
        refused again -- as it is for a site that was never shoehorned at
        all.
        """

        pytest.importorskip('openmm')

        folder, _ = self.template_folder(tmp_path / 'site')
        manager = self.manager()
        manager.load_template_from_folder(folder)

        impossible = {'active_site': {'bonds': {'rms': 0.0}}}

        # never shoehorned: nothing to fall back on
        manager.compare_active_site(self.active_site_builder(), mm_opt=False)
        manager.selection_criteria = impossible

        with pytest.raises(Exception, match='criteria'):
            manager.build_ff_from_template()

        # shoehorned, then handed a different site
        manager.selection_criteria = 'tight'
        builder = self.active_site_builder()
        builder.remove_residue('191')
        manager.compare_active_site(builder, mm_opt=False)

        assert manager.shoehorn('site') is True

        manager.compare_active_site(self.active_site_builder(), mm_opt=False)
        manager.selection_criteria = impossible

        with pytest.raises(Exception, match='criteria'):
            manager.build_ff_from_template()

    def test_shoehorn_puts_the_protonation_right(self, tmp_path):
        """
        A residue protonated unlike the template is what stops the atom
        mapping being solved at all, so it is put right before anything is
        measured on it.
        """

        pytest.importorskip('openmm')

        folder, _ = self.template_folder(tmp_path / 'site')
        manager = self.manager()
        manager.load_template_from_folder(folder)

        builder = self.active_site_builder()
        builder.update_protonation_state('130', 'ASH')

        manager.compare_active_site(builder, mm_opt=False)
        entry = manager._comparison['templates']['site']

        # the hydrogen is enough to make it a different site
        assert entry['status'] != 'measured'

        assert manager.shoehorn('site') is True
        assert 'ASH' not in builder.binding_modes['variants'].values()
        assert manager.compare_active_site(mm_opt=False) is True

    def test_shoehorn_wires_the_metals_the_way_the_template_does(
            self, tmp_path):
        """
        A template detected with a tighter cutoff grips its metals with fewer
        atoms than the same structure does at the default one. What the
        template does not coordinate is unbound rather than removed, and only
        a residue left coordinating nothing goes.
        """

        pytest.importorskip('openmm')

        folder, _ = self.template_folder(tmp_path / 'tight',
                                         metal_bond_cutoff=2.3)
        manager = self.manager()
        manager.load_template_from_folder(folder, name='tight')

        builder = self.active_site_builder()
        manager.compare_active_site(builder, mm_opt=False)

        bridging = [
            ligand for ligand in builder.binding_modes['ligands']
            if len(ligand['metals']) > 1
        ]
        assert bridging

        assert manager.shoehorn('tight') is True

        # the bridge is a grip in the template, and the residue stays
        assert not [
            ligand for ligand in builder.binding_modes['ligands']
            if len(ligand['metals']) > 1
        ]
        assert 'GLU73' in builder.active_site['residues']
        assert manager.compare_active_site(mm_opt=False) is True

    def test_shoehorn_moves_a_residue_between_metal_centers(self, tmp_path):
        """
        A residue the cutoffs put on the near metal is sometimes the one the
        template wants on the far one, and nothing but moving it will do.

        This is what filling each center's deficit one center at a time
        cannot express: on a real design it stole a residue that was already
        coordinating to patch one center, stranded the one that should have
        gone there, and failed. The assignment is solved over the whole site
        for that reason, and a residue moves where that is the cheapest
        answer.
        """

        pytest.importorskip('openmm')

        folder, _ = self.template_folder(tmp_path / 'site')
        manager = self.manager()
        manager.load_template_from_folder(folder)

        # ASP134 grips the first metal at 1.96 A and reaches the second with
        # its other oxygen at 3.96 A, so it is one the geometry allows to be
        # put on the wrong center
        builder = self.active_site_builder()
        first = self.metal_label(builder, 0)
        second = self.metal_label(builder, 1)

        builder.remove_metal_bond('134', metal=first)
        builder.add_metal_bond('134', second, atom='OD1')

        def metals_of(label):
            return {
                metal
                for ligand in builder.binding_modes['ligands']
                if ligand['residue'] == label for metal in ligand['metals']
            }

        assert metals_of('ASP134') == {second}

        manager.compare_active_site(builder, mm_opt=False)

        assert manager.shoehorn('site') is True

        # back where the template has it, rather than a second residue being
        # dragged in to cover for it
        assert metals_of('ASP134') == {first}
        assert manager.compare_active_site(mm_opt=False) is True

    def test_shoehorn_reports_a_residue_it_cannot_reach(self, tmp_path):
        """
        A residue the structure does not hold anywhere near the metals is
        what a shoehorning cannot do anything about, and it says so rather
        than building something else.
        """

        pytest.importorskip('openmm')

        folder, _ = self.template_folder(tmp_path / 'site')
        manager = self.manager()
        manager.load_template_from_folder(folder)

        builder = self.active_site_builder(
            self.structure_without('HIS', tmp_path / 'no_his.pdb'))
        manager.compare_active_site(builder, mm_opt=False)

        before = self.coordination(builder)

        assert manager.shoehorn('site') is False
        assert self.coordination(builder) == before

    def test_a_failed_shoehorn_puts_back_the_edits_it_made(self, tmp_path):
        """
        A run that gets part way and then cannot go on leaves the site as it
        found it: the record of the edits is what everything is derived from,
        so putting that back is the whole of undoing them.
        """

        pytest.importorskip('openmm')

        folder, _ = self.template_folder(tmp_path / 'site')
        manager = self.manager()
        manager.load_template_from_folder(folder)

        builder = self.active_site_builder()
        builder.remove_residue('187')
        builder.remove_residue('191')

        before = self.coordination(builder)
        assert 'HIS187' not in builder.active_site['residues']
        assert 'HIS191' not in builder.active_site['residues']

        # a radius that reaches the nearer histidine and not the further
        # one, so that the first is included and the second cannot be
        radius = sorted(
            self.metal_distance(builder, resid) for resid in ('187', '191'))
        assert radius[0] < radius[1]

        manager.compare_active_site(builder, mm_opt=False)

        assert manager.shoehorn('site',
                                max_include_radius=0.5 *
                                (radius[0] + radius[1])) is False

        # the histidine the first stage included is back out again
        assert 'HIS187' not in builder.active_site['residues']
        assert 'HIS191' not in builder.active_site['residues']
        assert self.coordination(builder) == before

    def metal_distance(self, builder, resid):
        """
        How far the sidechain of one residue reaches toward a metal center,
        which is what a shoehorning measures its radius by.
        """

        topology = builder.enzyme_topology
        positions = np.asarray(builder.enzyme_positions)
        residue = self.shells().sites._resolve_residue(topology, resid, None)
        metals = [
            positions[entry['index']]
            for entry in builder.binding_modes['metals']
        ]
        atoms = [
            atom for atom in residue.atoms()
            if atom.name not in util.BACKBONE_ATOM_NAMES and atom.name != 'CA'
            and atom.element is not None and atom.element.symbol != 'H'
        ]

        return min(
            float(np.linalg.norm(positions[atom.index] - metal))
            for atom in atoms for metal in metals)

    def test_a_pruned_template_matches_the_site_it_was_built_from(
            self, tmp_path):
        """
        A template is described by its force field's bonds and a query site by
        its connectivity matrix, so the pruning has to reach both or the two
        descriptions of one structure disagree: the bridging residue comes out
        with degree two on the query side and degree one on the template side,
        and template_connectivity then forces a bond back onto the very
        structure the template was built from.
        """

        pytest.importorskip('openmm')

        # an all-ones Hessian fits every metal term to zero, which leaves the
        # fixture's own GLU73 bridge to be decided on its geometry and pruned
        prepared = self.prepared()
        active_site = prepared['active_site']
        n_atoms = active_site['molecule'].number_of_atoms()
        hessian = np.ones((3 * n_atoms, 3 * n_atoms))

        forcefield = self.fitted_forcefield(active_site,
                                            hessian,
                                            prune_weak_bridge_bonds=True)

        detected = util.connectivity_bonds(
            active_site['connectivity_matrix'])

        assert len(forcefield.bonds) < len(detected)

        folder = tmp_path / 'site'
        folder.mkdir(parents=True, exist_ok=True)
        util.save_forcefield(folder / util.FORCEFIELD_FILE, forcefield)
        active_site['molecule'].write_xyz_file(
            str(folder / util.GEOMETRY_FILE))

        manager = self.manager()
        manager.load_template_from_folder(folder)
        template = manager.templates['site']

        builder = self.active_site_builder(prune_weak_bridge_bonds=True,
                                           do_qm_optimization=False)
        builder.build_forcefield(hessian=hessian)

        manager.compare_active_site(builder, mm_opt=False)
        entry = manager._comparison['templates']['site']

        # the same coordination on both sides, so nothing is forced
        assert entry['status'] == 'measured'
        assert entry['metal_bonds']['template_only'] == 0
        assert entry['metal_bonds']['query_only'] == 0

        _, changes, _ = manager._loader.build_forcefield_from_template(
            template, entry['mapping'], manager._comparison['active_site'],
            manager.builder._sites())

        assert changes == {'added': [], 'removed': []}

import h5py
import numpy as np
import pytest
from mpi4py import MPI

from veloxchem.molecularorbitals import MolecularOrbitals, molorb


class TestMolecularOrbitals:

    @staticmethod
    def make_restricted():

        orbitals = [np.array([[1.0, 0.0], [0.0, 1.0]])]
        energies = [np.array([-0.8, 0.3])]
        occupations = [np.array([2.0, 0.0])]

        return MolecularOrbitals(orbitals, energies, occupations, molorb.rest)

    @staticmethod
    def make_unrestricted():

        orbitals = [
            np.array([[1.0, 0.0], [0.0, 1.0]]),
            np.array([[0.0, 1.0], [1.0, 0.0]]),
        ]
        energies = [np.array([-0.7, 0.2]), np.array([-0.6, 0.4])]
        occupations = [np.array([1.0, 0.0]), np.array([1.0, 0.0])]

        return MolecularOrbitals(orbitals, energies, occupations, molorb.unrest)

    @staticmethod
    def make_restricted_openshell():

        orbitals = [np.array([[0.8, 0.6], [0.6, -0.8]])]
        energies = [np.array([-0.5, 0.1])]
        occupations = [np.array([1.0, 0.0]), np.array([0.0, 0.0])]

        return MolecularOrbitals(orbitals, energies, occupations,
                                 molorb.restopen)

    def test_empty_initialization(self):

        mol_orbs = MolecularOrbitals()

        assert mol_orbs.is_empty()
        assert mol_orbs.get_orbitals_type() is None

    def test_constructor_and_numpy_accessors_copy_inputs(self):

        orbitals = [np.array([[1.0, 0.0], [0.0, 1.0]])]
        energies = [np.array([-0.8, 0.3])]
        occupations = [np.array([2.0, 0.0])]

        mol_orbs = MolecularOrbitals(orbitals, energies, occupations,
                                     molorb.rest)

        orbitals[0][0, 0] = 9.0
        energies[0][0] = 9.0
        occupations[0][0] = 9.0

        alpha = mol_orbs.alpha_to_numpy()
        alpha[0, 0] = -1.0
        ea = mol_orbs.ea_to_numpy()
        ea[0] = -1.0
        occa = mol_orbs.occa_to_numpy()
        occa[0] = -1.0

        assert mol_orbs.number_aos() == 2
        assert mol_orbs.number_of_aos() == 2
        assert mol_orbs.number_mos() == 2
        assert mol_orbs.number_of_mos() == 2
        assert mol_orbs.alpha_to_numpy()[0, 0] == pytest.approx(1.0)
        assert mol_orbs.ea_to_numpy()[0] == pytest.approx(-0.8)
        assert mol_orbs.occa_to_numpy()[0] == pytest.approx(2.0)
        assert np.array_equal(mol_orbs.beta_to_numpy(),
                              mol_orbs.alpha_to_numpy())
        assert np.array_equal(mol_orbs.eb_to_numpy(), mol_orbs.ea_to_numpy())
        assert np.array_equal(mol_orbs.occb_to_numpy(),
                              mol_orbs.occa_to_numpy())

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    def test_constructor_rejects_invalid_inputs(self):

        with pytest.raises(AssertionError, match='Invalid orbitals type'):
            MolecularOrbitals([np.eye(2)], [np.zeros(2)], [np.zeros(2)], 'rest')

        with pytest.raises(
                AssertionError,
                match='Inconsistent orbitals, energies or occupation numbers'):
            MolecularOrbitals([np.eye(2)], [np.zeros(2)],
                              [np.zeros(2), np.zeros(2)], molorb.unrest)

    @pytest.mark.parametrize(
        'factory_name,scf_type,expected',
        [
            (
                'make_restricted',
                'restricted',
                (np.array([[4.0, 0.0], [0.0, 0.0]]),),
            ),
            (
                'make_unrestricted',
                'unrestricted',
                (
                    np.array([[1.0, 0.0], [0.0, 0.0]]),
                    np.array([[0.0, 0.0], [0.0, 1.0]]),
                ),
            ),
            (
                'make_restricted_openshell',
                'restricted_openshell',
                (
                    np.array([[0.64, 0.48], [0.48, 0.36]]),
                    np.zeros((2, 2)),
                ),
            ),
        ],
    )
    def test_get_density(self, factory_name, scf_type, expected):

        mol_orbs = getattr(self, factory_name)()
        density = mol_orbs.get_density(None, scf_type)

        assert len(density) == len(expected)
        for dens, ref in zip(density, expected):
            assert np.allclose(dens, ref)

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    def test_get_density_rejects_mismatched_scf_type(self):

        mol_orbs = self.make_restricted()

        with pytest.raises(AssertionError,
                           match='Invalid molecular orbitals type'):
            mol_orbs.get_density(None, 'unrestricted')

    @pytest.mark.parametrize(
        'factory_name,expected_type',
        [
            ('make_restricted', molorb.rest),
            ('make_unrestricted', molorb.unrest),
            ('make_restricted_openshell', molorb.restopen),
        ],
    )
    def test_hdf5_roundtrip_and_label_validity(self, tmp_path, factory_name,
                                               expected_type):

        filename = tmp_path / 'molecular_orbitals.h5'
        label = 'state1_'
        mol_orbs = getattr(self, factory_name)()

        mol_orbs.write_hdf5(str(filename), label=label)

        assert MolecularOrbitals.check_label_validity(str(filename),
                                                      label=label)
        assert not MolecularOrbitals.check_label_validity(str(filename),
                                                          label='missing_')

        loaded = MolecularOrbitals.read_hdf5(str(filename), label=label)

        assert loaded.get_orbitals_type() == expected_type
        assert np.array_equal(loaded.alpha_to_numpy(),
                              mol_orbs.alpha_to_numpy())
        assert np.array_equal(loaded.ea_to_numpy(), mol_orbs.ea_to_numpy())
        assert np.array_equal(loaded.occa_to_numpy(), mol_orbs.occa_to_numpy())

        if expected_type == molorb.unrest:
            assert np.array_equal(loaded.beta_to_numpy(),
                                  mol_orbs.beta_to_numpy())
            assert np.array_equal(loaded.eb_to_numpy(), mol_orbs.eb_to_numpy())

        if expected_type != molorb.rest:
            assert np.array_equal(loaded.occb_to_numpy(),
                                  mol_orbs.occb_to_numpy())

    def test_match_hdf5_uses_metadata_and_scf_type(self, tmp_path):

        filename = tmp_path / 'molecular_orbitals.h5'
        mol_orbs = self.make_unrestricted()
        mol_orbs.write_hdf5(str(filename), label='mo_test_')

        with h5py.File(filename, 'a') as handle:
            handle.create_dataset('nuclear_charges', data=np.array([8, 1, 1]))
            handle.create_dataset('basis_set', data=np.array([b'sto-3g']))

        assert MolecularOrbitals.match_hdf5(str(filename),
                                            np.array([8, 1, 1]),
                                            'STO-3G',
                                            'unrestricted',
                                            label='mo_test_')
        assert not MolecularOrbitals.match_hdf5(str(filename),
                                                np.array([8, 1]),
                                                'STO-3G',
                                                'unrestricted',
                                                label='mo_test_')
        assert not MolecularOrbitals.match_hdf5(str(filename),
                                                np.array([8, 1, 1]),
                                                'def2-svp',
                                                'unrestricted',
                                                label='mo_test_')
        assert not MolecularOrbitals.match_hdf5(str(filename),
                                                np.array([8, 1, 1]),
                                                'STO-3G',
                                                'restricted',
                                                label='mo_test_')

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    def test_read_hdf5_requires_complete_unrestricted_datasets(self, tmp_path):

        filename = tmp_path / 'broken_molecular_orbitals.h5'

        with h5py.File(filename, 'w') as handle:
            handle.create_dataset('C_alpha', data=np.eye(2))
            handle.create_dataset('E_alpha', data=np.zeros(2))
            handle.create_dataset('occ_alpha', data=np.zeros(2))
            handle.create_dataset('C_beta', data=np.eye(2))
            handle.create_dataset('scf_type', data=np.bytes_(['unrestricted']))

        with pytest.raises(AssertionError, match='E_beta not found'):
            MolecularOrbitals.read_hdf5(str(filename))

    def test_create_nto_and_is_nto(self):

        nto_orbitals = [np.eye(4)]
        nto_lambdas = [np.array([-0.8, -0.2, 0.2, 0.8])]

        nto_mol_orbs = MolecularOrbitals.create_nto(nto_orbitals, nto_lambdas,
                                                    molorb.rest)

        assert nto_mol_orbs.get_orbitals_type() == molorb.rest
        assert np.array_equal(nto_mol_orbs.ea_to_numpy(), np.zeros(4))
        assert np.array_equal(nto_mol_orbs.occa_to_numpy(), nto_lambdas[0])
        assert nto_mol_orbs.is_nto()

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    def test_create_nto_rejects_inconsistent_lambdas(self):

        with pytest.raises(AssertionError,
                           match='Inconsistent number of lambda values'):
            MolecularOrbitals.create_nto([np.eye(3)],
                                         [np.array([-0.4, 0.1, 0.4])],
                                         molorb.rest)

    def test_broadcast_with_comm_self(self):

        mol_orbs = self.make_restricted_openshell()
        broadcasted = mol_orbs.broadcast(MPI.COMM_SELF, root=0)

        assert broadcasted is not mol_orbs
        assert broadcasted.get_orbitals_type() == molorb.restopen
        assert np.array_equal(broadcasted.alpha_to_numpy(),
                              mol_orbs.alpha_to_numpy())
        assert np.array_equal(broadcasted.occb_to_numpy(),
                              mol_orbs.occb_to_numpy())

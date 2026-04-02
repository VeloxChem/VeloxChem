from mpi4py import MPI
from copy import deepcopy
import numpy as np
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.scfunrestdriver import ScfUnrestrictedDriver
from veloxchem.outputstream import OutputStream
from veloxchem.firstorderprop import FirstOrderProperties
from veloxchem.firstorderpropdriver import FirstOrderPropertyDriver


@pytest.mark.solvers
class TestFirstOrderPropertyDriver:

    def run_scf(self, mol, bas, xcfun_label, ref_dipole):

        if mol.get_multiplicity() == 1:
            scf_drv = ScfRestrictedDriver()
        else:
            scf_drv = ScfUnrestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.xcfun = xcfun_label
        scf_results = scf_drv.compute(mol, bas)

        prop_drv = FirstOrderPropertyDriver()
        prop_drv.property = 'electric dipole moment'
        prop = prop_drv.compute(mol, bas, scf_results)

        dipole_dict = prop['electric dipole moment']

        assert np.max(np.abs(dipole_dict['nuclear'])) < 1e-10
        assert np.max(np.abs(dipole_dict['total'] - ref_dipole)) < 1e-6

        return scf_results, dipole_dict

    def test_prop_rest(self):

        xyz_string = """3
        xyz
        O   -0.1858140  -1.1749469   0.7662596
        H   -0.1285513  -0.8984365   1.6808606
        H   -0.0582782  -0.3702550   0.2638279
        """
        basis_label = 'def2-svp'

        mol = Molecule.read_xyz_string(xyz_string)
        bas = MolecularBasis.read(mol, basis_label, ostream=None)

        ref_dipole = np.array([0.13230065, 0.77405262, 0.29507959])

        self.run_scf(mol, bas, 'hf', ref_dipole)

    def test_prop_unrest(self):

        xyz_string = """3
        xyz
        O   -0.1858140  -1.1749469   0.7662596
        H   -0.1285513  -0.8984365   1.6808606
        H   -0.0582782  -0.3702550   0.2638279
        """
        basis_label = 'def2-svp'

        mol = Molecule.read_xyz_string(xyz_string)
        mol.set_multiplicity(3)
        bas = MolecularBasis.read(mol, basis_label, ostream=None)

        ref_dipole = np.array([-0.0441419, -0.25826667, -0.09845354])

        self.run_scf(mol, bas, 'b3lyp', ref_dipole)

    def test_prop_rest_with_ecp(self):

        xyz_string = """2
        xyz
        Au 0 0 0
        Cl 0 0 2.3
        """
        basis_label = 'def2-svp'

        mol = Molecule.read_xyz_string(xyz_string)
        bas = MolecularBasis.read(mol, basis_label, ostream=None)

        ref_dipole = np.array([0.0, 0.0, -1.32557130])

        self.run_scf(mol, bas, 'b3lyp', ref_dipole)

    def test_prop_unrest_with_ecp(self):

        xyz_string = """2
        xyz
        Au 0 0 0
        H  0 0 1.55
        """
        basis_label = 'def2-svp'

        mol = Molecule.read_xyz_string(xyz_string)
        mol.set_charge(1)
        mol.set_multiplicity(2)
        bas = MolecularBasis.read(mol, basis_label, ostream=None)

        ref_dipole = np.array([0.0, 0.0, -0.315138231])

        self.run_scf(mol, bas, 'hf', ref_dipole)

    def test_prop_with_underscore_alias_and_multistate_density(self):

        xyz_string = """3
        xyz
        O   -0.1858140  -1.1749469   0.7662596
        H   -0.1285513  -0.8984365   1.6808606
        H   -0.0582782  -0.3702550   0.2638279
        """
        basis_label = 'def2-svp'

        mol = Molecule.read_xyz_string(xyz_string)
        bas = MolecularBasis.read(mol, basis_label, ostream=None)

        scf_results, ref_dipole = self.run_scf(
            mol, bas, 'hf', np.array([0.13230065, 0.77405262, 0.29507959]))

        prop_drv = FirstOrderPropertyDriver()
        prop_drv.property = 'electric_dipole_moment'

        prop = prop_drv.compute(mol, bas, scf_results)
        dipole_dict = prop['electric_dipole_moment']

        assert np.allclose(dipole_dict['total'], ref_dipole['total'])
        assert np.allclose(prop['electric dipole moment']['total'],
                           dipole_dict['total'])

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            total_density = np.stack([
                scf_results['D_alpha'] + scf_results['D_beta'],
                scf_results['D_alpha'] + scf_results['D_beta'],
            ])
        else:
            total_density = None

        multistate_dict = prop_drv.compute_electric_dipole_moment(
            mol, bas, total_density)

        assert multistate_dict['total'].shape == (2, 3)
        assert multistate_dict['electronic'].shape == (2, 3)
        assert multistate_dict['nuclear'].shape == (3,)
        assert np.allclose(multistate_dict['total'][0], dipole_dict['total'])
        assert np.allclose(multistate_dict['total'][1], dipole_dict['total'])

        prop_drv_copy = deepcopy(prop_drv)

        assert prop_drv_copy is not prop_drv
        assert prop_drv_copy.property == prop_drv.property
        assert prop_drv_copy.comm is prop_drv.comm
        assert prop_drv_copy.ostream is prop_drv.ostream

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    def test_prop_raises_for_unsupported_property(self):

        prop_drv = FirstOrderPropertyDriver()
        prop_drv.property = 'quadrupole moment'

        with pytest.raises(
                AssertionError,
                match='Property quadrupole moment not yet supported'):
            prop_drv.compute(None, None, None)

    def test_first_order_properties_compute_print_and_deepcopy(self, tmp_path):

        xyz_string = """3
        xyz
        O   -0.1858140  -1.1749469   0.7662596
        H   -0.1285513  -0.8984365   1.6808606
        H   -0.0582782  -0.3702550   0.2638279
        """
        basis_label = 'def2-svp'

        mol = Molecule.read_xyz_string(xyz_string)
        bas = MolecularBasis.read(mol, basis_label, ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.xcfun = 'hf'
        scf_results = scf_drv.compute(mol, bas)

        ostream = OutputStream(tmp_path / 'first_order_prop.out')
        scf_prop = FirstOrderProperties(ostream=ostream)
        scf_prop.compute_scf_prop(mol, bas, scf_results)

        dipole_moment = scf_prop.get_property('dipole moment')
        assert np.allclose(dipole_moment,
                           scf_prop.get_property('dipole_moment'))

        scf_prop.print_properties(mol, title='Ground State Test')

        excited_state_dipoles = np.vstack([dipole_moment, 2.0 * dipole_moment])
        scf_prop.properties['dipole moment'] = excited_state_dipoles
        scf_prop.properties['dipole_moment'] = excited_state_dipoles
        mol.set_charge(1)
        scf_prop.print_properties(mol,
                                  title='Excited State Test',
                                  states=[1, 2])
        ostream.close()

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            output_text = (tmp_path / 'first_order_prop.out').read_text()

            assert 'Ground State Test' in output_text
            assert 'Excited State 1' in output_text
            assert 'Excited State 2' in output_text
            assert 'Center of nuclear charge is chosen as the origin.' in output_text
            assert 'Debye' in output_text

        scf_prop_copy = deepcopy(scf_prop)

        assert scf_prop_copy is not scf_prop
        assert np.allclose(scf_prop_copy.get_property('dipole moment'),
                           scf_prop.get_property('dipole moment'))
        assert scf_prop_copy.comm is scf_prop.comm
        assert scf_prop_copy.ostream is scf_prop.ostream

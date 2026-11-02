from mpi4py import MPI
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.cubicresponsedriver import CubicResponseDriver
from veloxchem.thgdriver import ThgDriver
from veloxchem.outputstream import OutputStream


class TestthgFromCrf:

    def run_thg_from_crf(self, xcfun_label, *, ri):

        molecule_string = """
            O   0.0   0.0   0.0
            H   0.0   1.4   1.1
            H   0.0  -1.4   1.1
        """
        basis_set_label = 'def2-svp'

        molecule = Molecule.read_molecule_string(molecule_string, units='au')
        basis = MolecularBasis.read(molecule, basis_set_label, ostream=None)

        scf_settings = {}
        method_settings = {'xcfun': xcfun_label, 'grid_level': 1}

        scfdrv = ScfRestrictedDriver()
        scfdrv.ri = ri
        scfdrv.ostream.mute()
        scfdrv.update_settings(scf_settings, method_settings)
        scf_results = scfdrv.compute(molecule, basis)

        if scfdrv.rank == mpi_master():
            if ri:
                assert scf_results['ri_jk'] or scf_results['ri_coulomb']
                assert not scf_results['ri_fallback']

        w1 = 0.05
        w2 = w1
        w3 = w1
        damping = 0.1
        components = 'xyz'

        ref_thg_results = {
            'gamma': 0.0,
        }

        crf = CubicResponseDriver()
        crf.ostream.mute()
        rsp_settings = {}

        for a in components:
            for b in components:
                rsp_settings.update({
                    'b_frequencies': [w1],
                    'c_frequencies': [w2],
                    'd_frequencies': [w3],
                    'a_component': a,
                    'b_component': a,
                    'c_component': b,
                    'd_component': b,
                    'damping': damping,
                })
                crf.update_settings(rsp_settings, method_settings)
                crf_results = crf.compute(molecule, basis, scf_results)
                if MPI.COMM_WORLD.Get_rank() == mpi_master():
                    ref_thg_results['gamma'] += crf_results[('crf', w1, w2, w3)]

        for a in components:
            for b in components:
                rsp_settings.update({
                    'b_frequencies': [w1],
                    'c_frequencies': [w2],
                    'd_frequencies': [w3],
                    'a_component': a,
                    'b_component': b,
                    'c_component': a,
                    'd_component': b,
                    'damping': damping,
                })
                crf.update_settings(rsp_settings, method_settings)
                crf_results = crf.compute(molecule, basis, scf_results)
                if MPI.COMM_WORLD.Get_rank() == mpi_master():
                    ref_thg_results['gamma'] += crf_results[('crf', w1, w2, w3)]

        for a in components:
            for b in components:
                rsp_settings.update({
                    'b_frequencies': [w1],
                    'c_frequencies': [w2],
                    'd_frequencies': [w3],
                    'a_component': a,
                    'b_component': b,
                    'c_component': b,
                    'd_component': a,
                    'damping': damping,
                })
                crf.update_settings(rsp_settings, method_settings)
                crf_results = crf.compute(molecule, basis, scf_results)
                if MPI.COMM_WORLD.Get_rank() == mpi_master():
                    ref_thg_results['gamma'] += crf_results[('crf', w1, w2, w3)]

        rsp_settings = {
            'frequencies': [w1],
            'damping': damping,
        }

        thg = ThgDriver(MPI.COMM_WORLD,OutputStream())
        thg.ostream.mute()
        thg.update_settings(rsp_settings, method_settings)
        thg_results = thg.compute(molecule, basis, scf_results)

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            for key in ref_thg_results:
                ref_thg_results[key] /= 15.0
                ref_thg_results[key] *= -1.0  # rsp func. -> gamma

            tol = 1.0e-5

            for key in ref_thg_results:
                ref_val = ref_thg_results[key]
                calc_val = thg_results['THG'][(w1, w2, w3)]
                assert abs(abs(calc_val.real / ref_val.real) - 1.0) < tol

    @pytest.mark.timeconsuming
    @pytest.mark.parametrize('ri', [False, True])
    def test_thg_from_crf_lda(self, ri):

        self.run_thg_from_crf('slda', ri=ri)

    @pytest.mark.solvers
    @pytest.mark.parametrize('ri', [False, True])
    def test_thg_from_crf_gga(self, ri):

        self.run_thg_from_crf('pbe0', ri=ri)

    @pytest.mark.timeconsuming
    @pytest.mark.parametrize('ri', [False, True])
    def test_thg_from_crf_mgga(self, ri):

        self.run_thg_from_crf('m06-2x', ri=ri)

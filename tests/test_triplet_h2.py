from pathlib import Path
import numpy as np
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.mpitask import MpiTask
from veloxchem.scfunrestdriver import ScfUnrestrictedDriver
from veloxchem.firstorderprop import FirstOrderProperties


@pytest.mark.solvers
class TestTripletH2:

    def run_scf(self, inpfile, potfile, xcfun_label, ri, ref_e_scf, ref_dip):

        task = MpiTask([inpfile, None])
        task.input_dict['scf']['checkpoint_file'] = None

        if potfile is not None:
            task.input_dict['method_settings']['potfile'] = potfile

        if xcfun_label is not None:
            task.input_dict['method_settings']['xcfun'] = xcfun_label

        scf_drv = ScfUnrestrictedDriver(task.mpi_comm, task.ostream)
        scf_drv.ri = ri
        scf_drv.update_settings(task.input_dict['scf'],
                                task.input_dict['method_settings'])
        scf_results = scf_drv.compute(task.molecule, task.ao_basis,
                                      task.min_basis)

        scf_prop = FirstOrderProperties(task.mpi_comm, task.ostream)
        scf_prop.compute_scf_prop(task.molecule, task.ao_basis, scf_results)

        if task.mpi_rank == mpi_master():
            if ri:
                assert scf_results['ri_jk'] or scf_results['ri_coulomb']
                assert not scf_results['ri_fallback']

            e_scf = scf_drv.get_scf_energy()
            tol = 1.0e-5 if xcfun_label is not None else 1.0e-6
            assert np.max(np.abs(e_scf - ref_e_scf)) < tol

            dip = scf_prop.get_property('dipole moment')
            assert np.max(np.abs(dip - ref_dip)) < 1.0e-5

    @pytest.mark.parametrize('ri', [False, True])
    def test_uscf_hf(self, ri):

        here = Path(__file__).parent
        inpfile = str(here / 'data' / 'h2_t.inp')

        potfile = None
        xcfun_label = None

        ref_e_scf = {False: -0.77331601, True: -0.77330484}[ri]
        ref_dip = np.array([0.000000, 0.000000, 0.000000])

        self.run_scf(inpfile, potfile, xcfun_label, ri, ref_e_scf, ref_dip)

    @pytest.mark.parametrize('ri', [False, True])
    def test_uscf_lda(self, ri):

        here = Path(__file__).parent
        inpfile = str(here / 'data' / 'h2_t.inp')

        potfile = None
        xcfun_label = 'slater'

        ref_e_scf = {False: -0.71006575, True: -0.71502880}[ri]
        ref_dip = np.array([0.000000, 0.000000, 0.000000])

        self.run_scf(inpfile, potfile, xcfun_label, ri, ref_e_scf, ref_dip)

    @pytest.mark.parametrize('ri', [False, True])
    def test_uscf_gga(self, ri):

        here = Path(__file__).parent
        inpfile = str(here / 'data' / 'h2_t.inp')

        potfile = None
        xcfun_label = 'b3lyp'

        ref_e_scf = {False: -0.78739868, True: -0.78858776}[ri]
        ref_dip = np.array([0.000000, 0.000000, 0.000000])

        self.run_scf(inpfile, potfile, xcfun_label, ri, ref_e_scf, ref_dip)

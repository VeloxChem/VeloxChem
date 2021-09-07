from pathlib import Path
from pytest import approx
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.mpitask import MpiTask
from veloxchem.scfrestopendriver import ScfRestrictedOpenDriver
from veloxchem.molecularbasis import MolecularBasis


class TestScfRestrictedOpenShell:

    def run_scf(self, inpfile, potfile, xcfun_label, ref_e_scf):

        task = MpiTask([inpfile, None])
        task.input_dict['scf']['checkpoint_file'] = None
        task.ao_basis = MolecularBasis.read(task.molecule, 'def2-svp', '.',
                                            task.ostream)

        if potfile is not None:
            task.input_dict['method_settings']['potfile'] = potfile

        if xcfun_label is not None:
            task.input_dict['method_settings']['xcfun'] = xcfun_label

        scf_drv = ScfRestrictedOpenDriver(task.mpi_comm, task.ostream)
        scf_drv.update_settings(task.input_dict['scf'],
                                task.input_dict['method_settings'])

        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)

        if task.mpi_rank == mpi_master():
            e_scf = scf_drv.get_scf_energy()
            tol = 1.0e-5 if xcfun_label is not None else 1.0e-6
            assert e_scf == approx(ref_e_scf, abs=tol)

    @pytest.mark.parametrize(
        'inpfile, ref_e_scf',
        [
            ('heh.inp', -3.348471908695),
            ('li.inp', -7.425064044619),
            ('ts01.inp', -460.413199994),
            ('ts02.inp', -76.406503929),
            ('ts03.inp', -40.618499031),
        ],
        ids=['HeH', 'Li', 'ClH2', 'H3O', 'CH5'],
    )
    def test_scf_hf(self, inpfile, ref_e_scf):

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / inpfile)

        potfile = None

        xcfun_label = None

        self.run_scf(inpfile, potfile, xcfun_label, ref_e_scf)

    @pytest.mark.parametrize(
        'inpfile, xcfun_label, ref_e_scf',
        [('heh.inp', None, -3.34847190869),
         ('heh.inp', 'slater', -3.166481549682),
         ('heh.inp', 'b3lyp', -3.404225946805)],
        ids=['HeH-ROHF', 'HeH-Slater', 'HeH-B3LYP'],
    )
    def test_scf_dft(self, xcfun_label, inpfile, ref_e_scf):

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / inpfile)

        potfile = None

        self.run_scf(inpfile, potfile, xcfun_label, ref_e_scf)

from pathlib import Path
from pytest import approx
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.mpitask import MpiTask
from veloxchem.scfrestopendriver import ScfRestrictedOpenDriver
from veloxchem.molecularbasis import MolecularBasis


class TestScfRestrictedOpenShell:

    def run_scf(self, inpfile, potfile, xcfun_label, basis_label, ref_e_scf):

        task = MpiTask([inpfile, None])
        task.input_dict['scf']['checkpoint_file'] = None
        task.ao_basis = MolecularBasis.read(task.molecule, basis_label)

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
        'inpfile, xcfun_label, basis_label, ref_e_scf',
        [
            ('li.inp', None, 'def2-svp', -7.425064044619),
            ('ts01.inp', None, 'def2-svp', -460.413199994),
            ('ts02.inp', None, 'def2-svp', -76.406503929),
            ('ts03.inp', None, 'def2-svp', -40.618499031),
            ('heh.inp', None, 'def2-svp', -3.34847190869),
            ('heh.inp', 'slater', 'def2-svp', -3.166481549682),
            ('heh.inp', 'b3lyp', 'def2-svp', -3.404225946805),
            ('water_triplet.inp', None, 'def2-svp', -75.701357737185),
            ('water_triplet.inp', 'slater', 'def2-svp', -74.8660271496),
            ('water_triplet.inp', 'slda', 'def2-svp', -75.6935516775),
            ('water_triplet.inp', 'blyp', 'def2-svp', -76.055256325587),
            ('water_triplet.inp', 'b3lyp', 'def2-svp', -76.074465451578),
            ('water_triplet.inp', 'b3lyp', 'aug-cc-pvdz', -76.176630915242),
        ],
        ids=[
            'Li', 'ClH2', 'H3O', 'CH5', 'HeH', 'HeH-Slater', 'HeH-B3LYP',
            'H2O-HF', 'H2O-SLATER', 'H2O-SLDA', 'H2O-BLYP', 'H2O-B3LYP',
            'H2O-B3LYP-aDZ'
        ],
    )
    def test_scf(self, inpfile, xcfun_label, basis_label, ref_e_scf):

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / inpfile)

        potfile = None

        self.run_scf(inpfile, potfile, xcfun_label, basis_label, ref_e_scf)

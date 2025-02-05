from pathlib import Path
import numpy as np
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.veloxchemlib import bohr_in_angstrom
from veloxchem.mpitask import MpiTask
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.scfgradientdriver import ScfGradientDriver
from veloxchem.optimizationdriver import OptimizationDriver


@pytest.mark.filterwarnings(
    'ignore:.*tostring.*tobytes:DeprecationWarning:geometric')
@pytest.mark.solvers
class TestOptimizeSCF:

    def run_opt(self, inpfile, basis_label, ref_coords):

        task = MpiTask([inpfile, None])
        task.input_dict['scf']['checkpoint_file'] = None
        task.input_dict['method_settings']['basis'] = basis_label.upper()

        if task.mpi_rank == mpi_master():
            task.ao_basis = MolecularBasis.read(task.molecule,
                                                basis_label.upper(),
                                                ostream=None)
        else:
            task.ao_basis = None
        task.ao_basis = task.mpi_comm.bcast(task.ao_basis, root=mpi_master())

        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)
        scf_drv.update_settings(task.input_dict['scf'],
                                task.input_dict['method_settings'])
        scf_results = scf_drv.compute(task.molecule, task.ao_basis,
                                      task.min_basis)

        grad_drv = ScfGradientDriver(scf_drv)
        opt_drv = OptimizationDriver(grad_drv)
        opt_drv.update_settings({
            'coordsys': 'tric',
            'filename': task.input_dict['filename'],
        })
        opt_results = opt_drv.compute(task.molecule, task.ao_basis, scf_results)

        if task.mpi_rank == mpi_master():
            opt_mol = Molecule.read_xyz_string(opt_results['final_geometry'])
            opt_coords = opt_mol.get_coordinates_in_bohr()
            assert np.max(np.abs(opt_coords - ref_coords)) < 1.0e-5

            inpfile = Path(inpfile)

            optfile = Path(str(inpfile.with_name(inpfile.stem)) + '_optim.xyz')
            if optfile.is_file():
                optfile.unlink()

            logfile = inpfile.with_suffix('.log')
            if logfile.is_file():
                logfile.unlink()

            scffile = Path(str(inpfile.with_name(inpfile.stem)) + '_scf.h5')
            if scffile.is_file():
                scffile.unlink()

            scffile = Path(
                str(inpfile.with_name(inpfile.stem)) + '_scf.results.h5')
            if scffile.is_file():
                scffile.unlink()

    def test_nh3(self):

        here = Path(__file__).parent
        inpfile = str(here / 'data' / 'nh3.inp')

        basis_label = 'sto-3g'

        ref_coords = np.array([
            [-1.941114808233, 1.588257603470, -0.020632249491],
            [-1.965745778628, 2.619171893222, 0.031449914839],
            [-2.496768727202, 1.280798298868, 0.793505674044],
            [-2.536814090843, 1.362016810479, -0.833073376431],
        ]) / bohr_in_angstrom()

        self.run_opt(inpfile, basis_label, ref_coords)

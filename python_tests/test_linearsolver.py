from mpi4py import MPI
from datetime import datetime, timedelta
import numpy as np

from veloxchem.veloxchemlib import denmat
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.outputstream import OutputStream
from veloxchem.aodensitymatrix import AODensityMatrix
from veloxchem.aofockmatrix import AOFockMatrix
from veloxchem.linearsolver import LinearSolver


class TestLinearSolver:

    def test_update_settings(self):

        rsp_dict = {
            'eri_thresh': 1e-13,
            'qq_type': 'QQ',
            'batch_size': 99,
            'conv_thresh': 1e-7,
            'max_iter': 199,
            'lindep_thresh': 1e-10,
            #'restart': False,
            'checkpoint_file': 'mycheckpoint.h5',
            'timing': True,
            'profiling': True,
            'memory_profiling': True,
            'memory_tracing': True,
        }

        method_dict = {
            'grid_level': 5,
            'electric_field': (0, -0.002, 0.001),
            'use_split_comm': True,
        }

        lin_solver = LinearSolver(MPI.COMM_WORLD, OutputStream(None))

        for key, val in rsp_dict.items():
            assert getattr(lin_solver, key) != val
        for key, val in method_dict.items():
            assert getattr(lin_solver, key) != val

        lin_solver.update_settings(rsp_dict, method_dict)

        for key, val in rsp_dict.items():
            assert getattr(lin_solver, key) == val
        for key, val in method_dict.items():
            assert getattr(lin_solver, key) == val

    def get_molecule_and_basis(self):

        mol_str = """
            O   0.0   0.0   0.0
            H   0.0   1.4   1.1
            H   0.0  -1.4   1.1
        """
        mol = Molecule.read_molecule_string(mol_str, units='bohr')
        bas = MolecularBasis.read(mol, 'aug-cc-pvdz', ostream=None)

        return mol, bas

    def test_comp_fock_split_comm(self):

        comm = MPI.COMM_WORLD
        ostream = OutputStream(None)

        mol, bas = self.get_molecule_and_basis()
        nao = bas.get_dimension_of_basis(mol)

        dmat = np.diag(np.ones(nao))
        dens = AODensityMatrix([dmat], denmat.rest)
        fock = AOFockMatrix(dens)

        eri_dict = {'screening': None}
        dft_dict = {'molgrid': None, 'gs_density': None}
        pe_dict = {'V_es': None, 'pe_drv': None}

        solver = LinearSolver(comm, ostream)
        solver._comp_lr_fock_split_comm(fock, dens, mol, bas, eri_dict,
                                        dft_dict, pe_dict)

        assert fock.alpha_to_numpy(0).shape == dmat.shape

    def test_need_graceful_exit(self):

        comm = MPI.COMM_WORLD
        ostream = OutputStream(None)

        solver = LinearSolver(comm, ostream)
        solver.program_end_time = datetime.now() + timedelta(hours=1.0)

        assert solver._need_graceful_exit(1.5)
        assert not solver._need_graceful_exit(0.5)

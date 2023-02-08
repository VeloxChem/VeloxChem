from pathlib import Path

import numpy as np
from mpi4py import MPI
from veloxchem.aodensitymatrix import AODensityMatrix
from veloxchem.aofockmatrix import AOFockMatrix
from veloxchem.firstorderprop import FirstOrderProperties
from veloxchem.mpitask import MpiTask
from veloxchem.outputstream import OutputStream
from veloxchem.qqscheme import get_qq_scheme
from veloxchem.scfdriver import ScfDriver
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.veloxchemlib import (ElectronRepulsionIntegralsDriver, denmat,
                                    is_mpi_master)

from .addons import using_cppe


class TestScfRestricted:

    def run_scf(self, inpfile, potfile, xcfun_label, efield, ref_e_scf,
                ref_dip):

        task = MpiTask([inpfile, None])
        task.input_dict['scf']['checkpoint_file'] = None

        if potfile is not None:
            task.input_dict['method_settings']['potfile'] = potfile
        elif efield is not None:
            task.input_dict['method_settings']['electric_field'] = efield

        if xcfun_label is not None:
            task.input_dict['method_settings']['xcfun'] = xcfun_label

        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)
        scf_drv.update_settings(task.input_dict['scf'],
                                task.input_dict['method_settings'])
        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)

        scf_prop = FirstOrderProperties(task.mpi_comm, task.ostream)
        scf_prop.compute_scf_prop(task.molecule, task.ao_basis,
                                  scf_drv.scf_tensors)

        if is_mpi_master(task.mpi_comm):
            if xcfun_label is not None:
                ene_tol, dip_tol = 1.0e-4, 1.0e-4
            else:
                ene_tol, dip_tol = 1.0e-6, 1.0e-5

            e_scf = scf_drv.get_scf_energy()
            assert np.max(np.abs(e_scf - ref_e_scf)) < ene_tol

            dip = scf_prop.get_property('dipole moment')
            assert np.max(np.abs(dip - ref_dip)) < dip_tol

    def test_scf_hf(self):

        # vlxtag: RHF, Energy

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'water.inp')
        potfile = None

        xcfun_label = None
        electric_field = None
        ref_e_scf = -76.041697549811
        ref_dip = np.array([0.000000, 0.000000, 0.786770])

        self.run_scf(inpfile, potfile, xcfun_label, electric_field, ref_e_scf,
                     ref_dip)

        xcfun_label = None
        electric_field = '0, 0, 0.03'
        ref_e_scf = -76.0688447658
        ref_dip = np.array([0.000000, 0.000000, 1.023625])

        self.run_scf(inpfile, potfile, xcfun_label, electric_field, ref_e_scf,
                     ref_dip)

    def test_scf_dft(self):

        # vlxtag: RKS, Energy

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'water.inp')
        potfile = None

        xcfun_label = 'slda'
        electric_field = None
        ref_e_scf = -76.074208655329
        ref_dip = np.array([0.000000, 0.000000, 0.731288])

        self.run_scf(inpfile, potfile, xcfun_label, electric_field, ref_e_scf,
                     ref_dip)

        xcfun_label = 'slda'
        electric_field = '0.001, 0, 0'
        ref_e_scf = -76.0742132904
        ref_dip = np.array([0.009288, 0.000000, 0.731345])

        self.run_scf(inpfile, potfile, xcfun_label, electric_field, ref_e_scf,
                     ref_dip)

        xcfun_label = 'b3lyp'
        electric_field = None
        ref_e_scf = -76.443546270423
        ref_dip = np.array([0.000000, 0.000000, 0.731257])

        self.run_scf(inpfile, potfile, xcfun_label, electric_field, ref_e_scf,
                     ref_dip)

        xcfun_label = 'tpssh'
        electric_field = None
        ref_e_scf = -76.436760570286
        ref_dip = np.array([0.000000, 0.000000, 0.723071])

        self.run_scf(inpfile, potfile, xcfun_label, electric_field, ref_e_scf,
                     ref_dip)

        xcfun_label = 'm06'
        electric_field = None
        ref_e_scf = -76.405960204758
        ref_dip = np.array([0.000000, 0.000000, 0.748155])
        # note: reference dipole moment with finer grid is 0.748301

        self.run_scf(inpfile, potfile, xcfun_label, electric_field, ref_e_scf,
                     ref_dip)

    @using_cppe
    def test_scf_hf_pe(self):

        # vlxtag: RHF, Energy, PE

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'pe_water.inp')
        potfile = str(here / 'inputs' / 'pe_water.pot')

        xcfun_label = None
        electric_field = None

        ref_e_scf = -76.067159426565
        ref_dip = np.array([-0.039715, 0.098785, 0.945128])

        self.run_scf(inpfile, potfile, xcfun_label, electric_field, ref_e_scf,
                     ref_dip)

    @using_cppe
    def test_scf_dft_pe(self):

        # vlxtag: RKS, Energy, PE

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'pe_water.inp')
        potfile = str(here / 'inputs' / 'pe_water.pot')

        xcfun_label = 'b3lyp'
        electric_field = None

        ref_e_scf = -76.4687342577
        ref_dip = np.array([-0.048191, 0.098718, 0.902823])

        self.run_scf(inpfile, potfile, xcfun_label, electric_field, ref_e_scf,
                     ref_dip)

    def test_comp_fock_split_comm(self):

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'water.inp')
        task = MpiTask([inpfile, None])

        mol = task.molecule
        bas = task.ao_basis
        nao = bas.get_dimensions_of_basis(mol)

        dmat = np.diag(np.ones(nao))
        dens = AODensityMatrix([dmat], denmat.rest)
        fock = AOFockMatrix(dens)

        eri_drv = ElectronRepulsionIntegralsDriver(task.mpi_comm)
        screening = eri_drv.compute(get_qq_scheme('QQ_DEN'), 1.0e-12, mol, bas)

        solver = ScfDriver(task.mpi_comm, task.ostream)
        solver._comp_2e_fock_split_comm(fock, dens, mol, bas, screening)

        assert fock.alpha_to_numpy(0).shape == dmat.shape

    def test_update_settings(self):

        scf_dict = {
            'acc_type': 'DIIS',
            'max_iter': 199,
            'conv_thresh': 1e-7,
            'qq_type': 'QQ',
            'eri_thresh': 1e-13,
            'restart': False,
            'checkpoint_file': 'mycheckpoint.h5',
            'timing': True,
            'profiling': True,
            'memory_profiling': True,
            'memory_tracing': True,
        }

        method_dict = {
            'dispersion': True,
            'grid_level': 5,
            'electric_field': (0, -0.002, 0.001),
            'use_split_comm': True,
        }

        scf_drv = ScfDriver(MPI.COMM_WORLD, OutputStream(None))

        for key, val in scf_dict.items():
            assert getattr(scf_drv, key) != val
        for key, val in method_dict.items():
            assert getattr(scf_drv, key) != val

        scf_drv.update_settings(scf_dict, method_dict)

        for key, val in scf_dict.items():
            assert getattr(scf_drv, key) == val
        for key, val in method_dict.items():
            assert getattr(scf_drv, key) == val

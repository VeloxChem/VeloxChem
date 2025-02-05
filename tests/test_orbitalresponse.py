from mpi4py import MPI
from pathlib import Path
import numpy as np

from veloxchem.veloxchemlib import mpi_master
from veloxchem.mpitask import MpiTask
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.tdaeigensolver import TdaEigenSolver
from veloxchem.lreigensolver import LinearResponseEigenSolver
from veloxchem.tddftorbitalresponse import TddftOrbitalResponse
from veloxchem.checkpoint import read_rsp_hdf5


class TestOrbitalResponse:

    def run_orbitalresponse(self, inpfile, potfile, xcfun_label,
                            orbrsp_ref_file, is_tda):

        task = MpiTask([inpfile, None], MPI.COMM_WORLD)
        task.input_dict['scf']['checkpoint_file'] = None
        # use a smaller basis for lambda test

        rsp_dict = {'nstates': 3}
        orbrsp_dict = dict(rsp_dict)
        orbrsp_dict['state_deriv_index'] = '1'

        if potfile is not None:
            task.input_dict['method_settings']['potfile'] = potfile

        if xcfun_label is not None:
            task.input_dict['method_settings']['xcfun'] = xcfun_label

        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)
        scf_drv.update_settings(task.input_dict['scf'],
                                task.input_dict['method_settings'])
        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)

        # Our references: lambda and omega in AO basis

        if is_tda:
            tda_solver = TdaEigenSolver(task.mpi_comm, task.ostream)
            tda_solver.update_settings(rsp_dict,
                                       task.input_dict['method_settings'])
            rsp_results = tda_solver.compute(task.molecule, task.ao_basis,
                                             scf_drv.scf_tensors)

            orb_resp = TddftOrbitalResponse(task.mpi_comm, task.ostream)
            orbrsp_dict['tamm_dancoff'] = 'yes'
            lambda_ref = 'lambda_tda'
            omega_ref = 'omega_tda'
        else:
            rpa_solver = LinearResponseEigenSolver(task.mpi_comm, task.ostream)
            rpa_solver.update_settings({'nstates': 3},
                                       task.input_dict['method_settings'])
            rsp_results = rpa_solver.compute(task.molecule, task.ao_basis,
                                             scf_drv.scf_tensors)

            orb_resp = TddftOrbitalResponse(task.mpi_comm, task.ostream)
            lambda_ref = 'lambda_rpa'
            omega_ref = 'omega_rpa'

        orb_resp.update_settings(orbrsp_dict,
                                 task.input_dict['method_settings'])

        orb_resp.compute(task.molecule, task.ao_basis, scf_drv.scf_tensors,
                         rsp_results)
        orb_resp_results = orb_resp.cphf_results
        #omega_ao = orb_resp.compute_omega(task.molecule, task.ao_basis,
        #                                  scf_drv.scf_tensors)

        dft_dict = {'dft_func_label': 'HF'}
        pe_dict = {'potfile_text': ''}

        ref_lambda_ao, ref_omega_ao = read_rsp_hdf5(orbrsp_ref_file,
                                                    [lambda_ref, omega_ref],
                                                    task.molecule,
                                                    task.ao_basis, dft_dict,
                                                    pe_dict, task.ostream)

        dist_cphf_coefficients = orb_resp_results['dist_cphf_ov']
        dof = len(dist_cphf_coefficients)

        cphf_coefficients = []
        for x in range(dof):
            solution_vec = dist_cphf_coefficients[x].get_full_vector(0)
            if scf_drv.rank == mpi_master():
                cphf_coefficients.append(solution_vec)

        if task.mpi_rank == mpi_master():
            nocc = task.molecule.number_of_alpha_electrons()
            mo = scf_drv.scf_tensors['C_alpha']
            mo_occ = mo[:, :nocc]
            mo_vir = mo[:, nocc:]
            nvir = mo_vir.shape[1]
            nao = task.ao_basis.get_dimensions_of_basis()
            #lambda_ov = orb_resp_results['cphf_ov']

            lambda_ov = np.array(cphf_coefficients)
            dof = lambda_ov.shape[0]

            lambda_ov = lambda_ov.reshape(dof, nocc, nvir)
            lambda_ao = np.einsum('mi,sia,na->smn', mo_occ, lambda_ov, mo_vir)
            lambda_ao = lambda_ao.reshape(dof, nao, nao)

            assert np.max(np.abs(lambda_ao[0] - ref_lambda_ao)) < 5.0e-4
            # TODO: uncomment once TDDFT gradients are working
            #assert np.max(np.abs(omega_ao[0] - ref_omega_ao)) < 5.0e-4

    def test_tda_hf(self):

        here = Path(__file__).parent
        inpfile = str(here / 'data' / 'water_orbrsp.inp')
        orbrsp_ref_file = str(here / 'data' / 'orbital_response_hf_ref.h5')

        potfile = None

        xcfun_label = None

        self.run_orbitalresponse(inpfile, potfile, xcfun_label, orbrsp_ref_file,
                                 True)

    def test_rpa_hf(self):

        here = Path(__file__).parent
        inpfile = str(here / 'data' / 'water_orbrsp.inp')
        orbrsp_ref_file = str(here / 'data' / 'orbital_response_hf_ref.h5')

        potfile = None

        xcfun_label = None

        self.run_orbitalresponse(inpfile, potfile, xcfun_label, orbrsp_ref_file,
                                 False)

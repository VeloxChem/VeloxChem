from pathlib import Path
import numpy as np

from veloxchem.veloxchemlib import is_mpi_master
from veloxchem.veloxchemlib import rotatory_strength_in_cgs
from veloxchem.mpitask import MpiTask
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.tdaexcidriver import TDAExciDriver
from veloxchem.lreigensolver import LinearResponseEigenSolver


class TestECD:

    def run_ecd(self, inpfile, ref, flag):

        task = MpiTask([inpfile, None])
        task.input_dict['scf']['checkpoint_file'] = None

        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)
        scf_drv.update_settings(task.input_dict['scf'],
                                task.input_dict['method_settings'])
        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)

        if flag == 'tda':
            rsp_drv = TDAExciDriver(task.mpi_comm, task.ostream)
        elif flag == 'rpa':
            rsp_drv = LinearResponseEigenSolver(task.mpi_comm, task.ostream)
        rsp_drv.update_settings({'nstates': ref['eig'].shape[0]},
                                task.input_dict['method_settings'])
        rsp_results = rsp_drv.compute(task.molecule, task.ao_basis,
                                      scf_drv.scf_tensors)

        if is_mpi_master(task.mpi_comm):
            eigvals = rsp_results['eigenvalues']
            edip = np.array(rsp_results['electric_transition_dipoles'])
            vdip = np.array(rsp_results['velocity_transition_dipoles'])
            mdip = np.array(rsp_results['magnetic_transition_dipoles'])
            osc = rsp_results['oscillator_strengths']
            rot = rsp_results['rotatory_strengths'] / rotatory_strength_in_cgs()
            rot_len = np.sum(edip * mdip, axis=1)

            for i in range(edip.shape[0]):
                if np.vdot(edip[i, :], ref['edip'][i, :]) < 0.0:
                    edip[i, :] *= -1.0
                    vdip[i, :] *= -1.0
                    mdip[i, :] *= -1.0

            assert np.max(np.abs(eigvals - ref['eig'])) < 1.0e-6
            assert np.max(np.abs(edip - ref['edip'])) < 5.0e-4
            assert np.max(np.abs(vdip - ref['vdip'])) < 5.0e-4
            assert np.max(np.abs(mdip - ref['mdip'])) < 5.0e-4
            assert np.max(np.abs(osc - ref['osc'])) < 1.0e-4
            assert np.max(np.abs(rot - ref['rot'])) < 2.0e-4
            assert np.max(np.abs(rot_len - ref['rot_len'])) < 2.0e-4

    def gen_ref(self, eigvals, edip, vdip, mdip):

        eigvals = np.array(eigvals)

        edip = np.array(edip).reshape(-1, 3)
        vdip = -1.0 * np.array(vdip).reshape(-1, 3) / np.array([eigvals] * 3).T
        mdip = 0.5 * np.array(mdip).reshape(-1, 3)

        osc = (2.0 / 3.0) * np.sum(edip**2, axis=1) * eigvals
        rot = np.sum(vdip * mdip, axis=1)
        rot_len = np.sum(edip * mdip, axis=1)

        return {
            'eig': eigvals,
            'edip': edip,
            'vdip': vdip,
            'mdip': mdip,
            'osc': osc,
            'rot': rot,
            'rot_len': rot_len,
        }

    def test_ecd_tda(self):

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'glycol.inp')

        raw_vals = """
            0.44968592
            0.45609001
            0.51132434
            0.52111191
            0.57126077
            0.57312114
        """
        eigvals = [float(x) for x in raw_vals.split()]

        raw_vals = """
             2.14210711E-03  -0.14787623       1.89694436E-02
            -0.14472438       3.43098619E-02   5.32069410E-02
            -2.00790481E-04   6.71962747E-02   9.07739438E-02
             1.06830878E-02   3.41590700E-02  -3.55199456E-02
            -0.12788423      -0.12511955       0.19163530
            -0.15155190      -2.21837514E-02  -0.15241120
        """
        edip = [float(x) for x in raw_vals.split()]

        raw_vals = """
            -1.18119397E-02  -4.43814071E-02   7.85928311E-03
            -0.16837588       1.33309556E-02   2.04168528E-02
             9.53534722E-03   1.49679248E-03  -2.94521180E-02
             6.02394035E-02  -2.71178958E-02  -1.68953427E-03
             0.10503782       0.18427200      -0.14253944
             0.11414391      -7.80124105E-04   0.18522988
        """
        vdip = [float(x) for x in raw_vals.split()]

        raw_vals = """
            -0.35406650      -0.17360791      -0.50226426
             5.47343373E-02   0.55674013      -0.33409218
            -0.19733858       3.53211932E-02  -0.29769985
             5.08628175E-02   0.38406640      -0.27950090
              1.2969232       8.29998894E-02   0.32298048
            -0.81764222       -1.3916163       0.10611006
        """
        mdip = [float(x) for x in raw_vals.split()]

        self.run_ecd(inpfile, self.gen_ref(eigvals, edip, vdip, mdip), 'tda')

    def test_ecd_rpa(self):

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'glycol.inp')

        raw_vals = """
            0.44794104
            0.45442482
            0.51029299
            0.51999153
            0.56787602
            0.56968523
        """
        eigvals = [float(x) for x in raw_vals.split()]

        raw_vals = """
            -2.53540170E-03   0.13498859      -1.71295732E-02
            -0.13392821       3.21968518E-02   4.74412325E-02
             2.01155681E-03  -5.98212548E-02  -8.57406679E-02
            -5.37689074E-04  -3.42411727E-02   3.18913526E-02
             0.12031292       0.11600085      -0.17689970
            -0.14705050      -1.82117552E-02  -0.14025928
        """
        edip = [float(x) for x in raw_vals.split()]

        raw_vals = """
             1.17555187E-02   3.47857580E-02  -3.99991422E-03
            -0.17071140       1.26146771E-02   1.54008833E-02
            -7.69595477E-03   6.45660748E-03   3.28253166E-02
            -5.38950756E-02   2.88047828E-02  -2.69781501E-03
            -0.12318731      -0.23462833       0.16179921
             0.13485363      -7.44477974E-03   0.22645969
        """
        vdip = [float(x) for x in raw_vals.split()]

        raw_vals = """
             0.35006744       0.15407816       0.48287324
             6.54157274E-02   0.60117578      -0.32728041
             0.21230997      -4.24620884E-02   0.34282114
            -4.83514330E-02  -0.38934067       0.31384320
             -1.4301910      -0.11579874      -0.45509415
            -0.92108540       -1.6472000       7.15718386E-02
        """
        mdip = [float(x) for x in raw_vals.split()]

        self.run_ecd(inpfile, self.gen_ref(eigvals, edip, vdip, mdip), 'rpa')

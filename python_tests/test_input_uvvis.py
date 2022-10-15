from pathlib import Path
from unittest.mock import patch
import numpy as np
import tempfile
import sys

from veloxchem.veloxchemlib import is_mpi_master, mpi_barrier
from veloxchem.main import main


class TestInputUVVis:

    def create_input_file(self, lines, fname):

        if is_mpi_master():
            with open(fname, 'w') as f_inp:
                for line in lines.splitlines():
                    if line:
                        f_inp.write(f'{line}\n')
        mpi_barrier()

    def get_input_lines(self, xcfun, tamm_dancoff):

        input_lines = """
            @jobs
            task: response
            @end

            @method settings
            basis: aug-cc-pvdz
            @end

            @scf
            restart: no
            @end

            @response
            restart: no
            property: absorption
            tamm_dancoff: yes
            nstates: 5
            @end

            @molecule
            charge: 0
            multiplicity: 1
            units: au
            xyz:
            O   0.0   0.0   0.0
            H   0.0   1.4   1.1
            H   0.0  -1.4   1.1
            @end
        """

        if xcfun is not None:
            input_lines = input_lines.replace(
                '@method settings', f'@method settings\nxcfun:{xcfun}')

        if tamm_dancoff:
            return input_lines
        else:
            return input_lines.replace('tamm_dancoff: yes', 'tamm_dancoff: no')

    def run_input_uvvis(self, capsys, xcfun, tamm_dancoff, ref_data):

        with tempfile.TemporaryDirectory() as temp_dir:
            input_file = Path(temp_dir, 'vlx_uvvis')

            input_lines = self.get_input_lines(xcfun, tamm_dancoff)
            self.create_input_file(input_lines, input_file)

            with patch.object(sys, 'argv', ['vlx', str(input_file)]):
                main()
                captured = capsys.readouterr()

                if is_mpi_master():
                    exc_ene = []
                    osc_str = []

                    for line in captured.out.splitlines():
                        if 'Wave Function Model' in line:
                            wf_model = line.split(':')[1].strip()

                        if 'Total Energy' in line:
                            e_scf = float(line.split()[3])

                        if 'SCF converged in' in line:
                            n_scf_iterations = int(line.split()[4])

                        if 'Debye' in line and 'Total' in line:
                            dipole = float(line.split()[2])

                        if 'Number of States' in line:
                            nstates = int(line.split(':')[1].strip())

                        if 'converged in' in line and 'SCF' not in line:
                            n_uvvis_iterations = int(
                                line.split('converged in')[1].split()[0])

                        if 'Osc.Str.' in line:
                            content = line.split()
                            exc_ene.append(float(content[3]))
                            osc_str.append(float(content[-1]))

                    exc_ene = np.array(exc_ene)
                    osc_str = np.array(osc_str)

                    assert wf_model == ref_data['wf_model']
                    assert abs(e_scf - ref_data['e_scf']) < 1.0e-6
                    assert n_scf_iterations == ref_data['n_scf_iterations']
                    assert abs(dipole - ref_data['dipole']) < 1.0e-5

                    assert nstates == ref_data['nstates']
                    assert n_uvvis_iterations == ref_data['n_uvvis_iterations']
                    assert np.max(
                        np.abs(exc_ene - ref_data['exc_ene'])) < 1.0e-5
                    assert np.max(
                        np.abs(osc_str - ref_data['osc_str'])) < 1.0e-4

    def test_input_rhf_cis_uvvis(self, capsys):

        # vlxtag: RHF, Absorption, ECD, CIS

        xcfun = None
        tamm_dancoff = True

        ref_data = {
            'wf_model': 'Spin-Restricted Hartree-Fock',
            'e_scf': -76.0416975498,
            'n_scf_iterations': 9,
            'dipole': 0.786770,
            'nstates': 5,
            'n_uvvis_iterations': 9,
            'exc_ene': np.array(
                [0.32282999, 0.38491251, 0.41024825, 0.44793326, 0.47139996]),
            'osc_str': np.array([0.0530, 0.0000, 0.1070, 0.0045, 0.0248]),
        }

        self.run_input_uvvis(capsys, xcfun, tamm_dancoff, ref_data)

    def test_input_rhf_tdhf_uvvis(self, capsys):

        # vlxtag: RHF, Absorption, ECD, TDHF

        xcfun = None
        tamm_dancoff = False

        ref_data = {
            'wf_model': 'Spin-Restricted Hartree-Fock',
            'e_scf': -76.0416975498,
            'n_scf_iterations': 9,
            'dipole': 0.786770,
            'nstates': 5,
            'n_uvvis_iterations': 10,
            'exc_ene': np.array(
                [0.32139517, 0.38339622, 0.40934426, 0.44665234, 0.47004797]),
            'osc_str': np.array([0.0519, 0.0000, 0.1018, 0.0048, 0.0233]),
        }

        self.run_input_uvvis(capsys, xcfun, tamm_dancoff, ref_data)

    def test_input_rks_tda_uvvis(self, capsys):

        # vlxtag: RKS, Absorption, ECD, TDA

        xcfun = 'b3lyp'
        tamm_dancoff = True

        ref_data = {
            'wf_model': 'Spin-Restricted Kohn-Sham',
            'e_scf': -76.4435551612,
            'n_scf_iterations': 7,
            'dipole': 0.731256,
            'nstates': 5,
            'n_uvvis_iterations': 6,
            'exc_ene': np.array(
                [0.25719730, 0.30975371, 0.33858914, 0.37852621, 0.39052591]),
            'osc_str': np.array([0.0537, 0.0000, 0.0906, 0.0000, 0.0127]),
        }

        self.run_input_uvvis(capsys, xcfun, tamm_dancoff, ref_data)

    def test_input_rks_tddft_uvvis(self, capsys):

        # vlxtag: RKS, Absorption, ECD, TDDFT

        xcfun = 'b3lyp'
        tamm_dancoff = False

        ref_data = {
            'wf_model': 'Spin-Restricted Kohn-Sham',
            'e_scf': -76.4435551612,
            'n_scf_iterations': 7,
            'dipole': 0.731256,
            'nstates': 5,
            'n_uvvis_iterations': 7,
            'exc_ene': np.array(
                [0.25674877, 0.30969350, 0.33793671, 0.37842938, 0.39023439]),
            'osc_str': np.array([0.0520, 0.0000, 0.0848, 0.0000, 0.0116]),
        }

        self.run_input_uvvis(capsys, xcfun, tamm_dancoff, ref_data)

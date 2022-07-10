from pathlib import Path
from unittest.mock import patch
import tempfile
import sys

from veloxchem.veloxchemlib import is_mpi_master, mpi_barrier
from veloxchem.main import main


class TestInputEnergy:

    def create_input_file(self, lines, fname):

        if is_mpi_master():
            with open(fname, 'w') as f_inp:
                for line in lines.splitlines():
                    if line:
                        f_inp.write(f'{line}\n')
        mpi_barrier()

    def get_input_lines(self, xcfun, scf_type):

        input_lines = """
            @jobs
            task: scf
            @end

            @method settings
            basis: def2-svp
            @end

            @scf
            restart: no
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

        if scf_type == 'restricted':
            return input_lines

        elif scf_type == 'unrestricted':
            return input_lines.replace('task: scf', 'task: uscf').replace(
                'multiplicity: 1', 'multiplicity: 3')

        elif scf_type == 'restricted openshell':
            return input_lines.replace('task: scf', 'task: roscf').replace(
                'multiplicity: 1', 'multiplicity: 3')

    def run_input_energy(self, capsys, xcfun, scf_type, ref_data):

        with tempfile.TemporaryDirectory() as temp_dir:
            input_file = Path(temp_dir, 'vlx_scf')

            input_lines = self.get_input_lines(xcfun, scf_type)
            self.create_input_file(input_lines, input_file)

            with patch.object(sys, 'argv', ['vlx', str(input_file)]):
                main()
                captured = capsys.readouterr()

                if is_mpi_master():
                    for line in captured.out.splitlines():
                        if 'Wave Function Model' in line:
                            wf_model = line.split(':')[1].strip()

                        if 'Total Energy' in line:
                            e_scf = float(line.split()[3])

                        if 'SCF converged in' in line:
                            n_iterations = int(line.split()[4])

                        if 'Debye' in line and 'Total' in line:
                            dipole = float(line.split()[2])

                    assert wf_model == ref_data['wf_model']
                    assert abs(e_scf - ref_data['e_scf']) < 1.0e-6
                    assert n_iterations == ref_data['n_iterations']
                    assert abs(dipole - ref_data['dipole']) < 1.0e-5

    def test_input_rhf_energy(self, capsys):

        # vlxtag: RHF, Energy

        scf_type = 'restricted'
        xcfun = None

        ref_data = {
            'wf_model': 'Spin-Restricted Hartree-Fock',
            'e_scf': -75.9612569185,
            'n_iterations': 8,
            'dipole': 0.837406,
        }

        self.run_input_energy(capsys, xcfun, scf_type, ref_data)

    def test_input_uhf_energy(self, capsys):

        # vlxtag: UHF, Energy

        scf_type = 'unrestricted'
        xcfun = None

        ref_data = {
            'wf_model': 'Spin-Unrestricted Hartree-Fock',
            'e_scf': -75.7066294501,
            'n_iterations': 10,
            'dipole': 0.291780,
        }

        self.run_input_energy(capsys, xcfun, scf_type, ref_data)

    def test_input_rohf_energy(self, capsys):

        # vlxtag: ROHF, Energy

        scf_type = 'restricted openshell'
        xcfun = None

        ref_data = {
            'wf_model': 'Spin-Restricted Open-Shell Hartree-Fock',
            'e_scf': -75.7013577367,
            'n_iterations': 8,
            'dipole': 0.287662,
        }

        self.run_input_energy(capsys, xcfun, scf_type, ref_data)

    def test_input_rks_energy(self, capsys):

        # vlxtag: RKS, Energy

        scf_type = 'restricted'
        xcfun = 'b3lyp'

        ref_data = {
            'wf_model': 'Spin-Restricted Kohn-Sham',
            'e_scf': -76.3571498354,
            'n_iterations': 7,
            'dipole': 0.780568,
        }

        self.run_input_energy(capsys, xcfun, scf_type, ref_data)

    def test_input_uks_energy(self, capsys):

        # vlxtag: UKS, Energy

        scf_type = 'unrestricted'
        xcfun = 'b3lyp'

        ref_data = {
            'wf_model': 'Spin-Unrestricted Kohn-Sham',
            'e_scf': -76.0761283703,
            'n_iterations': 8,
            'dipole': 0.300485,
        }

        self.run_input_energy(capsys, xcfun, scf_type, ref_data)

    def test_input_roks_energy(self, capsys):

        # vlxtag: ROKS, Energy

        scf_type = 'restricted openshell'
        xcfun = 'b3lyp'

        ref_data = {
            'wf_model': 'Spin-Restricted Open-Shell Kohn-Sham',
            'e_scf': -76.0744709571,
            'n_iterations': 7,
            'dipole': 0.298955,
        }

        self.run_input_energy(capsys, xcfun, scf_type, ref_data)

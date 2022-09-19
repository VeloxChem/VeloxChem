from pathlib import Path
import tempfile
import pytest

from veloxchem.mpitask import MpiTask
from veloxchem.outputstream import OutputStream
from veloxchem.rspcdspec import CircularDichroismSpectrum
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.veloxchemlib import is_mpi_master


@pytest.mark.solvers
class TestCppEcd:

    def run_scf(self, task):

        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)
        scf_drv.update_settings(task.input_dict['scf'],
                                task.input_dict['method_settings'])
        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)

        return scf_drv.scf_tensors

    def run_cpp(self, inpfile, potfile, xcfun_label, ref_spectrum):

        task = MpiTask([inpfile, None])
        task.input_dict['scf']['checkpoint_file'] = None

        if potfile is not None:
            task.input_dict['method_settings']['potfile'] = potfile

        if xcfun_label is not None:
            task.input_dict['method_settings']['xcfun'] = xcfun_label

        scf_tensors = self.run_scf(task)

        ref_freqs = tuple([w for (w, Delta_epsilon) in ref_spectrum])

        cpp_prop = CircularDichroismSpectrum(
            {'frequencies': ref_freqs},
            task.input_dict['method_settings'],
        )
        cpp_prop.init_driver(task.mpi_comm, task.ostream)
        cpp_prop.compute(task.molecule, task.ao_basis, scf_tensors)

        assert cpp_prop.is_converged

        if is_mpi_master(task.mpi_comm):
            self.check_printout(cpp_prop)

            spectrum = cpp_prop.get_spectrum()
            for i, (w, Delta_epsilon) in enumerate(spectrum):
                ref_w, ref_Delta_epsilon = ref_spectrum[i]
                assert abs(w - ref_w) < 1.0e-6
                assert abs(Delta_epsilon / ref_Delta_epsilon - 1.0) < 1.0e-6

    def check_printout(self, cpp_prop):

        rsp_func = cpp_prop.rsp_property['response_functions']

        with tempfile.TemporaryDirectory() as temp_dir:
            fname = str(Path(temp_dir, 'cpp.out'))

            ostream = OutputStream(fname)
            cpp_prop.print_property(ostream)
            ostream.close()

            with open(fname, 'r') as f_out:
                lines = f_out.readlines()

            for key, val in rsp_func.items():
                key_found = False
                for line in lines:
                    if f'{key[0]}  ;  {key[1]}' in line:
                        content = line.split('>>')[1].split()
                        print_freq = float(content[0])
                        if abs(key[2] - print_freq) < 1e-4:
                            key_found = True
                            print_real = float(content[1])
                            print_imag = float(content[2].replace('j', ''))
                            assert abs(val.real - print_real) < 1.0e-6
                            assert abs(val.imag - print_imag) < 1.0e-6
                assert key_found

    def test_cpp_hf(self):

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'glycol.inp')

        potfile = None
        xcfun_label = None

        ref_spectrum = [
            (0.01, 6.570564e-05),
            (0.02, 0.0001317045),
            (0.03, 0.0001982910),
            (0.04, 0.0002657620),
            (0.05, 0.0003344180),
        ]

        self.run_cpp(inpfile, potfile, xcfun_label, ref_spectrum)

    def test_cpp_dft_slda(self):

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'glycol.inp')

        potfile = None
        xcfun_label = 'slda'

        ref_spectrum = [
            (0.01, 6.921420987263478e-05),
            (0.02, 0.00013859820885426657),
            (0.03, 0.0002083198761764199),
            (0.04, 0.00027854307847049683),
            (0.05, 0.00034590446460771457),
        ]

        self.run_cpp(inpfile, potfile, xcfun_label, ref_spectrum)

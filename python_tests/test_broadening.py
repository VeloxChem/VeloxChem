from pathlib import Path
import numpy as np
import pytest

from veloxchem.veloxchemlib import is_single_node
from veloxchem import LinearResponseEigenSolver


class TestBroadening:

    def read_tddft_output(self, output_file):

        nstates = 0

        with output_file.open('r') as fh:
            for line in fh:
                if 'Osc.Str.' in line:
                    content = line.split()
                    state_str = content[2].replace(':', '')
                    state_id = int(state_str.replace('S', '')) - 1
                    nstates = state_id + 1

        exc_ene = np.zeros(nstates)
        osc_str = np.zeros(nstates)
        rot_str = np.zeros(nstates)

        with output_file.open('r') as fh:
            for line in fh:
                if 'Osc.Str.' in line:
                    content = line.split()
                    state_str = content[2].replace(':', '')
                    state_id = int(state_str.replace('S', '')) - 1
                    # excitation energy in a.u.
                    exc_ene[state_id] = float(content[3])
                    # oscillator strength
                    osc_str[state_id] = float(content[8])
                if 'Rot.Str.' in line:
                    content = line.split()
                    state_str = content[2].replace(':', '')
                    state_id = int(state_str.replace('S', '')) - 1
                    # rotatory strength in 10**(-40) cgs
                    rot_str[state_id] = float(content[6])

        return exc_ene, osc_str, rot_str

    @pytest.mark.skipif(not is_single_node(), reason="single node only")
    def test_ecd_broadening(self):

        tddft_file = Path(__file__).parent / 'inputs' / 'noradr_tddft.txt'

        exc_ene, osc_str, rot_str = self.read_tddft_output(tddft_file)
        xplot, yplot = LinearResponseEigenSolver.lorentzian_ecd_spectrum(
            exc_ene, rot_str, 0.1825, 0.5, 0.0025)

        ecd_file = Path(__file__).parent / 'inputs' / 'noradr_ecd.txt'

        xref = []
        yref = []
        with ecd_file.open('r') as fh:
            for line in fh:
                content = line.split()
                xref.append(float(content[1]))
                yref.append(float(content[2]))
        xref = np.array(xref)
        yref = np.array(yref)

        assert xplot.size == xref.size
        assert np.max(np.abs(xplot - xref)) < 5.0e-4

        assert yplot.size == yref.size
        for i in range(32):
            if abs(yref[i]) > 2.0:
                assert np.abs(yplot[i] / yref[i] - 1.0) < 1.0e-2
            else:
                assert np.abs(yplot[i] / yref[i] - 1.0) < 5.0e-2

        try:
            import matplotlib.pyplot as plt

            plt.xlabel(r'Photon energy [eV]')
            plt.ylabel(r'$\Delta\epsilon$($\omega$) [L mol$^{-1}$cm$^{-1}$]')

            plt.plot(xref,
                     yref,
                     linestyle='-',
                     linewidth=1.5,
                     color='C0',
                     label='ECD (CPP)')

            plt.plot(xplot,
                     yplot,
                     linestyle='--',
                     linewidth=1.5,
                     color='C1',
                     label='ECD (TDDFT)')

            plt.legend()

            fig = plt.gcf()
            fig.set_size_inches(8.0, 5.0)
            pdf_fname = '_vlx_test_ecd_broadening_.pdf'
            fig.savefig(pdf_fname, dpi=300, bbox_inches='tight')
            fig.clear()

            pdf_file = Path(pdf_fname)
            if pdf_file.is_file():
                pdf_file.unlink()

        except ImportError:
            pass

    @pytest.mark.skipif(not is_single_node(), reason="single node only")
    def test_absorption_broadening(self):

        tddft_file = Path(__file__).parent / 'inputs' / 'noradr_tddft.txt'

        exc_ene, osc_str, rot_str = self.read_tddft_output(tddft_file)
        xplot, yplot = LinearResponseEigenSolver.lorentzian_absorption_spectrum(
            exc_ene, osc_str, 0.0025, 0.5, 0.0025)

        abs_file = Path(__file__).parent / 'inputs' / 'noradr_absorption.txt'

        xref = []
        yref = []
        with abs_file.open('r') as fh:
            for line in fh:
                content = line.split()
                xref.append(float(content[1]))
                yref.append(float(content[2]))
        xref = np.array(xref)
        yref = np.array(yref)

        assert xplot.size == xref.size
        assert np.max(np.abs(xplot - xref)) < 5.0e-4

        assert yplot.size == yref.size
        for i in range(100):
            if abs(yref[i]) < 0.5:
                assert np.abs(yplot[i] - yref[i]) < 3.0e-2
            else:
                assert np.abs(yplot[i] / yref[i] - 1.0) < 5.0e-2

        try:
            import matplotlib.pyplot as plt

            plt.xlabel(r'Photon energy [eV]')
            plt.ylabel(r'$\sigma$($\omega$) [a.u.]')

            plt.plot(xref,
                     yref,
                     linestyle='-',
                     linewidth=1.5,
                     color='C0',
                     label='Absorption (CPP)')

            plt.plot(xplot,
                     yplot,
                     linestyle='--',
                     linewidth=1.5,
                     color='C1',
                     label='Absorption (TDDFT)')

            plt.legend()

            fig = plt.gcf()
            fig.set_size_inches(8.0, 5.0)
            pdf_fname = '_vlx_test_absorption_broadening_.pdf'
            fig.savefig(pdf_fname, dpi=300, bbox_inches='tight')
            fig.clear()

            pdf_file = Path(pdf_fname)
            if pdf_file.is_file():
                pdf_file.unlink()

        except ImportError:
            pass

from pathlib import Path
import numpy as np

from broaden_ecd_spectrum import lorentzian_ecd_ev
from broaden_ecd_spectrum import read_energies_and_rotatory_strengths


class TestBroadening:

    def test_ecd_broadening(self):

        tddft_file = Path(__file__).parent / 'test_data' / 'noradr_tddft.txt'

        exc_ene, rot_str = read_energies_and_rotatory_strengths(tddft_file)
        xplot, yplot = lorentzian_ecd_ev(exc_ene, rot_str, 0.1825, 0.5, 0.0025)

        ecd_file = Path(__file__).parent / 'test_data' / 'noradr_ecd.txt'

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

            pdf_file = Path(pdf_fname)
            if pdf_file.is_file():
                pdf_file.unlink()

        except ImportError:
            pass

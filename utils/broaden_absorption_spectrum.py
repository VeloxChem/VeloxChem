from pathlib import Path
import numpy as np
import math
import sys


def lorentzian_absorption_ev(x,
                             y,
                             xmin=None,
                             xmax=None,
                             xstep=0.0002,
                             gamma=0.0045563353):
    """
    Broadens absorption stick spectrum.

    :param x:
        Excitation energies in a.u.
    :param y:
        Oscillator strengths.
    :param xmin:
        Minimal excitation energy in a.u.
    :param xmax:
        Maximum excitation energy in a.u.
    :param xstep:
        Step size of excitation energy in a.u.
    :param gamma:
        The broadening parameter in a.u.

    :return:
        The excitation energies in eV and sigma(w) in a.u.
    """

    if xmin is None:
        xmin = max(np.min(x) - 0.01, xstep)

    if xmax is None:
        xmax = np.max(x) + 0.01

    x_i = np.arange(xmin, xmax + xstep / 100.0, xstep, dtype=np.float64)
    y_i = np.zeros_like(x_i)

    fine_structure_constant = 7.2973525693e-3
    factor = 2.0 * math.pi * fine_structure_constant

    for i in range(x_i.size):
        for j in range(y.size):
            y_i[i] += factor * gamma / ((x_i[i] - x[j])**2 + gamma**2) * y[j]

    hartree_in_ev = 27.211386246
    x_i *= hartree_in_ev

    return x_i, y_i


def lorentzian_absorption_nm(x,
                             y,
                             wmin=None,
                             wmax=None,
                             wstep=0.5,
                             gamma=0.0045563353):
    """
    Broadens absorption stick spectrum.

    :param x:
        Excitation energies in a.u.
    :param y:
        Oscillator strengths.
    :param xmin:
        Minimal wavelength in nm.
    :param xmax:
        Maximum wavelength in nm.
    :param xstep:
        Step size of wavelength in nm.
    :param gamma:
        The broadening parameter in a.u.

    :return:
        The wavelengths in nm and sigma(w) in a.u.
    """

    hartree_in_inverse_nm = 0.0219474631
    hartree_nm = 1.0 / hartree_in_inverse_nm

    if wmin is None:
        wmin = max(int(hartree_nm / np.max(x)) - 50, wstep)

    if wmax is None:
        wmax = int(hartree_nm / np.min(x)) + 75

    w_i = np.arange(wmin, wmax + wstep / 100.0, wstep, dtype=np.float64)
    x_i = hartree_nm / w_i
    y_i = np.zeros_like(x_i)

    fine_structure_constant = 7.2973525693e-3
    factor = 2.0 * math.pi * fine_structure_constant

    for i in range(x_i.size):
        for j in range(y.size):
            y_i[i] += factor * gamma / ((x_i[i] - x[j])**2 + gamma**2) * y[j]

    return w_i, y_i


def read_energies_and_oscillator_strengths(output_file):
    """
    Reads excitation energies and oscillator strengths from output file.

    :param output_file:
        The output file (Path object).

    :return:
        The excitation energies in a.u. and oscillator strengths.
    """

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

    with output_file.open('r') as fh:
        for line in fh:
            if 'Osc.Str.' in line:
                content = line.split()
                state_str = content[2].replace(':', '')
                state_id = int(state_str.replace('S', '')) - 1
                exc_ene[state_id] = float(content[3])
                osc_str[state_id] = float(content[8])

    return exc_ene, osc_str


if __name__ == '__main__':

    if len(sys.argv) <= 1:
        fname = Path(__file__).name
        print(f'Usage: python3 {fname} veloxchem_output_file.out [ev|nm]')
        sys.exit(1)

    output_file = Path(sys.argv[1])
    x_unit = sys.argv[2].lower() if len(sys.argv) > 2 else 'ev'

    exc_ene, osc_str = read_energies_and_oscillator_strengths(output_file)

    if x_unit == 'ev':
        xplot, yplot = lorentzian_absorption_ev(exc_ene, osc_str)
    elif x_unit == 'nm':
        xplot, yplot = lorentzian_absorption_nm(exc_ene, osc_str)

    for i in range(xplot.size):
        print(xplot[i], yplot[i])

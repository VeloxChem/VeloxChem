from pathlib import Path
import numpy as np
import sys


def lorentzian_ecd_ev(exc_ene,
                      rot_str,
                      e_min=None,
                      e_max=None,
                      e_step=0.0002,
                      gamma=0.0045563353):
    """
    Broadens ECD stick spectrum.

    :param exc_ene:
        Excitation energies in a.u.
    :param rot_str:
        Rotatory strengths in a.u.
    :param e_min:
        Minimal excitation energy in a.u.
    :param e_max:
        Maximum excitation energy in a.u.
    :param e_step:
        Step size of excitation energy in a.u.
    :param gamma:
        The broadening parameter in a.u.

    :return:
        The excitation energies in eV and Delta_epsilon in L mol^-1 cm^-1
    """

    if e_min is None:
        e_min = max(np.min(exc_ene) - 0.01, e_step)

    if e_max is None:
        e_max = np.max(exc_ene) + 0.01

    x_i = np.arange(e_min, e_max + e_step / 100.0, e_step, dtype=np.float64)
    y_i = np.zeros_like(x_i)

    extinction_coefficient_from_beta = 19.6036975758
    factor = extinction_coefficient_from_beta / 3

    for i in range(x_i.size):
        for s in range(exc_ene.size):
            y_i[i] += factor * gamma / (
                (x_i[i] - exc_ene[s])**2 + gamma**2) * exc_ene[s] * rot_str[s]

    hartree_in_ev = 27.211386246
    x_i *= hartree_in_ev

    return x_i, y_i


def lorentzian_ecd_nm(exc_ene,
                      rot_str,
                      w_min=None,
                      w_max=None,
                      w_step=0.5,
                      gamma=0.0045563353):
    """
    Broadens ECD stick spectrum.

    :param exc_ene:
        Excitation energies in a.u.
    :param rot_str:
        Rotatory strengths in a.u.
    :param xmin:
        Minimal wavelength in nm.
    :param xmax:
        Maximum wavelength in nm.
    :param xstep:
        Step size of wavelength in nm.
    :param gamma:
        The broadening parameter in a.u.

    :return:
        The wavelengths in nm and Delta_epsilon in L mol^-1 cm^-1
    """

    hartree_in_inverse_nm = 0.0219474631
    hartree_nm = 1.0 / hartree_in_inverse_nm

    if w_min is None:
        w_min = max(int(hartree_nm / np.max(exc_ene)) - 50, w_step)

    if w_max is None:
        w_max = int(hartree_nm / np.min(exc_ene)) + 75

    w_i = np.arange(w_min, w_max + w_step / 100.0, w_step, dtype=np.float64)
    x_i = hartree_nm / w_i
    y_i = np.zeros_like(x_i)

    extinction_coefficient_from_beta = 19.6036975758
    factor = extinction_coefficient_from_beta / 3

    for i in range(x_i.size):
        for s in range(exc_ene.size):
            y_i[i] += factor * gamma / (
                (x_i[i] - exc_ene[s])**2 + gamma**2) * exc_ene[s] * rot_str[s]

    return w_i, y_i


def read_energies_and_rotatory_strengths(output_file):
    """
    Reads excitation energies and rotatory strengths from output file.

    :param output_file:
        The output file (Path object).

    :return:
        The excitation energies in a.u. and rotatory strengths in a.u.
    """

    nstates = 0

    with output_file.open('r') as fh:
        for line in fh:
            if 'Rot.Str.' in line:
                content = line.split()
                state_str = content[2].replace(':', '')
                state_id = int(state_str.replace('S', '')) - 1
                nstates = state_id + 1

    exc_ene = np.zeros(nstates)
    rot_str = np.zeros(nstates)

    with output_file.open('r') as fh:
        for line in fh:
            if 'Osc.Str.' in line:
                content = line.split()
                state_str = content[2].replace(':', '')
                state_id = int(state_str.replace('S', '')) - 1
                exc_ene[state_id] = float(content[3])
            if 'Rot.Str.' in line:
                content = line.split()
                state_str = content[2].replace(':', '')
                state_id = int(state_str.replace('S', '')) - 1
                rot_str[state_id] = float(content[4])

    return exc_ene, rot_str


if __name__ == '__main__':

    if len(sys.argv) <= 1:
        fname = Path(__file__).name
        print(f'Usage: python3 {fname} veloxchem_output_file.out [ev|nm]')
        sys.exit(1)

    output_file = Path(sys.argv[1])
    x_unit = sys.argv[2].lower() if len(sys.argv) > 2 else 'ev'

    exc_ene, rot_str = read_energies_and_rotatory_strengths(output_file)

    if x_unit == 'ev':
        xplot, yplot = lorentzian_ecd_ev(exc_ene, rot_str)
    elif x_unit == 'nm':
        xplot, yplot = lorentzian_ecd_nm(exc_ene, rot_str)

    for i in range(xplot.size):
        print(xplot[i], yplot[i])

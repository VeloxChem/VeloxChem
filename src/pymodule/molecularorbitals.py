import numpy as np
import h5py
import math
import sys

from .veloxchemlib import MolecularOrbitals
from .veloxchemlib import molorb
from .veloxchemlib import to_angular_momentum
from .outputstream import OutputStream
from .errorhandler import assert_msg_critical


def _print_orbitals(self,
                    molecule,
                    ao_basis,
                    all_orbs=False,
                    ostream=OutputStream(sys.stdout)):

    norb = self.number_mos()

    ao_map = ao_basis.get_ao_basis_map(molecule)

    if self.get_orbitals_type() == molorb.rest:

        ostream.print_blank()

        ostream.print_header("Spin Restricted Orbitals")
        ostream.print_header("------------------------")

        nocc = molecule.number_of_electrons() // 2

        if all_orbs:
            nstart, nend = 0, norb
        else:
            nstart, nend = max(0, nocc - 5), min(norb, nocc + 5)

        rvecs = self.alpha_to_numpy()
        reigs = self.ea_to_numpy()
        rnocc = [2.0 if x < nocc else 0.0 for x in range(norb)]

        for i in range(nstart, nend):
            _print_coefficients(reigs[i], rnocc[i], i, rvecs[:, i], ao_map,
                                0.15, ostream)

        ostream.print_blank()

    elif self.get_orbitals_type() == molorb.unrest:

        ostream.print_blank()

        ostream.print_header("Spin Unrestricted Alpha Orbitals")
        ostream.print_header("--------------------------------")

        # FIX ME

        ostream.print_blank()

        ostream.print_header("Spin Unrestricted Beta Orbitals")
        ostream.print_header("-------------------------------")

        # FIX ME

    else:

        errmsg = "MolecularOrbitals.get_density:"
        errmsg += " Invalid molecular orbitals type"
        assert_msg_critical(False, errmsg)


def _print_coefficients(eigval, focc, iorb, coeffs, ao_map, thresh, ostream):
    ostream.print_blank()

    valstr = "Molecular Orbital No.{:4d}:".format(iorb + 1)
    ostream.print_header(valstr.ljust(92))
    valstr = 26 * "-"
    ostream.print_header(valstr.ljust(92))

    valstr = "Occupation: {:.1f} Energy: {:10.5f} au".format(focc, eigval)
    ostream.print_header(valstr.ljust(92))

    tuplist = []

    for i in range(coeffs.shape[0]):

        if math.fabs(coeffs[i]) > thresh:
            atomidx = int(ao_map[i].split()[0])
            anglmom = to_angular_momentum(ao_map[i][-3].upper())
            valstr = "(" + ao_map[i] + ": {:8.2f}".format(coeffs[i]) + ") "
            tuplist.append((atomidx, anglmom, valstr))

    valstr = ""
    curidx = 0

    for t in sorted(tuplist):
        valstr += t[-1]
        curidx += 1

        if curidx == 3:
            ostream.print_header(valstr.ljust(92))
            valstr = ""
            curidx = 0

    if curidx > 0:
        ostream.print_header(valstr.ljust(92))


def _get_density(self, molecule):

    if self.get_orbitals_type() == molorb.rest:

        nelec = molecule.number_of_electrons()
        return self.get_ao_density(nelec)

    elif self.get_orbitals_type() == molorb.unrest:

        nalpha = molecule.number_of_alpha_electrons()
        nbeta = molecule.number_of_beta_electrons()
        return self.get_ao_density(nalpha, nbeta)

    else:

        errmsg = "MolecularOrbitals.get_density:"
        errmsg += " Invalid molecular orbitals type"
        assert_msg_critical(False, errmsg)


def _write_hdf5(self, fname):

    hf = h5py.File(fname, 'w')

    name = "alpha_orbitals"
    array = self.alpha_to_numpy()
    hf.create_dataset(name, data=array, compression="gzip")

    name = "alpha_energies"
    array = self.ea_to_numpy()
    hf.create_dataset(name, data=array, compression="gzip")

    if self.get_orbitals_type() == molorb.unrest:

        name = "beta_orbitals"
        array = self.beta_to_numpy()
        hf.create_dataset(name, data=array, compression="gzip")

        name = "beta_energies"
        array = self.eb_to_numpy()
        hf.create_dataset(name, data=array, compression="gzip")

    hf.close()


@staticmethod
def _read_hdf5(fname):

    hf = h5py.File(fname, 'r')

    assert_msg_critical(
        len(list(hf.keys())) == 2 or len(list(hf.keys())) == 4,
        "MolecularOrbitals.read_hdf5: Incorrect number of datasets!")

    if len(list(hf.keys())) == 2:
        orbs_type = molorb.rest
    elif len(list(hf.keys())) == 4:
        orbs_type = molorb.unrest

    orbs = []
    enes = []

    orbs.append(np.array(hf.get("alpha_orbitals")))
    enes.append(np.array(hf.get("alpha_energies")))

    if orbs_type == molorb.unrest:
        orbs.append(np.array(hf.get("beta_orbitals")))
        enes.append(np.array(hf.get("beta_energies")))

    hf.close()

    return MolecularOrbitals(orbs, enes, orbs_type)


MolecularOrbitals.get_density = _get_density
MolecularOrbitals.write_hdf5 = _write_hdf5
MolecularOrbitals.read_hdf5 = _read_hdf5
MolecularOrbitals.print_orbitals = _print_orbitals

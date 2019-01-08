from .veloxchemlib import MolecularOrbitals
from .veloxchemlib import molorb
from .veloxchemlib import assert_msg_critical
import h5py
import numpy as np


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


MolecularOrbitals.write_hdf5 = _write_hdf5
MolecularOrbitals.read_hdf5 = _read_hdf5

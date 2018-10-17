from .VeloxChemLib import MolecularOrbitals
from .VeloxChemLib import molorb
from .VeloxChemLib import assert_msg_critical
import h5py
import numpy as np


def _write_hdf5(self, fname):

    hf = h5py.File(fname, 'w')

    count = 0

    for index in range(self.get_number_of_orbitals_matrices()):

        if self.get_orbitals_type() == molorb.rest:
            name = str(count) + "_total_" + str(index)
            array = self.total_to_numpy(index)
            hf.create_dataset(name, data=array, compression="gzip")
            count += 1

        else:
            name = str(count) + "_alpha_" + str(index)
            array = self.alpha_to_numpy(index)
            hf.create_dataset(name, data=array, compression="gzip")
            count += 1

            name = str(count) + "_beta_" + str(index)
            array = self.beta_to_numpy(index)
            hf.create_dataset(name, data=array, compression="gzip")
            count += 1

    hf.close()


@staticmethod
def _read_hdf5(fname):

    orbtype = {
        "total": molorb.rest,
        "alpha": molorb.unrest,
        "beta": molorb.unrest
    }

    hf = h5py.File(fname, 'r')

    dens = []
    types = []

    ordered_keys = []
    for key in list(hf.keys()):
        i = int(key.split("_")[0])
        ordered_keys.append((i, key))
    ordered_keys.sort()

    for i, key in ordered_keys:
        type_str = key.split("_")[1]
        dens.append(np.array(hf.get(key)))
        types.append(orbtype[type_str])

    hf.close()

    assert_msg_critical(
        len(set(types)) == 1,
        "MolecularOrbitals.read_hdf5: inconsistent orbital type!")

    assert_msg_critical(
        types[0] in list(orbtype.values()),
        "MolecularOrbitals.read_hdf5: invalid orbital type!")

    return MolecularOrbitals.from_numpy_list(dens, types[0])


MolecularOrbitals.write_hdf5 = _write_hdf5
MolecularOrbitals.read_hdf5 = _read_hdf5

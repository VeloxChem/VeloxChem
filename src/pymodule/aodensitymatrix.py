from .veloxchemlib import AODensityMatrix
from .veloxchemlib import denmat
from .veloxchemlib import assert_msg_critical
import h5py
import numpy as np


def _write_hdf5(self, fname):

    hf = h5py.File(fname, 'w')

    count = 0

    for index in range(self.number_of_density_matrices()):

        if self.get_density_type() == denmat.rest:
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

    dentype = {
        "total": denmat.rest,
        "alpha": denmat.unrest,
        "beta": denmat.unrest
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
        types.append(dentype[type_str])

    hf.close()

    assert_msg_critical(
        len(set(types)) == 1,
        "AODensityMatrix.read_hdf5: inconsistent density type!")

    assert_msg_critical(
        types[0] in list(dentype.values()),
        "AODensityMatrix.read_hdf5: invalid density type!")

    return AODensityMatrix(dens, types[0])


AODensityMatrix.write_hdf5 = _write_hdf5
AODensityMatrix.read_hdf5 = _read_hdf5

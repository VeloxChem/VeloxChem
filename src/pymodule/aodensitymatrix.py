from .VeloxChemLib import AODensityMatrix
from .VeloxChemLib import denmat
from .VeloxChemLib import assert_msg_critical
import h5py
import numpy as np


def _write_hdf5(self, fname):

    hf = h5py.File(fname, 'w')

    for i in range(self.get_number_of_density_matrices()):

        if self.get_density_type() == denmat.rest:
            name = "total_" + str(i)
            array = self.total_to_numpy(i)
            hf.create_dataset(name, data=array, compression="gzip")

        else:
            name = "alpha_" + str(i)
            array = self.alpha_to_numpy(i)
            hf.create_dataset(name, data=array, compression="gzip")

            name = "beta_" + str(i)
            array = self.beta_to_numpy(i)
            hf.create_dataset(name, data=array, compression="gzip")

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

    for key in list(hf.keys()):
        type_str, id_str = key.split("_")
        dens.append(np.array(hf.get(key)))
        types.append(dentype[type_str])

    hf.close()

    assert_msg_critical(
        len(set(types)) == 1,
        "AODensityMatrix.read_hdf5: inconsistent density type!")

    assert_msg_critical(
        types[0] in list(dentype.values()),
        "AODensityMatrix.read_hdf5: invalid density type!")

    return AODensityMatrix.from_numpy_list(dens, types[0])


AODensityMatrix.write_hdf5 = _write_hdf5
AODensityMatrix.read_hdf5 = _read_hdf5

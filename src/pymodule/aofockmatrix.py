from .VeloxChemLib import AOFockMatrix
from .VeloxChemLib import fockmat
from .VeloxChemLib import assert_msg_critical
import h5py
import numpy as np


def _write_hdf5(self, fname):

    focktype = {
        fockmat.restjk: "restjk",
        fockmat.restjkx: "restjkx",
        fockmat.restj: "restj",
        fockmat.restk: "restk",
        fockmat.restkx: "restkx"
    }

    hf = h5py.File(fname, 'w')

    factors = []
    for i in range(self.get_number_of_fock_matrices()):
        factors.append(self.get_scale_factor(i))
    hf.create_dataset("factors", data=factors, compression="gzip")

    for i in range(self.get_number_of_fock_matrices()):
        index = self.get_density_identifier(i)
        name = focktype[self.get_fock_type(i)] + "_" + str(index)
        array = self.to_numpy(i)
        hf.create_dataset(name, data=array, compression="gzip")

    hf.close()


@staticmethod
def _read_hdf5(fname):

    focktype = {
        "restjk": fockmat.restjk,
        "restjkx": fockmat.restjkx,
        "restj": fockmat.restj,
        "restk": fockmat.restk,
        "restkx": fockmat.restkx
    }

    hf = h5py.File(fname, 'r')

    focks = []
    types = []
    factors = list(hf.get("factors"))
    indices = []

    keys = list(hf.keys())
    for key in keys:
        if key != "factors":
            type_str, id_str = key.split("_")
            focks.append(np.array(hf.get(key)))
            types.append(focktype[type_str])
            indices.append(int(id_str))

    hf.close()

    for ftype in set(types):
        assert_msg_critical(
            ftype in list(focktype.values()),
            "AOFockMatrix.read_hdf5: invalid Fock types!")

    return AOFockMatrix.from_numpy_list(focks, types, factors, indices)


AOFockMatrix.write_hdf5 = _write_hdf5
AOFockMatrix.read_hdf5 = _read_hdf5

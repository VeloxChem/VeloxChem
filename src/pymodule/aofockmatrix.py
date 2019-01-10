from .veloxchemlib import AOFockMatrix
from .veloxchemlib import fockmat
from .veloxchemlib import assert_msg_critical
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
    for i in range(self.number_of_fock_matrices()):
        factors.append(self.get_scale_factor(i))
    hf.create_dataset("factors", data=factors, compression="gzip")

    for i in range(self.number_of_fock_matrices()):
        index = self.get_density_identifier(i)
        name = str(i) + "_" + focktype[self.get_fock_type(i)] + "_" + str(index)
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

    ordered_keys = []
    for key in list(hf.keys()):
        if key == "factors":
            continue
        i = int(key.split("_")[0])
        ordered_keys.append((i, key))
    ordered_keys.sort()

    for i, key in ordered_keys:
        type_str, index_str = key.split("_")[1:]
        focks.append(np.array(hf.get(key)))
        types.append(focktype[type_str])
        indices.append(int(index_str))

    hf.close()

    for ftype in set(types):
        assert_msg_critical(
            ftype in list(focktype.values()),
            "AOFockMatrix.read_hdf5: invalid Fock types!")

    return AOFockMatrix(focks, types, factors, indices)


AOFockMatrix.write_hdf5 = _write_hdf5
AOFockMatrix.read_hdf5 = _read_hdf5

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "moldata/Molecule.hpp"
#include "orbdata/MolecularBasis.hpp"
#include "t4c_func/FockDriver.hpp"
#include "t4c_func/T4CScreener.hpp"

namespace py = pybind11;

void export_aoeris(py::module_ &m) {
    m.def("compute_ao_eris", [](const CMolecule &mol, const CMolecularBasis &basis) {
        // Construct screener and Fock driver using partition method
        CT4CScreener screener;
        screener.partition(basis, mol, "eri");
        CFockDriver fock_drv;
        int nao = basis.basis_functions().size();
        int ithreshold = 0; // Use default threshold (tightest)
        auto eri_tensor = fock_drv.compute_eri(screener, nao, ithreshold);
        // Convert to numpy array (shape: nao, nao, nao, nao)
        py::array_t<double> result({nao, nao, nao, nao});
        auto buf = result.mutable_unchecked<4>();
        const double* data = eri_tensor.values();
        // CDense4DTensor is stored in (i,j,k,l) order, contiguous
        size_t idx = 0;
        for (int i = 0; i < nao; ++i)
            for (int j = 0; j < nao; ++j)
                for (int k = 0; k < nao; ++k)
                    for (int l = 0; l < nao; ++l, ++idx)
                        buf(i, j, k, l) = data[idx];
        return result;
    }, R"pbdoc(Returns AO four-index electron repulsion integrals (μν|λσ) as a numpy array.)pbdoc");
}

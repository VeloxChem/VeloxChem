#include "ExportOneElecInts.hpp"

#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "AngularMomentumIntegrals.hpp"
#include "ExportGeneral.hpp"
#include "ErrorHandler.hpp"
#include "LinearMomentumIntegrals.hpp"

namespace py = pybind11;
using namespace py::literals;

namespace vlx_oneeints {

// Exports classes/functions in src/onee_ints to python

void
export_oneeints(py::module& m)
{
    m.def("compute_linear_momentum_integrals",
            [](const CMolecule&           molecule,
               const CMolecularBasis&     basis) -> py::list {
                auto linmom = onee::computeLinearMomentumIntegrals(molecule, basis);
                py::list ret;
                ret.append(vlx_general::pointer_to_numpy(linmom[0].values(), {linmom[0].getNumberOfRows(), linmom[0].getNumberOfColumns()}));
                ret.append(vlx_general::pointer_to_numpy(linmom[1].values(), {linmom[1].getNumberOfRows(), linmom[1].getNumberOfColumns()}));
                ret.append(vlx_general::pointer_to_numpy(linmom[2].values(), {linmom[2].getNumberOfRows(), linmom[2].getNumberOfColumns()}));
                return ret;
            },
            "Computes linear momentum integrals.");

    m.def("compute_angular_momentum_integrals",
            [](const CMolecule&           molecule,
               const CMolecularBasis&     basis,
               const std::vector<double>& origin) -> py::list {
                auto angmom = onee::computeAngularMomentumIntegrals(molecule, basis, origin);
                py::list ret;
                ret.append(vlx_general::pointer_to_numpy(angmom[0].values(), {angmom[0].getNumberOfRows(), angmom[0].getNumberOfColumns()}));
                ret.append(vlx_general::pointer_to_numpy(angmom[1].values(), {angmom[1].getNumberOfRows(), angmom[1].getNumberOfColumns()}));
                ret.append(vlx_general::pointer_to_numpy(angmom[2].values(), {angmom[2].getNumberOfRows(), angmom[2].getNumberOfColumns()}));
                return ret;
            },
            "Computes angular momentum integrals.");
}

}  // namespace vlx_oneeints

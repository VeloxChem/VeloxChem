//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include <pybind11/pybind11.h>

#include "ExportResponse.hpp"

#include "ExcitationVector.hpp"

namespace py = pybind11;

namespace vlx_rsp { // vlx_rsp namespace
    
// Exports classes/functions in src/response to python
    
// Helper function for printing CExcitationVector
    
static std::string
CExcitationVector_str(const CExcitationVector& self)
{
    return self.getString();
}

void export_response(py::module& m)
{
    
    // CExcitationVector class
    
    py::class_< CExcitationVector, std::shared_ptr<CExcitationVector> >
        (
            m, "ExcitationVector"
        )
        .def(py::init<>())
        .def(py::init<const szblock,
                      const std::vector<int32_t>&,
                      const std::vector<int32_t>&,
                      const std::vector<double>&,
                      const std::vector<double>&>())
        .def(py::init<const szblock,
                      const int32_t,
                      const int32_t,
                      const int32_t,
                      const int32_t,
                      const bool>())
        .def(py::init<const szblock,
                      const int32_t,
                      const int32_t,
                      const int32_t,
                      const int32_t>())
        .def(py::init<const CExcitationVector&>())
        .def("__str__", &CExcitationVector_str)
        .def("number_excitations", &CExcitationVector::getNumberOfExcitations)
        .def("bra_unique_indexes", &CExcitationVector::getBraUniqueIndexes)
        .def("ket_unique_indexes", &CExcitationVector::getKetUniqueIndexes)
        .def("get_zmatrix", &CExcitationVector::getMatrixZ)
        .def("get_ymatrix", &CExcitationVector::getMatrixY)
        .def("get_zdensity", &CExcitationVector::getDensityZ)
        .def("get_ydensity", &CExcitationVector::getDensityY)
    ;
}
    
} // vlx_rsp namespace

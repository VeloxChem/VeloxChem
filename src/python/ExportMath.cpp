//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

#include "DenseMatrix.hpp"

#include "ExportGeneral.hpp"
#include "ExportMath.hpp"

namespace bp = boost::python;

namespace np = boost::python::numpy;

namespace bp_math { // bp_math namespace

// Helper function for printing CDenseMatrix

std::string
CDenseMatrix_str(const CDenseMatrix& self)
{
    return self.getString();
}

// Helper function for converting CDenseMatrix to numpy array

np::ndarray
CDenseMatrix_to_numpy(const CDenseMatrix& self)
{
    return bp_general::pointer_to_numpy(self.values(),
                                        self.getNumberOfRows(),
                                        self.getNumberOfColumns());
}

// Helper function for converting numpy array to CDenseMatrix

CDenseMatrix
CDenseMatrix_from_numpy(const np::ndarray& arr)
{
    const double* data = reinterpret_cast<double*>(arr.get_data());

    if (data == nullptr) return CDenseMatrix();

    auto size = static_cast<int32_t> (arr.shape(0) * arr.shape(1));

    if (size == 0) return CDenseMatrix();

    std::vector<double> vec (data, data + size);

    int32_t nrows = static_cast<int32_t>(arr.shape(0));

    int32_t ncols = static_cast<int32_t>(arr.shape(1));

    return CDenseMatrix(vec, nrows, ncols);
}

// Exports classes/functions in src/math to python

void export_math()
{
    // initialize numpy

    Py_Initialize();

    np::initialize();

    // CDenseMatrix class
    // Note: CDenseMatrix has several constructors

    bp::class_< CDenseMatrix, std::shared_ptr<CDenseMatrix> >
        (
            "DenseMatrix",
            bp::init<
                const std::vector<double>&,
                const int32_t,
                const int32_t
                >()
        )
        .def(bp::init<>())
        .def(bp::init<const int32_t>())
        .def(bp::init<const int32_t, const int32_t>())
        .def(bp::init<const CDenseMatrix&>())
        .def("__str__", &CDenseMatrix_str)
        .def("to_numpy", &CDenseMatrix_to_numpy)
        .def("from_numpy", &CDenseMatrix_from_numpy)
        .staticmethod("from_numpy")
        .def("zero", &CDenseMatrix::zero)
        .def("symmetrize", &CDenseMatrix::symmetrize)
        .def("slice", &CDenseMatrix::slice)
        .def(bp::self == bp::other<CDenseMatrix>())
    ;
}

} // bp_math namespace

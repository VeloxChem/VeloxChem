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
CDenseMatrix_str (const CDenseMatrix& self)
{
    return self.getString();
}

// Helper function for converting CDenseMatrix to numpy array

np::ndarray
dense_to_numpy(const CDenseMatrix& mat)
{
    return bp_general::pointer_to_numpy(mat.values(),
                                        mat.getNumberOfRows(),
                                        mat.getNumberOfColumns());
}

// Exports classes/functions in src/math to python

void export_math()
{
    // initialize numpy

    Py_Initialize();

    np::initialize();

    // exposing functions

    bp::def("to_numpy", &dense_to_numpy);

    // CDenseMatrix class
    // Note: CDenseMatrix has several constructors

    bp::class_< CDenseMatrix, std::shared_ptr<CDenseMatrix> >
        (
            "CDenseMatrix",
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
        .def(bp::self == bp::other<CDenseMatrix>())
    ;
}

} // bp_math namespace

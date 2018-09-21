//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef ExportMath_hpp
#define ExportMath_hpp

#include <boost/python/numpy.hpp>

#include "DenseMatrix.hpp"

namespace np = boost::python::numpy;

namespace bp_math { // bp_math namespace

/**
 Converts numpy array to CDenseMatrix

 @param arr the numpy array.
 @return a CDenseMatrix object.
 */
CDenseMatrix
CDenseMatrix_from_numpy(const np::ndarray& arr);

/**
 Exports classes/functions in src/math to python.
 */
void export_math();

} // bp_math namespace

#endif /* ExportMath_hpp */

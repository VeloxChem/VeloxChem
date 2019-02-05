//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef AODensityMatrixSetter_hpp
#define AODensityMatrixSetter_hpp

#include "DenseMatrix.hpp"
#include "AODensityMatrix.hpp"

namespace vlxden { // vlxden namespace
    
CAODensityMatrix getRestDensityMatrixForH2O();

CDenseMatrix getCoccMatrixForH2O();

CDenseMatrix getCvirMatrixForH2O();

CAODensityMatrix getGenOODensityMatrixForH2O();

} // vlxden namespace

#endif /* AODensityMatrixSetter_hpp */

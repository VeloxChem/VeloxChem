//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#ifndef AODensityMatrixSetter_hpp
#define AODensityMatrixSetter_hpp

#include "DenseMatrix.hpp"
#include "AODensityMatrix.hpp"

namespace vlxden { // vlxden namespace
    
CAODensityMatrix getRestDensityMatrixForH2O();

CDenseMatrix getCoccMatrixForH2O();

CDenseMatrix getCvirMatrixForH2O();

CAODensityMatrix getGenOODensityMatrixForH2O();

CAODensityMatrix getUnrestDensityMatrixForH2O();

} // vlxden namespace

#endif /* AODensityMatrixSetter_hpp */

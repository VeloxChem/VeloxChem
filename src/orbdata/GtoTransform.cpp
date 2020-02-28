//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "GtoTransform.hpp"

namespace gtotra {  // gtotra namespace
    
CDenseMatrix
to_veloxchem(const CDenseMatrix&    matrix,
             const CMolecularBasis& basis,
             const CMolecule&       molecule)
{
    auto vlxmat = matrix;
    
    // set up pointers to matrices data
    
    auto vlxdat = vlxmat.values();

    auto daldat = matrix.values();
    
    // set up ordering matrix
    
    auto idsmap = basis.getIndexMapForDalton(molecule);
    
    auto vlxidx = idsmap.data(0);
    
    auto dalidx = idsmap.data(1);
    
    // reorder square matrix
    
    auto ndim = matrix.getNumberOfColumns();
    
    for (int32_t i = 0; i < ndim; i++)
    {
        for (int32_t j = 0; j < ndim; j++)
        {
            vlxdat[vlxidx[i] * ndim + vlxidx[j]] = daldat[dalidx[i] * ndim + dalidx[j]];
        }
    }
    
    return vlxmat;
}
    
CDenseMatrix
to_dalton(const CDenseMatrix&    matrix,
          const CMolecularBasis& basis,
          const CMolecule&       molecule)
{
    auto dalmat = matrix;
    
    // set up pointers to matrices data
    
    auto vlxdat = matrix.values();
    
    auto daldat = dalmat.values();
    
    // set up ordering matrix
    
    auto idsmap = basis.getIndexMapForDalton(molecule);
    
    auto vlxidx = idsmap.data(0);
    
    auto dalidx = idsmap.data(1);
    
    // reorder square matrix
    
    auto ndim = matrix.getNumberOfColumns();
    
    for (int32_t i = 0; i < ndim; i++)
    {
        for (int32_t j = 0; j < ndim; j++)
        {
            daldat[dalidx[i] * ndim + dalidx[j]] = vlxdat[vlxidx[i] * ndim + vlxidx[j]];
        }
    }
    
    return dalmat;
}
    
}  // namespace gtotra

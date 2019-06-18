//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef GtoRecFunc_hpp
#define GtoRecFunc_hpp

#include <cstdint>

#include "GtoContainer.hpp"
#include "MemBlock2D.hpp"

namespace gtorec {  // gtorec namespace

/**
 Computes values of primitive S-type Gaussian functions at grid point (LDA).

 @param gtoValues the vector of primitive Gaussian function values at grid
 point.
 @param distances the vector of distances between Gaussian function centers
 and grid point.
 @param gtoContainer the GTO's container.
 @param redDimensions the vector of reduced dimensions.
 @param iGtoBlock the index of GTOs block.
 */
void compGtoTypeSForLDA(CMemBlock2D<double>&        gtoValues,
                        const CMemBlock2D<double>&  distances,
                        const CGtoContainer&        gtoContainer,
                        const CMemBlock2D<int32_t>& redDimensions,
                        const int32_t               iGtoBlock);

/**
 Computes values of primitive P-type Gaussian functions at grid point (LDA).

 @param gtoValues the vector of primitive Gaussian function values at grid
 point.
 @param distances the vector of distances between Gaussian function centers
 and grid point.
 @param redDimensions the vector of reduced dimensions.
 @param iGtoBlock the index of GTOs block.
 */
void compGtoTypePForLDA(CMemBlock2D<double>&        gtoValues,
                        const CMemBlock2D<double>&  distances,
                        const CMemBlock2D<int32_t>& redDimensions,
                        const int32_t               iGtoBlock);

/**
 Computes values of primitive D-type Gaussian functions at grid point (LDA).

 @param gtoValues the vector of primitive Gaussian function values at grid
 point.
 @param distances the vector of distances between Gaussian function centers
 and grid point.
 @param redDimensions the vector of reduced dimensions.
 @param iGtoBlock the index of GTOs block.
 */
void compGtoTypeDForLDA(CMemBlock2D<double>&        gtoValues,
                        const CMemBlock2D<double>&  distances,
                        const CMemBlock2D<int32_t>& redDimensions,
                        const int32_t               iGtoBlock);

/**
 Computes values of primitive F-type Gaussian functions at grid point (LDA).

 @param gtoValues the vector of primitive Gaussian function values at grid
 point.
 @param distances the vector of distances between Gaussian function centers
 and grid point.
 @param redDimensions the vector of reduced dimensions.
 @param iGtoBlock the index of GTOs block.
 */
void compGtoTypeFForLDA(CMemBlock2D<double>&        gtoValues,
                        const CMemBlock2D<double>&  distances,
                        const CMemBlock2D<int32_t>& redDimensions,
                        const int32_t               iGtoBlock);

/**
 Computes values of primitive G-type Gaussian functions at grid point (LDA).

 @param gtoValues the vector of primitive Gaussian function values at grid
 point.
 @param distances the vector of distances between Gaussian function centers
 and grid point.
 @param redDimensions the vector of reduced dimensions.
 @param iGtoBlock the index of GTOs block.
 */
void compGtoTypeGForLDA(CMemBlock2D<double>&        gtoValues,
                        const CMemBlock2D<double>&  distances,
                        const CMemBlock2D<int32_t>& redDimensions,
                        const int32_t               iGtoBlock);

/**
 Computes values of primitive S-type Gaussian functions at grid point (GGA).

 @param gtoValues the vector of primitive Gaussian function values at grid
 point.
 @param distances the vector of distances between Gaussian function centers
 and grid point.
 @param gtoContainer the GTO's container.
 @param redDimensions the vector of reduced dimensions.
 @param iGtoBlock the index of GTOs block.
 */
void compGtoTypeSForGGA(CMemBlock2D<double>&        gtoValues,
                        const CMemBlock2D<double>&  distances,
                        const CGtoContainer&        gtoContainer,
                        const CMemBlock2D<int32_t>& redDimensions,
                        const int32_t               iGtoBlock);

/**
 Computes values of primitive P-type Gaussian functions at grid point (GGA).

 @param gtoValues the vector of primitive Gaussian function values at grid
 point.
 @param distances the vector of distances between Gaussian function centers
 and grid point.
 @param redDimensions the vector of reduced dimensions.
 @param iGtoBlock the index of GTOs block.
 */
void compGtoTypePForGGA(CMemBlock2D<double>&        gtoValues,
                        const CMemBlock2D<double>&  distances,
                        const CMemBlock2D<int32_t>& redDimensions,
                        const int32_t               iGtoBlock);

/**
 Computes values of primitive D-type Gaussian functions at grid point (GGA).

 @param gtoValues the vector of primitive Gaussian function values at grid
 point.
 @param distances the vector of distances between Gaussian function centers
 and grid point.
 @param redDimensions the vector of reduced dimensions.
 @param iGtoBlock the index of GTOs block.
 */
void compGtoTypeDForGGA(CMemBlock2D<double>&        gtoValues,
                        const CMemBlock2D<double>&  distances,
                        const CMemBlock2D<int32_t>& redDimensions,
                        const int32_t               iGtoBlock);

/**
 Computes values of primitive F-type Gaussian functions at grid point (GGA).

 @param gtoValues the vector of primitive Gaussian function values at grid
 point.
 @param distances the vector of distances between Gaussian function centers
 and grid point.
 @param redDimensions the vector of reduced dimensions.
 @param iGtoBlock the index of GTOs block.
 */
void compGtoTypeFForGGA(CMemBlock2D<double>&        gtoValues,
                        const CMemBlock2D<double>&  distances,
                        const CMemBlock2D<int32_t>& redDimensions,
                        const int32_t               iGtoBlock);

/**
 Computes values of primitive G-type Gaussian functions at grid point (GGA).

 @param gtoValues the vector of primitive Gaussian function values at grid
 point.
 @param distances the vector of distances between Gaussian function centers
 and grid point.
 @param redDimensions the vector of reduced dimensions.
 @param iGtoBlock the index of GTOs block.
 */
void compGtoTypeGForGGA(CMemBlock2D<double>&        gtoValues,
                        const CMemBlock2D<double>&  distances,
                        const CMemBlock2D<int32_t>& redDimensions,
                        const int32_t               iGtoBlock);

/**
 Computes values of primitive S-type Gaussian functions at grid point (MGGA).

 @param gtoValues the vector of primitive Gaussian function values at grid
 point.
 @param distances the vector of distances between Gaussian function centers
 and grid point.
 @param gtoContainer the GTO's container.
 @param redDimensions the vector of reduced dimensions.
 @param iGtoBlock the index of GTOs block.
 */
void compGtoTypeSForMGGA(CMemBlock2D<double>&        gtoValues,
                         const CMemBlock2D<double>&  distances,
                         const CGtoContainer&        gtoContainer,
                         const CMemBlock2D<int32_t>& redDimensions,
                         const int32_t               iGtoBlock);

/**
 Computes values of primitive P-type Gaussian functions at grid point (MGGA).

 @param gtoValues the vector of primitive Gaussian function values at grid
 point.
 @param distances the vector of distances between Gaussian function centers
 and grid point.
 @param redDimensions the vector of reduced dimensions.
 @param iGtoBlock the index of GTOs block.
 */
void compGtoTypePForMGGA(CMemBlock2D<double>&        gtoValues,
                         const CMemBlock2D<double>&  distances,
                         const CMemBlock2D<int32_t>& redDimensions,
                         const int32_t               iGtoBlock);

/**
 Computes values of primitive D-type Gaussian functions at grid point (MGGA).

 @param gtoValues the vector of primitive Gaussian function values at grid
 point.
 @param distances the vector of distances between Gaussian function centers
 and grid point.
 @param redDimensions the vector of reduced dimensions.
 @param iGtoBlock the index of GTOs block.
 */
void compGtoTypeDForMGGA(CMemBlock2D<double>&        gtoValues,
                         const CMemBlock2D<double>&  distances,
                         const CMemBlock2D<int32_t>& redDimensions,
                         const int32_t               iGtoBlock);

/**
 Computes values of primitive F-type Gaussian functions at grid point (MGGA).

 @param gtoValues the vector of primitive Gaussian function values at grid
 point.
 @param distances the vector of distances between Gaussian function centers
 and grid point.
 @param redDimensions the vector of reduced dimensions.
 @param iGtoBlock the index of GTOs block.
 */
void compGtoTypeFForMGGA(CMemBlock2D<double>&        gtoValues,
                         const CMemBlock2D<double>&  distances,
                         const CMemBlock2D<int32_t>& redDimensions,
                         const int32_t               iGtoBlock);

/**
 Computes values of primitive G-type Gaussian functions at grid point (MGGA).

 @param gtoValues the vector of primitive Gaussian function values at grid
 point.
 @param distances the vector of distances between Gaussian function centers
 and grid point.
 @param redDimensions the vector of reduced dimensions.
 @param iGtoBlock the index of GTOs block.
 */
void compGtoTypeGForMGGA(CMemBlock2D<double>&        gtoValues,
                         const CMemBlock2D<double>&  distances,
                         const CMemBlock2D<int32_t>& redDimensions,
                         const int32_t               iGtoBlock);

}  // namespace gtorec

#endif /* GtoRecFunc_hpp */

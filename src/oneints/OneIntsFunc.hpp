//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef OneIntsFunc_hpp
#define OneIntsFunc_hpp

#include <cstdint>
#include <cmath>

#include "MemBlock2D.hpp"
#include "GtoBlock.hpp"
#include "GtoContainer.hpp"
#include "DenseMatrix.hpp"

namespace intsfunc { // intsfunc namespace
    
    /**
     Computes distances between specific contracted GTO on bra side with all
     contracted GTOs on ket side.

     @param abDistances the vector of Cartesian R(AB) = A - B distances.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compDistancesAB(      CMemBlock2D<double>& abDistances,
                         const CGtoBlock&           braGtoBlock,
                         const CGtoBlock&           ketGtoBlock,
                         const int32_t              iContrGto);
    
    /**
     Computes Obara-Saika factors for overlap integrals.

     @param osFactors the vector of Obara-Saika factors.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compFactorsForOverlap(      CMemBlock2D<double>& osFactors,
                               const CGtoBlock&           braGtoBlock,
                               const CGtoBlock&           ketGtoBlock,
                               const int32_t              iContrGto);
    
    /**
     Computes Obara-Saika factors for kinetic energy integrals.
     
     @param osFactors the vector of Obara-Saika factors.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compFactorsForKineticEnergy(      CMemBlock2D<double>& osFactors,
                                     const CGtoBlock&           braGtoBlock,
                                     const CGtoBlock&           ketGtoBlock,
                                     const int32_t              iContrGto);
    
    /**
     Computes Obara-Saika factors for nuclear potential integrals.
     
     @param osFactors the vector of Obara-Saika factors.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compFactorsForNuclearPotential(      CMemBlock2D<double>& osFactors,
                                        const CGtoBlock&           braGtoBlock,
                                        const CGtoBlock&           ketGtoBlock,
                                        const int32_t              iContrGto);
    
    /**
     Computes Obara-Saika factors for electronic potential integrals.
     
     @param osFactors the vector of Obara-Saika factors.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compFactorsForElectronicPotential(      CMemBlock2D<double>& osFactors,
                                           const CGtoBlock&           braGtoBlock,
                                           const CGtoBlock&           ketGtoBlock,
                                           const int32_t              iContrGto);
    
    /**
     Computes Obara-Saika factors for linear momentum integrals.
     
     @param osFactors the vector of Obara-Saika factors.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compFactorsForLinearMomentum(      CMemBlock2D<double>& osFactors,
                                      const CGtoBlock&           braGtoBlock,
                                      const CGtoBlock&           ketGtoBlock,
                                      const int32_t              iContrGto);
    
    /**
     Computes vector of distances between center P of combined primitive GTOs
     and center A of primitive GTO on bra side.
     
     @param paDistances the vector of Cartesian R(PA) = P - A distances.
     @param abDistances the vector of Cartesian R(AB) = A - B distances.
     @param osFactors the vector of Obara-Saika factors.
     @param nFactors the fundamental dimensions of Obara-Saika factors vectors.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compDistancesPA(      CMemBlock2D<double>& paDistances,
                         const CMemBlock2D<double>& abDistances,
                         const CMemBlock2D<double>& osFactors,
                         const int32_t              nFactors,
                         const CGtoBlock&           braGtoBlock,
                         const CGtoBlock&           ketGtoBlock,
                         const int32_t              iContrGto);
    
    /**
    Computes set of tensors of distances between center P of combined primitive
    GTOs and center A of primitive GTO on bra side.
    
    @param paDistances the set of tensors of Cartesian R(PA) = P - A distances.
    @param abDistances the vector of Cartesian R(AB) = A - B distances.
    @param osFactors the vector of Obara-Saika factors.
    @param nFactors the fundamental dimensions of Obara-Saika factors vectors.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compTensorsPA(      CMemBlock2D<double>& paDistances,
                       const CMemBlock2D<double>& abDistances,
                       const CMemBlock2D<double>& osFactors,
                       const int32_t              nFactors,
                       const CGtoBlock&           braGtoBlock,
                       const CGtoBlock&           ketGtoBlock,
                       const int32_t              iContrGto);
    
    /**
     Computes vector of distances between center P of combined primitive GTOs
     and center B of primitive GTO on ket side.
     
     @param pbDistances the vector of Cartesian R(PB) = P - B distances.
     @param abDistances the vector of Cartesian R(AB) = A - B distances.
     @param osFactors the vector of Obara-Saika factors.
     @param nFactors the fundamental dimension of Obara-Saika factors vectors.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compDistancesPB(      CMemBlock2D<double>& pbDistances,
                         const CMemBlock2D<double>& abDistances,
                         const CMemBlock2D<double>& osFactors,
                         const int32_t              nFactors,
                         const CGtoBlock&           braGtoBlock,
                         const CGtoBlock&           ketGtoBlock,
                         const int32_t              iContrGto);
    
    /**
    Computes set of tensors of distances between center P of combined primitive
    GTOs and center B of primitive GTO on ket side.
    
    @param pbDistances the set of tensors of Cartesian R(PB) = P - B distances.
    @param abDistances the vector of Cartesian R(AB) = A - B distances.
    @param osFactors the vector of Obara-Saika factors.
    @param nFactors the fundamental dimension of Obara-Saika factors vectors.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoBlock the GTOs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void compTensorsPB(      CMemBlock2D<double>& pbDistances,
                       const CMemBlock2D<double>& abDistances,
                       const CMemBlock2D<double>& osFactors,
                       const int32_t              nFactors,
                       const CGtoBlock&           braGtoBlock,
                       const CGtoBlock&           ketGtoBlock,
                       const int32_t              iContrGto);
    
    /**
     Computes coordinates of combined Gaussian function, which is obtained by
     applying Gaussian product rule to two Gaussian functions.

     @param pCoordinates the vector of coordinates for combined Gaussian
            functions.
     @param osFactors the vector of Obara-Saika factors.
     @param nFactors the fundamental dimension of Obara-Saika factors vectors.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compCoordinatesForP(      CMemBlock2D<double>& pCoordinates,
                             const CMemBlock2D<double>& osFactors,
                             const int32_t              nFactors,
                             const CGtoBlock&           braGtoBlock,
                             const CGtoBlock&           ketGtoBlock,
                             const int32_t              iContrGto);
    
    /**
     Computes vector of distances between center P of combined primitive GTOs
     and center A of primitive GTO on bra side.

     @param paDistances the vector of Cartesian R(PA) = P - A distances.
     @param pCoordinates the vector of coordinates for combined Gaussian
            functions.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compDistancesPA(      CMemBlock2D<double>& paDistances,
                         const CMemBlock2D<double>& pCoordinates,
                         const CGtoBlock&           braGtoBlock,
                         const CGtoBlock&           ketGtoBlock,
                         const int32_t              iContrGto);
    
    /**
     Computes set of tensors of distances between center P of combined primitive
     GTOs and center A of primitive GTO on bra side.
     
     @param paDistances the set of tensors of Cartesian R(PA) = P - A distances.
     @param pCoordinates the vector of coordinates for combined Gaussian
     functions.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compTensorsPA(      CMemBlock2D<double>& paDistances,
                       const CMemBlock2D<double>& pCoordinates,
                       const CGtoBlock&           braGtoBlock,
                       const CGtoBlock&           ketGtoBlock,
                       const int32_t              iContrGto);
    
    /**
     Computes vector of distances between center P of combined primitive GTOs
     and center B of primitive GTO on bra side.
     
     @param pbDistances the vector of Cartesian R(PB) = P - B distances.
     @param pCoordinates the vector of coordinates for combined Gaussian
     functions.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compDistancesPB(      CMemBlock2D<double>& pbDistances,
                         const CMemBlock2D<double>& pCoordinates,
                         const CGtoBlock&           braGtoBlock,
                         const CGtoBlock&           ketGtoBlock,
                         const int32_t              iContrGto);
    
    /**
     Computes set of tensors of distances between center P of combined primitive
     GTOs and center B of primitive GTO on bra side.
     
     @param pbDistances the set of tensors of Cartesian R(PB) = P - B distances.
     @param pCoordinates the vector of coordinates for combined Gaussian
     functions.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compTensorsPB(      CMemBlock2D<double>& pbDistances,
                       const CMemBlock2D<double>& pCoordinates,
                       const CGtoBlock&           braGtoBlock,
                       const CGtoBlock&           ketGtoBlock,
                       const int32_t              iContrGto);
    
    /**
     Computes set of tensors of distances between center P of combined primitive
     GTOs and center C of specific point charge.
     
     @param pcDistances the set of tensors of Cartesian R(PC) = P - C distances.
     @param pCoordinates the vector of coordinates for combined Gaussian
     functions.
     @param cCoordinates the vector of coordinates of point charges.
     @param orderOfTensor the max. order of R(PC) tensors.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     @param iPointCharge the index of point charge in vector point charges
     coordinates.
     */
    void compTensorsPC(      CMemBlock2D<double>& pcDistances,
                       const CMemBlock2D<double>& pCoordinates,
                       const CMemBlock2D<double>& cCoordinates,
                       const int32_t              orderOfTensor,
                       const CGtoBlock&           braGtoBlock,
                       const CGtoBlock&           ketGtoBlock,
                       const int32_t              iContrGto,
                       const int32_t              iPointCharge);
    
    /**
     Computes vector of distances between center P of combined primitive GTOs
     and center C of specific point charge.
     
     @param pcDistances the vector of Cartesian R(PB) = P - C distances.
     @param pCoordinates the vector of coordinates for combined Gaussian
            functions.
     @param cCoordinates the vector of coordinates of point charges.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     @param iPointCharge the index of point charge in vector point charges
            coordinates.
     */
    void compDistancesPC(      CMemBlock2D<double>& pcDistances,
                         const CMemBlock2D<double>& pCoordinates,
                         const CMemBlock2D<double>& cCoordinates,
                         const CGtoBlock&           braGtoBlock,
                         const CGtoBlock&           ketGtoBlock,
                         const int32_t              iContrGto,
                         const int32_t              iPointCharge);
    
    
    /**
     Computes set of tensors of distances between center P of combined primitive
     GTOs and center C.

     @param pcDistances the set of tensors of Cartesian R(PC) = P - C distances.
     @param pCoordinates the vector of coordinates for combined Gaussian
     functions.
     @param xCoordinateC the Cartesian X coordinate of center C.
     @param yCoordinateC the Cartesian Y coordinate of center C.
     @param zCoordinateC the Cartesian Z coordinate of center C.
     @param orderOfTensor the max. order of R(PC) tensors.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compTensorsPC(      CMemBlock2D<double>& pcDistances,
                       const CMemBlock2D<double>& pCoordinates,
                       const double               xCoordinateC,
                       const double               yCoordinateC,
                       const double               zCoordinateC,
                       const int32_t              orderOfTensor,
                       const CGtoBlock&           braGtoBlock,
                       const CGtoBlock&           ketGtoBlock,
                       const int32_t              iContrGto);
    
    /**
     Computes vector of distances between center P of combined primitive GTOs
     and center C.
     
     @param pcDistances the vector of Cartesian R(PB) = P - C distances.
     @param pCoordinates the vector of coordinates for combined Gaussian
     functions.
     @param xCoordinateC the Cartesian X coordinate of center C.
     @param yCoordinateC the Cartesian Y coordinate of center C.
     @param zCoordinateC the Cartesian Z coordinate of center C.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compDistancesPC(      CMemBlock2D<double>& pcDistances,
                         const CMemBlock2D<double>& pCoordinates,
                         const double               xCoordinateC,
                         const double               yCoordinateC,
                         const double               zCoordinateC,
                         const CGtoBlock&           braGtoBlock,
                         const CGtoBlock&           ketGtoBlock,
                         const int32_t              iContrGto);
    
    /**
     Gets number of components in distances tensor.

     @param tensorOrder the order of distances tensor.
     @return the number of components in distances tensor.
     */
    int32_t getNumberOfComponentsInDistancesTensor(const int32_t tensorOrder);
    
    /**
     Computes tensor product of two tensors.

     @param aTensor the tensor A = B * C.
     @param bTensor the tensor B.
     @param cTensor the tensor C.
     @param bDimensions the tensorial dimensions of tensor B.
     @param cDimensions the tensorial dimensions of tensor C.
     @param nBlocks the number of tensorial blocks.
     */
    void compTensorsProduct(      CMemBlock2D<double>& aTensor,
                            const CMemBlock2D<double>& bTensor,
                            const CMemBlock2D<double>& cTensor,
                            const int32_t              bDimensions,
                            const int32_t              cDimensions,
                            const int32_t              nBlocks);
    
} // intsfunc namespace

#endif /* OneIntsFunc_hpp */

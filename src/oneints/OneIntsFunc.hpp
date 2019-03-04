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
     coordinates.
     */
    void compDistancesPC(      CMemBlock2D<double>& pcDistances,
                         const CMemBlock2D<double>& pCoordinates,
                         const double               xCoordinateC,
                         const double               yCoordinateC,
                         const double               zCoordinateC,
                         const CGtoBlock&           braGtoBlock,
                         const CGtoBlock&           ketGtoBlock,
                         const int32_t              iContrGto);
    
} // intsfunc namespace

#endif /* OneIntsFunc_hpp */

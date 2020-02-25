//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include <cstdint>

#include "GtoBlock.hpp"
#include "MemBlock2D.hpp"
#include "RecursionMap.hpp"

namespace kinrecfunc {  // kinrecfunc namespace

/**
Computes batch of primitive (P|T|P) kinetic energy integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compKineticEnergyForPP(CMemBlock2D<double>&       primBuffer,
                            const CRecursionMap&       recursionMap,
                            const CMemBlock2D<double>& osFactors,
                            const CMemBlock2D<double>& paDistances,
                            const CGtoBlock&           braGtoBlock,
                            const CGtoBlock&           ketGtoBlock,
                            const int32_t              iContrGto);

/**
Computes batch of primitive (P|T|D) kinetic energy integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compKineticEnergyForPD(CMemBlock2D<double>&       primBuffer,
                            const CRecursionMap&       recursionMap,
                            const CMemBlock2D<double>& osFactors,
                            const CMemBlock2D<double>& paDistances,
                            const CGtoBlock&           braGtoBlock,
                            const CGtoBlock&           ketGtoBlock,
                            const int32_t              iContrGto);

/**
Computes batch of primitive (D|T|P) kinetic energy integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compKineticEnergyForDP(CMemBlock2D<double>&       primBuffer,
                            const CRecursionMap&       recursionMap,
                            const CMemBlock2D<double>& osFactors,
                            const CMemBlock2D<double>& paDistances,
                            const CGtoBlock&           braGtoBlock,
                            const CGtoBlock&           ketGtoBlock,
                            const int32_t              iContrGto);

/**
Computes batch of primitive (P|T|F) kinetic energy integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compKineticEnergyForPF(CMemBlock2D<double>&       primBuffer,
                            const CRecursionMap&       recursionMap,
                            const CMemBlock2D<double>& osFactors,
                            const CMemBlock2D<double>& paDistances,
                            const CGtoBlock&           braGtoBlock,
                            const CGtoBlock&           ketGtoBlock,
                            const int32_t              iContrGto);

/**
Computes batch of primitive (F|T|P) kinetic energy integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compKineticEnergyForFP(CMemBlock2D<double>&       primBuffer,
                            const CRecursionMap&       recursionMap,
                            const CMemBlock2D<double>& osFactors,
                            const CMemBlock2D<double>& paDistances,
                            const CGtoBlock&           braGtoBlock,
                            const CGtoBlock&           ketGtoBlock,
                            const int32_t              iContrGto);

/**
Computes batch of primitive (P|T|G) kinetic energy integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compKineticEnergyForPG(CMemBlock2D<double>&       primBuffer,
                            const CRecursionMap&       recursionMap,
                            const CMemBlock2D<double>& osFactors,
                            const CMemBlock2D<double>& paDistances,
                            const CGtoBlock&           braGtoBlock,
                            const CGtoBlock&           ketGtoBlock,
                            const int32_t              iContrGto);

/**
Computes batch of primitive (G|T|P) kinetic energy integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compKineticEnergyForGP(CMemBlock2D<double>&       primBuffer,
                            const CRecursionMap&       recursionMap,
                            const CMemBlock2D<double>& osFactors,
                            const CMemBlock2D<double>& paDistances,
                            const CGtoBlock&           braGtoBlock,
                            const CGtoBlock&           ketGtoBlock,
                            const int32_t              iContrGto);

}  // namespace kinrecfunc

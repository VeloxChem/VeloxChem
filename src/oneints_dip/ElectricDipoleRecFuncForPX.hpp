//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include <cstdint>

#include "GtoBlock.hpp"
#include "MemBlock2D.hpp"
#include "RecursionMap.hpp"

namespace ediprecfunc {  // ediprecfunc namespace

/**
Computes batch of primitive (P|D|P) electric dipole integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param nOSFactors the number of Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compElectricDipoleForPP(CMemBlock2D<double>&       primBuffer,
                             const CRecursionMap&       recursionMap,
                             const CMemBlock2D<double>& osFactors,
                             const int32_t              nOSFactors,
                             const CMemBlock2D<double>& paDistances,
                             const CGtoBlock&           braGtoBlock,
                             const CGtoBlock&           ketGtoBlock,
                             const int32_t              iContrGto);

/**
Computes batch of primitive (P|D|D) electric dipole integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param nOSFactors the number of Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compElectricDipoleForPD(CMemBlock2D<double>&       primBuffer,
                             const CRecursionMap&       recursionMap,
                             const CMemBlock2D<double>& osFactors,
                             const int32_t              nOSFactors,
                             const CMemBlock2D<double>& paDistances,
                             const CGtoBlock&           braGtoBlock,
                             const CGtoBlock&           ketGtoBlock,
                             const int32_t              iContrGto);

/**
Computes sub-batch (0,27) of primitive (P|D|D) electric dipole integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param nOSFactors the number of Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compElectricDipoleForPD_0_27(CMemBlock2D<double>&       primBuffer,
                                  const CRecursionMap&       recursionMap,
                                  const CMemBlock2D<double>& osFactors,
                                  const int32_t              nOSFactors,
                                  const CMemBlock2D<double>& paDistances,
                                  const CGtoBlock&           braGtoBlock,
                                  const CGtoBlock&           ketGtoBlock,
                                  const int32_t              iContrGto);

/**
Computes sub-batch (27,54) of primitive (P|D|D) electric dipole integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param nOSFactors the number of Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compElectricDipoleForPD_27_54(CMemBlock2D<double>&       primBuffer,
                                   const CRecursionMap&       recursionMap,
                                   const CMemBlock2D<double>& osFactors,
                                   const int32_t              nOSFactors,
                                   const CMemBlock2D<double>& paDistances,
                                   const CGtoBlock&           braGtoBlock,
                                   const CGtoBlock&           ketGtoBlock,
                                   const int32_t              iContrGto);

/**
Computes batch of primitive (D|D|P) electric dipole integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param nOSFactors the number of Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compElectricDipoleForDP(CMemBlock2D<double>&       primBuffer,
                             const CRecursionMap&       recursionMap,
                             const CMemBlock2D<double>& osFactors,
                             const int32_t              nOSFactors,
                             const CMemBlock2D<double>& paDistances,
                             const CGtoBlock&           braGtoBlock,
                             const CGtoBlock&           ketGtoBlock,
                             const int32_t              iContrGto);

/**
Computes sub-batch (0,27) of primitive (D|D|P) electric dipole integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param nOSFactors the number of Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compElectricDipoleForDP_0_27(CMemBlock2D<double>&       primBuffer,
                                  const CRecursionMap&       recursionMap,
                                  const CMemBlock2D<double>& osFactors,
                                  const int32_t              nOSFactors,
                                  const CMemBlock2D<double>& paDistances,
                                  const CGtoBlock&           braGtoBlock,
                                  const CGtoBlock&           ketGtoBlock,
                                  const int32_t              iContrGto);

/**
Computes sub-batch (27,54) of primitive (D|D|P) electric dipole integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param nOSFactors the number of Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compElectricDipoleForDP_27_54(CMemBlock2D<double>&       primBuffer,
                                   const CRecursionMap&       recursionMap,
                                   const CMemBlock2D<double>& osFactors,
                                   const int32_t              nOSFactors,
                                   const CMemBlock2D<double>& paDistances,
                                   const CGtoBlock&           braGtoBlock,
                                   const CGtoBlock&           ketGtoBlock,
                                   const int32_t              iContrGto);

/**
Computes batch of primitive (P|D|F) electric dipole integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param nOSFactors the number of Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compElectricDipoleForPF(CMemBlock2D<double>&       primBuffer,
                             const CRecursionMap&       recursionMap,
                             const CMemBlock2D<double>& osFactors,
                             const int32_t              nOSFactors,
                             const CMemBlock2D<double>& paDistances,
                             const CGtoBlock&           braGtoBlock,
                             const CGtoBlock&           ketGtoBlock,
                             const int32_t              iContrGto);

/**
Computes sub-batch (0,45) of primitive (P|D|F) electric dipole integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param nOSFactors the number of Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compElectricDipoleForPF_0_45(CMemBlock2D<double>&       primBuffer,
                                  const CRecursionMap&       recursionMap,
                                  const CMemBlock2D<double>& osFactors,
                                  const int32_t              nOSFactors,
                                  const CMemBlock2D<double>& paDistances,
                                  const CGtoBlock&           braGtoBlock,
                                  const CGtoBlock&           ketGtoBlock,
                                  const int32_t              iContrGto);

/**
Computes sub-batch (45,90) of primitive (P|D|F) electric dipole integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param nOSFactors the number of Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compElectricDipoleForPF_45_90(CMemBlock2D<double>&       primBuffer,
                                   const CRecursionMap&       recursionMap,
                                   const CMemBlock2D<double>& osFactors,
                                   const int32_t              nOSFactors,
                                   const CMemBlock2D<double>& paDistances,
                                   const CGtoBlock&           braGtoBlock,
                                   const CGtoBlock&           ketGtoBlock,
                                   const int32_t              iContrGto);

/**
Computes batch of primitive (F|D|P) electric dipole integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param nOSFactors the number of Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compElectricDipoleForFP(CMemBlock2D<double>&       primBuffer,
                             const CRecursionMap&       recursionMap,
                             const CMemBlock2D<double>& osFactors,
                             const int32_t              nOSFactors,
                             const CMemBlock2D<double>& paDistances,
                             const CGtoBlock&           braGtoBlock,
                             const CGtoBlock&           ketGtoBlock,
                             const int32_t              iContrGto);

/**
Computes sub-batch (0,45) of primitive (F|D|P) electric dipole integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param nOSFactors the number of Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compElectricDipoleForFP_0_45(CMemBlock2D<double>&       primBuffer,
                                  const CRecursionMap&       recursionMap,
                                  const CMemBlock2D<double>& osFactors,
                                  const int32_t              nOSFactors,
                                  const CMemBlock2D<double>& paDistances,
                                  const CGtoBlock&           braGtoBlock,
                                  const CGtoBlock&           ketGtoBlock,
                                  const int32_t              iContrGto);

/**
Computes sub-batch (45,90) of primitive (F|D|P) electric dipole integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param nOSFactors the number of Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compElectricDipoleForFP_45_90(CMemBlock2D<double>&       primBuffer,
                                   const CRecursionMap&       recursionMap,
                                   const CMemBlock2D<double>& osFactors,
                                   const int32_t              nOSFactors,
                                   const CMemBlock2D<double>& paDistances,
                                   const CGtoBlock&           braGtoBlock,
                                   const CGtoBlock&           ketGtoBlock,
                                   const int32_t              iContrGto);

/**
Computes batch of primitive (P|D|G) electric dipole integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param nOSFactors the number of Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compElectricDipoleForPG(CMemBlock2D<double>&       primBuffer,
                             const CRecursionMap&       recursionMap,
                             const CMemBlock2D<double>& osFactors,
                             const int32_t              nOSFactors,
                             const CMemBlock2D<double>& paDistances,
                             const CGtoBlock&           braGtoBlock,
                             const CGtoBlock&           ketGtoBlock,
                             const int32_t              iContrGto);

/**
Computes sub-batch (0,45) of primitive (P|D|G) electric dipole integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param nOSFactors the number of Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compElectricDipoleForPG_0_45(CMemBlock2D<double>&       primBuffer,
                                  const CRecursionMap&       recursionMap,
                                  const CMemBlock2D<double>& osFactors,
                                  const int32_t              nOSFactors,
                                  const CMemBlock2D<double>& paDistances,
                                  const CGtoBlock&           braGtoBlock,
                                  const CGtoBlock&           ketGtoBlock,
                                  const int32_t              iContrGto);

/**
Computes sub-batch (45,90) of primitive (P|D|G) electric dipole integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param nOSFactors the number of Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compElectricDipoleForPG_45_90(CMemBlock2D<double>&       primBuffer,
                                   const CRecursionMap&       recursionMap,
                                   const CMemBlock2D<double>& osFactors,
                                   const int32_t              nOSFactors,
                                   const CMemBlock2D<double>& paDistances,
                                   const CGtoBlock&           braGtoBlock,
                                   const CGtoBlock&           ketGtoBlock,
                                   const int32_t              iContrGto);

/**
Computes sub-batch (90,135) of primitive (P|D|G) electric dipole integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param nOSFactors the number of Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compElectricDipoleForPG_90_135(CMemBlock2D<double>&       primBuffer,
                                    const CRecursionMap&       recursionMap,
                                    const CMemBlock2D<double>& osFactors,
                                    const int32_t              nOSFactors,
                                    const CMemBlock2D<double>& paDistances,
                                    const CGtoBlock&           braGtoBlock,
                                    const CGtoBlock&           ketGtoBlock,
                                    const int32_t              iContrGto);

/**
Computes batch of primitive (G|D|P) electric dipole integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param nOSFactors the number of Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compElectricDipoleForGP(CMemBlock2D<double>&       primBuffer,
                             const CRecursionMap&       recursionMap,
                             const CMemBlock2D<double>& osFactors,
                             const int32_t              nOSFactors,
                             const CMemBlock2D<double>& paDistances,
                             const CGtoBlock&           braGtoBlock,
                             const CGtoBlock&           ketGtoBlock,
                             const int32_t              iContrGto);

/**
Computes sub-batch (0,45) of primitive (G|D|P) electric dipole integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param nOSFactors the number of Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compElectricDipoleForGP_0_45(CMemBlock2D<double>&       primBuffer,
                                  const CRecursionMap&       recursionMap,
                                  const CMemBlock2D<double>& osFactors,
                                  const int32_t              nOSFactors,
                                  const CMemBlock2D<double>& paDistances,
                                  const CGtoBlock&           braGtoBlock,
                                  const CGtoBlock&           ketGtoBlock,
                                  const int32_t              iContrGto);

/**
Computes sub-batch (45,90) of primitive (G|D|P) electric dipole integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param nOSFactors the number of Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compElectricDipoleForGP_45_90(CMemBlock2D<double>&       primBuffer,
                                   const CRecursionMap&       recursionMap,
                                   const CMemBlock2D<double>& osFactors,
                                   const int32_t              nOSFactors,
                                   const CMemBlock2D<double>& paDistances,
                                   const CGtoBlock&           braGtoBlock,
                                   const CGtoBlock&           ketGtoBlock,
                                   const int32_t              iContrGto);

/**
Computes sub-batch (90,135) of primitive (G|D|P) electric dipole integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param nOSFactors the number of Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compElectricDipoleForGP_90_135(CMemBlock2D<double>&       primBuffer,
                                    const CRecursionMap&       recursionMap,
                                    const CMemBlock2D<double>& osFactors,
                                    const int32_t              nOSFactors,
                                    const CMemBlock2D<double>& paDistances,
                                    const CGtoBlock&           braGtoBlock,
                                    const CGtoBlock&           ketGtoBlock,
                                    const int32_t              iContrGto);

}  // namespace ediprecfunc

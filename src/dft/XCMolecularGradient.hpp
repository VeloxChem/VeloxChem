//
//                           VELOXCHEM 1.0-RC2
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact
//
//  SPDX-License-Identifier: LGPL-3.0-or-later
//
//  This file is part of VeloxChem.
//
//  VeloxChem is free software: you can redistribute it and/or modify it under
//  the terms of the GNU Lesser General Public License as published by the Free
//  Software Foundation, either version 3 of the License, or (at your option)
//  any later version.
//
//  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
//  License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.

#ifndef XCMolecularGradient_hpp
#define XCMolecularGradient_hpp

#include <mpi.h>

#include <cstdint>
#include <string>
#include <vector>

#include "AODensityMatrix.hpp"
#include "DenseMatrix.hpp"
#include "DensityGrid.hpp"
#include "DensityGridQuad.hpp"
#include "GtoContainer.hpp"
#include "MolecularBasis.hpp"
#include "MolecularGrid.hpp"
#include "Molecule.hpp"
#include "XCCubicHessianGrid.hpp"
#include "XCGradientGrid.hpp"
#include "XCHessianGrid.hpp"

/**
 Class CXCMolecularGradient implements integration of exchange-correlation contribution to molecular
 gradient.

 @author Z. Rinkevicius
 */
class CXCMolecularGradient
{
    /**
     The rank of associated local MPI process.
     */
    int32_t _locRank;

    /**
     The total number of local MPI processes.
     */
    int32_t _locNodes;

    /**
     The MPI communicator.
     */
    MPI_Comm _locComm;

    /**
     The threshold of density screening.
     */
    double _thresholdOfDensity;

    /**
     Integrates first-order exchnage-correlation functional contribution to
     molecular gradient.

     @param molecularGradient the molecular gradient object.
     @param molecule the molecule.
     @param basis the molecular basis.
     @param xcFuncType the type of exchange-correlation functional.
     @param densityMatrix the AO density matrix (to be contracted with GTO
            gradient).
     @param molecularGrid the molecular grid.
     @param gsDenistyGrid the ground state density grid.
     @param xcGradientGrid the exchange-correlation gradient grid.
     */
    void _compVxcContrib(CDenseMatrix&           molecularGradient,
                         const CMolecule&        molecule,
                         const CMolecularBasis&  basis,
                         const xcfun             xcFuncType,
                         const CAODensityMatrix& densityMatrix,
                         const CMolecularGrid&   molecularGrid,
                         const CDensityGrid&     gsDensityGrid,
                         const CXCGradientGrid&  xcGradientGrid) const;

    /**
     Integrates second-order exchnage-correlation functional contribution to
     molecular gradient.

     @param molecularGradient the molecular gradient object.
     @param molecule the molecule.
     @param basis the molecular basis.
     @param xcFuncType the type of exchange-correlation functional.
     @param densityMatrix the AO density matrix (to be contracted with GTO
            gradient).
     @param molecularGrid the molecular grid.
     @param gsDenistyGrid the ground state density grid.
     @param rwDenistyGrid the perturbed density grid.
     @param xcGradientGrid the exchange-correlation gradient grid.
     @param xcHessianGrid the exchange-correlation hessian grid.
     */
    void _compFxcContrib(CDenseMatrix&           molecularGradient,
                         const CMolecule&        molecule,
                         const CMolecularBasis&  basis,
                         const xcfun             xcFuncType,
                         const CAODensityMatrix& densityMatrix,
                         const CMolecularGrid&   molecularGrid,
                         const CDensityGrid&     gsDensityGrid,
                         const CDensityGrid&     rwDensityGrid,
                         const CXCGradientGrid&  xcGradientGrid,
                         const CXCHessianGrid&   xcHessianGrid) const;

    /**
     Integrates third-order exchnage-correlation functional contribution to
     molecular gradient.

     @param molecularGradient the molecular gradient object.
     @param molecule the molecule.
     @param basis the molecular basis.
     @param xcFuncType the type of exchange-correlation functional.
     @param densityMatrix the AO density matrix (to be contracted with GTO
            gradient).
     @param molecularGrid the molecular grid.
     @param gsDenistyGrid the ground state density grid.
     @param rwDenistyGridQuad the perturbed density grid for quadratic
            response.
     @param xcGradientGrid the exchange-correlation gradient grid.
     @param xcHessianGrid the exchange-correlation hessian grid.
     @param xcCubicHessianGrid the exchange-correlation cubic hessian grid.
     */
    void _compGxcContrib(CDenseMatrix&              molecularGradient,
                         const CMolecule&           molecule,
                         const CMolecularBasis&     basis,
                         const xcfun                xcFuncType,
                         const CAODensityMatrix&    densityMatrix,
                         const CMolecularGrid&      molecularGrid,
                         const CDensityGrid&        gsDensityGrid,
                         const CDensityGridQuad&    rwDensityGridQuad,
                         const CXCGradientGrid&     xcGradientGrid,
                         const CXCHessianGrid&      xcHessianGrid,
                         const CXCCubicHessianGrid& xcCubicHessianGrid) const;

    /**
     Integrates TDDFT exchnage-correlation functional contribution to molecular
     gradient.

     @param molecularGradient the molecular gradient object.
     @param molecule the molecule.
     @param basis the molecular basis.
     @param xcFuncType the type of exchange-correlation functional.
     @param gsDensityMatrix the ground state AO density matrix.
     @param rwDensityMatrixOne the perturbed AO density matrix (relaxed_density_ao).
     @param rwDensityMatrixTwo the perturbed AO density matrix (x_minus_y_ao).
     @param molecularGrid the molecular grid.
     @param gsDenistyGrid the ground state density grid.
     @param rwDenistyGridOne the perturbed density grid from rwDensityMatrixOne.
     @param rwDenistyGridTwo the perturbed density grid from rwDensityMatrixTwo.
     @param rwDenistyGridQuad the perturbed density grid for quadratic
            response (QRF).
     @param xcGradientGrid the exchange-correlation gradient grid.
     @param xcHessianGrid the exchange-correlation hessian grid.
     @param xcCubicHessianGrid the exchange-correlation cubic hessian grid.
     */
    void _compTddftContrib(CDenseMatrix&              molecularGradient,
                           const CMolecule&           molecule,
                           const CMolecularBasis&     basis,
                           const xcfun                xcFuncType,
                           const CAODensityMatrix&    gsDensityMatrix,
                           const CAODensityMatrix&    rwdensityMatrixOne,
                           const CAODensityMatrix&    rwdensityMatrixTwo,
                           const CMolecularGrid&      molecularGrid,
                           const CDensityGrid&        gsDensityGrid,
                           const CDensityGrid&        rwDensityGridOne,
                           const CDensityGrid&        rwDensityGridTwo,
                           const CDensityGridQuad&    rwDensityGridQuad,
                           const CXCGradientGrid&     xcGradientGrid,
                           const CXCHessianGrid&      xcHessianGrid,
                           const CXCCubicHessianGrid& xcCubicHessianGrid) const;

    /**
     Integrates batch of grid points for LDA first-order exchnage-correlation
     functional contribution to molecular gradient.

     @param molecularGradient the molecular gradient object.
     @param densityMatrix the AO density matrix (to be contracted with GTO
            gradient).
     @param molecule the molecule.
     @param basis the molecular basis.
     @param molecularGrid the molecular grid.
     @param gsDenistyGrid the ground state density grid.
     @param xcGradientGrid the exchange-correlation gradient grid.
     @param gridOffset the batch offset in vector grid points.
     @param nGridPoints the number of grid points in batch.
     */
    void _compVxcBatchForLDA(CDenseMatrix&           molecularGradient,
                             const CAODensityMatrix& densityMatrix,
                             const CMolecule&        molecule,
                             const CMolecularBasis&  basis,
                             const CMolecularGrid&   molecularGrid,
                             const CDensityGrid&     gsDensityGrid,
                             const CXCGradientGrid&  xcGradientGrid,
                             const int32_t           gridOffset,
                             const int32_t           nGridPoints) const;

    /**
     Integrates batch of grid points for GGA first-order exchnage-correlation
     functional contribution to molecular gradient.

     @param molecularGradient the molecular gradient object.
     @param densityMatrix the AO density matrix (to be contracted with GTO
            gradient).
     @param molecule the molecule.
     @param basis the molecular basis.
     @param molecularGrid the molecular grid.
     @param gsDenistyGrid the ground state density grid.
     @param xcGradientGrid the exchange-correlation gradient grid.
     @param gridOffset the batch offset in vector grid points.
     @param nGridPoints the number of grid points in batch.
     */
    void _compVxcBatchForGGA(CDenseMatrix&           molecularGradient,
                             const CAODensityMatrix& densityMatrix,
                             const CMolecule&        molecule,
                             const CMolecularBasis&  basis,
                             const CMolecularGrid&   molecularGrid,
                             const CDensityGrid&     gsDensityGrid,
                             const CXCGradientGrid&  xcGradientGrid,
                             const int32_t           gridOffset,
                             const int32_t           nGridPoints) const;

    /**
     Integrates batch of grid points for LDA second-order exchnage-correlation
     functional contribution to molecular gradient.

     @param molecularGradient the molecular gradient object.
     @param densityMatrix the AO density matrix (to be contracted with GTO
            gradient).
     @param molecule the molecule.
     @param basis the molecular basis.
     @param molecularGrid the molecular grid.
     @param gsDenistyGrid the ground state density grid.
     @param rwDenistyGrid the perturbed density grid.
     @param xcGradientGrid the exchange-correlation gradient grid.
     @param xcHessianGrid the exchange-correlation hessian grid.
     @param gridOffset the batch offset in vector grid points.
     @param nGridPoints the number of grid points in batch.
     */
    void _compFxcBatchForLDA(CDenseMatrix&           molecularGradient,
                             const CAODensityMatrix& densityMatrix,
                             const CMolecule&        molecule,
                             const CMolecularBasis&  basis,
                             const CMolecularGrid&   molecularGrid,
                             const CDensityGrid&     gsDensityGrid,
                             const CDensityGrid&     rwDensityGrid,
                             const CXCGradientGrid&  xcGradientGrid,
                             const CXCHessianGrid&   xcHessianGrid,
                             const int32_t           gridOffset,
                             const int32_t           nGridPoints) const;

    /**
     Integrates batch of grid points for GGA second-order exchnage-correlation
     functional contribution to molecular gradient.

     @param molecularGradient the molecular gradient object.
     @param densityMatrix the AO density matrix (to be contracted with GTO
            gradient).
     @param molecule the molecule.
     @param basis the molecular basis.
     @param molecularGrid the molecular grid.
     @param gsDenistyGrid the ground state density grid.
     @param rwDenistyGrid the perturbed density grid.
     @param xcGradientGrid the exchange-correlation gradient grid.
     @param xcHessianGrid the exchange-correlation hessian grid.
     @param gridOffset the batch offset in vector grid points.
     @param nGridPoints the number of grid points in batch.
     */
    void _compFxcBatchForGGA(CDenseMatrix&           molecularGradient,
                             const CAODensityMatrix& densityMatrix,
                             const CMolecule&        molecule,
                             const CMolecularBasis&  basis,
                             const CMolecularGrid&   molecularGrid,
                             const CDensityGrid&     gsDensityGrid,
                             const CDensityGrid&     rwDensityGrid,
                             const CXCGradientGrid&  xcGradientGrid,
                             const CXCHessianGrid&   xcHessianGrid,
                             const int32_t           gridOffset,
                             const int32_t           nGridPoints) const;

    /**
     Integrates batch of grid points for LDA third-order exchnage-correlation
     functional contribution to molecular gradient.

     @param molecularGradient the molecular gradient object.
     @param densityMatrix the AO density matrix (to be contracted with GTO
            gradient).
     @param molecule the molecule.
     @param basis the molecular basis.
     @param molecularGrid the molecular grid.
     @param gsDenistyGrid the ground state density grid.
     @param rwDenistyGridQuad the perturbed density grid for quadratic
            response.
     @param xcGradientGrid the exchange-correlation gradient grid.
     @param xcHessianGrid the exchange-correlation hessian grid.
     @param xcCubicHessianGrid the exchange-correlation cubic hessian grid.
     @param gridOffset the batch offset in vector grid points.
     @param nGridPoints the number of grid points in batch.
     */
    void _compGxcBatchForLDA(CDenseMatrix&              molecularGradient,
                             const CAODensityMatrix&    densityMatrix,
                             const CMolecule&           molecule,
                             const CMolecularBasis&     basis,
                             const CMolecularGrid&      molecularGrid,
                             const CDensityGrid&        gsDensityGrid,
                             const CDensityGridQuad&    rwDensityGridQuad,
                             const CXCGradientGrid&     xcGradientGrid,
                             const CXCHessianGrid&      xcHessianGrid,
                             const CXCCubicHessianGrid& xcCubicHessianGrid,
                             const int32_t              gridOffset,
                             const int32_t              nGridPoints) const;

    /**
     Integrates batch of grid points for GGA third-order exchnage-correlation
     functional contribution to molecular gradient.

     @param molecularGradient the molecular gradient object.
     @param densityMatrix the AO density matrix (to be contracted with GTO
            gradient).
     @param molecule the molecule.
     @param basis the molecular basis.
     @param molecularGrid the molecular grid.
     @param gsDenistyGrid the ground state density grid.
     @param rwDenistyGridQuad the perturbed density grid for quadratic
            response.
     @param xcGradientGrid the exchange-correlation gradient grid.
     @param xcHessianGrid the exchange-correlation hessian grid.
     @param xcCubicHessianGrid the exchange-correlation cubic hessian grid.
     @param gridOffset the batch offset in vector grid points.
     @param nGridPoints the number of grid points in batch.
     */
    void _compGxcBatchForGGA(CDenseMatrix&              molecularGradient,
                             const CAODensityMatrix&    densityMatrix,
                             const CMolecule&           molecule,
                             const CMolecularBasis&     basis,
                             const CMolecularGrid&      molecularGrid,
                             const CDensityGrid&        gsDensityGrid,
                             const CDensityGridQuad&    rwDensityGridQuad,
                             const CXCGradientGrid&     xcGradientGrid,
                             const CXCHessianGrid&      xcHessianGrid,
                             const CXCCubicHessianGrid& xcCubicHessianGrid,
                             const int32_t              gridOffset,
                             const int32_t              nGridPoints) const;

    /**
     Integrates batch of grid points for LDA TDDFT exchnage-correlation
     functional contribution to molecular gradient.

     @param molecularGradient the molecular gradient object.
     @param gsdensityMatrix the ground state AO density matrix.
     @param rwDensityMatrixOne the perturbed AO density matrix (relaxed_density_ao).
     @param rwDensityMatrixTwo the perturbed AO density matrix (x_minus_y_ao).
     @param molecule the molecule.
     @param basis the molecular basis.
     @param molecularGrid the molecular grid.
     @param gsDenistyGrid the ground state density grid.
     @param rwDenistyGridOne the perturbed density grid from rwDensityMatrixOne.
     @param rwDenistyGridTwo the perturbed density grid from rwDensityMatrixTwo.
     @param rwDenistyGridQuad the perturbed density grid for quadratic
            response (QRF).
     @param xcGradientGrid the exchange-correlation gradient grid.
     @param xcHessianGrid the exchange-correlation hessian grid.
     @param xcCubicHessianGrid the exchange-correlation cubic hessian grid.
     @param gridOffset the batch offset in vector grid points.
     @param nGridPoints the number of grid points in batch.
     */
    void _compTddftBatchForLDA(CDenseMatrix&              molecularGradient,
                               const CAODensityMatrix&    gsDensityMatrix,
                               const CAODensityMatrix&    rwDensityMatrixOne,
                               const CAODensityMatrix&    rwDensityMatrixTwo,
                               const CMolecule&           molecule,
                               const CMolecularBasis&     basis,
                               const CMolecularGrid&      molecularGrid,
                               const CDensityGrid&        gsDensityGrid,
                               const CDensityGrid&        rwDensityGridOne,
                               const CDensityGrid&        rwDensityGridTwo,
                               const CDensityGridQuad&    rwDensityGridQuad,
                               const CXCGradientGrid&     xcGradientGrid,
                               const CXCHessianGrid&      xcHessianGrid,
                               const CXCCubicHessianGrid& xcCubicHessianGrid,
                               const int32_t              gridOffset,
                               const int32_t              nGridPoints) const;

    /**
     Integrates batch of grid points for GGA TDDFT exchnage-correlation
     functional contribution to molecular gradient.

     @param molecularGradient the molecular gradient object.
     @param gsdensityMatrix the ground state AO density matrix.
     @param rwDensityMatrixOne the perturbed AO density matrix (relaxed_density_ao).
     @param rwDensityMatrixTwo the perturbed AO density matrix (x_minus_y_ao).
     @param molecule the molecule.
     @param basis the molecular basis.
     @param molecularGrid the molecular grid.
     @param gsDenistyGrid the ground state density grid.
     @param rwDenistyGridOne the perturbed density grid from rwDensityMatrixOne.
     @param rwDenistyGridTwo the perturbed density grid from rwDensityMatrixTwo.
     @param rwDenistyGridQuad the perturbed density grid for quadratic
            response (QRF).
     @param xcGradientGrid the exchange-correlation gradient grid.
     @param xcHessianGrid the exchange-correlation hessian grid.
     @param xcCubicHessianGrid the exchange-correlation cubic hessian grid.
     @param gridOffset the batch offset in vector grid points.
     @param nGridPoints the number of grid points in batch.
     */
    void _compTddftBatchForGGA(CDenseMatrix&              molecularGradient,
                               const CAODensityMatrix&    gsDensityMatrix,
                               const CAODensityMatrix&    rwDensityMatrixOne,
                               const CAODensityMatrix&    rwDensityMatrixTwo,
                               const CMolecule&           molecule,
                               const CMolecularBasis&     basis,
                               const CMolecularGrid&      molecularGrid,
                               const CDensityGrid&        gsDensityGrid,
                               const CDensityGrid&        rwDensityGridOne,
                               const CDensityGrid&        rwDensityGridTwo,
                               const CDensityGridQuad&    rwDensityGridQuad,
                               const CXCGradientGrid&     xcGradientGrid,
                               const CXCHessianGrid&      xcHessianGrid,
                               const CXCCubicHessianGrid& xcCubicHessianGrid,
                               const int32_t              gridOffset,
                               const int32_t              nGridPoints) const;

    /**
     Distributes density values into density grid for LDA.

     @param densityGrid the pointer to density grid object.
     @param densityMatrix the AO density matrix.
     @param aoIdentifiers the identifiers of AOs.
     @param braGtoValues the GTOs values buffer on bra side.
     @param ketGtoValuesX the GTOs gradient X values buffer on ket side.
     @param ketGtoValuesY the GTOs gradient Y values buffer on ket side.
     @param ketGtoValuesZ the GTOs gradient Z values buffer on ket side.
     @param nGridPoints the number of grid points in grid points batch.
     */
    void _distGradientDensityValuesForLDA(CDensityGrid&              densityGrid,
                                          const CAODensityMatrix&    densityMatrix,
                                          const CMemBlock<int32_t>&  aoIdentifiers,
                                          const CMemBlock2D<double>& braGtoValues,
                                          const CMemBlock2D<double>& ketGtoValuesX,
                                          const CMemBlock2D<double>& ketGtoValuesY,
                                          const CMemBlock2D<double>& ketGtoValuesZ,
                                          const int32_t              nGridPoints) const;

    /**
     Distributes density values into density grid for GGA.

     @param densityGrid the pointer to density grid object.
     @param densityMatrix the AO density matrix.
     @param aoIdentifiers the identifiers of AOs.
     @param braGtoValues the GTOs values buffer on bra side.
     @param braGtoValuesX the GTOs gradient X values buffer on bra side.
     @param braGtoValuesY the GTOs gradient Y values buffer on bra side.
     @param braGtoValuesZ the GTOs gradient Z values buffer on bra side.
     @param ketGtoValuesX the GTOs gradient X values buffer on ket side.
     @param ketGtoValuesY the GTOs gradient Y values buffer on ket side.
     @param ketGtoValuesZ the GTOs gradient Z values buffer on ket side.
     @param ketGtoValuesXX the GTOs hessian XX values buffer on ket side.
     @param ketGtoValuesXY the GTOs hessian XY values buffer on ket side.
     @param ketGtoValuesXZ the GTOs hessian XZ values buffer on ket side.
     @param ketGtoValuesYY the GTOs hessian YY values buffer on ket side.
     @param ketGtoValuesYZ the GTOs hessian YZ values buffer on ket side.
     @param ketGtoValuesZZ the GTOs hessian ZZ values buffer on ket side.
     @param nGridPoints the number of grid points in grid points batch.
     */
    void _distGradientDensityValuesForGGA(CDensityGrid&              densityGrid,
                                          const CAODensityMatrix&    densityMatrix,
                                          const CMemBlock<int32_t>&  aoIdentifiers,
                                          const CMemBlock2D<double>& braGtoValues,
                                          const CMemBlock2D<double>& braGtoValuesX,
                                          const CMemBlock2D<double>& braGtoValuesY,
                                          const CMemBlock2D<double>& braGtoValuesZ,
                                          const CMemBlock2D<double>& ketGtoValuesX,
                                          const CMemBlock2D<double>& ketGtoValuesY,
                                          const CMemBlock2D<double>& ketGtoValuesZ,
                                          const CMemBlock2D<double>& ketGtoValuesXX,
                                          const CMemBlock2D<double>& ketGtoValuesXY,
                                          const CMemBlock2D<double>& ketGtoValuesXZ,
                                          const CMemBlock2D<double>& ketGtoValuesYY,
                                          const CMemBlock2D<double>& ketGtoValuesYZ,
                                          const CMemBlock2D<double>& ketGtoValuesZZ,
                                          const int32_t              nGridPoints) const;

    /**
     Accumulates LDA first-order exchnage-correlation functional contribution
     to molecular gradient of an atom.

     @param molecularGradient the molecular gradient object.
     @param iAtom the index of the atom.
     @param gradientDensityGrid the nuclear gradient grid.
     @param molecularGrid the molecular grid.
     @param gsDenistyGrid the ground state density grid.
     @param xcGradientGrid the exchange-correlation gradient grid.
     @param gridOffset the batch offset in vector grid points.
     @param gridBlockPosition the position of grid block.
     @param nGridPoints the number of grid points in batch.
     */
    void _accumulateVxcContribForLDA(CDenseMatrix&          molecularGradient,
                                     const int32_t          iAtom,
                                     const CDensityGrid&    gradientDensityGrid,
                                     const CMolecularGrid&  molecularGrid,
                                     const CDensityGrid&    gsDensityGrid,
                                     const CXCGradientGrid& xcGradientGrid,
                                     const int32_t          gridOffset,
                                     const int32_t          gridBlockPosition,
                                     const int32_t          nGridPoints) const;

    /**
     Accumulates GGA first-order exchnage-correlation functional contribution
     to molecular gradient of an atom.

     @param molecularGradient the molecular gradient object.
     @param iAtom the index of the atom.
     @param gradientDensityGrid the nuclear gradient grid.
     @param molecularGrid the molecular grid.
     @param gsDenistyGrid the ground state density grid.
     @param xcGradientGrid the exchange-correlation gradient grid.
     @param gridOffset the batch offset in vector grid points.
     @param gridBlockPosition the position of grid block.
     @param nGridPoints the number of grid points in batch.
     */
    void _accumulateVxcContribForGGA(CDenseMatrix&          molecularGradient,
                                     const int32_t          iAtom,
                                     const CDensityGrid&    gradientDensityGrid,
                                     const CMolecularGrid&  molecularGrid,
                                     const CDensityGrid&    gsDensityGrid,
                                     const CXCGradientGrid& xcGradientGrid,
                                     const int32_t          gridOffset,
                                     const int32_t          gridBlockPosition,
                                     const int32_t          nGridPoints) const;

    /**
     Accumulates LDA second-order exchnage-correlation functional contribution
     to molecular gradient of an atom.

     @param molecularGradient the molecular gradient object.
     @param iAtom the index of the atom.
     @param gradientDensityGrid the nuclear gradient grid.
     @param molecularGrid the molecular grid.
     @param gsDenistyGrid the ground state density grid.
     @param rwDenistyGrid the perturbed density grid.
     @param xcGradientGrid the exchange-correlation gradient grid.
     @param xcHessianGrid the exchange-correlation hessian grid.
     @param gridOffset the batch offset in vector grid points.
     @param gridBlockPosition the position of grid block.
     @param nGridPoints the number of grid points in batch.
     */
    void _accumulateFxcContribForLDA(CDenseMatrix&          molecularGradient,
                                     const int32_t          iAtom,
                                     const CDensityGrid&    gradientDensityGrid,
                                     const CMolecularGrid&  molecularGrid,
                                     const CDensityGrid&    gsDensityGrid,
                                     const CDensityGrid&    rwDensityGrid,
                                     const CXCGradientGrid& xcGradientGrid,
                                     const CXCHessianGrid&  xcHessianGrid,
                                     const int32_t          gridOffset,
                                     const int32_t          gridBlockPosition,
                                     const int32_t          nGridPoints) const;

    /**
     Accumulates GGA second-order exchnage-correlation functional contribution
     to molecular gradient of an atom.

     @param molecularGradient the molecular gradient object.
     @param iAtom the index of the atom.
     @param gradientDensityGrid the nuclear gradient grid.
     @param molecularGrid the molecular grid.
     @param gsDenistyGrid the ground state density grid.
     @param rwDenistyGrid the perturbed density grid.
     @param xcGradientGrid the exchange-correlation gradient grid.
     @param xcHessianGrid the exchange-correlation hessian grid.
     @param gridOffset the batch offset in vector grid points.
     @param gridBlockPosition the position of grid block.
     @param nGridPoints the number of grid points in batch.
     */
    void _accumulateFxcContribForGGA(CDenseMatrix&          molecularGradient,
                                     const int32_t          iAtom,
                                     const CDensityGrid&    gradientDensityGrid,
                                     const CMolecularGrid&  molecularGrid,
                                     const CDensityGrid&    gsDensityGrid,
                                     const CDensityGrid&    rwDensityGrid,
                                     const CXCGradientGrid& xcGradientGrid,
                                     const CXCHessianGrid&  xcHessianGrid,
                                     const int32_t          gridOffset,
                                     const int32_t          gridBlockPosition,
                                     const int32_t          nGridPoints) const;

    /**
     Accumulates LDA third-order exchnage-correlation functional contribution
     to molecular gradient of an atom.

     @param molecularGradient the molecular gradient object.
     @param iAtom the index of the atom.
     @param gradientDensityGrid the nuclear gradient grid.
     @param molecularGrid the molecular grid.
     @param gsDenistyGrid the ground state density grid.
     @param rwDenistyGridQuad the perturbed density grid for quadratic
            response.
     @param xcGradientGrid the exchange-correlation gradient grid.
     @param xcHessianGrid the exchange-correlation hessian grid.
     @param xcCubicHessianGrid the exchange-correlation cubic hessian grid.
     @param gridOffset the batch offset in vector grid points.
     @param gridBlockPosition the position of grid block.
     @param nGridPoints the number of grid points in batch.
     */
    void _accumulateGxcContribForLDA(CDenseMatrix&              molecularGradient,
                                     const int32_t              iAtom,
                                     const CDensityGrid&        gradientDensityGrid,
                                     const CMolecularGrid&      molecularGrid,
                                     const CDensityGrid&        gsDensityGrid,
                                     const CDensityGridQuad&    rwDensityGridQuad,
                                     const CXCGradientGrid&     xcGradientGrid,
                                     const CXCHessianGrid&      xcHessianGrid,
                                     const CXCCubicHessianGrid& xcCubicHessianGrid,
                                     const int32_t              gridOffset,
                                     const int32_t              gridBlockPosition,
                                     const int32_t              nGridPoints) const;

    /**
     Accumulates GGA third-order exchnage-correlation functional contribution
     to molecular gradient of an atom.

     @param molecularGradient the molecular gradient object.
     @param iAtom the index of the atom.
     @param gradientDensityGrid the nuclear gradient grid.
     @param molecularGrid the molecular grid.
     @param gsDenistyGrid the ground state density grid.
     @param rwDenistyGridQuad the perturbed density grid for quadratic
            response.
     @param xcGradientGrid the exchange-correlation gradient grid.
     @param xcHessianGrid the exchange-correlation hessian grid.
     @param xcCubicHessianGrid the exchange-correlation cubic hessian grid.
     @param gridOffset the batch offset in vector grid points.
     @param gridBlockPosition the position of grid block.
     @param nGridPoints the number of grid points in batch.
     */
    void _accumulateGxcContribForGGA(CDenseMatrix&              molecularGradient,
                                     const int32_t              iAtom,
                                     const CDensityGrid&        gradientDensityGrid,
                                     const CMolecularGrid&      molecularGrid,
                                     const CDensityGrid&        gsDensityGrid,
                                     const CDensityGridQuad&    rwDensityGridQuad,
                                     const CXCGradientGrid&     xcGradientGrid,
                                     const CXCHessianGrid&      xcHessianGrid,
                                     const CXCCubicHessianGrid& xcCubicHessianGrid,
                                     const int32_t              gridOffset,
                                     const int32_t              gridBlockPosition,
                                     const int32_t              nGridPoints) const;

    /**
     Gets size of block in grid batch.

     @return the size of block in grid batch.
     */
    int32_t _getSizeOfBlock() const;

   public:
    /**
     Creates a XC molecular gradient integrator object using MPI info.

     @param comm the MPI communicator.
     */
    CXCMolecularGradient(MPI_Comm comm);

    /**
     Destroys a XC molecular gradient integrator object.
     */
    ~CXCMolecularGradient();

    /**
     Integrates 1st-order exchnage-correlation functional contribution to
     molecular gradient.

     @param gsDensityMatrix the ground state AO density matrix object.
     @param molecule the molecule.
     @param basis the molecular basis.
     @param molecularGrid the molecular grid.
     @param xcFuncLabel the label of exchange-correlation functional.
     @return the molecular gradient.
     */
    CDenseMatrix integrateVxcGradient(const CAODensityMatrix& gsDensityMatrix,
                                      const CMolecule&        molecule,
                                      const CMolecularBasis&  basis,
                                      const CMolecularGrid&   molecularGrid,
                                      const std::string&      xcFuncLabel) const;

    /**
     Integrates 1st-order exchnage-correlation functional contribution to
     molecular gradient.

     @param rwDensityMatrix the perturbed AO density matrix object (to be
            contracted with GTO gradient).
     @param gsDensityMatrix the ground state AO density matrix object.
     @param molecule the molecule.
     @param basis the molecular basis.
     @param molecularGrid the molecular grid.
     @param xcFuncLabel the label of exchange-correlation functional.
     @return the molecular gradient.
     */
    CDenseMatrix integrateVxcGradient(const CAODensityMatrix& rwDensityMatrix,
                                      const CAODensityMatrix& gsDensityMatrix,
                                      const CMolecule&        molecule,
                                      const CMolecularBasis&  basis,
                                      const CMolecularGrid&   molecularGrid,
                                      const std::string&      xcFuncLabel) const;

    /**
     Integrates 2nd-order exchnage-correlation functional contribution to
     molecular gradient.

     @param rwDensityMatrixOne the perturbed AO density matrix object.
     @param rwDensityMatrixTwo the perturbed AO density matrix object (to be
            contracted with GTO gradient).
     @param gsDensityMatrix the ground state AO density matrix object.
     @param molecule the molecule.
     @param basis the molecular basis.
     @param molecularGrid the molecular grid.
     @param xcFuncLabel the label of exchange-correlation functional.
     @return the molecular gradient.
     */
    CDenseMatrix integrateFxcGradient(const CAODensityMatrix& rwDensityMatrixOne,
                                      const CAODensityMatrix& rwDensityMatrixTwo,
                                      const CAODensityMatrix& gsDensityMatrix,
                                      const CMolecule&        molecule,
                                      const CMolecularBasis&  basis,
                                      const CMolecularGrid&   molecularGrid,
                                      const std::string&      xcFuncLabel) const;

    /**
     Integrates 3rd-order exchnage-correlation functional contribution to
     molecular gradient.

     @param rwDensityMatrixOne the perturbed AO density matrix object.
     @param rwDensityMatrixTwo the perturbed AO density matrix object.
     @param gsDensityMatrix the ground state AO density matrix object (to be
            contracted with GTO gradient).
     @param molecule the molecule.
     @param basis the molecular basis.
     @param molecularGrid the molecular grid.
     @param xcFuncLabel the label of exchange-correlation functional.
     @return the molecular gradient.
     */
    CDenseMatrix integrateGxcGradient(const CAODensityMatrix& rwDensityMatrixOne,
                                      const CAODensityMatrix& rwDensityMatrixTwo,
                                      const CAODensityMatrix& gsDensityMatrix,
                                      const CMolecule&        molecule,
                                      const CMolecularBasis&  basis,
                                      const CMolecularGrid&   molecularGrid,
                                      const std::string&      xcFuncLabel) const;

    /**
     Integrates TDDFT exchnage-correlation functional contribution to molecular
     gradient.

     @param rwDensityMatrixOne the perturbed AO density matrix (relaxed_density_ao).
     @param rwDensityMatrixTwo the perturbed AO density matrix (x_minus_y_ao).
     @param gsDensityMatrix the ground state AO density matrix.
     @param molecule the molecule.
     @param basis the molecular basis.
     @param molecularGrid the molecular grid.
     @param xcFuncLabel the label of exchange-correlation functional.
     @return the molecular gradient.
     */
    CDenseMatrix integrateTddftGradient(const CAODensityMatrix& rwDensityMatrixOne,
                                        const CAODensityMatrix& rwDensityMatrixTwo,
                                        const CAODensityMatrix& gsDensityMatrix,
                                        const CMolecule&        molecule,
                                        const CMolecularBasis&  basis,
                                        const CMolecularGrid&   molecularGrid,
                                        const std::string&      xcFuncLabel) const;
};

#endif /* XCMolecularGradient_hpp */

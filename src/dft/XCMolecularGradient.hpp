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

    void _compVxcContribForLDA(CDenseMatrix&           molecularGradient,
                               const CAODensityMatrix& densityMatrix,
                               const CMolecule&        molecule,
                               const CMolecularBasis&  basis,
                               const CMolecularGrid&   molecularGrid,
                               const CDensityGrid&     gsDensityGrid,
                               const CXCGradientGrid&  xcGradientGrid) const;

    void _compVxcBatchForLDA(CDenseMatrix&           molecularGradient,
                             const CAODensityMatrix& densityMatrix,
                             const CMolecule&        molecule,
                             const CMolecularBasis&  basis,
                             const CMolecularGrid&   molecularGrid,
                             const CDensityGrid&     gsDensityGrid,
                             const CXCGradientGrid&  xcGradientGrid,
                             const int32_t           gridOffset,
                             const int32_t           nGridPoints) const;

    void _distGradientDensityValuesForLda(CDensityGrid&              densityGrid,
                                          const CAODensityMatrix&    densityMatrix,
                                          const CMemBlock<int32_t>&  aoIdentifiers,
                                          const CMemBlock2D<double>& braGtoValues,
                                          const CMemBlock2D<double>& ketGtoValuesX,
                                          const CMemBlock2D<double>& ketGtoValuesY,
                                          const CMemBlock2D<double>& ketGtoValuesZ,
                                          const int32_t              nGridPoints) const;

    void _accumulateVxcContribForLDA(CDenseMatrix&          molecularGradient,
                                     const int32_t          iAtom,
                                     const CDensityGrid&    gradientDensityGrid,
                                     const CMolecularGrid&  molecularGrid,
                                     const CDensityGrid&    gsDensityGrid,
                                     const CXCGradientGrid& xcGradientGrid,
                                     const int32_t          gridOffset,
                                     const int32_t          gridBlockPosition,
                                     const int32_t          nGridPoints) const;

    int32_t _getSizeOfBlock() const;

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
     Integrates exchnage-correlation functional contribution to molecular gradient.

     @param aoDensityMatrix the AO density matrix object.
     @param molecule the molecule.
     @param basis the molecular basis.
     @param molecularGrid the molecular grid.
     @param xcFuncLabel the label of exchange-correlation functional.
     @return the molecular gradient.
     */
    CDenseMatrix integrateVxcGradient(const CAODensityMatrix& aoDensityMatrix,
                                      const CMolecule&        molecule,
                                      const CMolecularBasis&  basis,
                                      const CMolecularGrid&   molecularGrid,
                                      const std::string&      xcFuncLabel) const;

    /**
     Integrates exchnage-correlation functional contribution to molecular gradient.

     @param rwDensityMatrix the perturbed AO density matrix object.
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
     Integrates exchnage-correlation functional contribution to molecular gradient.

     @param rwDensityMatrixOne the perturbed AO density matrix object.
     @param rwDensityMatrixTwo the perturbed AO density matrix object.
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
     Integrates exchnage-correlation functional contribution to molecular gradient.

     @param rwDensityMatrixOne the perturbed AO density matrix object.
     @param rwDensityMatrixTwo the perturbed AO density matrix object.
     @param gsDensityMatrix the ground state AO density matrix object.
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
     Integrates exchnage-correlation functional contribution to molecular gradient.

     @param rwDensityMatrixOne the perturbed AO density matrix object (relaxed_dens).
     @param rwDensityMatrixTwo the perturbed AO density matrix object (x_minus_y_ao).
     @param gsDensityMatrix the ground state AO density matrix object.
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

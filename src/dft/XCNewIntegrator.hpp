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

#ifndef XCNewIntegrator_hpp
#define XCNewIntegrator_hpp

#include <list>
#include <array>
#include <string>

#include "AODensityMatrix.hpp"
#include "AOFockMatrix.hpp"
#include "AOKohnShamMatrix.hpp"
#include "DenseMatrix.hpp"
#include "DensityGridQuad.hpp"
#include "GridBox.hpp"
#include "GtoContainer.hpp"
#include "MemBlock2D.hpp"
#include "MolecularGrid.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "MultiTimer.hpp"
#include "XCFunctional.hpp"
#include "XCGradientGrid.hpp"
#include "XCHessianGrid.hpp"
#include "XCCubicHessianGrid.hpp"

/**
 Class CXCNewIntegrator implements XC integrator.

 @author X. Li
 */
class CXCNewIntegrator
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
     Screening threshold for GTO values on grid points.
     */
    double _screeningThresholdForGTOValues;

    /**
     Screening threshold for density values on grid points.
     */
    double _screeningThresholdForDensityValues;

    /**
     Integrates first-order LDA exchnage-correlation functional contribution to
     AO Kohn-Sham matrix.

     @param molecule the molecule.
     @param basis the molecular basis.
     @param densityMatrix the AO density matrix.
     @param molecularGrid the molecular grid.
     @param xcFunctional the exchange-correlation functional.
     @return the AO Kohn-Sham matrix.
     */
    CAOKohnShamMatrix _integrateVxcFockForLDA(const CMolecule&        molecule,
                                              const CMolecularBasis&  basis,
                                              const CAODensityMatrix& densityMatrix,
                                              const CMolecularGrid&   molecularGrid,
                                              const CXCFunctional&    xcFunctional) const;

    /**
     Integrates first-order GGA exchnage-correlation functional contribution to
     AO Kohn-Sham matrix.

     @param molecule the molecule.
     @param basis the molecular basis.
     @param densityMatrix the AO density matrix.
     @param molecularGrid the molecular grid.
     @param xcFunctional the exchange-correlation functional.
     @return the AO Kohn-Sham matrix.
     */
    CAOKohnShamMatrix _integrateVxcFockForGGA(const CMolecule&        molecule,
                                              const CMolecularBasis&  basis,
                                              const CAODensityMatrix& densityMatrix,
                                              const CMolecularGrid&   molecularGrid,
                                              const CXCFunctional&    xcFunctional) const;

    /**
     Integrates second-order LDA exchnage-correlation functional contribution
     to AO Fock matrix.

     @param aoFockMatrix the AO Fock matrix.
     @param molecule the molecule.
     @param basis the molecular basis.
     @param rwDensityMatrix the perturbed density matrix.
     @param gsDensityMatrix the ground state density matrix.
     @param molecularGrid the molecular grid.
     @param xcFunctional the exchange-correlation functional.
     */
    void _integrateFxcFockForLDA(CAOFockMatrix&          aoFockMatrix,
                                 const CMolecule&        molecule,
                                 const CMolecularBasis&  basis,
                                 const CAODensityMatrix& rwDensityMatrix,
                                 const CAODensityMatrix& gsDensityMatrix,
                                 const CMolecularGrid&   molecularGrid,
                                 const CXCFunctional&    xcFunctional) const;

    /**
     Integrates second-order GGA exchnage-correlation functional contribution
     to AO Fock matrix.

     @param aoFockMatrix the AO Fock matrix.
     @param molecule the molecule.
     @param basis the molecular basis.
     @param rwDensityMatrix the perturbed density matrix.
     @param gsDensityMatrix the ground state density matrix.
     @param molecularGrid the molecular grid.
     @param xcFunctional the exchange-correlation functional.
     */
    void _integrateFxcFockForGGA(CAOFockMatrix&          aoFockMatrix,
                                 const CMolecule&        molecule,
                                 const CMolecularBasis&  basis,
                                 const CAODensityMatrix& rwDensityMatrix,
                                 const CAODensityMatrix& gsDensityMatrix,
                                 const CMolecularGrid&   molecularGrid,
                                 const CXCFunctional&    xcFunctional) const;

    /**
     Integrates third-order LDA exchnage-correlation functional contribution
     to AO Fock matrix.

     @param aoFockMatrix the AO Fock matrix.
     @param molecule the molecule.
     @param basis the molecular basis.
     @param rwDensityMatrix the perturbed one-time transformed densities.
     @param rw2DensityMatrix the two-time transformed densities.
     @param gsDensityMatrix the ground state density matrix.
     @param molecularGrid the molecular grid.
     @param xcFunctional the exchange-correlation functional.
     @param quadMode a string that specifies which densities should be combined.
     */
    void _integrateKxcFockForLDA(CAOFockMatrix&          aoFockMatrix,
                                 const CMolecule&        molecule,
                                 const CMolecularBasis&  basis,
                                 const CAODensityMatrix& rwDensityMatrix,
                                 const CAODensityMatrix& rw2DensityMatrix,
                                 const CAODensityMatrix& gsDensityMatrix,
                                 const CMolecularGrid&   molecularGrid,
                                 const CXCFunctional&    xcFunctional,
                                 const std::string&      quadMode) const;

    /**
     Integrates third-order GGA exchnage-correlation functional contribution
     to AO Fock matrix.

     @param aoFockMatrix the AO Fock matrix.
     @param molecule the molecule.
     @param basis the molecular basis.
     @param rwDensityMatrix the perturbed one-time transformed densities.
     @param rw2DensityMatrix the two-time transformed densities.
     @param gsDensityMatrix the ground state density matrix.
     @param molecularGrid the molecular grid.
     @param xcFunctional the exchange-correlation functional.
     @param quadMode a string that specifies which densities should be combined.
     */
    void _integrateKxcFockForGGA(CAOFockMatrix&          aoFockMatrix,
                                 const CMolecule&        molecule,
                                 const CMolecularBasis&  basis,
                                 const CAODensityMatrix& rwDensityMatrix,
                                 const CAODensityMatrix& rw2DensityMatrix,
                                 const CAODensityMatrix& gsDensityMatrix,
                                 const CMolecularGrid&   molecularGrid,
                                 const CXCFunctional&    xcFunctional,
                                 const std::string&      quadMode) const;

    /**
     Integrates LDA contribution to (first-order) Vxc matrix.

     @param npoints the number of grid points.
     @param weights the weights of grid points.
     @param gtoValues the GTO values on grid points.
     @param xcGradientGrid the exchange-correlation gradient grid.
     @param timer the timer.
     @return the contribution as a CDenseMatrix object.
     */
    CDenseMatrix _integratePartialVxcFockForLDA(const int32_t          npoints,
                                                const double*          weights,
                                                const CDenseMatrix&    gtoValues,
                                                const CXCGradientGrid& xcGradientGrid,
                                                CMultiTimer&           timer) const;

    /**
     Integrates GGA contribution to AO Kohn-Sham matrix.

     @param npoints the number of grid points.
     @param weights the weights of grid points.
     @param gtoValues the GTO values on grid points.
     @param gtoValuesX the GTO gradient X values on grid points.
     @param gtoValuesY the GTO gradient Y values on grid points.
     @param gtoValuesZ the GTO gradient Z values on grid points.
     @param xcGradientGrid the exchange-correlation gradient grid.
     @param densityGrid the density grid.
     @param timer the timer.
     @return the contribution as a CDenseMatrix object.
     */
    CDenseMatrix _integratePartialVxcFockForGGA(const int32_t          npoints,
                                                const double*          weights,
                                                const CDenseMatrix&    gtoValues,
                                                const CDenseMatrix&    gtoValuesX,
                                                const CDenseMatrix&    gtoValuesY,
                                                const CDenseMatrix&    gtoValuesZ,
                                                const CXCGradientGrid& xcGradientGrid,
                                                const CDensityGrid&    densityGrid,
                                                CMultiTimer&           timer) const;

    /**
     Integrates LDA contribution to (second-order) Fxc matrix.

     @param npoints the number of grid points.
     @param weights the weights of grid points.
     @param gtoValues the GTO values on grid points.
     @param xcHessianGrid the exchange-correlation hessian grid.
     @param rwDensityGrid the perturbed density grid.
     @param timer the timer.
     @return the contribution as a CDenseMatrix object.
     */
    CDenseMatrix _integratePartialFxcFockForLDA(const int32_t          npoints,
                                                const double*          weights,
                                                const CDenseMatrix&    gtoValues,
                                                const CXCHessianGrid&  xcHessianGrid,
                                                const CDensityGrid&    rwDensityGrid,
                                                CMultiTimer&           timer) const;

    /**
     Integrates GGA contribution to (second-order) Fxc matrix.

     @param npoints the number of grid points.
     @param weights the weights of grid points.
     @param gtoValues the GTO values on grid points.
     @param gtoValuesX the GTO gradient X values on grid points.
     @param gtoValuesY the GTO gradient Y values on grid points.
     @param gtoValuesZ the GTO gradient Z values on grid points.
     @param xcGradientGrid the exchange-correlation gradient grid.
     @param xcHessianGrid the exchange-correlation hessian grid.
     @param rwDensityGrid the perturbed density grid.
     @param gsDensityGrid the ground-state density grid.
     @param timer the timer.
     @return the contribution as a CDenseMatrix object.
     */
    CDenseMatrix _integratePartialFxcFockForGGA(const int32_t          npoints,
                                                const double*          weights,
                                                const CDenseMatrix&    gtoValues,
                                                const CDenseMatrix&    gtoValuesX,
                                                const CDenseMatrix&    gtoValuesY,
                                                const CDenseMatrix&    gtoValuesZ,
                                                const CXCGradientGrid& xcGradientGrid,
                                                const CXCHessianGrid&  xcHessianGrid,
                                                const CDensityGrid&    rwDensityGrid,
                                                const CDensityGrid&    gsDensityGrid,
                                                CMultiTimer&           timer) const;

    /**
     Integrates LDA contribution to (third-order) Kxc matrix.

     @param gridBlockPosition the starting position of the grid box.
     @param npoints the number of grid points.
     @param weights the weights of grid points.
     @param gtoValues the GTO values on grid points.
     @param xcHessianGrid the exchange-correlation hessian grid.
     @param xcCubicHessianGrid the exchange-correlation cubic hessian grid.
     @param rwDensityGridQuad the products of one-time transformed densities on grid points.
     @param rw2DensityMatrix the two-time transformed densities on grid points.
     @param iFock the index of the AO Fock matrix.
     @param timer the timer.
     @return the contribution as a CDenseMatrix object.
     */
    CDenseMatrix _integratePartialKxcFockForLDA(const int32_t              gridBlockPosition,
                                                const int32_t              npoints,
                                                const double*              weights,
                                                const CDenseMatrix&        gtoValues,
                                                const CXCHessianGrid&      xcHessianGrid,
                                                const CXCCubicHessianGrid& xcCubicHessianGrid,
                                                const CDensityGridQuad&    rwDensityGridQuad,
                                                const CDensityGrid&        rw2DensityGrid,
                                                const int32_t              iFock,
                                                CMultiTimer&               timer) const;

    /**
     Integrates GGA contribution to (third-order) Kxc matrix.

     @param gridblockpos the starting position of the grid box.
     @param npoints the number of grid points.
     @param weights the weights of grid points.
     @param gtoValues the GTO values on grid points.
     @param gtoValuesX the GTO gradient X values on grid points.
     @param gtoValuesY the GTO gradient Y values on grid points.
     @param gtoValuesZ the GTO gradient Z values on grid points.
     @param xcGradientGrid the exchange-correlation gradient grid.
     @param xcHessianGrid the exchange-correlation hessian grid.
     @param xcCubicHessianGrid the exchange-correlation cubic hessian grid.
     @param rwDensityGridQuad the products of one-time transformed densities on grid points.
     @param rw2DensityMatrix the two-time transformed densities on grid points.
     @param gsDensityGrid the ground-state density grid.
     @param iFock the index of the AO Fock matrix.
     @param timer the timer.
     @return the contribution as a CDenseMatrix object.
     */
    CDenseMatrix _integratePartialKxcFockForGGA(const int32_t              gridblockpos,
                                                const int32_t              npoints,
                                                const double*              weights,
                                                const CDenseMatrix&        gtoValues,
                                                const CDenseMatrix&        gtoValuesX,
                                                const CDenseMatrix&        gtoValuesY,
                                                const CDenseMatrix&        gtoValuesZ,
                                                const CXCGradientGrid&     xcGradientGrid,
                                                const CXCHessianGrid&      xcHessianGrid,
                                                const CXCCubicHessianGrid& xcCubicHessianGrid,
                                                const CDensityGridQuad&    rwDensityGridQuad,
                                                const CDensityGrid&        rw2DensityGrid,
                                                const CDensityGrid&        gsDensityGrid,
                                                const int32_t              iFock,
                                                CMultiTimer&               timer) const;

   public:
    /**
     Creates an XC integrator object using MPI info.

     @param comm the MPI communicator.
     */
    CXCNewIntegrator(MPI_Comm comm);

    /**
     Destroys an XC integrator object.
     */
    ~CXCNewIntegrator();

    /**
     Integrates first-order exchnage-correlation functional contribution to AO
     Kohn-Sham matrix.

     @param molecule the molecule.
     @param basis the molecular basis.
     @param densityMatrix the AO density matrix object.
     @param molecularGrid the molecular grid.
     @param xcFuncLabel the label of exchange-correlation functional.
     @return the AO Kohn-Sham matrix.
     */
    CAOKohnShamMatrix integrateVxcFock(const CMolecule&        molecule,
                                       const CMolecularBasis&  basis,
                                       const CAODensityMatrix& densityMatrix,
                                       const CMolecularGrid&   molecularGrid,
                                       const std::string&      xcFuncLabel) const;

    /**
     Integrates second-order exchnage-correlation functional contribution to AO
     Fock matrix.

     @param aoFockMatrix the AO Fock matrix.
     @param molecule the molecule.
     @param basis the molecular basis.
     @param rwDensityMatrix the perturbed density matrix.
     @param gsDensityMatrix the ground state density matrix.
     @param molecularGrid the molecular grid.
     @param xcFuncLabel the label of exchange-correlation functional.
     */
    void integrateFxcFock(CAOFockMatrix&          aoFockMatrix,
                          const CMolecule&        molecule,
                          const CMolecularBasis&  basis,
                          const CAODensityMatrix& rwDensityMatrix,
                          const CAODensityMatrix& gsDensityMatrix,
                          const CMolecularGrid&   molecularGrid,
                          const std::string&      xcFuncLabel) const;

    /**
     Integrates third-order exchnage-correlation functional contribution to AO
     Fock matrix in quadratic response.

     @param aoFockMatrix the AO Fock matrix.
     @param molecule the molecule.
     @param basis the molecular basis.
     @param rwDensityMatrix the one-time transformed densities.
     @param rw2DensityMatrix the two-time transformed densities.
     @param gsDensityMatrix the ground state density matrix.
     @param molecularGrid the molecular grid.
     @param xcFuncLabel the label of exchange-correlation functional.
     @param quadMode a string that specifies which densities should be combined.
     */
    void integrateKxcFock(CAOFockMatrix&          aoFockMatrix,
                          const CMolecule&        molecule,
                          const CMolecularBasis&  basis,
                          const CAODensityMatrix& rwDensityMatrix,
                          const CAODensityMatrix& rw2DensityMatrix,
                          const CAODensityMatrix& gsDensityMatrix,
                          const CMolecularGrid&   molecularGrid,
                          const std::string&      xcFuncLabel,
                          const std::string&      quadMode) const;

    /**
     Computes GTOs values on grid points.

     @param molecule the molecule.
     @param basis the molecular basis.
     @param molecularGrid the molecular grid.
     @return the GTO values on grid points.
     */
    CDenseMatrix computeGtoValuesOnGridPoints(const CMolecule&        molecule,
                                              const CMolecularBasis&  basis,
                                              const CMolecularGrid&   molecularGrid) const;
};

#endif /* XCNewIntegrator_hpp */

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

#include <array>
#include <string>

#include "AODensityMatrix.hpp"
#include "AOFockMatrix.hpp"
#include "AOKohnShamMatrix.hpp"
#include "Dense4DTensor.hpp"
#include "DenseMatrix.hpp"
#include "DensityGridQuad.hpp"
#include "GridBox.hpp"
#include "GtoContainer.hpp"
#include "MemBlock2D.hpp"
#include "MolecularBasis.hpp"
#include "MolecularGrid.hpp"
#include "Molecule.hpp"
#include "MultiTimer.hpp"
#include "XCCubicHessianGrid.hpp"
#include "XCFunctional.hpp"
#include "XCGradientGrid.hpp"
#include "XCHessianGrid.hpp"
#include "XCNewFunctional.hpp"
#include "XCPairDensityFunctional.hpp"

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
     Integrates first-order LDA exchange-correlation functional contribution to
     AO Kohn-Sham matrix.

     @param molecule the molecule.
     @param basis the molecular basis.
     @param densityMatrix the AO density matrix.
     @param molecularGrid the molecular grid.
     @param xcFunctional the exchange-correlation functional.
     @param flag the flag for closed/open shell.
     @return the AO Kohn-Sham matrix.
     */
    CAOKohnShamMatrix _integrateVxcFockForLDA(const CMolecule&        molecule,
                                              const CMolecularBasis&  basis,
                                              const CAODensityMatrix& densityMatrix,
                                              const CMolecularGrid&   molecularGrid,
                                              const CXCNewFunctional& xcFunctional,
                                              const std::string&      flag=std::string("closedshell")) const;

    /**
     Integrates first-order GGA exchange-correlation functional contribution to
     AO Kohn-Sham matrix.

     @param molecule the molecule.
     @param basis the molecular basis.
     @param densityMatrix the AO density matrix.
     @param molecularGrid the molecular grid.
     @param xcFunctional the exchange-correlation functional.
     @param flag the flag for closed/open shell.
     @return the AO Kohn-Sham matrix.
     */
    CAOKohnShamMatrix _integrateVxcFockForGGA(const CMolecule&        molecule,
                                              const CMolecularBasis&  basis,
                                              const CAODensityMatrix& densityMatrix,
                                              const CMolecularGrid&   molecularGrid,
                                              const CXCNewFunctional& xcFunctional,
                                              const std::string&      flag=std::string("closedshell")) const;

    /**
     Integrates second-order LDA exchange-correlation functional contribution
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
                                 const CXCNewFunctional& xcFunctional) const;

    /**
     Integrates second-order GGA exchange-correlation functional contribution
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
                                 const CXCNewFunctional& xcFunctional) const;

    /**
     Integrates third-order LDA exchange-correlation functional contribution
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
     Integrates third-order GGA exchange-correlation functional contribution
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
     Integrates first-order LDA pair-density functional contribution to
     AO Kohn-Sham matrix and MO "Q-matrix".

     @param aoFockMatrix the AO Fock matrix.
     @param moTwoBodyGradient the MO Two-body energy gradient term.
     @param molecule the molecule.
     @param basis the molecular basis.
     @param DensityMatrix the AO density matrix object.
     @param TwoBodyDensityMatrix the MO two-body active density matrix.
     @param ActiveMOs the active molecular orbitals.
     @param molecularGrid the molecular grid.
     @param fvxc the exchange-correlation functional.
     */
     void _integrateVxcPDFTForLDA(CAOKohnShamMatrix&              aoFockMatrix,
                                  CDense4DTensor&                 moTwoBodyGradient,
                                  const CMolecule&                molecule,
                                  const CMolecularBasis&          basis,
                                  const CAODensityMatrix&         DensityMatrix,
                                  const CDense4DTensor&           TwoBodyDensityMatrix,
                                  const CDenseMatrix&             ActiveMOs,
                                  const CMolecularGrid&           molecularGrid,
                                  const CXCPairDensityFunctional& fvxc) const;

    /**
     Integrates first-order GGA pair-density functional contribution to
     AO Kohn-Sham matrix and MO "Q-matrix".

     @param aoFockMatrix the AO Fock matrix.
     @param moTwoBodyGradient the MO Two-body energy gradient term.
     @param molecule the molecule.
     @param basis the molecular basis.
     @param DensityMatrix the AO density matrix object.
     @param TwoBodyDensityMatrix the MO two-body active density matrix.
     @param ActiveMOs the active molecular orbitals.
     @param molecularGrid the molecular grid.
     @param fvxc the exchange-correlation functional.
     */
     void _integrateVxcPDFTForGGA(CAOKohnShamMatrix&      aoFockMatrix,
                                  CDense4DTensor&         moTwoBodyGradient,
                                  const CMolecule&        molecule,
                                  const CMolecularBasis&  basis,
                                  const CAODensityMatrix& DensityMatrix,
                                  const CDense4DTensor&   TwoBodyDensityMatrix,
                                  const CDenseMatrix&     ActiveMOs,
                                  const CMolecularGrid&   molecularGrid,
                                  const CXCNewFunctional& fvxc) const;

    /**
     Integrates LDA contribution to (first-order) Vxc matrix.

     @param npoints the number of grid points.
     @param weights the weights of grid points.
     @param gtoValues the GTO values on grid points.
     @param vrho the 1st-order functional derivative wrt density.
     @param timer the timer.
     @return the contribution as a CDenseMatrix object.
     */
    CDenseMatrix _integratePartialVxcFockForLDA(const int32_t       npoints,
                                                const double*       weights,
                                                const CDenseMatrix& gtoValues,
                                                const double*       vrho,
                                                CMultiTimer&        timer) const;

    /**
     Integrates LDA contribution to (first-order) Vxc matrix.

     @param npoints the number of grid points.
     @param weights the weights of grid points.
     @param gtoValues the GTO values on grid points.
     @param vrho the 1st-order functional derivative wrt density.
     @param timer the timer.
     @return the alpha and beta contribution as a list of CDenseMatrix objects.
     */
    std::vector<CDenseMatrix> _integratePartialVxcFockForLDAOpenShell(const int32_t          npoints,
                                                                      const double*          weights,
                                                                      const CDenseMatrix&    gtoValues,
                                                                      const double*          vrho,
                                                                      CMultiTimer&           timer) const;

    /**
     Integrates GGA contribution to AO Kohn-Sham matrix.

     @param npoints the number of grid points.
     @param weights the weights of grid points.
     @param gtoValues the GTO values on grid points.
     @param gtoValuesX the GTO gradient X values on grid points.
     @param gtoValuesY the GTO gradient Y values on grid points.
     @param gtoValuesZ the GTO gradient Z values on grid points.
     @param rhograd the gradient density.
     @param vrho the 1st-order functional derivative wrt rho.
     @param vsigma the 1st-order functional derivative wrt sigma.
     @param timer the timer.
     @return the contribution as a CDenseMatrix object.
     */
    CDenseMatrix _integratePartialVxcFockForGGA(const int32_t          npoints,
                                                const double*          weights,
                                                const CDenseMatrix&    gtoValues,
                                                const CDenseMatrix&    gtoValuesX,
                                                const CDenseMatrix&    gtoValuesY,
                                                const CDenseMatrix&    gtoValuesZ,
                                                const double*          rhograd,
                                                const double*          vrho,
                                                const double*          vsigma,
                                                CMultiTimer&           timer) const;

    /**
     Integrates GGA contribution to AO Kohn-Sham matrix.

     @param npoints the number of grid points.
     @param weights the weights of grid points.
     @param gtoValues the GTO values on grid points.
     @param gtoValuesX the GTO gradient X values on grid points.
     @param gtoValuesY the GTO gradient Y values on grid points.
     @param gtoValuesZ the GTO gradient Z values on grid points.
     @param rhograd the gradient density.
     @param vrho the 1st-order functional derivative wrt rho.
     @param vsigma the 1st-order functional derivative wrt sigma.
     @param timer the timer.
     @return the alpha and beta contribution as a list of CDenseMatrix objects.
     */
    std::vector<CDenseMatrix> _integratePartialVxcFockForGGAOpenShell(const int32_t          npoints,
                                                                      const double*          weights,
                                                                      const CDenseMatrix&    gtoValues,
                                                                      const CDenseMatrix&    gtoValuesX,
                                                                      const CDenseMatrix&    gtoValuesY,
                                                                      const CDenseMatrix&    gtoValuesZ,
                                                                      const double*          rhograd,
                                                                      const double*          vrho,
                                                                      const double*          vsigma,
                                                                      CMultiTimer&           timer) const;

    /**
     Integrates LDA contribution to (second-order) Fxc matrix.

     @param npoints the number of grid points.
     @param weights the weights of grid points.
     @param gtoValues the GTO values on grid points.
     @param rhow the pointer to perturbed density.
     @param v2rho2 the 2nd-order functional derivative wrt density.
     @param timer the timer.
     @return the contribution as a CDenseMatrix object.
     */
    CDenseMatrix _integratePartialFxcFockForLDA(const int32_t         npoints,
                                                const double*         weights,
                                                const CDenseMatrix&   gtoValues,
                                                const double*         rhow,
                                                const double*         v2rho2,
                                                CMultiTimer&          timer) const;

    /**
     Integrates GGA contribution to (second-order) Fxc matrix.

     @param npoints the number of grid points.
     @param weights the weights of grid points.
     @param gtoValues the GTO values on grid points.
     @param gtoValuesX the GTO gradient X values on grid points.
     @param gtoValuesY the GTO gradient Y values on grid points.
     @param gtoValuesZ the GTO gradient Z values on grid points.
     @param rhow the pointer to perturbed density.
     @param rhograd the pointer to density gradient.
     @param rhowgrad the pointer to perturbed density gradient.
     @param v2rho2 the 2nd-order functional derivative wrt density.
     @param v2rhosigma the 2nd-order functional derivative wrt density and
            density gradient.
     @param v2sigma2 the 2nd-order functional derivative wrt density gradient.
     @param timer the timer.
     @return the contribution as a CDenseMatrix object.
     */
    CDenseMatrix _integratePartialFxcFockForGGA(const int32_t          npoints,
                                                const double*          weights,
                                                const CDenseMatrix&    gtoValues,
                                                const CDenseMatrix&    gtoValuesX,
                                                const CDenseMatrix&    gtoValuesY,
                                                const CDenseMatrix&    gtoValuesZ,
                                                const double*          rhow,
                                                const double*          rhograd,
                                                const double*          rhowgrad,
                                                const double*          vsigma,
                                                const double*          v2rho2,
                                                const double*          v2rhosigma,
                                                const double*          v2sigma2,
                                                CMultiTimer&           timer) const;

    /**
     Integrates LDA contribution to (third-order) Kxc matrix.

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
    CDenseMatrix _integratePartialKxcFockForLDA(const int32_t              npoints,
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
    CDenseMatrix _integratePartialKxcFockForGGA(const int32_t              npoints,
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
     Integrates first-order exchange-correlation functional contribution to AO
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
                                       CMolecularGrid&         molecularGrid,
                                       const std::string&      xcFuncLabel) const;

    /**
     Integrates second-order exchange-correlation functional contribution to AO
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
                          CMolecularGrid&         molecularGrid,
                          const std::string&      xcFuncLabel) const;

    /**
     Integrates third-order exchange-correlation functional contribution to AO
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
                          CMolecularGrid&         molecularGrid,
                          const std::string&      xcFuncLabel,
                          const std::string&      quadMode) const;

    /**
     Integrates first-order pair-density functional contribution to AO
     Fock matrix and MO "Q-matrix".

     @param aoFockMatrix the AO Fock matrix.
     @param moTwoBodyGradient the MO Two-body energy gradient term.
     @param molecule the molecule.
     @param basis the molecular basis.
     @param densityMatrix the AO density matrix object.
     @param TwoBodyDensityMatrix the MO two-body active density matrix.
     @param ActiveMOs the active molecular orbitals.
     @param molecularGrid the molecular grid.
     @param xcFuncLabel the label of exchange-correlation functional.
     */
    void integrateVxcPDFT(CAOKohnShamMatrix&      aoFockMatrix,
                          CDense4DTensor&         moTwoBodyGradient,
                          const CMolecule&        molecule,
                          const CMolecularBasis&  basis,
                          const CAODensityMatrix& DensityMatrix,
                          const CDense4DTensor&   TwoBodyDensityMatrix,
                          const CDenseMatrix&     ActiveMOs,
                          CMolecularGrid&         molecularGrid,
                          const std::string&      xcFuncLabel) const;

    /**
     Computes GTOs values on grid points.

     @param molecule the molecule.
     @param basis the molecular basis.
     @param molecularGrid the molecular grid.
     @return the GTO values on grid points.
     */
    CDenseMatrix computeGtoValuesOnGridPoints(const CMolecule& molecule, const CMolecularBasis& basis, CMolecularGrid& molecularGrid) const;

    /**
     Computes GTOs values and derivatives on grid points.

     @param molecule the molecule.
     @param basis the molecular basis.
     @param molecularGrid the molecular grid.
     @return the GTO values and derivatives on grid points.
     */
    std::vector<CDenseMatrix> computeGtoValuesAndDerivativesOnGridPoints(const CMolecule&       molecule,
                                                                         const CMolecularBasis& basis,
                                                                         CMolecularGrid&        molecularGrid) const;

    /**
     Computes fucntional derivatives for LDA.

     @param xcFuncLabel the label of exchange-correlation functional.
     @param npoints the number of grid points.
     @param rho the constant pointer to densities.
     @param exc the pointer to functional values.
     @param vrho the pointer to functional derivative w.r.t. densities.
     */
    void computeExcVxcForLDA(const std::string& xcFuncLabel, const int32_t npoints, const double* rho, double* exc, double* vrho) const;

    /**
     Computes fucntional derivatives for GGA.

     @param xcFuncLabel the label of exchange-correlation functional.
     @param npoints the number of grid points.
     @param rho the constant pointer to densities.
     @param sigma the constant pointer to density gradients.
     @param exc the pointer to functional values.
     @param vrho the pointer to functional derivative w.r.t. densities.
     @param vsigma the pointer to functional derivative w.r.t. density
            gradients.
     */
    void computeExcVxcForGGA(const std::string& xcFuncLabel,
                             const int32_t      npoints,
                             const double*      rho,
                             const double*      sigma,
                             double*            exc,
                             double*            vrho,
                             double*            vsigma) const;

    /**
     Computes 2nd-order fucntional derivatives for LDA.

     @param xcFuncLabel the label of exchange-correlation functional.
     @param npoints the number of grid points.
     @param rho the constant pointer to densities.
     @param v2rho2 the pointer to 2nd-order functional derivative w.r.t.
            densities.
     */
    void computeFxcForLDA(const std::string& xcFuncLabel, const int32_t npoints, const double* rho, double* v2rho2) const;

    /**
     Computes 2nd-order fucntional derivatives for GGA.

     @param xcFuncLabel the label of exchange-correlation functional.
     @param npoints the number of grid points.
     @param rho the constant pointer to densities.
     @param sigma the constant pointer to density gradients.
     @param v2rho2 the pointer to 2nd-order functional derivative w.r.t.
            densities.
     @param v2rhosigma the pointer to 2nd-order functional derivative w.r.t.
            densities and density gradients.
     @param v2sigma2 the pointer to 2nd-order functional derivative w.r.t.
            density gradients.
     */
    void computeFxcForGGA(const std::string& xcFuncLabel,
                          const int32_t      npoints,
                          const double*      rho,
                          const double*      sigma,
                          double*            v2rho2,
                          double*            v2rhosigma,
                          double*            v2sigma2) const;
};

#endif /* XCNewIntegrator_hpp */

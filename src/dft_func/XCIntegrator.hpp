//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//
//  Redistribution and use in source and binary forms, with or without modification,
//  are permitted provided that the following conditions are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this
//     list of conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//  3. Neither the name of the copyright holder nor the names of its contributors
//     may be used to endorse or promote products derived from this software without
//     specific prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
//  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
//  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
//  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
//  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
//  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
//  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef XCIntegrator_hpp
#define XCIntegrator_hpp

#include <array>
#include <string>

#include "AODensityMatrix.hpp"
#include "AOKohnShamMatrix.hpp"
#include "Dense4DTensor.hpp"
#include "DenseMatrix.hpp"
#include "GridBox.hpp"
#include "GtoBlock.hpp"
#include "MolecularBasis.hpp"
#include "MolecularGrid.hpp"
#include "Molecule.hpp"
#include "MultiTimer.hpp"
#include "XCFunctional.hpp"
#include "XCPairDensityFunctional.hpp"

/**
 Class CXCIntegrator implements XC integrator.
 */
class CXCIntegrator
{
   private:
    /**
     Screening threshold for GTO values on grid points.
     */
    double _screeningThresholdForGTOValues;

   public:
    /**
     Creates an XC integrator object.
     */
    CXCIntegrator();

    /**
     Integrates first-order exchange-correlation functional contribution to AO
     Kohn-Sham matrix.

     @param molecule the molecule.
     @param basis the molecular basis.
     @param gsDensityPointers the pointers to AO density matrices.
     @param molecularGrid the molecular grid.
     @param xcFuncLabel the label of exchange-correlation functional.
     @return the AO Kohn-Sham matrix.
     */
    auto integrateVxcFock(const CMolecule&                  molecule,
                          const CMolecularBasis&            basis,
                          const std::vector<const double*>& gsDensityPointers,
                          const CMolecularGrid&             molecularGrid,
                          const std::string&                xcFuncLabel) const -> CAOKohnShamMatrix;

    /**
     Integrates first-order exchange-correlation functional contribution to AO
     Kohn-Sham matrix.

     @param molecule the molecule.
     @param basis the molecular basis.
     @param gsDensityPointers the pointers to AO density matrices.
     @param molecularGrid the molecular grid.
     @param fvxc the exchange-correlation functional.
     @return the AO Kohn-Sham matrix.
     */
    auto integrateVxcFock(const CMolecule&                  molecule,
                          const CMolecularBasis&            basis,
                          const std::vector<const double*>& gsDensityPointers,
                          const CMolecularGrid&             molecularGrid,
                          const CXCFunctional&              fvxc) const -> CAOKohnShamMatrix;

    /**
     Integrates second-order exchange-correlation functional contribution to AO
     Fock matrix.

     @param aoFockPointers the pointers to AO Fock matrices.
     @param molecule the molecule.
     @param basis the molecular basis.
     @param rwDensityPointers the pointers to perturbed density matrices.
     @param gsDensityPointers the pointers to ground state density matrices.
     @param molecularGrid the molecular grid.
     @param xcFuncLabel the label of exchange-correlation functional.
     */
    auto integrateFxcFock(const std::vector<double*>&       aoFockPointers,
                          const CMolecule&                  molecule,
                          const CMolecularBasis&            basis,
                          const std::vector<const double*>& rwDensityPointers,
                          const std::vector<const double*>& gsDensityPointers,
                          const CMolecularGrid&             molecularGrid,
                          const std::string&                xcFuncLabel) const -> void;

    /**
     Integrates second-order exchange-correlation functional contribution to AO
     Fock matrix.

     @param aoFockPointers the pointers to AO Fock matrices.
     @param molecule the molecule.
     @param basis the molecular basis.
     @param rwDensityPointers the pointers to perturbed density matrices.
     @param gsDensityPointers the pointers to ground state density matrices.
     @param molecularGrid the molecular grid.
     @param fvxc the exchange-correlation functional.
     */
    auto integrateFxcFock(const std::vector<double*>&       aoFockPointers,
                          const CMolecule&                  molecule,
                          const CMolecularBasis&            basis,
                          const std::vector<const double*>& rwDensityPointers,
                          const std::vector<const double*>& gsDensityPointers,
                          const CMolecularGrid&             molecularGrid,
                          const CXCFunctional&              fvxc) const -> void;

    /**
     Integrates third-order exchange-correlation functional contribution to AO
     Fock matrix in quadratic response.

     @param aoFockPointers the pointers to AO Fock matrices.
     @param molecule the molecule.
     @param basis the molecular basis.
     @param rwDensityPointers the pointers to one-time transformed densities.
     @param rw2DensityPointers the pointers to two-time transformed densities.
     @param gsDensityPointers the pointers to ground state density matrix.
     @param molecularGrid the molecular grid.
     @param xcFuncLabel the label of exchange-correlation functional.
     @param quadMode a string that specifies which densities should be combined.
     */
    auto integrateKxcFock(const std::vector<double*>& aoFockPointers,
                          const CMolecule&        molecule,
                          const CMolecularBasis&  basis,
                          const std::vector<const double*>& rwDensityPointers,
                          const std::vector<const double*>& rw2DensityPointers,
                          const std::vector<const double*>& gsDensityPointers,
                          const CMolecularGrid&   molecularGrid,
                          const std::string&      xcFuncLabel,
                          const std::string&      quadMode) const -> void;

    /**
     Integrates third-order exchange-correlation functional contribution to AO
     Fock matrix in quadratic response.

     @param aoFockPointers the pointers to AO Fock matrices.
     @param molecule the molecule.
     @param basis the molecular basis.
     @param rwDensityPointers the pointers to one-time transformed densities.
     @param rw2DensityPointers the pointers to two-time transformed densities.
     @param gsDensityPointers the pointers to ground state density matrix.
     @param molecularGrid the molecular grid.
     @param fvxc the exchange-correlation functional.
     @param quadMode a string that specifies which densities should be combined.
     */
    auto integrateKxcFock(const std::vector<double*>& aoFockPointers,
                          const CMolecule&        molecule,
                          const CMolecularBasis&  basis,
                          const std::vector<const double*>& rwDensityPointers,
                          const std::vector<const double*>& rw2DensityPointers,
                          const std::vector<const double*>& gsDensityPointers,
                          const CMolecularGrid&   molecularGrid,
                          const CXCFunctional&    fvxc,
                          const std::string&      quadMode) const -> void;

    /**
     Integrates fourth-order exchnage-correlation functional contribution to AO
     Fock matrix in cubic response.

     @param aoFockPointers the pointers to AO Fock matrix.
     @param molecule the molecule.
     @param basis the molecular basis.
     @param rwDensityPointers the pointers to one-time transformed densities.
     @param rw2DensityPointers the pointers to two-time transformed densities.
     @param rw3DensityPointers the pointers to three-time transformed densities.
     @param gsDensityPointers the pointers to ground state density matrix.
     @param molecularGrid the molecular grid.
     @param xcFuncLabel the label of exchange-correlation functional.
     @param cubeMode a string that specifies which densities should be combined.
     */
    auto integrateKxcLxcFock(const std::vector<double*>& aoFockPointers,
                             const CMolecule&        molecule,
                             const CMolecularBasis&  basis,
                             const std::vector<const double*>& rwDensityPointers,
                             const std::vector<const double*>& rw2DensityPointers,
                             const std::vector<const double*>& rw3DensityPointers,
                             const std::vector<const double*>& gsDensityPointers,
                             const CMolecularGrid&   molecularGrid,
                             const std::string&      xcFuncLabel,
                             const std::string&      cubeMode) const -> void;

    /**
     Integrates fourth-order exchnage-correlation functional contribution to AO
     Fock matrix in cubic response.

     @param aoFockPointers the pointers to AO Fock matrix.
     @param molecule the molecule.
     @param basis the molecular basis.
     @param rwDensityPointers the pointers to one-time transformed densities.
     @param rw2DensityPointers the pointers to two-time transformed densities.
     @param rw3DensityPointers the pointers to three-time transformed densities.
     @param gsDensityPointers the pointers to ground state density matrix.
     @param molecularGrid the molecular grid.
     @param fvxc the exchange-correlation functional.
     @param cubeMode a string that specifies which densities should be combined.
     */
    auto integrateKxcLxcFock(const std::vector<double*>& aoFockPointers,
                             const CMolecule&        molecule,
                             const CMolecularBasis&  basis,
                             const std::vector<const double*>& rwDensityPointers,
                             const std::vector<const double*>& rw2DensityPointers,
                             const std::vector<const double*>& rw3DensityPointers,
                             const std::vector<const double*>& gsDensityPointers,
                             const CMolecularGrid&   molecularGrid,
                             const CXCFunctional&    fvxc,
                             const std::string&      cubeMode) const -> void;

    /**
     Integrates first-order pair-density functional contribution to AO
     Fock matrix and MO "Q-matrix".

     @param aoFockMatrix the AO Fock matrix.
     @param tensorWxc the MO Two-body energy gradient term.
     @param molecule the molecule.
     @param basis the molecular basis.
     @param densityMatrixPointer the pointer to AO density matrix.
     @param twoBodyDensityMatrix the MO two-body active density matrix.
     @param activeMOs the active molecular orbitals.
     @param molecularGrid the molecular grid.
     @param fvxc exchange-correlation functional.
     @param rs_omega range-separation parameter.
     */
    auto integrateVxcPDFT(CAOKohnShamMatrix&                  aoFockMatrix,
                          CDense4DTensor&                     tensorWxc,
                          const CMolecule&                    molecule,
                          const CMolecularBasis&              basis,
                          const double*                       densityMatrixPointer,
                          const CDenseMatrix&                 twoBodyDensityMatrix,
                          const CDenseMatrix&                 activeMOs,
                          const CMolecularGrid&               molecularGrid,
                          const CXCPairDensityFunctional&     fvxc,
                          const double                        rs_omega) const -> void;

    /**
     Computes GTOs values on grid points.

     @param molecule the molecule.
     @param basis the molecular basis.
     @param molecularGrid the molecular grid.
     @return the GTO values on grid points.
     */
    auto computeGtoValuesOnGridPoints(const CMolecule&       molecule,
                                      const CMolecularBasis& basis,
                                      const CMolecularGrid&  molecularGrid) const -> CDenseMatrix;

    auto computeGtoValuesAndDerivativesOnGridPoints(const CMolecule&       molecule,
                                                    const CMolecularBasis& basis,
                                                    const CMolecularGrid&  molecularGrid) const -> std::vector<CDenseMatrix>;

    auto computeGtoValuesAndSecondOrderDerivativesOnGridPoints(const CMolecule&       molecule,
                                                               const CMolecularBasis& basis,
                                                               const CMolecularGrid&  molecularGrid) const -> std::vector<CDenseMatrix>;
};

#endif /* XCIntegrator_hpp */

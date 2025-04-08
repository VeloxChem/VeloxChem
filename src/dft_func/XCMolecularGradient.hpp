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

#ifndef XCMolecularGradient_hpp
#define XCMolecularGradient_hpp

#include <array>
#include <list>
#include <string>

#include "DenseMatrix.hpp"
#include "MolecularBasis.hpp"
#include "MolecularGrid.hpp"
#include "Molecule.hpp"
#include "XCFunctional.hpp"
#include "XCPairDensityFunctional.hpp"

/**
 Class CXCMolecularGradient implements XC molecular gradient.
 */
class CXCMolecularGradient
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
    CXCMolecularGradient();

    /**
     Integrates first-order exchnage-correlation functional contribution to
     molecular gradient.

     @param molecule the molecule.
     @param basis the molecular basis.
     @param gsDensityPointers the pointers to ground state AO density matrix.
     @param molecularGrid the molecular grid.
     @param xcFuncLabel the label of exchange-correlation functional.
     @return the molecular gradient.
     */
    CDenseMatrix integrateVxcGradient(const CMolecule&                  molecule,
                                      const CMolecularBasis&            basis,
                                      const std::vector<const double*>& gsDensityPointers,
                                      const CMolecularGrid&             molecularGrid,
                                      const std::string&                xcFuncLabel) const;

    /**
     Integrates first-order exchnage-correlation functional contribution to
     molecular gradient.

     @param molecule the molecule.
     @param basis the molecular basis.
     @param rwDensityPointers the pointers to perturbed AO density matrix (to
            be contracted with GTO gradient).
     @param gsDensityPointers the pointers to ground state AO density matrix.
     @param molecularGrid the molecular grid.
     @param xcFuncLabel the label of exchange-correlation functional.
     @return the molecular gradient.
     */
    CDenseMatrix integrateVxcGradient(const CMolecule&                  molecule,
                                      const CMolecularBasis&            basis,
                                      const std::vector<const double*>& rwDensityPointers,
                                      const std::vector<const double*>& gsDensityPointers,
                                      const CMolecularGrid&             molecularGrid,
                                      const std::string&                xcFuncLabel) const;

    /**
     Integrates second-order exchnage-correlation functional contribution to
     molecular gradient.

     @param molecule the molecule.
     @param basis the molecular basis.
     @param rwDensityPointersOne the pointers to perturbed AO density matrix.
     @param rwDensityPointersTwo the pointers to perturbed AO density matrix (to be
            contracted with GTO gradient).
     @param gsDensityPointers the pointers to ground state AO density matrix.
     @param molecularGrid the molecular grid.
     @param xcFunctional the exchange-correlation functional.
     @return the molecular gradient.
     */
    CDenseMatrix integrateFxcGradient(const CMolecule&                  molecule,
                                      const CMolecularBasis&            basis,
                                      const std::vector<const double*>& rwDensityPointersOne,
                                      const std::vector<const double*>& rwDensityPointersTwo,
                                      const std::vector<const double*>& gsDensityPointers,
                                      const CMolecularGrid&             molecularGrid,
                                      const std::string&                xcFuncLabel) const;

    /**
     Integrates third-order exchnage-correlation functional contribution to
     molecular gradient.

     @param molecule the molecule.
     @param basis the molecular basis.
     @param rwDensityPointersOne the pointers to perturbed AO density matrix.
     @param rwDensityPointersTwo the pointers to perturbed AO density matrix.
     @param gsDensityPointers the pointers to ground state AO density matrix (to be
            contracted with GTO gradient).
     @param molecularGrid the molecular grid.
     @param xcFuncLabel the label of exchange-correlation functional.
     @return the molecular gradient.
     */
    CDenseMatrix integrateKxcGradient(const CMolecule&                  molecule,
                                      const CMolecularBasis&            basis,
                                      const std::vector<const double*>& rwDensityPointersOne,
                                      const std::vector<const double*>& rwDensityPointersTwo,
                                      const std::vector<const double*>& gsDensityPointers,
                                      const CMolecularGrid&             molecularGrid,
                                      const std::string&                xcFuncLabel) const;

    /**
     Integrates first-order exchange-correlation functional contribution to
     PDFT molecular gradient.

     @param molecule the molecule.
     @param basis the molecular basis.
     @param densityMatrix the AO density matrix object.
     @param twoBodyDensityMatrix the MO two-body active density matrix.
     @param activeMOs the active molecular orbitals.
     @param molecularGrid the molecular grid.
     @param fvxc exchange-correlation functional.
     @param rs_omega range-separation parameter.
     @return the molecular gradient.
     */
    CDenseMatrix integrateVxcPDFTGradient(const CMolecule&                    molecule,
                                          const CMolecularBasis&              basis,
                                          const double*                       densityMatrixPointer,
                                          const CDenseMatrix&                 twoBodyDensityMatrix,
                                          const CDenseMatrix&                 activeMOs,
                                          const CMolecularGrid&               molecularGrid,
                                          const CXCPairDensityFunctional&     fvxc,
                                          const double                        rs_omega) const;
};

#endif /* XCMolecularGradient_hpp */

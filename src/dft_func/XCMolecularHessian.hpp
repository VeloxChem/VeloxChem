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

#ifndef XCMolecularHessian_hpp
#define XCMolecularHessian_hpp

#include <vector>
#include <string>

#include "DenseMatrix.hpp"
#include "MolecularBasis.hpp"
#include "MolecularGrid.hpp"
#include "Molecule.hpp"
#include "XCFunctional.hpp"

/**
 Class CXCMolecularHessian implements XC molecular Hessian.
 */
class CXCMolecularHessian
{
    /**
     Screening threshold for GTO values on grid points.
     */
    double _screeningThresholdForGTOValues;

   public:
    /**
     Creates an XC integrator object.
     */
    CXCMolecularHessian();

    /**
     Integrates exchnage-correlation functional contribution to molecular
     Hessian.

     @param molecule the molecule.
     @param basis the molecular basis.
     @param gsDensityPointers the pointers to ground state AO density matrix.
     @param molecularGrid the molecular grid.
     @param xcFuncLabel the label of exchange-correlation functional.
     @return the molecular Hessian.
     */
    CDenseMatrix integrateExcHessian(const CMolecule&                  molecule,
                                     const CMolecularBasis&            basis,
                                     const std::vector<const double*>& gsDensityPointers,
                                     const CMolecularGrid&             molecularGrid,
                                     const std::string&                xcFuncLabel) const;

    /**
     Integrates exchnage-correlation functional contribution to molecular
     gradient of Vxc matrix element.

     @param molecule the molecule.
     @param basis the molecular basis.
     @param gsDensityPointers the pointers to ground state AO density matrix.
     @param molecularGrid the molecular grid.
     @param xcFuncLabel the label of exchange-correlation functional.
     @param atomIdxVec the indices of the atoms with respect to which gradient is
     computed.
     @return the Vxc gradient.
     */
    std::vector<CDenseMatrix> integrateVxcFockGradient(const CMolecule&                  molecule,
                                                       const CMolecularBasis&            basis,
                                                       const std::vector<const double*>& gsDensityPointers,
                                                       const CMolecularGrid&             molecularGrid,
                                                       const std::string&                xcFuncLabel,
                                                       const std::vector<int>&           atomIdxVec) const;
};

#endif /* XCMolecularHessian_hpp */

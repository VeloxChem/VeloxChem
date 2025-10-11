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

#include "XCMolecularHessian.hpp"

#include "ErrorHandler.hpp"
#include "FunctionalParser.hpp"
#include "XCMolecularHessianForLDA.hpp"
#include "XCMolecularHessianForGGA.hpp"
#include "XCMolecularHessianForMGGA.hpp"

CXCMolecularHessian::CXCMolecularHessian()

    : _screeningThresholdForGTOValues(1.0e-12)
{
}

auto
CXCMolecularHessian::integrateExcHessian(const CMolecule&        molecule,
                                         const CMolecularBasis&  basis,
                                         const std::vector<const double*>& gsDensityPointers,
                                         const CMolecularGrid&   molecularGrid,
                                         const std::string&      xcFuncLabel) const -> CDenseMatrix
{
    auto fvxc = vxcfuncs::getExchangeCorrelationFunctional(xcFuncLabel);

    auto xcfuntype = fvxc.getFunctionalType();

    if (gsDensityPointers.size() == 1)
    {
        if (xcfuntype == xcfun::lda)
        {
            return xchesslda::integrateExcHessianForLdaClosedShell(molecule, basis, gsDensityPointers, molecularGrid, _screeningThresholdForGTOValues, fvxc);
        }
        else if (xcfuntype == xcfun::gga)
        {
            return xchessgga::integrateExcHessianForGgaClosedShell(molecule, basis, gsDensityPointers, molecularGrid, _screeningThresholdForGTOValues, fvxc);
        }
        else
        {
            return xchessmgga::integrateExcHessianForMetaGgaClosedShell(molecule, basis, gsDensityPointers, molecularGrid, _screeningThresholdForGTOValues, fvxc);
        }
    }
    else
    {
        if (xcfuntype == xcfun::lda)
        {
            return xchesslda::integrateExcHessianForLdaOpenShell(molecule, basis, gsDensityPointers, molecularGrid, _screeningThresholdForGTOValues, fvxc);
        }
        else
        {
            std::string erropenshell("XCMolecularHessian.integrateExcHessian: Not implemented for open-shell GGA/meta-GGA");

            errors::assertMsgCritical(false, erropenshell);
        }
    }

    return CDenseMatrix();
}

auto
CXCMolecularHessian::integrateVxcFockGradient(const CMolecule&        molecule,
                                              const CMolecularBasis&  basis,
                                              const std::vector<const double*>& gsDensityPointers,
                                              const CMolecularGrid&   molecularGrid,
                                              const std::string&      xcFuncLabel,
                                              const std::vector<int>& atomIdxVec) const -> std::vector<CDenseMatrix>
{
    auto fvxc = vxcfuncs::getExchangeCorrelationFunctional(xcFuncLabel);

    auto xcfuntype = fvxc.getFunctionalType();

    if (gsDensityPointers.size() == 1)
    {
        if (xcfuntype == xcfun::lda)
        {
            return xchesslda::integrateVxcFockGradientForLdaClosedShell(molecule, basis, gsDensityPointers, molecularGrid, _screeningThresholdForGTOValues, fvxc, atomIdxVec);
        }
        else if (xcfuntype == xcfun::gga)
        {
            return xchessgga::integrateVxcFockGradientForGgaClosedShell(molecule, basis, gsDensityPointers, molecularGrid, _screeningThresholdForGTOValues, fvxc, atomIdxVec);
        }
        else
        {
            return xchessmgga::integrateVxcFockGradientForMetaGgaClosedShell(molecule, basis, gsDensityPointers, molecularGrid, _screeningThresholdForGTOValues, fvxc, atomIdxVec);
        }
    }
    else
    {
        if (xcfuntype == xcfun::lda)
        {
            return xchesslda::integrateVxcFockGradientForLdaOpenShell(molecule, basis, gsDensityPointers, molecularGrid, _screeningThresholdForGTOValues, fvxc, atomIdxVec);
        }
        else if (xcfuntype == xcfun::gga)
        {
            return xchessgga::integrateVxcFockGradientForGgaOpenShell(molecule, basis, gsDensityPointers, molecularGrid, _screeningThresholdForGTOValues, fvxc, atomIdxVec);
        }
        else
        {
            std::string erropenshell("XCMolecularHessian.integrateVxcFockGradient: Not implemented for open-shell meta-GGA");

            errors::assertMsgCritical(false, erropenshell);
        }
    }

    return std::vector<CDenseMatrix>();
}

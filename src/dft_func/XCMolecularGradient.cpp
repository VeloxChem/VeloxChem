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

#include "XCMolecularGradient.hpp"

#include <omp.h>

#include <cmath>
#include <cstring>

#include "DenseMatrix.hpp"
#include "ErrorHandler.hpp"
#include "FunctionalParser.hpp"
#include "XCFunctional.hpp"
#include "XCMolecularGradientForLDA.hpp"
#include "XCMolecularGradientForGGA.hpp"
#include "XCMolecularGradientForMGGA.hpp"
#include "XCMolecularGradientForPDFT.hpp"

CXCMolecularGradient::CXCMolecularGradient()

    : _screeningThresholdForGTOValues(1.0e-12)
{
}

CDenseMatrix
CXCMolecularGradient::integrateVxcGradient(const CMolecule&        molecule,
                                           const CMolecularBasis&  basis,
                                           const std::vector<const double*>& gsDensityPointers,
                                           const CMolecularGrid&   molecularGrid,
                                           const std::string&      xcFuncLabel) const
{
    return integrateVxcGradient(molecule, basis, gsDensityPointers, gsDensityPointers, molecularGrid, xcFuncLabel);
}

CDenseMatrix
CXCMolecularGradient::integrateVxcGradient(const CMolecule&        molecule,
                                           const CMolecularBasis&  basis,
                                           const std::vector<const double*>& rwDensityPointers,
                                           const std::vector<const double*>& gsDensityPointers,
                                           const CMolecularGrid&   molecularGrid,
                                           const std::string&      xcFuncLabel) const
{
    auto fvxc = vxcfuncs::getExchangeCorrelationFunctional(xcFuncLabel);

    auto xcfuntype = fvxc.getFunctionalType();

    if (gsDensityPointers.size() == 1)
    {
        if (xcfuntype == xcfun::lda)
        {
            return xcgradlda::integrateVxcGradientForLdaClosedShell(molecule, basis, rwDensityPointers, gsDensityPointers, molecularGrid, _screeningThresholdForGTOValues, fvxc);
        }
        else if (xcfuntype == xcfun::gga)
        {
            return xcgradgga::integrateVxcGradientForGgaClosedShell(molecule, basis, rwDensityPointers, gsDensityPointers, molecularGrid, _screeningThresholdForGTOValues, fvxc);
        }
        else if (xcfuntype == xcfun::mgga)
        {
            return xcgradmgga::integrateVxcGradientForMetaGgaClosedShell(molecule, basis, rwDensityPointers, gsDensityPointers, molecularGrid, _screeningThresholdForGTOValues, fvxc);
        }
        else
        {
            std::string errxcfuntype("XCMolecularGradient.integrateVxcGradient: Only implemented for LDA/GGA/meta-GGA");

            errors::assertMsgCritical(false, errxcfuntype);
        }
    }
    else
    {
        if (xcfuntype == xcfun::lda)
        {
            return xcgradlda::integrateVxcGradientForLdaOpenShell(molecule, basis, rwDensityPointers, gsDensityPointers, molecularGrid, _screeningThresholdForGTOValues, fvxc);
        }
        else if (xcfuntype == xcfun::gga)
        {
            return xcgradgga::integrateVxcGradientForGgaOpenShell(molecule, basis, rwDensityPointers, gsDensityPointers, molecularGrid, _screeningThresholdForGTOValues, fvxc);
        }
        else if (xcfuntype == xcfun::mgga)
        {
            return xcgradmgga::integrateVxcGradientForMetaGgaOpenShell(molecule, basis, rwDensityPointers, gsDensityPointers, molecularGrid, _screeningThresholdForGTOValues, fvxc);
        }
        else
        {
            std::string errxcfuntype("XCMolecularGradient.integrateVxcGradient: Only implemented for LDA/GGA/meta-GGA");

            errors::assertMsgCritical(false, errxcfuntype);
        }
    }

    return CDenseMatrix();
}

CDenseMatrix
CXCMolecularGradient::integrateFxcGradient(const CMolecule&        molecule,
                                           const CMolecularBasis&  basis,
                                           const std::vector<const double*>& rwDensityPointersOne,
                                           const std::vector<const double*>& rwDensityPointersTwo,
                                           const std::vector<const double*>& gsDensityPointers,
                                           const CMolecularGrid&   molecularGrid,
                                           const std::string&      xcFuncLabel) const
{
    auto fvxc = vxcfuncs::getExchangeCorrelationFunctional(xcFuncLabel);

    auto xcfuntype = fvxc.getFunctionalType();

    if (gsDensityPointers.size() == 1)
    {
        if (xcfuntype == xcfun::lda)
        {
            return xcgradlda::integrateFxcGradientForLdaClosedShell(molecule, basis, rwDensityPointersOne, rwDensityPointersTwo, gsDensityPointers, molecularGrid, _screeningThresholdForGTOValues, fvxc);
        }
        else if (xcfuntype == xcfun::gga)
        {
            return xcgradgga::integrateFxcGradientForGgaClosedShell(molecule, basis, rwDensityPointersOne, rwDensityPointersTwo, gsDensityPointers, molecularGrid, _screeningThresholdForGTOValues, fvxc);
        }
        else
        {
            std::string errxcfuntype("XCMolecularGradient.integrateFxcGradient: Only implemented for LDA/GGA");

            errors::assertMsgCritical(false, errxcfuntype);
        }
    }
    else
    {
        std::string erropenshell("XCMolecularGradient.integrateFxcGradient: Not implemented for open-shell");

        errors::assertMsgCritical(false, erropenshell);
    }

    return CDenseMatrix();
}

CDenseMatrix
CXCMolecularGradient::integrateKxcGradient(const CMolecule&        molecule,
                                           const CMolecularBasis&  basis,
                                           const std::vector<const double*>& rwDensityPointersOne,
                                           const std::vector<const double*>& rwDensityPointersTwo,
                                           const std::vector<const double*>& gsDensityPointers,
                                           const CMolecularGrid&   molecularGrid,
                                           const std::string&      xcFuncLabel) const
{
    auto fvxc = vxcfuncs::getExchangeCorrelationFunctional(xcFuncLabel);

    auto xcfuntype = fvxc.getFunctionalType();

    if (gsDensityPointers.size() == 1)
    {
        if (xcfuntype == xcfun::lda)
        {
            return xcgradlda::integrateKxcGradientForLdaClosedShell(molecule, basis, rwDensityPointersOne, rwDensityPointersTwo, gsDensityPointers, molecularGrid, _screeningThresholdForGTOValues, fvxc);
        }
        else if (xcfuntype == xcfun::gga)
        {
            return xcgradgga::integrateKxcGradientForGgaClosedShell(molecule, basis, rwDensityPointersOne, rwDensityPointersTwo, gsDensityPointers, molecularGrid, _screeningThresholdForGTOValues, fvxc);
        }
        else
        {
            std::string errxcfuntype("XCMolecularGradient.integrateKxcGradient: Only implemented for LDA/GGA");

            errors::assertMsgCritical(false, errxcfuntype);
        }
    }
    else
    {
        std::string erropenshell("XCMolecularGradient.integrateKxcGradient: Not implemented for open-shell");

        errors::assertMsgCritical(false, erropenshell);
    }

    return CDenseMatrix();
}

auto
CXCMolecularGradient::integrateVxcPDFTGradient(const CMolecule&                    molecule,
                                        const CMolecularBasis&              basis,
                                        const double*                       densityMatrixPointer,
                                        const CDenseMatrix&                 twoBodyDensityMatrix,
                                        const CDenseMatrix&                 activeMOs,
                                        const CMolecularGrid&               molecularGrid,
                                        const CXCPairDensityFunctional&     fvxc,
                                        const double                        rs_omega) const -> CDenseMatrix
{
    auto xcfuntype = fvxc.getFunctionalType();

    if (xcfuntype == "PLDA")
    {
        return xcgradpdft::integrateVxcPDFTGradientForLDA(molecule,
                                                  basis,
                                                  densityMatrixPointer,
                                                  twoBodyDensityMatrix,
                                                  activeMOs,
                                                  molecularGrid,
                                                 _screeningThresholdForGTOValues,
                                                  fvxc, rs_omega);
    }
    else if (xcfuntype == "PGGA")
    {
        return xcgradpdft::integrateVxcPDFTGradientForGGA(molecule,
                                                  basis,
                                                  densityMatrixPointer,
                                                  twoBodyDensityMatrix,
                                                  activeMOs,
                                                  molecularGrid,
                                                  _screeningThresholdForGTOValues,
                                                  fvxc, rs_omega);
    }
    else
    {
        std::string errxcfuntype("XCMolecularGradient.integrateVxcPDFTGradient: Only implemented for PLDA/PGGA");

        errors::assertMsgCritical(false, errxcfuntype);
    }
    return CDenseMatrix();
}

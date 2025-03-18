//
//                              VELOXCHEM
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2024 by VeloxChem developers. All rights reserved.
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

#include "XCMolecularGradient.hpp"

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "AODensityMatrix.hpp"
#include "DensityGridGenerator.hpp"
#include "DftSubMatrix.hpp"
#include "ErrorHandler.hpp"
#include "FunctionalParser.hpp"
#include "GridScreener.hpp"
#include "GtoFunc.hpp"
#include "GtoValues.hpp"
#include "Prescreener.hpp"
#include "XCFunctional.hpp"
#include "XCMolecularGradientForLDA.hpp"
#include "XCMolecularGradientForGGA.hpp"
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
        else
        {
            std::string errxcfuntype("XCMolecularGradient.integrateVxcGradient: Only implemented for LDA/GGA");

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
        else
        {
            std::string errxcfuntype("XCMolecularGradient.integrateVxcGradient: Only implemented for LDA/GGA");

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

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
            std::string errxcfuntype("XCMolecularHessian.integrateExcHessian: Only implemented for LDA/GGA");

            errors::assertMsgCritical(false, errxcfuntype);
        }
    }
    else
    {
        std::string erropenshell("XCMolecularHessian.integrateExcHessian: Not implemented for open-shell");

        errors::assertMsgCritical(false, erropenshell);
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
            return xchesslda::integrateVxcFockGradientForLDA(molecule, basis, gsDensityPointers, molecularGrid, _screeningThresholdForGTOValues, fvxc, atomIdxVec);
        }
        else if (xcfuntype == xcfun::gga)
        {
            return xchessgga::integrateVxcFockGradientForGGA(molecule, basis, gsDensityPointers, molecularGrid, _screeningThresholdForGTOValues, fvxc, atomIdxVec);
        }
        else
        {
            return xchessmgga::integrateVxcFockGradientForMetaGGA(molecule, basis, gsDensityPointers, molecularGrid, _screeningThresholdForGTOValues, fvxc, atomIdxVec);
        }
    }
    else
    {
        std::string erropenshell("XCMolecularHessian.integrateVxcFockGradient: Not implemented for open-shell");

        errors::assertMsgCritical(false, erropenshell);
    }

    return std::vector<CDenseMatrix>();
}

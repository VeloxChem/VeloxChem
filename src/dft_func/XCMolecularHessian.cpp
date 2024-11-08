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
            return xchesslda::integrateExcHessianForLDA(molecule, basis, gsDensityPointers, molecularGrid, _screeningThresholdForGTOValues, fvxc);
        }
        else if (xcfuntype == xcfun::gga)
        {
            //return _integrateExcHessianForGGA(molecule, basis, gsDensityMatrix, molecularGrid, fvxc);
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
                                              const int               atomIdx) const -> std::vector<CDenseMatrix>
{
    auto fvxc = vxcfuncs::getExchangeCorrelationFunctional(xcFuncLabel);

    auto xcfuntype = fvxc.getFunctionalType();

    if (gsDensityPointers.size() == 1)
    {
        if (xcfuntype == xcfun::lda)
        {
            return xchesslda::integrateVxcFockGradientForLDA(molecule, basis, gsDensityPointers, molecularGrid, _screeningThresholdForGTOValues, fvxc, atomIdx);
        }
        else if (xcfuntype == xcfun::gga)
        {
            //return _integrateVxcFockGradientForGGA(molecule, basis, gsDensityMatrix, molecularGrid, fvxc, atomIdx);
        }
        else
        {
            std::string errxcfuntype("XCMolecularHessian.integrateVxcFockGradient: Only implemented for LDA/GGA");

            errors::assertMsgCritical(false, errxcfuntype);
        }
    }
    else
    {
        std::string erropenshell("XCMolecularHessian.integrateVxcFockGradient: Not implemented for open-shell");

        errors::assertMsgCritical(false, erropenshell);
    }

    return std::vector<CDenseMatrix>();
}

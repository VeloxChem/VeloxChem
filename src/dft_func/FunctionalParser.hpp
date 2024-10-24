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

#ifndef FunctionalParser_hpp
#define FunctionalParser_hpp

#include <string>
#include <vector>

#include "XCFunctional.hpp"
#include "XCPairDensityFunctional.hpp"

namespace vxcfuncs {  // vxcfuncs namespace

/**
 Gets labels of available exchange-correlation functional.

 @return a vector of labels of available exchange-correlation functionals.
 */
std::vector<std::string> getAvailableFunctionals();

/**
 Converts exchange-correlation functional label to exchange-correlation
 functional object.

 @param xcLabel the label of exchange-correlation functional.
 @return the exchange-correlation functional object.
 */
CXCFunctional getExchangeCorrelationFunctional(const std::string &xcLabel);

/**
 Gets labels of available pair-density exchange-correlation functional components.

 @return a vector of labels of available exchange-correlation functional components.
 */
std::vector<std::string> getAvailablePairDensityFunctionals();

}  // namespace vxcfuncs

#endif /* FunctionalParser_hpp */

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

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


#ifndef SimdKineticEnergyFunc_hpp
#define SimdKineticEnergyFunc_hpp

#include <cstddef>
#include <string>
#include <vector>

#include "BasisFunction.hpp"
#include "ErrorHandler.hpp"
#include "SimdMatrix.hpp"
#include "SimdKineticEnergyRecSS.hpp"
#include "SimdKineticEnergyRecSP.hpp"
#include "SimdKineticEnergyRecSD.hpp"
#include "SimdKineticEnergyRecSF.hpp"
#include "SimdKineticEnergyRecSG.hpp"
#include "SimdKineticEnergyRecSH.hpp"
#include "SimdKineticEnergyRecSI.hpp"
#include "SimdKineticEnergyRecPS.hpp"
#include "SimdKineticEnergyRecDS.hpp"
#include "SimdKineticEnergyRecFS.hpp"
#include "SimdKineticEnergyRecGS.hpp"
#include "SimdKineticEnergyRecHS.hpp"
#include "SimdKineticEnergyRecIS.hpp"
#include "SimdKineticEnergyRecPP.hpp"
#include "SimdKineticEnergyRecPD.hpp"
#include "SimdKineticEnergyRecPF.hpp"
#include "SimdKineticEnergyRecPG.hpp"
#include "SimdKineticEnergyRecPH.hpp"
#include "SimdKineticEnergyRecPI.hpp"
#include "SimdKineticEnergyRecDP.hpp"
#include "SimdKineticEnergyRecDD.hpp"
#include "SimdKineticEnergyRecDF.hpp"
#include "SimdKineticEnergyRecDG.hpp"
#include "SimdKineticEnergyRecDH.hpp"
#include "SimdKineticEnergyRecDI.hpp"
#include "SimdKineticEnergyRecFP.hpp"
#include "SimdKineticEnergyRecFD.hpp"
#include "SimdKineticEnergyRecFF.hpp"
#include "SimdKineticEnergyRecFG.hpp"
#include "SimdKineticEnergyRecFH.hpp"
#include "SimdKineticEnergyRecFI.hpp"
#include "SimdKineticEnergyRecGP.hpp"
#include "SimdKineticEnergyRecGD.hpp"
#include "SimdKineticEnergyRecGF.hpp"
#include "SimdKineticEnergyRecGG.hpp"
#include "SimdKineticEnergyRecGH.hpp"
#include "SimdKineticEnergyRecGI.hpp"
#include "SimdKineticEnergyRecHP.hpp"
#include "SimdKineticEnergyRecHD.hpp"
#include "SimdKineticEnergyRecHF.hpp"
#include "SimdKineticEnergyRecHG.hpp"
#include "SimdKineticEnergyRecHH.hpp"
#include "SimdKineticEnergyRecHI.hpp"
#include "SimdKineticEnergyRecIP.hpp"
#include "SimdKineticEnergyRecID.hpp"
#include "SimdKineticEnergyRecIF.hpp"
#include "SimdKineticEnergyRecIG.hpp"
#include "SimdKineticEnergyRecIH.hpp"
#include "SimdKineticEnergyRecII.hpp"

namespace simdkin {  // simdkin namespace

/// @brief Computes the kinetic energy integral of two basis functions centered on
/// the same atom.
/// @param bra The basis function on bra side.
/// @param ket The basis function on ket side.
/// @return The kinetic energy integral of the basis functions.
/// @note The kinetic energy of two solid harmonic primitives on the same center
/// is diagonal in the angular components and does not depend on the position of
/// the atom, as the overlap is, so the integral is a single value. Basis
/// functions of different angular momenta have no integral and give zero.
auto one_center_kinetic_energy(const CBasisFunction &bra, const CBasisFunction &ket) -> double;

/// @brief Computes the kinetic energy integrals of a combination of basis
/// functions over the atom pairs of a block.
/// @param values The values of the combination of basis functions to compute.
/// @param nvalues The number of atom pairs the combination reaches.
/// @param bra The basis function on bra side.
/// @param ket The basis function on ket side.
/// @param harmonics The solid harmonics of the vectors between the atoms.
/// @param coordinates The coordinates of the atom pairs.
/// @param threshold The screening threshold of the integrals.
/// @note The geometry of the kinetic energy reaches the sum of the angular
/// momenta of the basis functions, as the overlap does, so the harmonics are
/// created to that sum. The squared distance enters to higher powers, but that
/// is carried by the coordinates and not by the harmonics.
auto compute_kinetic_energy(double                         *values,
                            const size_t                    nvalues,
                            const CBasisFunction           &bra,
                            const CBasisFunction           &ket,
                            const std::vector<CSimdMatrix> &harmonics,
                            const CSimdMatrix              &coordinates,
                            const double                    threshold) -> void;

}  // namespace simdkin

#endif /* SimdKineticEnergyFunc_hpp */

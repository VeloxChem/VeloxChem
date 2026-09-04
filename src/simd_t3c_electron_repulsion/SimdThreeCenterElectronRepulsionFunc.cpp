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



#include "SimdThreeCenterElectronRepulsionFunc.hpp"

#include <string>

#include "ErrorHandler.hpp"
#include "SimdThreeCenterElectronRepulsionRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionRecDSS.hpp"
#include "SimdThreeCenterElectronRepulsionRecFSS.hpp"
#include "SimdThreeCenterElectronRepulsionRecGSS.hpp"
#include "SimdThreeCenterElectronRepulsionRecHSS.hpp"
#include "SimdThreeCenterElectronRepulsionRecISS.hpp"
#include "SimdThreeCenterElectronRepulsionRecSPS.hpp"
#include "SimdThreeCenterElectronRepulsionRecSDS.hpp"
#include "SimdThreeCenterElectronRepulsionRecSFS.hpp"
#include "SimdThreeCenterElectronRepulsionRecSGS.hpp"
#include "SimdThreeCenterElectronRepulsionRecSHS.hpp"
#include "SimdThreeCenterElectronRepulsionRecSIS.hpp"
#include "SimdThreeCenterElectronRepulsionRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionRecSSG.hpp"
#include "SimdThreeCenterElectronRepulsionRecSSH.hpp"
#include "SimdThreeCenterElectronRepulsionRecSSI.hpp"
#include "SimdThreeCenterElectronRepulsionRecSSP.hpp"
#include "SimdThreeCenterElectronRepulsionRecSSS.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_electron_repulsion(double                         *values,
                           const size_t                    npairs,
                           const size_t                    natoms,
                           const size_t                    iatom,
                           const CBasisFunction           &a_function,
                           const CBasisFunction           &b_function,
                           const CBasisFunction           &c_function,
                           const std::vector<CSimdMatrix> &ab_harmonics,
                           const std::vector<CSimdMatrix> &bc_harmonics,
                           const CSimdMatrix              &ab_coordinates,
                           const CSimdMatrix              &bc_coordinates,
                           const double                    threshold) -> void
{
    const auto la = a_function.get_angular_momentum();

    const auto lb = b_function.get_angular_momentum();

    const auto lc = c_function.get_angular_momentum();

    // NOTE: the kernel of three S type functions needs no solid harmonics, as the
    // harmonics of angular momentum zero are one for every atom pair.

    if ((la == 0) && (lb == 0) && (lc == 0))
    {
        compute_sss_electron_repulsion(
            values, npairs, natoms, iatom, a_function, b_function, c_function, ab_coordinates, bc_coordinates, threshold);

        return;
    }

    // NOTE: the kernel of two S type functions and one of angular momentum one
    // takes the solid harmonics of angular momentum one of both the atom pairs
    // and the atoms on c side, as the addition theorem which splits the vector
    // from the product center of a pair to the atom on c side carries the
    // angular momentum onto both of them.

    if ((la == 0) && (lb == 0) && (lc == 1))
    {
        compute_ssp_electron_repulsion(values,
                                       npairs,
                                       natoms,
                                       iatom,
                                       a_function,
                                       b_function,
                                       c_function,
                                       ab_harmonics[0],
                                       bc_harmonics[0],
                                       ab_coordinates,
                                       bc_coordinates,
                                       threshold);

        return;
    }

    // NOTE: the kernel of two S type functions and one of angular momentum two
    // takes the solid harmonics of angular momenta one and two of both sides, as
    // the addition theorem splits the vector to the atom on c side into three
    // bidegrees and the middle one couples the orders of the two.

    if ((la == 0) && (lb == 0) && (lc == 2))
    {
        compute_ssd_electron_repulsion(values,
                                       npairs,
                                       natoms,
                                       iatom,
                                       a_function,
                                       b_function,
                                       c_function,
                                       ab_harmonics,
                                       bc_harmonics,
                                       ab_coordinates,
                                       bc_coordinates,
                                       threshold);

        return;
    }

    if ((la == 0) && (lb == 0) && (lc == 3))
    {
        compute_ssf_electron_repulsion(values,
                                       npairs,
                                       natoms,
                                       iatom,
                                       a_function,
                                       b_function,
                                       c_function,
                                       ab_harmonics,
                                       bc_harmonics,
                                       ab_coordinates,
                                       bc_coordinates,
                                       threshold);

        return;
    }

    if ((la == 0) && (lb == 0) && (lc == 4))
    {
        compute_ssg_electron_repulsion(values,
                                       npairs,
                                       natoms,
                                       iatom,
                                       a_function,
                                       b_function,
                                       c_function,
                                       ab_harmonics,
                                       bc_harmonics,
                                       ab_coordinates,
                                       bc_coordinates,
                                       threshold);

        return;
    }

    if ((la == 0) && (lb == 0) && (lc == 5))
    {
        compute_ssh_electron_repulsion(values,
                                       npairs,
                                       natoms,
                                       iatom,
                                       a_function,
                                       b_function,
                                       c_function,
                                       ab_harmonics,
                                       bc_harmonics,
                                       ab_coordinates,
                                       bc_coordinates,
                                       threshold);

        return;
    }

    if ((la == 0) && (lb == 0) && (lc == 6))
    {
        compute_ssi_electron_repulsion(values,
                                       npairs,
                                       natoms,
                                       iatom,
                                       a_function,
                                       b_function,
                                       c_function,
                                       ab_harmonics,
                                       bc_harmonics,
                                       ab_coordinates,
                                       bc_coordinates,
                                       threshold);

        return;
    }

    // NOTE: the combinations which carry the angular momentum on the a or b side
    // expand the harmonic of the vector to the atom on c side back onto the two
    // the block already carries, so they read the same harmonics as the ones
    // above but weigh them with a polynomial in the Boys function.

    if ((la == 1) && (lb == 0) && (lc == 0))
    {
        compute_pss_electron_repulsion(values,
                                       npairs,
                                       natoms,
                                       iatom,
                                       a_function,
                                       b_function,
                                       c_function,
                                       ab_harmonics,
                                       bc_harmonics,
                                       ab_coordinates,
                                       bc_coordinates,
                                       threshold);

        return;
    }

    if ((la == 0) && (lb == 1) && (lc == 0))
    {
        compute_sps_electron_repulsion(values,
                                       npairs,
                                       natoms,
                                       iatom,
                                       a_function,
                                       b_function,
                                       c_function,
                                       ab_harmonics,
                                       bc_harmonics,
                                       ab_coordinates,
                                       bc_coordinates,
                                       threshold);

        return;
    }

    if ((la == 2) && (lb == 0) && (lc == 0))
    {
        compute_dss_electron_repulsion(values,
                                       npairs,
                                       natoms,
                                       iatom,
                                       a_function,
                                       b_function,
                                       c_function,
                                       ab_harmonics,
                                       bc_harmonics,
                                       ab_coordinates,
                                       bc_coordinates,
                                       threshold);

        return;
    }

    if ((la == 0) && (lb == 2) && (lc == 0))
    {
        compute_sds_electron_repulsion(values,
                                       npairs,
                                       natoms,
                                       iatom,
                                       a_function,
                                       b_function,
                                       c_function,
                                       ab_harmonics,
                                       bc_harmonics,
                                       ab_coordinates,
                                       bc_coordinates,
                                       threshold);

        return;
    }

    if ((la == 3) && (lb == 0) && (lc == 0))
    {
        compute_fss_electron_repulsion(values,
                                       npairs,
                                       natoms,
                                       iatom,
                                       a_function,
                                       b_function,
                                       c_function,
                                       ab_harmonics,
                                       bc_harmonics,
                                       ab_coordinates,
                                       bc_coordinates,
                                       threshold);

        return;
    }

    if ((la == 0) && (lb == 3) && (lc == 0))
    {
        compute_sfs_electron_repulsion(values,
                                       npairs,
                                       natoms,
                                       iatom,
                                       a_function,
                                       b_function,
                                       c_function,
                                       ab_harmonics,
                                       bc_harmonics,
                                       ab_coordinates,
                                       bc_coordinates,
                                       threshold);

        return;
    }

    if ((la == 4) && (lb == 0) && (lc == 0))
    {
        compute_gss_electron_repulsion(values,
                                       npairs,
                                       natoms,
                                       iatom,
                                       a_function,
                                       b_function,
                                       c_function,
                                       ab_harmonics,
                                       bc_harmonics,
                                       ab_coordinates,
                                       bc_coordinates,
                                       threshold);

        return;
    }

    if ((la == 0) && (lb == 4) && (lc == 0))
    {
        compute_sgs_electron_repulsion(values,
                                       npairs,
                                       natoms,
                                       iatom,
                                       a_function,
                                       b_function,
                                       c_function,
                                       ab_harmonics,
                                       bc_harmonics,
                                       ab_coordinates,
                                       bc_coordinates,
                                       threshold);

        return;
    }

    if ((la == 5) && (lb == 0) && (lc == 0))
    {
        compute_hss_electron_repulsion(values,
                                       npairs,
                                       natoms,
                                       iatom,
                                       a_function,
                                       b_function,
                                       c_function,
                                       ab_harmonics,
                                       bc_harmonics,
                                       ab_coordinates,
                                       bc_coordinates,
                                       threshold);

        return;
    }

    if ((la == 0) && (lb == 5) && (lc == 0))
    {
        compute_shs_electron_repulsion(values,
                                       npairs,
                                       natoms,
                                       iatom,
                                       a_function,
                                       b_function,
                                       c_function,
                                       ab_harmonics,
                                       bc_harmonics,
                                       ab_coordinates,
                                       bc_coordinates,
                                       threshold);

        return;
    }

    if ((la == 6) && (lb == 0) && (lc == 0))
    {
        compute_iss_electron_repulsion(values,
                                       npairs,
                                       natoms,
                                       iatom,
                                       a_function,
                                       b_function,
                                       c_function,
                                       ab_harmonics,
                                       bc_harmonics,
                                       ab_coordinates,
                                       bc_coordinates,
                                       threshold);

        return;
    }

    if ((la == 0) && (lb == 6) && (lc == 0))
    {
        compute_sis_electron_repulsion(values,
                                       npairs,
                                       natoms,
                                       iatom,
                                       a_function,
                                       b_function,
                                       c_function,
                                       ab_harmonics,
                                       bc_harmonics,
                                       ab_coordinates,
                                       bc_coordinates,
                                       threshold);

        return;
    }

    // NOTE: the remaining combinations of angular momenta have no kernel yet. The
    // driver aborts on them rather than leaving their values undefined or filling
    // them with something which is not the integral: a tensor whose blocks were
    // partly computed and partly not would carry no sign of which of the two any
    // one of its values is.

    errors::assertMsgCritical(false,
                              std::string("SimdThreeCenterElectronRepulsionFunc.compute_electron_repulsion: Combination of angular momenta ") +
                                  std::to_string(la) + std::string(", ") + std::to_string(lb) + std::string(" and ") + std::to_string(lc) +
                                  std::string(" is not implemented"));
}

}  // namespace simdt3ceri

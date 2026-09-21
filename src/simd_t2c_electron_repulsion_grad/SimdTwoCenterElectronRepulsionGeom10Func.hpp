//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//

#ifndef SimdTwoCenterElectronRepulsionGeom10Func_hpp
#define SimdTwoCenterElectronRepulsionGeom10Func_hpp

#include <cstddef>
#include <string>

#include "BasisFunction.hpp"
#include "ErrorHandler.hpp"
#include "SimdMatrix.hpp"

#include "SimdElectronRepulsionGeom10RecSS.hpp"
#include "SimdElectronRepulsionGeom10RecSP.hpp"
#include "SimdElectronRepulsionGeom10RecSD.hpp"
#include "SimdElectronRepulsionGeom10RecSF.hpp"
#include "SimdElectronRepulsionGeom10RecSG.hpp"
#include "SimdElectronRepulsionGeom10RecSH.hpp"
#include "SimdElectronRepulsionGeom10RecSI.hpp"
#include "SimdElectronRepulsionGeom10RecSK.hpp"
#include "SimdElectronRepulsionGeom10RecSL.hpp"
#include "SimdElectronRepulsionGeom10RecPS.hpp"
#include "SimdElectronRepulsionGeom10RecPP.hpp"
#include "SimdElectronRepulsionGeom10RecPD.hpp"
#include "SimdElectronRepulsionGeom10RecPF.hpp"
#include "SimdElectronRepulsionGeom10RecPG.hpp"
#include "SimdElectronRepulsionGeom10RecPH.hpp"
#include "SimdElectronRepulsionGeom10RecPI.hpp"
#include "SimdElectronRepulsionGeom10RecPK.hpp"
#include "SimdElectronRepulsionGeom10RecPL.hpp"
#include "SimdElectronRepulsionGeom10RecDS.hpp"
#include "SimdElectronRepulsionGeom10RecDP.hpp"
#include "SimdElectronRepulsionGeom10RecDD.hpp"
#include "SimdElectronRepulsionGeom10RecDF.hpp"
#include "SimdElectronRepulsionGeom10RecDG.hpp"
#include "SimdElectronRepulsionGeom10RecDH.hpp"
#include "SimdElectronRepulsionGeom10RecDI.hpp"
#include "SimdElectronRepulsionGeom10RecDK.hpp"
#include "SimdElectronRepulsionGeom10RecDL.hpp"
#include "SimdElectronRepulsionGeom10RecFS.hpp"
#include "SimdElectronRepulsionGeom10RecFP.hpp"
#include "SimdElectronRepulsionGeom10RecFD.hpp"
#include "SimdElectronRepulsionGeom10RecFF.hpp"
#include "SimdElectronRepulsionGeom10RecFG.hpp"
#include "SimdElectronRepulsionGeom10RecFH.hpp"
#include "SimdElectronRepulsionGeom10RecFI.hpp"
#include "SimdElectronRepulsionGeom10RecFK.hpp"
#include "SimdElectronRepulsionGeom10RecFL.hpp"
#include "SimdElectronRepulsionGeom10RecGS.hpp"
#include "SimdElectronRepulsionGeom10RecGP.hpp"
#include "SimdElectronRepulsionGeom10RecGD.hpp"
#include "SimdElectronRepulsionGeom10RecGF.hpp"
#include "SimdElectronRepulsionGeom10RecGG.hpp"
#include "SimdElectronRepulsionGeom10RecGH.hpp"
#include "SimdElectronRepulsionGeom10RecGI.hpp"
#include "SimdElectronRepulsionGeom10RecGK.hpp"
#include "SimdElectronRepulsionGeom10RecGL.hpp"
#include "SimdElectronRepulsionGeom10RecHS.hpp"
#include "SimdElectronRepulsionGeom10RecHP.hpp"
#include "SimdElectronRepulsionGeom10RecHD.hpp"
#include "SimdElectronRepulsionGeom10RecHF.hpp"
#include "SimdElectronRepulsionGeom10RecHG.hpp"
#include "SimdElectronRepulsionGeom10RecHH.hpp"
#include "SimdElectronRepulsionGeom10RecHI.hpp"
#include "SimdElectronRepulsionGeom10RecHK.hpp"
#include "SimdElectronRepulsionGeom10RecHL.hpp"
#include "SimdElectronRepulsionGeom10RecIS.hpp"
#include "SimdElectronRepulsionGeom10RecIP.hpp"
#include "SimdElectronRepulsionGeom10RecID.hpp"
#include "SimdElectronRepulsionGeom10RecIF.hpp"
#include "SimdElectronRepulsionGeom10RecIG.hpp"
#include "SimdElectronRepulsionGeom10RecIH.hpp"
#include "SimdElectronRepulsionGeom10RecII.hpp"
#include "SimdElectronRepulsionGeom10RecIK.hpp"
#include "SimdElectronRepulsionGeom10RecIL.hpp"
#include "SimdElectronRepulsionGeom10RecKS.hpp"
#include "SimdElectronRepulsionGeom10RecKP.hpp"
#include "SimdElectronRepulsionGeom10RecKD.hpp"
#include "SimdElectronRepulsionGeom10RecKF.hpp"
#include "SimdElectronRepulsionGeom10RecKG.hpp"
#include "SimdElectronRepulsionGeom10RecKH.hpp"
#include "SimdElectronRepulsionGeom10RecKI.hpp"
#include "SimdElectronRepulsionGeom10RecKK.hpp"
#include "SimdElectronRepulsionGeom10RecKL.hpp"
#include "SimdElectronRepulsionGeom10RecLS.hpp"
#include "SimdElectronRepulsionGeom10RecLP.hpp"
#include "SimdElectronRepulsionGeom10RecLD.hpp"
#include "SimdElectronRepulsionGeom10RecLF.hpp"
#include "SimdElectronRepulsionGeom10RecLG.hpp"
#include "SimdElectronRepulsionGeom10RecLH.hpp"
#include "SimdElectronRepulsionGeom10RecLI.hpp"
#include "SimdElectronRepulsionGeom10RecLK.hpp"
#include "SimdElectronRepulsionGeom10RecLL.hpp"

namespace simdt2cerigrad {  // simdt2cerigrad namespace

/// @brief Computes the derivative of the two-center electron repulsion integrals
/// with respect to the position of the atom on bra side.
/// @param values The values of the derivative, three components of one row per
/// pair of angular components, each of one value per atom pair.
/// @param nvalues The number of atom pairs.
/// @param bra The basis function on bra side.
/// @param ket The basis function on ket side.
/// @param coordinates The coordinates of the atom pairs.
/// @param buffer The scratch of the combination of basis functions.
/// @note The derivative of the atom on ket side is not computed. The two centers
/// of a two-center integral move against each other, so their derivatives sum to
/// zero and the caller takes the second as the negative of the first.
/// @note A combination of angular momenta which has no kernel stops rather than
/// returning zeros: a gradient of zeros is a gradient which looks converged
/// everywhere, and nothing above would notice.
inline auto
compute_electron_repulsion_geom_10(double               *values,
                                   const size_t          nvalues,
                                   const CBasisFunction &bra,
                                   const CBasisFunction &ket,
                                   const CSimdMatrix    &coordinates,
                                   CSimdMatrix          &buffer) -> void
{
    const auto lbra = bra.get_angular_momentum();

    const auto lket = ket.get_angular_momentum();

    errors::assertMsgCritical((lbra >= 0) && (lbra < 9) && (lket >= 0) && (lket < 9),
                              std::string("SimdTwoCenterElectronRepulsionGeom10Func: Angular momentum is out of range"));

    switch (lbra * 9 + lket)
    {
        case 0:
            simdt2ceri::compute_geom_10_ss_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 1:
            simdt2ceri::compute_geom_10_sp_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 2:
            simdt2ceri::compute_geom_10_sd_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 3:
            simdt2ceri::compute_geom_10_sf_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 4:
            simdt2ceri::compute_geom_10_sg_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 5:
            simdt2ceri::compute_geom_10_sh_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 6:
            simdt2ceri::compute_geom_10_si_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 7:
            simdt2ceri::compute_geom_10_sk_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 8:
            simdt2ceri::compute_geom_10_sl_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 9:
            simdt2ceri::compute_geom_10_ps_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 10:
            simdt2ceri::compute_geom_10_pp_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 11:
            simdt2ceri::compute_geom_10_pd_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 12:
            simdt2ceri::compute_geom_10_pf_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 13:
            simdt2ceri::compute_geom_10_pg_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 14:
            simdt2ceri::compute_geom_10_ph_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 15:
            simdt2ceri::compute_geom_10_pi_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 16:
            simdt2ceri::compute_geom_10_pk_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 17:
            simdt2ceri::compute_geom_10_pl_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 18:
            simdt2ceri::compute_geom_10_ds_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 19:
            simdt2ceri::compute_geom_10_dp_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 20:
            simdt2ceri::compute_geom_10_dd_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 21:
            simdt2ceri::compute_geom_10_df_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 22:
            simdt2ceri::compute_geom_10_dg_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 23:
            simdt2ceri::compute_geom_10_dh_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 24:
            simdt2ceri::compute_geom_10_di_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 25:
            simdt2ceri::compute_geom_10_dk_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 26:
            simdt2ceri::compute_geom_10_dl_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 27:
            simdt2ceri::compute_geom_10_fs_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 28:
            simdt2ceri::compute_geom_10_fp_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 29:
            simdt2ceri::compute_geom_10_fd_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 30:
            simdt2ceri::compute_geom_10_ff_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 31:
            simdt2ceri::compute_geom_10_fg_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 32:
            simdt2ceri::compute_geom_10_fh_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 33:
            simdt2ceri::compute_geom_10_fi_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 34:
            simdt2ceri::compute_geom_10_fk_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 35:
            simdt2ceri::compute_geom_10_fl_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 36:
            simdt2ceri::compute_geom_10_gs_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 37:
            simdt2ceri::compute_geom_10_gp_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 38:
            simdt2ceri::compute_geom_10_gd_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 39:
            simdt2ceri::compute_geom_10_gf_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 40:
            simdt2ceri::compute_geom_10_gg_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 41:
            simdt2ceri::compute_geom_10_gh_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 42:
            simdt2ceri::compute_geom_10_gi_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 43:
            simdt2ceri::compute_geom_10_gk_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 44:
            simdt2ceri::compute_geom_10_gl_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 45:
            simdt2ceri::compute_geom_10_hs_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 46:
            simdt2ceri::compute_geom_10_hp_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 47:
            simdt2ceri::compute_geom_10_hd_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 48:
            simdt2ceri::compute_geom_10_hf_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 49:
            simdt2ceri::compute_geom_10_hg_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 50:
            simdt2ceri::compute_geom_10_hh_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 51:
            simdt2ceri::compute_geom_10_hi_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 52:
            simdt2ceri::compute_geom_10_hk_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 53:
            simdt2ceri::compute_geom_10_hl_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 54:
            simdt2ceri::compute_geom_10_is_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 55:
            simdt2ceri::compute_geom_10_ip_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 56:
            simdt2ceri::compute_geom_10_id_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 57:
            simdt2ceri::compute_geom_10_if_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 58:
            simdt2ceri::compute_geom_10_ig_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 59:
            simdt2ceri::compute_geom_10_ih_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 60:
            simdt2ceri::compute_geom_10_ii_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 61:
            simdt2ceri::compute_geom_10_ik_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 62:
            simdt2ceri::compute_geom_10_il_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 63:
            simdt2ceri::compute_geom_10_ks_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 64:
            simdt2ceri::compute_geom_10_kp_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 65:
            simdt2ceri::compute_geom_10_kd_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 66:
            simdt2ceri::compute_geom_10_kf_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 67:
            simdt2ceri::compute_geom_10_kg_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 68:
            simdt2ceri::compute_geom_10_kh_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 69:
            simdt2ceri::compute_geom_10_ki_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 70:
            simdt2ceri::compute_geom_10_kk_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 71:
            simdt2ceri::compute_geom_10_kl_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 72:
            simdt2ceri::compute_geom_10_ls_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 73:
            simdt2ceri::compute_geom_10_lp_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 74:
            simdt2ceri::compute_geom_10_ld_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 75:
            simdt2ceri::compute_geom_10_lf_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 76:
            simdt2ceri::compute_geom_10_lg_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 77:
            simdt2ceri::compute_geom_10_lh_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 78:
            simdt2ceri::compute_geom_10_li_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 79:
            simdt2ceri::compute_geom_10_lk_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        case 80:
            simdt2ceri::compute_geom_10_ll_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer);
            break;
        default:
            errors::assertMsgCritical(
                false, std::string("SimdTwoCenterElectronRepulsionGeom10Func: No kernel for the combination of angular momenta"));
    }
}

}  // namespace simdt2cerigrad

#endif /* SimdTwoCenterElectronRepulsionGeom10Func_hpp */

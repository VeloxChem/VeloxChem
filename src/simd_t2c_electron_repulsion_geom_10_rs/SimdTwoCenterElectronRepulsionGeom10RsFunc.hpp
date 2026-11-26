#ifndef SimdTwoCenterElectronRepulsionGeom10RsFunc_hpp
#define SimdTwoCenterElectronRepulsionGeom10RsFunc_hpp

#include <cstddef>
#include <string>

#include "BasisFunction.hpp"
#include "ErrorHandler.hpp"
#include "SimdMatrix.hpp"

#include "SimdElectronRepulsionGeom10RsRecSS.hpp"
#include "SimdElectronRepulsionGeom10RsRecSP.hpp"
#include "SimdElectronRepulsionGeom10RsRecSD.hpp"
#include "SimdElectronRepulsionGeom10RsRecSF.hpp"
#include "SimdElectronRepulsionGeom10RsRecSG.hpp"
#include "SimdElectronRepulsionGeom10RsRecSH.hpp"
#include "SimdElectronRepulsionGeom10RsRecSI.hpp"
#include "SimdElectronRepulsionGeom10RsRecSK.hpp"
#include "SimdElectronRepulsionGeom10RsRecSL.hpp"
#include "SimdElectronRepulsionGeom10RsRecPS.hpp"
#include "SimdElectronRepulsionGeom10RsRecPP.hpp"
#include "SimdElectronRepulsionGeom10RsRecPD.hpp"
#include "SimdElectronRepulsionGeom10RsRecPF.hpp"
#include "SimdElectronRepulsionGeom10RsRecPG.hpp"
#include "SimdElectronRepulsionGeom10RsRecPH.hpp"
#include "SimdElectronRepulsionGeom10RsRecPI.hpp"
#include "SimdElectronRepulsionGeom10RsRecPK.hpp"
#include "SimdElectronRepulsionGeom10RsRecPL.hpp"
#include "SimdElectronRepulsionGeom10RsRecDS.hpp"
#include "SimdElectronRepulsionGeom10RsRecDP.hpp"
#include "SimdElectronRepulsionGeom10RsRecDD.hpp"
#include "SimdElectronRepulsionGeom10RsRecDF.hpp"
#include "SimdElectronRepulsionGeom10RsRecDG.hpp"
#include "SimdElectronRepulsionGeom10RsRecDH.hpp"
#include "SimdElectronRepulsionGeom10RsRecDI.hpp"
#include "SimdElectronRepulsionGeom10RsRecDK.hpp"
#include "SimdElectronRepulsionGeom10RsRecDL.hpp"
#include "SimdElectronRepulsionGeom10RsRecFS.hpp"
#include "SimdElectronRepulsionGeom10RsRecFP.hpp"
#include "SimdElectronRepulsionGeom10RsRecFD.hpp"
#include "SimdElectronRepulsionGeom10RsRecFF.hpp"
#include "SimdElectronRepulsionGeom10RsRecFG.hpp"
#include "SimdElectronRepulsionGeom10RsRecFH.hpp"
#include "SimdElectronRepulsionGeom10RsRecFI.hpp"
#include "SimdElectronRepulsionGeom10RsRecFK.hpp"
#include "SimdElectronRepulsionGeom10RsRecFL.hpp"
#include "SimdElectronRepulsionGeom10RsRecGS.hpp"
#include "SimdElectronRepulsionGeom10RsRecGP.hpp"
#include "SimdElectronRepulsionGeom10RsRecGD.hpp"
#include "SimdElectronRepulsionGeom10RsRecGF.hpp"
#include "SimdElectronRepulsionGeom10RsRecGG.hpp"
#include "SimdElectronRepulsionGeom10RsRecGH.hpp"
#include "SimdElectronRepulsionGeom10RsRecGI.hpp"
#include "SimdElectronRepulsionGeom10RsRecGK.hpp"
#include "SimdElectronRepulsionGeom10RsRecGL.hpp"
#include "SimdElectronRepulsionGeom10RsRecHS.hpp"
#include "SimdElectronRepulsionGeom10RsRecHP.hpp"
#include "SimdElectronRepulsionGeom10RsRecHD.hpp"
#include "SimdElectronRepulsionGeom10RsRecHF.hpp"
#include "SimdElectronRepulsionGeom10RsRecHG.hpp"
#include "SimdElectronRepulsionGeom10RsRecHH.hpp"
#include "SimdElectronRepulsionGeom10RsRecHI.hpp"
#include "SimdElectronRepulsionGeom10RsRecHK.hpp"
#include "SimdElectronRepulsionGeom10RsRecHL.hpp"
#include "SimdElectronRepulsionGeom10RsRecIS.hpp"
#include "SimdElectronRepulsionGeom10RsRecIP.hpp"
#include "SimdElectronRepulsionGeom10RsRecID.hpp"
#include "SimdElectronRepulsionGeom10RsRecIF.hpp"
#include "SimdElectronRepulsionGeom10RsRecIG.hpp"
#include "SimdElectronRepulsionGeom10RsRecIH.hpp"
#include "SimdElectronRepulsionGeom10RsRecII.hpp"
#include "SimdElectronRepulsionGeom10RsRecIK.hpp"
#include "SimdElectronRepulsionGeom10RsRecIL.hpp"
#include "SimdElectronRepulsionGeom10RsRecKS.hpp"
#include "SimdElectronRepulsionGeom10RsRecKP.hpp"
#include "SimdElectronRepulsionGeom10RsRecKD.hpp"
#include "SimdElectronRepulsionGeom10RsRecKF.hpp"
#include "SimdElectronRepulsionGeom10RsRecKG.hpp"
#include "SimdElectronRepulsionGeom10RsRecKH.hpp"
#include "SimdElectronRepulsionGeom10RsRecKI.hpp"
#include "SimdElectronRepulsionGeom10RsRecKK.hpp"
#include "SimdElectronRepulsionGeom10RsRecKL.hpp"
#include "SimdElectronRepulsionGeom10RsRecLS.hpp"
#include "SimdElectronRepulsionGeom10RsRecLP.hpp"
#include "SimdElectronRepulsionGeom10RsRecLD.hpp"
#include "SimdElectronRepulsionGeom10RsRecLF.hpp"
#include "SimdElectronRepulsionGeom10RsRecLG.hpp"
#include "SimdElectronRepulsionGeom10RsRecLH.hpp"
#include "SimdElectronRepulsionGeom10RsRecLI.hpp"
#include "SimdElectronRepulsionGeom10RsRecLK.hpp"
#include "SimdElectronRepulsionGeom10RsRecLL.hpp"

namespace simdt2cerigrad {  // simdt2cerigrad namespace

/// @brief Computes the derivative with respect to the center of the bra of the
/// two-center integrals of the Coulomb operator and of the attenuated one, for one
/// combination of angular momenta.
/// @param values The buffer the values are written into, which holds **six** blocks
/// of nvalues: three Cartesian components of one operator and then three of the
/// other.
/// @param nvalues The number of values of one block.
/// @param bra The basis function on the bra side.
/// @param ket The basis function on the ket side.
/// @param coordinates The coordinates of the pairs of centers.
/// @param buffer The scratch the recursion is carried out in.
/// @param omega The range separation parameter.
/// @note The momenta are checked rather than trusted. A combination outside the
/// range would otherwise fall through the switch and leave the buffer as it was,
/// and a derivative of zeros is a gradient which looks converged everywhere.
inline auto
compute_rs_electron_repulsion_geom_10(double               *values,
                                      const size_t          nvalues,
                                      const CBasisFunction &bra,
                                      const CBasisFunction &ket,
                                      const CSimdMatrix    &coordinates,
                                      CSimdMatrix          &buffer,
                                      const double          omega) -> void
{
    const auto lbra = bra.get_angular_momentum();

    const auto lket = ket.get_angular_momentum();

    errors::assertMsgCritical(
        (lbra >= 0) && (lbra < 9) && (lket >= 0) && (lket < 9),
        std::string("SimdTwoCenterElectronRepulsionGeom10RsFunc: Angular momentum is out of range"));

    switch (lbra * 9 + lket)
    {
        case 0:
            simdt2ceri::compute_rs_geom_10_ss_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 1:
            simdt2ceri::compute_rs_geom_10_sp_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 2:
            simdt2ceri::compute_rs_geom_10_sd_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 3:
            simdt2ceri::compute_rs_geom_10_sf_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 4:
            simdt2ceri::compute_rs_geom_10_sg_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 5:
            simdt2ceri::compute_rs_geom_10_sh_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 6:
            simdt2ceri::compute_rs_geom_10_si_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 7:
            simdt2ceri::compute_rs_geom_10_sk_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 8:
            simdt2ceri::compute_rs_geom_10_sl_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 9:
            simdt2ceri::compute_rs_geom_10_ps_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 10:
            simdt2ceri::compute_rs_geom_10_pp_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 11:
            simdt2ceri::compute_rs_geom_10_pd_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 12:
            simdt2ceri::compute_rs_geom_10_pf_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 13:
            simdt2ceri::compute_rs_geom_10_pg_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 14:
            simdt2ceri::compute_rs_geom_10_ph_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 15:
            simdt2ceri::compute_rs_geom_10_pi_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 16:
            simdt2ceri::compute_rs_geom_10_pk_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 17:
            simdt2ceri::compute_rs_geom_10_pl_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 18:
            simdt2ceri::compute_rs_geom_10_ds_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 19:
            simdt2ceri::compute_rs_geom_10_dp_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 20:
            simdt2ceri::compute_rs_geom_10_dd_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 21:
            simdt2ceri::compute_rs_geom_10_df_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 22:
            simdt2ceri::compute_rs_geom_10_dg_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 23:
            simdt2ceri::compute_rs_geom_10_dh_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 24:
            simdt2ceri::compute_rs_geom_10_di_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 25:
            simdt2ceri::compute_rs_geom_10_dk_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 26:
            simdt2ceri::compute_rs_geom_10_dl_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 27:
            simdt2ceri::compute_rs_geom_10_fs_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 28:
            simdt2ceri::compute_rs_geom_10_fp_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 29:
            simdt2ceri::compute_rs_geom_10_fd_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 30:
            simdt2ceri::compute_rs_geom_10_ff_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 31:
            simdt2ceri::compute_rs_geom_10_fg_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 32:
            simdt2ceri::compute_rs_geom_10_fh_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 33:
            simdt2ceri::compute_rs_geom_10_fi_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 34:
            simdt2ceri::compute_rs_geom_10_fk_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 35:
            simdt2ceri::compute_rs_geom_10_fl_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 36:
            simdt2ceri::compute_rs_geom_10_gs_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 37:
            simdt2ceri::compute_rs_geom_10_gp_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 38:
            simdt2ceri::compute_rs_geom_10_gd_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 39:
            simdt2ceri::compute_rs_geom_10_gf_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 40:
            simdt2ceri::compute_rs_geom_10_gg_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 41:
            simdt2ceri::compute_rs_geom_10_gh_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 42:
            simdt2ceri::compute_rs_geom_10_gi_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 43:
            simdt2ceri::compute_rs_geom_10_gk_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 44:
            simdt2ceri::compute_rs_geom_10_gl_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 45:
            simdt2ceri::compute_rs_geom_10_hs_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 46:
            simdt2ceri::compute_rs_geom_10_hp_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 47:
            simdt2ceri::compute_rs_geom_10_hd_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 48:
            simdt2ceri::compute_rs_geom_10_hf_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 49:
            simdt2ceri::compute_rs_geom_10_hg_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 50:
            simdt2ceri::compute_rs_geom_10_hh_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 51:
            simdt2ceri::compute_rs_geom_10_hi_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 52:
            simdt2ceri::compute_rs_geom_10_hk_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 53:
            simdt2ceri::compute_rs_geom_10_hl_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 54:
            simdt2ceri::compute_rs_geom_10_is_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 55:
            simdt2ceri::compute_rs_geom_10_ip_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 56:
            simdt2ceri::compute_rs_geom_10_id_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 57:
            simdt2ceri::compute_rs_geom_10_if_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 58:
            simdt2ceri::compute_rs_geom_10_ig_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 59:
            simdt2ceri::compute_rs_geom_10_ih_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 60:
            simdt2ceri::compute_rs_geom_10_ii_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 61:
            simdt2ceri::compute_rs_geom_10_ik_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 62:
            simdt2ceri::compute_rs_geom_10_il_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 63:
            simdt2ceri::compute_rs_geom_10_ks_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 64:
            simdt2ceri::compute_rs_geom_10_kp_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 65:
            simdt2ceri::compute_rs_geom_10_kd_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 66:
            simdt2ceri::compute_rs_geom_10_kf_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 67:
            simdt2ceri::compute_rs_geom_10_kg_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 68:
            simdt2ceri::compute_rs_geom_10_kh_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 69:
            simdt2ceri::compute_rs_geom_10_ki_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 70:
            simdt2ceri::compute_rs_geom_10_kk_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 71:
            simdt2ceri::compute_rs_geom_10_kl_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 72:
            simdt2ceri::compute_rs_geom_10_ls_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 73:
            simdt2ceri::compute_rs_geom_10_lp_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 74:
            simdt2ceri::compute_rs_geom_10_ld_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 75:
            simdt2ceri::compute_rs_geom_10_lf_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 76:
            simdt2ceri::compute_rs_geom_10_lg_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 77:
            simdt2ceri::compute_rs_geom_10_lh_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 78:
            simdt2ceri::compute_rs_geom_10_li_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 79:
            simdt2ceri::compute_rs_geom_10_lk_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        case 80:
            simdt2ceri::compute_rs_geom_10_ll_electron_repulsion(
                values, nvalues, bra, ket, coordinates, buffer, omega);
            break;
        default:
            break;
    }
}

}  // namespace simdt2cerigrad

#endif /* SimdTwoCenterElectronRepulsionGeom10RsFunc_hpp */

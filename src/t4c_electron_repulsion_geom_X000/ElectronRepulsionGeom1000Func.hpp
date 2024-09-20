#ifndef ElectronRepulsionGeom1000Func_hpp
#define ElectronRepulsionGeom1000Func_hpp

#include <cstddef>
#include <iostream>
#include <utility>

#include "ElectronRepulsionRecSSSS.hpp"
#include "GtoPairBlock.hpp"

namespace erifunc {  // erifunc namespace

/// Computes vector of integrals of requested four center integral.
/// @param distributor  The pointer to distributor of integrals.
/// @param bra_gto_pair_block The basis function pairs block on bra side.
/// @param ket_gto_pair_block The basis function pairs block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
/// @param bra_eq_ket True if basis function pairs blocks on bra and ket are the same, False otherwise.
template <class T>
inline auto
compute_geom_1000(T&                               distributor,
                  const CGtoPairBlock&             bra_gto_pair_block,
                  const CGtoPairBlock&             ket_gto_pair_block,
                  const std::pair<size_t, size_t>& bra_indices,
                  const std::pair<size_t, size_t>& ket_indices -> void
{
    const auto bra_angmoms = bra_gto_pair_block.angular_momentums();

    const auto ket_angmoms = ket_gto_pair_block.angular_momentums();

    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_ge_ssss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_sssp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_sssd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_sssf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({0, 4})))
    {
        erirec::comp_electron_repulsion_sssg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_sspp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_sspd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_sspf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({1, 4})))
    {
        erirec::comp_electron_repulsion_sspg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_ssdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_ssdf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({2, 4})))
    {
        erirec::comp_electron_repulsion_ssdg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_ssff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({3, 4})))
    {
        erirec::comp_electron_repulsion_ssfg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({4, 4})))
    {
        erirec::comp_electron_repulsion_ssgg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_spss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_spsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_spsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_spsf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({0, 4})))
    {
        erirec::comp_electron_repulsion_spsg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_sppp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_sppd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_sppf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({1, 4})))
    {
        erirec::comp_electron_repulsion_sppg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_spdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_spdf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({2, 4})))
    {
        erirec::comp_electron_repulsion_spdg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_spff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({3, 4})))
    {
        erirec::comp_electron_repulsion_spfg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({4, 4})))
    {
        erirec::comp_electron_repulsion_spgg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_sdss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_sdsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_sdsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_sdsf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({0, 4})))
    {
        erirec::comp_electron_repulsion_sdsg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_sdpp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_sdpd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_sdpf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({1, 4})))
    {
        erirec::comp_electron_repulsion_sdpg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_sddd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_sddf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({2, 4})))
    {
        erirec::comp_electron_repulsion_sddg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_sdff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({3, 4})))
    {
        erirec::comp_electron_repulsion_sdfg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({4, 4})))
    {
        erirec::comp_electron_repulsion_sdgg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_sfss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_sfsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_sfsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_sfsf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({0, 4})))
    {
        erirec::comp_electron_repulsion_sfsg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_sfpp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_sfpd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_sfpf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({1, 4})))
    {
        erirec::comp_electron_repulsion_sfpg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_sfdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_sfdf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({2, 4})))
    {
        erirec::comp_electron_repulsion_sfdg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_sfff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({3, 4})))
    {
        erirec::comp_electron_repulsion_sffg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({4, 4})))
    {
        erirec::comp_electron_repulsion_sfgg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 4})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_sgss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 4})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_sgsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 4})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_sgsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 4})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_sgsf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 4})) && (ket_angmoms == std::pair<int, int>({0, 4})))
    {
        erirec::comp_electron_repulsion_sgsg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 4})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_sgpp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 4})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_sgpd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 4})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_sgpf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 4})) && (ket_angmoms == std::pair<int, int>({1, 4})))
    {
        erirec::comp_electron_repulsion_sgpg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 4})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_sgdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 4})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_sgdf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 4})) && (ket_angmoms == std::pair<int, int>({2, 4})))
    {
        erirec::comp_electron_repulsion_sgdg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 4})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_sgff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 4})) && (ket_angmoms == std::pair<int, int>({3, 4})))
    {
        erirec::comp_electron_repulsion_sgfg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 4})) && (ket_angmoms == std::pair<int, int>({4, 4})))
    {
        erirec::comp_electron_repulsion_sggg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_ppss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_ppsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_ppsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_ppsf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({0, 4})))
    {
        erirec::comp_electron_repulsion_ppsg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_pppp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_pppd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_pppf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({1, 4})))
    {
        erirec::comp_electron_repulsion_pppg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_ppdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_ppdf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({2, 4})))
    {
        erirec::comp_electron_repulsion_ppdg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_ppff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({3, 4})))
    {
        erirec::comp_electron_repulsion_ppfg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({4, 4})))
    {
        erirec::comp_electron_repulsion_ppgg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({1, 2})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_pdss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({1, 2})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_pdsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({1, 2})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_pdsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({1, 2})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_pdsf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({1, 2})) && (ket_angmoms == std::pair<int, int>({0, 4})))
    {
        erirec::comp_electron_repulsion_pdsg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({1, 2})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_pdpp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({1, 2})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_pdpd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({1, 2})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_pdpf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({1, 2})) && (ket_angmoms == std::pair<int, int>({1, 4})))
    {
        erirec::comp_electron_repulsion_pdpg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({1, 2})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_pddd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({1, 2})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_pddf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({1, 2})) && (ket_angmoms == std::pair<int, int>({2, 4})))
    {
        erirec::comp_electron_repulsion_pddg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({1, 2})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_pdff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({1, 2})) && (ket_angmoms == std::pair<int, int>({3, 4})))
    {
        erirec::comp_electron_repulsion_pdfg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({1, 2})) && (ket_angmoms == std::pair<int, int>({4, 4})))
    {
        erirec::comp_electron_repulsion_pdgg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({1, 3})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_pfss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({1, 3})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_pfsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({1, 3})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_pfsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({1, 3})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_pfsf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({1, 3})) && (ket_angmoms == std::pair<int, int>({0, 4})))
    {
        erirec::comp_electron_repulsion_pfsg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({1, 3})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_pfpp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({1, 3})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_pfpd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({1, 3})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_pfpf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({1, 3})) && (ket_angmoms == std::pair<int, int>({1, 4})))
    {
        erirec::comp_electron_repulsion_pfpg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({1, 3})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_pfdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({1, 3})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_pfdf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({1, 3})) && (ket_angmoms == std::pair<int, int>({2, 4})))
    {
        erirec::comp_electron_repulsion_pfdg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({1, 3})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_pfff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({1, 3})) && (ket_angmoms == std::pair<int, int>({3, 4})))
    {
        erirec::comp_electron_repulsion_pffg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({1, 3})) && (ket_angmoms == std::pair<int, int>({4, 4})))
    {
        erirec::comp_electron_repulsion_pfgg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({1, 4})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_pgss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({1, 4})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_pgsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({1, 4})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_pgsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({1, 4})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_pgsf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({1, 4})) && (ket_angmoms == std::pair<int, int>({0, 4})))
    {
        erirec::comp_electron_repulsion_pgsg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({1, 4})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_pgpp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({1, 4})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_pgpd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({1, 4})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_pgpf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({1, 4})) && (ket_angmoms == std::pair<int, int>({1, 4})))
    {
        erirec::comp_electron_repulsion_pgpg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({1, 4})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_pgdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({1, 4})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_pgdf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({1, 4})) && (ket_angmoms == std::pair<int, int>({2, 4})))
    {
        erirec::comp_electron_repulsion_pgdg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({1, 4})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_pgff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({1, 4})) && (ket_angmoms == std::pair<int, int>({3, 4})))
    {
        erirec::comp_electron_repulsion_pgfg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({1, 4})) && (ket_angmoms == std::pair<int, int>({4, 4})))
    {
        erirec::comp_electron_repulsion_pggg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({2, 2})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_ddss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({2, 2})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_ddsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({2, 2})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_ddsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({2, 2})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_ddsf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({2, 2})) && (ket_angmoms == std::pair<int, int>({0, 4})))
    {
        erirec::comp_electron_repulsion_ddsg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({2, 2})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_ddpp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({2, 2})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_ddpd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({2, 2})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_ddpf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({2, 2})) && (ket_angmoms == std::pair<int, int>({1, 4})))
    {
        erirec::comp_electron_repulsion_ddpg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({2, 2})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_dddd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({2, 2})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_dddf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({2, 2})) && (ket_angmoms == std::pair<int, int>({2, 4})))
    {
        erirec::comp_electron_repulsion_dddg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({2, 2})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_ddff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({2, 2})) && (ket_angmoms == std::pair<int, int>({3, 4})))
    {
        erirec::comp_electron_repulsion_ddfg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({2, 2})) && (ket_angmoms == std::pair<int, int>({4, 4})))
    {
        erirec::comp_electron_repulsion_ddgg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({2, 3})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_dfss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({2, 3})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_dfsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({2, 3})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_dfsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({2, 3})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_dfsf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({2, 3})) && (ket_angmoms == std::pair<int, int>({0, 4})))
    {
        erirec::comp_electron_repulsion_dfsg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({2, 3})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_dfpp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({2, 3})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_dfpd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({2, 3})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_dfpf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({2, 3})) && (ket_angmoms == std::pair<int, int>({1, 4})))
    {
        erirec::comp_electron_repulsion_dfpg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({2, 3})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_dfdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({2, 3})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_dfdf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({2, 3})) && (ket_angmoms == std::pair<int, int>({2, 4})))
    {
        erirec::comp_electron_repulsion_dfdg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({2, 3})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_dfff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({2, 3})) && (ket_angmoms == std::pair<int, int>({3, 4})))
    {
        erirec::comp_electron_repulsion_dffg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({2, 3})) && (ket_angmoms == std::pair<int, int>({4, 4})))
    {
        erirec::comp_electron_repulsion_dfgg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({2, 4})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_dgss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({2, 4})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_dgsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({2, 4})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_dgsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({2, 4})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_dgsf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({2, 4})) && (ket_angmoms == std::pair<int, int>({0, 4})))
    {
        erirec::comp_electron_repulsion_dgsg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({2, 4})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_dgpp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({2, 4})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_dgpd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({2, 4})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_dgpf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({2, 4})) && (ket_angmoms == std::pair<int, int>({1, 4})))
    {
        erirec::comp_electron_repulsion_dgpg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({2, 4})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_dgdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({2, 4})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_dgdf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({2, 4})) && (ket_angmoms == std::pair<int, int>({2, 4})))
    {
        erirec::comp_electron_repulsion_dgdg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({2, 4})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_dgff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({2, 4})) && (ket_angmoms == std::pair<int, int>({3, 4})))
    {
        erirec::comp_electron_repulsion_dgfg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({2, 4})) && (ket_angmoms == std::pair<int, int>({4, 4})))
    {
        erirec::comp_electron_repulsion_dggg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({3, 3})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_ffss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({3, 3})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_ffsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({3, 3})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_ffsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({3, 3})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_ffsf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({3, 3})) && (ket_angmoms == std::pair<int, int>({0, 4})))
    {
        erirec::comp_electron_repulsion_ffsg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({3, 3})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_ffpp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({3, 3})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_ffpd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({3, 3})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_ffpf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({3, 3})) && (ket_angmoms == std::pair<int, int>({1, 4})))
    {
        erirec::comp_electron_repulsion_ffpg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({3, 3})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_ffdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({3, 3})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_ffdf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({3, 3})) && (ket_angmoms == std::pair<int, int>({2, 4})))
    {
        erirec::comp_electron_repulsion_ffdg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({3, 3})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_ffff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({3, 3})) && (ket_angmoms == std::pair<int, int>({3, 4})))
    {
        erirec::comp_electron_repulsion_fffg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({3, 3})) && (ket_angmoms == std::pair<int, int>({4, 4})))
    {
        erirec::comp_electron_repulsion_ffgg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({3, 4})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_fgss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({3, 4})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_fgsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({3, 4})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_fgsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({3, 4})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_fgsf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({3, 4})) && (ket_angmoms == std::pair<int, int>({0, 4})))
    {
        erirec::comp_electron_repulsion_fgsg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({3, 4})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_fgpp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({3, 4})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_fgpd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({3, 4})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_fgpf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({3, 4})) && (ket_angmoms == std::pair<int, int>({1, 4})))
    {
        erirec::comp_electron_repulsion_fgpg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({3, 4})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_fgdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({3, 4})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_fgdf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({3, 4})) && (ket_angmoms == std::pair<int, int>({2, 4})))
    {
        erirec::comp_electron_repulsion_fgdg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({3, 4})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_fgff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({3, 4})) && (ket_angmoms == std::pair<int, int>({3, 4})))
    {
        erirec::comp_electron_repulsion_fgfg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({3, 4})) && (ket_angmoms == std::pair<int, int>({4, 4})))
    {
        erirec::comp_electron_repulsion_fggg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({4, 4})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_ggss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({4, 4})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_ggsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({4, 4})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_ggsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({4, 4})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_ggsf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({4, 4})) && (ket_angmoms == std::pair<int, int>({0, 4})))
    {
        erirec::comp_electron_repulsion_ggsg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({4, 4})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_ggpp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({4, 4})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_ggpd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({4, 4})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_ggpf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({4, 4})) && (ket_angmoms == std::pair<int, int>({1, 4})))
    {
        erirec::comp_electron_repulsion_ggpg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({4, 4})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_ggdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({4, 4})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_ggdf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({4, 4})) && (ket_angmoms == std::pair<int, int>({2, 4})))
    {
        erirec::comp_electron_repulsion_ggdg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({4, 4})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_ggff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({4, 4})) && (ket_angmoms == std::pair<int, int>({3, 4})))
    {
        erirec::comp_electron_repulsion_ggfg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({4, 4})) && (ket_angmoms == std::pair<int, int>({4, 4})))
    {
        erirec::comp_electron_repulsion_gggg(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
        return;
    }

    std::cout << " *** Integral not found in call tree :" << bra_angmoms.first << " , " << bra_angmoms.second << " , " << ket_angmoms.first << " , "
              << ket_angmoms.second << std::endl;
}

}  // namespace erifunc

#endif /* ElectronRepulsionGeom1000Func_hpp */

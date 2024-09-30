#ifndef ElectronRepulsionGeom1000Func_hpp
#define ElectronRepulsionGeom1000Func_hpp

#include <cstddef>
#include <iostream>
#include <utility>

#include "ElectronRepulsionGeom1000RecSSSS.hpp"
#include "ElectronRepulsionGeom1000RecSSSP.hpp"
#include "ElectronRepulsionGeom1000RecSSSD.hpp"
#include "ElectronRepulsionGeom1000RecSSSF.hpp"
#include "ElectronRepulsionGeom1000RecSSSG.hpp"
#include "ElectronRepulsionGeom1000RecSSPP.hpp"
#include "ElectronRepulsionGeom1000RecSSPD.hpp"
#include "ElectronRepulsionGeom1000RecSSPF.hpp"
#include "ElectronRepulsionGeom1000RecSSDD.hpp"
#include "ElectronRepulsionGeom1000RecSSDF.hpp"
#include "ElectronRepulsionGeom1000RecSSFF.hpp"
#include "ElectronRepulsionGeom1000RecSPSS.hpp"
#include "ElectronRepulsionGeom1000RecSPSP.hpp"
#include "ElectronRepulsionGeom1000RecSPSD.hpp"
#include "ElectronRepulsionGeom1000RecSPSF.hpp"
#include "ElectronRepulsionGeom1000RecSPPP.hpp"
#include "ElectronRepulsionGeom1000RecSPPD.hpp"
#include "ElectronRepulsionGeom1000RecSPPF.hpp"
#include "ElectronRepulsionGeom1000RecSPDD.hpp"
#include "ElectronRepulsionGeom1000RecSPDF.hpp"
#include "ElectronRepulsionGeom1000RecSPFF.hpp"
#include "ElectronRepulsionGeom1000RecSDSS.hpp"
#include "ElectronRepulsionGeom1000RecSDSP.hpp"
#include "ElectronRepulsionGeom1000RecSDSD.hpp"
#include "ElectronRepulsionGeom1000RecSDSF.hpp"
#include "ElectronRepulsionGeom1000RecSDPP.hpp"
#include "ElectronRepulsionGeom1000RecSDPD.hpp"
#include "ElectronRepulsionGeom1000RecSDPF.hpp"
#include "ElectronRepulsionGeom1000RecSDDD.hpp"
#include "ElectronRepulsionGeom1000RecSDDF.hpp"
#include "ElectronRepulsionGeom1000RecSDFF.hpp"
#include "ElectronRepulsionGeom1000RecSFSS.hpp"
#include "ElectronRepulsionGeom1000RecSFSP.hpp"
#include "ElectronRepulsionGeom1000RecSFSD.hpp"
#include "ElectronRepulsionGeom1000RecSFSF.hpp"
#include "ElectronRepulsionGeom1000RecSFPP.hpp"
#include "ElectronRepulsionGeom1000RecSFPD.hpp"
#include "ElectronRepulsionGeom1000RecSFPF.hpp"
#include "ElectronRepulsionGeom1000RecSFDD.hpp"
#include "ElectronRepulsionGeom1000RecSFDF.hpp"
#include "ElectronRepulsionGeom1000RecSFFF.hpp"
#include "ElectronRepulsionGeom1000RecPPSS.hpp"
#include "ElectronRepulsionGeom1000RecPPSP.hpp"
#include "ElectronRepulsionGeom1000RecPPPP.hpp"
#include "ElectronRepulsionGeom1000RecPSSS.hpp"
#include "ElectronRepulsionGeom1000RecPSSP.hpp"
#include "ElectronRepulsionGeom1000RecPSSD.hpp"
#include "ElectronRepulsionGeom1000RecPSSF.hpp"
#include "ElectronRepulsionGeom1000RecPSPP.hpp"
#include "ElectronRepulsionGeom1000RecPSPD.hpp"
#include "ElectronRepulsionGeom1000RecPSPF.hpp"
#include "ElectronRepulsionGeom1000RecPSDD.hpp"
#include "ElectronRepulsionGeom1000RecPSDF.hpp"
#include "ElectronRepulsionGeom1000RecPSFF.hpp"
#include "GtoPairBlock.hpp"

namespace erifunc {  // erifunc namespace

/// Computes vector of integrals of requested four center integral.
/// @param distributor  The pointer to distributor of integrals.
/// @param bra_gto_pair_block The basis function pairs block on bra side.
/// @param ket_gto_pair_block The basis function pairs block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
compute_geom_1000(T&                               distributor,
                  const CGtoPairBlock&             bra_gto_pair_block,
                  const CGtoPairBlock&             ket_gto_pair_block,
                  const std::pair<size_t, size_t>& bra_indices,
                  const std::pair<size_t, size_t>& ket_indices) -> void
{
    const auto bra_angmoms = bra_gto_pair_block.angular_momentums();

    const auto ket_angmoms = ket_gto_pair_block.angular_momentums();
    
    // leading [SS|XX] terms

    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom_1000_ssss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom_1000_sssp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom_1000_sssd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_geom_1000_sssf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom_1000_sspp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom_1000_sspd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_geom_1000_sspf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom_1000_ssdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_geom_1000_ssdf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_geom_1000_ssff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    // leading [SP|XX] terms
    
    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom_1000_spss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom_1000_spsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom_1000_spsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_geom_1000_spsf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom_1000_sppp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom_1000_sppd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_geom_1000_sppf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom_1000_spdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_geom_1000_spdf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_geom_1000_spff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }

    // leading [PS|XX] terms
    
    if ((bra_angmoms == std::pair<int, int>({1, 0})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom_1000_psss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 0})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom_1000_pssp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 0})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom_1000_pssd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 0})) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        erirec::comp_electron_repulsion_geom_1000_pssf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 0})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom_1000_pspp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 0})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom_1000_pspd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 0})) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        erirec::comp_electron_repulsion_geom_1000_pspf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 0})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom_1000_psdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 0})) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        erirec::comp_electron_repulsion_geom_1000_psdf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 0})) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        erirec::comp_electron_repulsion_geom_1000_psff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    // leading [SD|XX] terms
    
   if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom_1000_sdss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({0, 1})))
     {
         erirec::comp_electron_repulsion_geom_1000_sdsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
         return;
     }
    
    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({0, 2})))
     {
         erirec::comp_electron_repulsion_geom_1000_sdsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
         return;
     }
    
    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({0, 3})))
     {
         erirec::comp_electron_repulsion_geom_1000_sdsf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
         return;
     }
    
    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({1, 1})))
     {
         erirec::comp_electron_repulsion_geom_1000_sdpp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
         return;
     }
    
    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({1, 2})))
     {
         erirec::comp_electron_repulsion_geom_1000_sdpd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
         return;
     }
    
    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({1, 3})))
     {
         erirec::comp_electron_repulsion_geom_1000_sdpf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
         return;
     }
    
    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({2, 2})))
     {
         erirec::comp_electron_repulsion_geom_1000_sddd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
         return;
     }
    
    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({2, 3})))
     {
         erirec::comp_electron_repulsion_geom_1000_sddf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
         return;
     }
    
    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({3, 3})))
     {
         erirec::comp_electron_repulsion_geom_1000_sdff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
         return;
     }
    
    // leading [SF|XX] terms
    
   if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom_1000_sfss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({0, 1})))
     {
         erirec::comp_electron_repulsion_geom_1000_sfsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
         return;
     }
    
    if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({0, 2})))
     {
         erirec::comp_electron_repulsion_geom_1000_sfsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
         return;
     }
    
    if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({0, 3})))
     {
         erirec::comp_electron_repulsion_geom_1000_sfsf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
         return;
     }
    
    if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({1, 1})))
     {
         erirec::comp_electron_repulsion_geom_1000_sfpp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
         return;
     }
    
    if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({1, 2})))
     {
         erirec::comp_electron_repulsion_geom_1000_sfpd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
         return;
     }
    
    if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({1, 3})))
     {
         erirec::comp_electron_repulsion_geom_1000_sfpf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
         return;
     }
    
    if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({2, 2})))
     {
         erirec::comp_electron_repulsion_geom_1000_sfdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
         return;
     }
    
    if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({2, 3})))
     {
         erirec::comp_electron_repulsion_geom_1000_sfdf(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
         return;
     }
    
    if ((bra_angmoms == std::pair<int, int>({0, 3})) && (ket_angmoms == std::pair<int, int>({3, 3})))
     {
         erirec::comp_electron_repulsion_geom_1000_sfff(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
         return;
     }
    
    // leading [PP|XX] terms 
    
    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom_1000_ppss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom_1000_ppsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom_1000_pppp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }

    //std::cout << " *** Integral not found in call tree :" << bra_angmoms.first << " , " << bra_angmoms.second << " , " << ket_angmoms.first << " , "
    //          << ket_angmoms.second << std::endl;
}

}  // namespace erifunc

#endif /* ElectronRepulsionGeom1000Func_hpp */

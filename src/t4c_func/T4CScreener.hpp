#ifndef T4CScreener_hpp
#define T4CScreener_hpp

#include <cstdint>
#include <vector>
#include <string>

#include "GtoPairBlock.hpp"
#include "BlockedGtoPairBlock.hpp"


/// Class CT4CScreener provides methods for storing blocked GTOs pair blocks
/// partitioned according to Cauchy–Schwarz relationship.
class CT4CScreener
{
    /// Vectot of partitioned GTOs pair blocks.
    std::vector<CBlockedGtoPairBlock> _gto_pair_blocks;
    
    public:
    
    /// Creates an four center integrals screener.
    CT4CScreener() = default;
    
    
    /// Partitions GTOs pair blocks for given molecule and molecular basis for given
    /// type of four center integral.
    /// - Parameter basis : the molecular basis.
    /// - Parameter molecule : the molecule.
    /// - Parameter label : the label of four center integral.
    auto partition(const CMolecularBasis& basis,
                   const CMolecule&       molecule,
                   const std::string&     label) -> void;
    
    /// Gets vector of blocked GTOs pair blocks.
    auto gto_pair_blocks() const -> std::vector<CBlockedGtoPairBlock>;
};


#endif /* T4CScreener_hpp */

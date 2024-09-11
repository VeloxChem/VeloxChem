#ifndef T4CScreener_hpp
#define T4CScreener_hpp

#include <cstdint>
#include <string>
#include <vector>

#include "BlockedGtoPairBlock.hpp"
#include "GtoPairBlock.hpp"

/// @brief Class CT4CScreener provides methods for storing blocked GTOs pair blocks
/// partitioned according to Cauchy–Schwarz relationship.
class CT4CScreener
{
   public:
    /// Creates an four center integrals screener.
    CT4CScreener() = default;

    /// @brief The default copy constructor.
    /// @param other The integrals screener to be copied.
    CT4CScreener(const CT4CScreener &other) = delete;

    /// @brief The default move constructor.
    /// @param other The integrals screener  to be moved.
    CT4CScreener(CT4CScreener &&other) noexcept = delete;

    /// @brief The default destructor.
    ~CT4CScreener() = default;

    /// @brief The default copy assignment operator.
    /// @param other The integrals screener to be copy assigned.
    /// @return The assigned integrals screener.
    auto operator=(const CT4CScreener &other) -> CT4CScreener & = delete;

    /// @brief The default move assignment operator.
    /// @param other The integrals screener to be move assigned.
    /// @return The assigned integrals screener .
    auto operator=(CT4CScreener &&other) noexcept -> CT4CScreener & = delete;

    /// @brief The equality operator.
    /// @param other The integrals screener  to be compared.
    /// @return True if integrals screeners  are equal, False otherwise.
    auto operator==(const CT4CScreener &other) const -> bool = delete;

    /// @brief The non-equality operator.
    /// @param other The integrals screener to be compared.
    /// @return True if integrals screeners  are not equal, False otherwise.
    auto operator!=(const CT4CScreener &other) const -> bool = delete;

    /// @brief Partitions basis function pair blocks for given molecule and molecular basis for given
    /// type of four center integral.
    /// @param basis  The molecular basis.
    /// @param molecule  The molecule.
    /// @param label  The label of four center integral.
    auto partition(const CMolecularBasis &basis, const CMolecule &molecule, const std::string &label) -> void;

    /// @brief Gets vector of blocked basis function pair blocks.
    /// @return The vector of  blocked basis function pair blocks.
    auto gto_pair_blocks() const -> std::vector<CBlockedGtoPairBlock>;

   private:
    /// @brief Vector of partitioned basis function pair blocks.
    std::vector<CBlockedGtoPairBlock> _gto_pair_blocks;
};

#endif /* T4CScreener_hpp */

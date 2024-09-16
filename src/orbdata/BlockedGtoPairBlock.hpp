#ifndef BlockedGtoPairBlock_hpp
#define BlockedGtoPairBlock_hpp

#include <array>
#include <vector>

#include "GtoPairBlock.hpp"

/// @brief Class CBlockedGtoPairBlock stores data about packed basis function pairs and provides set of
/// methods for manipulating with packed basis function pairs.
class CBlockedGtoPairBlock
{
   public:
    /// Creates an empty  blocked basis function pairs block.
    CBlockedGtoPairBlock() = default;
    
    /// @brief Creates a blocked basis function pairs block.
    /// @param gto_pair_blocks  The array of  basis function pairs blocks.
    CBlockedGtoPairBlock(const std::array<CGtoPairBlock, 16>& gto_pair_blocks);

    /// @brief Creates a blocked basis function pairs block.
    /// @param gto_pair_blocks  The vector of  basis function pairs blocks.
    /// @param block_indices  The vector of indexes to distribute the vector of blocked basis function pairs blocks.
    CBlockedGtoPairBlock(const std::vector<CGtoPairBlock> &gto_pair_blocks, const std::vector<int> &block_indices);

    /// Creates a blocked basis function pairs block.
    /// @param gto_pair_block  The basis function pairs block.
    /// @param integrals The vector of integrals to partition  basis function pairs block.
    CBlockedGtoPairBlock(const CGtoPairBlock &gto_pair_block, const std::vector<double> &integrals);

    /// @brief The default copy constructor.
    /// @param other The blocked basis function pairs block to be copied.
    CBlockedGtoPairBlock(const CBlockedGtoPairBlock &other);

    /// @brief The default move constructor.
    /// @param other The blocked basis functions pairs block to be moved.
    CBlockedGtoPairBlock(CBlockedGtoPairBlock &&other) noexcept;

    /// @brief The default destructor.
    ~CBlockedGtoPairBlock() = default;

    /// @brief The default copy assignment operator.
    /// @param other The blocked basis function pairs block to be copy assigned.
    /// @return The assigned basis function pairs block.
    auto operator=(const CBlockedGtoPairBlock &other) -> CBlockedGtoPairBlock &;

    /// @brief The default move assignment operator.
    /// @param other The blocked basis function pairs block to be move assigned.
    /// @return The assigned basis function pairs block.
    auto operator=(CBlockedGtoPairBlock &&other) noexcept -> CBlockedGtoPairBlock &;

    /// @brief The equality operator.
    /// @param other The blocked basis function pairs block to be compared.
    /// @return True if blocked basis function pairs blocks are equal, False otherwise.
    auto operator==(const CBlockedGtoPairBlock &other) const -> bool;

    /// @brief The equality operator.
    /// @param other The blocked basis function pairs  block to be compared.
    /// @return True if blocked basis function pairs blocks are not equal, False
    /// otherwise.
    auto operator!=(const CBlockedGtoPairBlock &other) const -> bool;

    /// @brief Gets selected basis function pairs block.
    /// @param index  The index of basis function pairs block.
    /// @return The selected basis function pairs blocks.
    auto gto_pair_block(const int index) const -> CGtoPairBlock;

    /// @brief Checks if selected basis function pairs block is empty.
    /// @param index  The index of basis function pairs block.
    /// @return True if selected basis function pairs block is empty, False otherwise.
    auto is_empty_gto_pair_block(const int index) const -> bool;
    
    /// @brief Gets array of basis function pairs blocks.
    /// @return The array of basis function pairs blocks.
    auto gto_pair_blocks() const -> std::array<CGtoPairBlock, 16>;

   private:
    /// @brief The array of basis function pairs blocks.
    std::array<CGtoPairBlock, 16> _gto_pair_blocks;
};

#endif /* BlockedGtoPairBlock_hpp */

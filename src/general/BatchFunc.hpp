#ifndef BatchFunc_hpp
#define BatchFunc_hpp

#include <utility>

#include "CustomConstrains.hpp"

namespace batch {  // batch

/// @brief Gets number of batches to partition vector.
/// @param nelements The number of elements in vector.
/// @param bsize The batch size used to partition vector.
/// @return The number of batches.
template <Integral T>
inline auto
number_of_batches(const T nelements, const T bsize) -> T
{
    const auto nbatches = nelements / bsize;

    return ((nelements % bsize) != 0) ? nbatches + 1 : nbatches;
}

/// @brief Gets range of specific batch in partitioned vector.
/// @param ibatch The index of batch.
/// @param nelements The number of elements in vector.
/// @param bsize The batch size used to partition vector.
/// @return The  [first, last) pair  for requested batch.
template <Integral T>
inline auto
batch_range(const T ibatch, const T nelements, const T bsize) -> std::pair<T, T>
{
    const auto first = ibatch * bsize;

    if (first > nelements) return {nelements, nelements};

    const auto last = first + bsize;

    if (last > nelements)
    {
        return {first, nelements};
    }
    else
    {
        return {first, last};
    }
}

/// @brief Gets range of specific batch in partitioned vector.
/// @param ibatch The index of batch.
/// @param nelements The number of elements in vector.
/// @param bsize The batch size used to partition vector.
/// @param position The shift of pair positions.
/// @return The  [first, last) pair  for requested batch.
template <Integral T>
inline auto
batch_range(const T ibatch, const T nelements, const T bsize, const T position) -> std::pair<T, T>
{
    auto range = batch::batch_range(ibatch, nelements, bsize);
    
    return {range.first + position, range.second + position};
}

}  // namespace batch

#endif /* BatchFunc_hpp */

#ifndef BatchFunc_hpp
#define BatchFunc_hpp

#include <cstdint>
#include <utility>

namespace batch {  // batch namespace

/**
 Gets starting index of requested batch.

 @param ibatch the index of batch.
 @param nelements the number of elements in vector.
 @param nbatches the number of batches to partition vector.
 @return the starting index of batch.
 */
inline auto
getBatchIndex(const int64_t ibatch, const int64_t nelements, const int64_t nbatches) -> int64_t
{
    const auto bdim = nelements / nbatches;

    const auto brem = nelements % nbatches;

    if (const auto boff = bdim * ibatch; boff > nelements)
    {
        return nelements;
    }
    else
    {
        return boff + ((ibatch <= brem) ? ibatch : brem);
    }
};

/**
 Gets number of batches to partition vector.

 @param nelements the number of elements in vector.
 @param bsize the batch size used to partition vector.
 @return the number of batches.
 */
inline auto
getNumberOfBatches(const int64_t nelements, const int64_t bsize) -> int64_t
{
    const auto nbatches = nelements / bsize;

    if ((nelements % bsize) != 0)
    {
        return nbatches + 1;
    }
    else
    {
        return nbatches;
    }
};

/**
 Gets range of specific batch in partitioned vector.

 @param ibatch the index of batch.
 @param nelements the number of elements in vector.
 @param bsize the batch size used to partition vector.
 @return the number of loop passes.
 */
inline auto
getBatchRange(const int64_t ibatch, const int64_t nelements, const int64_t bsize) -> std::pair<int64_t, int64_t>
{
    const auto first = ibatch * bsize;

    if (first > nelements)
    {
        return {nelements, nelements};
    }

    const auto last = first + bsize;

    if (last > nelements)
    {
        return {first, nelements};
    }
    else
    {
        return {first, last};
    }
};

}  // namespace batch

#endif /* BatchFunc_hpp */

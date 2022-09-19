#ifndef CBufferTest_hpp
#define CBufferTest_hpp

#include <gtest/gtest.h>

#include <type_traits>

#include "Buffer.hpp"

template <typename T, typename B, auto M = Dynamic, auto N = Dynamic>
class CBufferParameters
{
    using size_type = std::decay_t<decltype(Dynamic)>;

   public:
    using value_type   = T;
    using backend_type = B;

    static constexpr size_type NRows = M;
    static constexpr size_type NCols = N;
};

template <typename T>
class CBufferTest : public testing::Test
{
};

namespace detail {
using implementations = ::testing::Types<
    /* @{ BufferX tests */
    CBufferParameters<int32_t, mem::Host, 1, Dynamic>,
    CBufferParameters<float, mem::Host, 1, Dynamic>,
    CBufferParameters<double, mem::Host, 1, Dynamic>,
    /* @} */
    /* @{ BufferN tests */
    CBufferParameters<int32_t, mem::Host, 1, 10>,
    CBufferParameters<float, mem::Host, 1, 10>,
    CBufferParameters<double, mem::Host, 1, 10>,
    /* @} */
    /* @{ BufferXY tests */
    CBufferParameters<int32_t, mem::Host, Dynamic, Dynamic>,
    CBufferParameters<float, mem::Host, Dynamic, Dynamic>,
    CBufferParameters<double, mem::Host, Dynamic, Dynamic>,
    /* @} */
    /* @{ BufferXN tests */
    CBufferParameters<int32_t, mem::Host, Dynamic, 10>,
    CBufferParameters<float, mem::Host, Dynamic, 10>,
    CBufferParameters<double, mem::Host, Dynamic, 10>,
    /* @} */
    /* @{ BufferMY tests */
    CBufferParameters<int32_t, mem::Host, 10, Dynamic>,
    CBufferParameters<float, mem::Host, 10, Dynamic>,
    CBufferParameters<double, mem::Host, 10, Dynamic>,
    /* @} */
    /* @{ BufferMN tests */
    CBufferParameters<int32_t, mem::Host, 10, 5>,
    CBufferParameters<float, mem::Host, 10, 5>,
    CBufferParameters<double, mem::Host, 10, 5>
    /* @} */
    >;
}  // namespace detail

#endif  // CBufferTest_hpp

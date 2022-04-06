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

#endif  // CBufferTest_hpp

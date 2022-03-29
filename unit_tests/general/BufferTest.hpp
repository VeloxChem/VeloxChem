#ifndef BufferTest_hpp
#define BufferTest_hpp

#include <gtest/gtest.h>

#include "Buffer.hpp"

// how to parametrize over multiple template type paramater
// https://stackoverflow.com/a/52207657/2528668

template <typename T>
class BufferXTest : public testing::Test
{
};

#endif  // BufferTest_hpp

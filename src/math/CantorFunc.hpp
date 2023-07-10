#ifndef CantorFunc_hpp
#define CantorFunc_hpp

#include <cmath>

#include "T2Pair.hpp"

namespace mathconst {  // mathfunc namespace

/**
 Computes Cantor pairing index for two integer non-negative numbers.

 @return the PI constant.
 */
inline auto
getCantorPair(const T2Pair& pair) -> int64_t
{
    const auto x = pair.first;
    
    const auto y = pair.second;
    
    return (x + y + 1) * (x + y) / 2 + y;
}

/**
 Reduces Cantor pairing index to two integer non-negative numbers.

 @return the PI constant.
 */
inline auto
getCantorPair(const int64_t index) -> T2Pair
{
    const auto w = static_cast<int64_t>(std::floor(0.5 * (std::sqrt(8 * index + 1) - 1)));
        
    const auto y = index - (w * w + w) / 2;
    
    return {w - y, y};
}

}  // namespace mathfunc

#endif /* CantorFunc_hpp */

#ifndef CustomConstrains_hpp
#define CustomConstrains_hpp

#include <concepts>

template <typename T>
concept Integral = std::is_integral_v<T>;

template <typename T>
concept FloatingPoint = std::is_floating_point_v<T>;

#endif /* CustomConstrains_hpp */

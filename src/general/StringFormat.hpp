#ifndef StringFormat_hpp
#define StringFormat_hpp

#include <string>

namespace format {  /// format

/// @brief Creates upper case copy of given string.
/// @param source The string to be trasnformed to upper case variant.
/// @return The upper cased string.
auto upper_case(const std::string &source) -> std::string;

/// @brief Creates lower case copy of given string.
/// @param source The string to be trasnformed to lower case variant.
/// @return The lower cased string.
auto lower_case(const std::string &source) -> std::string;

}  // namespace format

#endif /* StringFormat_hpp */

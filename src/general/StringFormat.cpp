#include "StringFormat.hpp"

#include <algorithm>
#include <cctype>
#include <iterator>

namespace format {  /// format

auto
upper_case(const std::string &source) -> std::string
{
    auto destination = std::string();

    destination.reserve(source.size());

    std::ranges::transform(source, std::back_inserter(destination), [](auto c) { return std::toupper(c); });

    return destination;
}

auto
lower_case(const std::string &source) -> std::string
{
    auto destination = std::string();

    destination.reserve(source.size());

    std::ranges::transform(source, std::back_inserter(destination), [](auto c) { return std::tolower(c); });

    return destination;
}

}  // namespace format

#include "ChemicalElement.hpp"

#include <algorithm>
#include <iterator>

namespace chem_elem {

auto
name(const int id) -> std::string
{
    return _names.at(id);
}

auto
label(const int id) -> std::string
{
    auto label = _names.at(id);

    if (label.size() == 2)
    {
        label[1] = std::tolower(label[1]);
    }
    return label;
}

auto
identifier(const std::string &name) -> int
{
    if (auto it = std::ranges::find(_names, name); it != _names.end())
    {
        return static_cast<int>(std::distance(_names.begin(), it));
    }
    else
    {
        return -1;
    }
}

auto
mass(const int id) -> double
{
    return _masses.at(id);
}

auto
max_angular_momentum(const int id) -> int
{
    if ((id > 0) && (id < 5)) return 0;

    if ((id > 4) && (id < 21)) return 1;

    if ((id > 20) && (id < 57)) return 2;

    if ((id > 56) && (id < 87)) return 3;

    return -1;
}

auto
max_identifier() -> int
{
    return static_cast<int>(_names.size()) - 1;
}

}  // namespace chem_elem

#ifndef GeomKeys_hpp
#define GeomKeys_hpp

#include <cstdint>
#include <utility>
#include <array>

using TGeomPair = std::pair<int64_t, char>;

using T2GeomKey = std::array<TGeomPair, 2>;

using T3GeomKey = std::array<TGeomPair, 3>;

using T4GeomKey = std::array<TGeomPair, 4>;

#endif /* GeomKeys_hpp */

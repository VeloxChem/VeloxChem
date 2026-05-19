//
//  Tabula — custom-recursion molecular-integral machinery.
//  Per-thread load balance of a two-center driver's task loop.
//

#ifndef TabulaThreadBalance_hpp
#define TabulaThreadBalance_hpp

#include <vector>

namespace tabula {  // tabula namespace

/// @brief Per-thread load balance of the task parallel loop in a two-center
/// driver's most recent compute run.
struct ThreadBalance
{
    /// @brief Wall time of the whole parallel region.
    double wall{0.0};
    /// @brief Busy wall seconds of each OpenMP thread.
    std::vector<double> busy;
    /// @brief Number of tasks handled by each OpenMP thread.
    std::vector<long> pairs;
};

}  // namespace tabula

#endif /* TabulaThreadBalance_hpp */

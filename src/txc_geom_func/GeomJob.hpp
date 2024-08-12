#ifndef GeomJob_hpp
#define GeomJob_hpp

#include <cmath>
#include <cstdint>
#include <map>
#include <vector>

#include "GeomArchetypeTask.hpp"
#include "MolecularBasis.hpp"

template <class T, class U, class W>
class CGeomJob
{
   public:
    CGeomJob(){};

    auto
    set_result_holders(const std::vector<T>& results) -> void
    {
        _results = results;
    }

    auto
    set_operators(const std::vector<std::string>& operators) -> void
    {
        _operators = operators;
    }

    auto
    set_job_type(const std::string& jobtype) -> void
    {
        _jobtype = jobtype;
    }

    auto
    setTasks(const std::vector<CGeomArchetypeTask<U, W>>& tasks) -> void
    {
        _tasks = tasks;
    }

    auto
    allocate(const CMolecularBasis& basis) -> void
    {
        if (!_tasks.empty())
        {
            const auto geom_order = _tasks[0].getArchetype();

            for (auto& result : _results)
            {
                result.allocate(basis, geom_comps);
            }
        }
    }

    auto
    get_raw_results() -> T*
    {
        return _results.data();
    }

    auto
    get_operators() -> const std::vector<std::string>
    {
        return _operators;
    }

    auto
    get_tasks() -> const std::vector<CGeomArchetypeTask<U, W>>
    {
        return _tasks;
    }

    auto
    get_results() -> std::vector<T>
    {
        return _results;
    }

   private:
    /// @brief Vector of results from set of geometrical derivatives computations.
    std::vector<T> _results;

    /// @brief Vector of requested operator labels.
    std::vector<std::string> _operators;

    /// @brief Type of geometrical derivatives evaluation.
    std::string _jobtype;

    /// @brief Vector of tasks.
    std::vector<CGeomArchetypeTask<U, W>> _tasks;
};

#endif /* GeomJob_hpp */

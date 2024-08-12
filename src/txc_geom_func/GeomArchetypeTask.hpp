#ifndef GeomArchetypeTask_hpp
#define GeomArchetypeTask_hpp

#include <cstdint>
#include <map>
#include <vector>

template <class T, class U>
class CGeomArchetypeTask
{
   public:
    CGeoIntArchetypeTask(){};

    auto
    set_archetype(const T& archetype) -> void
    {
        _archetype = archetype;
    }

    auto
    get_archetype() const -> T
    {
        return _archetype;
    }

    auto
    set_contraction_patterns(const std::map<size_t, std::vector<T>>& cpatterns) -> void
    {
        _cpatterns = cpatterns;
    }

    auto
    get_contraction_patterns_for_perm(const size_t perm_key) const -> std::vector<T>
    {
        return _cpatterns.at(perm_key);
    }

    auto
    add_integral_task(const T& atoms, const std::map<U, std::map<size_t, size_t>>& cpatterns) -> void
    {
        _integral_tasks[atoms] = cpatterns;
    }

    auto
    getIntegralTask(const T& atoms) const -> std::map<U, std::map<size_t, size_t>>
    {
        return _integral_tasks.at(atoms);
    }

   private:
    /// @brief Geometrical derivative archetype.
    T _archetype;

    /// @brief Contraction patterns map.
    std::map<size_t, std::vector<T>> _cpatterns;

    /// Integral tasks {Atom index: {T4Index (unique atoms only) integral: {contraction pattern permutation:
    /// coefficient}}}
    std::map<T, std::map<U, std::map<size_t, size_t>>> _integral_tasks;
};

#endif /* GeomArchetypeTask_hpp */

#ifndef TwoElectronExpectation_hpp
#define TwoElectronExpectation_hpp

#include <cstdint>
#include <map>
#include <vector>

#include "MolecularBasis.hpp"

template <class T>
class CTwoElectronExpectation
{
   public:
    CTwoElectronExpectation() = default;

    auto
    add_contr_pairs_for_atoms_tuple(const T& atoms, const std::map<std::size_t, std::vector<std::size_t>>& keys) -> void
    {
        std::map<std::size_t, std::map<std::size_t, std::map<std::size_t, double>>> storage;

        for (const auto& [contr_1, contr_2] : keys)
        {
            std::map<std::size_t, std::map<std::size_t, double>> new_entry;

            for (const auto& inner : contr_2)
            {
                new_entry.insert({inner, std::map<std::size_t, double>()});
            }

            storage.insert({contr_1, new_entry});
        }

        _results.insert({atoms, storage});
    }

    auto
    allocate(const CMolecularBasis& basis, const size_t nvals) -> void
    {
        std::map<std::size_t, double> tensor;

        for (size_t i = 0; i < nvals; i++)
        {
            tensor[i] = 0.0;
        }

        for (auto& [atoms, mmvalues] : _results)
        {
            for (auto& [key, mvalues] : mmvalues)
            {
                for (auto& [key, values] : mvalues)
                {
                    mvalues[key] = tensor;
                }
            }
        }
    }

    auto
    get_map() -> std::map<T, std::map<std::size_t, std::map<std::size_t, std::map<std::size_t, double>>>>&
    {
        return _results;
    }

    auto
    get_value(const T& atoms, const std::size_t key1, const std::size_t key2, const std::size_t component) const
        -> double
    {
        return _results.at(atoms).at(key1).at(key2).at(component);
    }

   private:
    /// @brief The map of two-electron expectation values.
    std::map<T, std::map<std::size_t, std::map<std::size_t, std::map<std::size_t, double>>>> _results;
};

#endif /* TwoElectronExpectation_hpp */

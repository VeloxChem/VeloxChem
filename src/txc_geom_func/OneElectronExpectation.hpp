#ifndef OneElectronExpectation_hpp
#define OneElectronExpectation_hpp

#include <cstdint>
#include <map>
#include <vector>

#include "MolecularBasis.hpp"

template <class T>
class COneElectronExpectation
{
   public:
    COneElectronExpectation() = default;

    auto
    add_contractions_for_atoms_tuple(const T& atoms, const std::vector<std::size_t>& keys) -> void
    {
        std::map<std::size_t, std::map<std::size_t, double>> storage;

        for (const auto key : keys)
        {
            storage.insert({key, std::map<std::size_t, double>()});
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

        for (auto& [atoms, mvalues] : _results)
        {
            for (auto& [key, values] : mvalues)
            {
                mvalues[key] = tensor;
            }
        }
    }

    auto
    get_map() -> std::map<T, std::map<std::size_t, std::map<std::size_t, double>>>&
    {
        return _results;
    }

    auto
    get_map_at_atoms(const T& atoms) -> std::map<std::size_t, std::map<std::size_t, double>>&
    {
        return _results.at(atoms);
    }

    auto
    get_ptr_to_map_at_atoms(const T& atoms) -> std::map<std::size_t, std::map<std::size_t, double>>*
    {
        return &_results.at(atoms);
    }

    auto
    get_value(const T& atoms, const std::size_t key, const std::size_t component) const -> double
    {
        return _results.at(atoms).at(key).at(component);
    }

   private:
    /// @brief Map to store expectation values results.
    std::map<T, std::map<std::size_t, std::map<std::size_t, double>>> _results;
};

#endif /* OneElectronExpectation_hpp */

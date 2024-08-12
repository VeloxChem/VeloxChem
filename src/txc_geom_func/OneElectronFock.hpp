#ifndef OneElectronFock_hpp
#define OneElectronFock_hpp

#include <cstdint>
#include <map>
#include <vector>

#include "Matrices.hpp"
#include "MolecularBasis.hpp"

template <class T, class U>
class COneElectronFock
{
   public:
    COneElectronFock(){};

    ~COneElectronFock()
    {
        for (auto& [atoms, mat] : _results)
        {
            delete mat;
        }
    }

    COneElectronFock(const COneElectronFock& other)
    {
        for (const auto& [atoms, mats] : other._results)
        {
            _results.insert({atoms, new CMatrices(*mats)});
        }
    }

    auto
    allocate(const CMolecularBasis& basis, const U& garch) -> void
    {
        CMatrices matrices;

        for (size_t i = 0; i < to_components(garch); i++)
        {
            matrices.add(matfunc::makeMatrix(basis, mat_t::gen), i);
        }

        matrices.zero();

        for (const auto& [atoms, mat] : _results)
        {
            _results[atoms] = new CMatrices(matrices);
        }
    }

    auto
    add_atoms_tuple(const T& atoms) -> void
    {
        _results.insert({atoms, new CMatrices()});
    }

    auto
    get_map() -> std::map<T, CMatrices*>&
    {
        return _results;
    }

    auto
    get_map_at_atoms(const T& atoms) -> CMatrices*
    {
        return _results.at(atoms);
    }

    auto
    get_keys() -> std::vector<T>
    {
        std::vector<T> keys;

        for (const auto& [atoms, mats] : _results)
        {
            keys.push_back(atoms);
        }

        return keys;
    }

    auto
    get_value(const T& atoms, const std::size_t component) -> CMatrix*
    {
        return _results.at(atoms)->getMatrix(component);
    }

    auto
    get_value(const T& atoms, const std::size_t component) const -> const CMatrix*
    {
        return _results.at(atoms)->getMatrix(component);
    }

   private:
    /// @brief Map for storing one-electron Fock matrices.
    std::map<T, CMatrices*> _results;
};

#endif /* OneElectronFock_hpp */

#ifndef TwoElectronFock_hpp
#define TwoElectronFock_hpp

template <class T, class U>
class CTwoElectronFock
{
   public:
    CTwoElectronFock(){};

    ~CTwoElectronFock()
    {
        for (auto& [atoms, mlist] : _results)
        {
            for (auto& mvalue : mlist)
            {
                delete mvalue.second;
            }
        }
    }

    CTwoElectronFock(const CTwoElFock& other)
    {
        for (const auto& [atoms, mlist] : other._results)
        {
            std::map<size_t, CMatrices*> storage;

            for (const auto& [key, mats] : mlist)
            {
                storage.insert({key, new CMatrices(*mats)});
            }

            _results[atoms] = storage;
        }
    }

    auto
    allocate(const CMolecularBasis& basis, const const U& garch) -> void
    {
        CMatrices matrices;

        for (size_t i = 0; i < to_components(garch); i++)
        {
            matrices.add(matfunc::makeMatrix(basis, mat_t::gen), i);
        }

        matrices.zero();

        for (auto& [atoms, mlist] : _results)
        {
            std::map<std::size_t, CMatrices*> storage;

            for (const auto& [key, mats] : mlist)
            {
                storage.insert({key, new CMatrices(matrices)});
            }

            _results[atoms] = storage;
        }
    }

    auto
    add_contractions_for_atoms_tuple(const T& atoms, const std::vector<std::size_t>& keys) -> void
    {
        std::map<size_t, CMatrices*> storage;

        for (const auto key : keys)
        {
            storage.insert({key, new CMatrices()});
        }

        _results.insert({atoms, storage});
    }

    auto
    get_map() -> std::map<T, std::map<std::size_t, CMatrices*>>&
    {
        return _results;
    }

    auto
    getValue(const T& atoms, const std::size_t key, const std::size_t component) -> CMatrix*
    {
        return _results.at(atoms).at(key)->getMatrix(component);
    }

    auto
    getValue(const T4Index& atoms, const std::size_t key, const std::size_t component) const -> const CMatrix*
    {
        return _results.at(atoms).at(key)->getMatrix(component);
    }

   private:
    /// @brief The map of matrices.
    std::map<T, std::map<std::size_t, CMatrices*>> _results;
};

#endif /* TwoElectronFock_hpp */

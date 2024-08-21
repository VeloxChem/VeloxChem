#include "MolecularBasis.hpp"

#include <algorithm>
#include <numeric>
#include <ranges>

#include "ChemicalElement.hpp"
#include "TensorComponents.hpp"

CMolecularBasis::CMolecularBasis()

    : _basis_sets{}

    , _indices{}
{
}

CMolecularBasis::CMolecularBasis(const std::vector<CAtomBasis> &basis_sets, const std::vector<int> &indices)

    : _basis_sets(basis_sets)

    , _indices(indices)
{
}

CMolecularBasis::CMolecularBasis(const CMolecularBasis &other)

    : _basis_sets(other._basis_sets)

    , _indices(other._indices)
{
}

CMolecularBasis::CMolecularBasis(CMolecularBasis &&other) noexcept

    : _basis_sets{}

    , _indices{}
{
    std::swap(_basis_sets, other._basis_sets);

    std::swap(_indices, other._indices);
}

auto
CMolecularBasis::operator=(const CMolecularBasis &other) -> CMolecularBasis &
{
    _basis_sets = other._basis_sets;

    _indices = other._indices;

    return *this;
}

auto
CMolecularBasis::operator=(CMolecularBasis &&other) noexcept -> CMolecularBasis &
{
    std::swap(_basis_sets, other._basis_sets);

    std::swap(_indices, other._indices);

    return *this;
}

auto
CMolecularBasis::operator==(const CMolecularBasis &other) const -> bool
{
    if (_indices != other._indices)
    {
        return false;
    }
    else
    {
        return _basis_sets == other._basis_sets;
    }
}

auto
CMolecularBasis::operator!=(const CMolecularBasis &other) const -> bool
{
    return !(*this == other);
}

auto
CMolecularBasis::add(const CAtomBasis &basis) -> void
{
    auto pos = std::ranges::find_if(
        _basis_sets, [&](const auto &abas) { return (abas.get_name() == basis.get_name()) && (abas.get_identifier() == basis.get_identifier()); });

    if (pos == _basis_sets.end())
    {
        _indices.push_back(static_cast<int>(_basis_sets.size()));

        _basis_sets.push_back(basis);
    }
    else
    {
        _indices.push_back(static_cast<int>(std::distance(_basis_sets.begin(), pos)));
    }
}

auto
CMolecularBasis::reduce_to_valence_basis() const -> CMolecularBasis
{
    std::vector<CAtomBasis> rbasis_sets;

    rbasis_sets.reserve(_basis_sets.size());

    std::ranges::transform(_basis_sets, std::back_inserter(rbasis_sets), [](const auto &abas) { return abas.reduce_to_valence_basis(); });

    return CMolecularBasis(rbasis_sets, _indices);
}

auto
CMolecularBasis::slice(const std::vector<int> &atoms) const -> CMolecularBasis
{
    CMolecularBasis mbasis;

    std::ranges::for_each(atoms, [&](const int i) { mbasis.add(_basis_sets[_indices.at(i)]); });

    return mbasis;
}

auto
CMolecularBasis::basis_sets() const -> std::vector<CAtomBasis>
{
    return _basis_sets;
}

auto
CMolecularBasis::basis_sets_indices() const -> std::vector<int>
{
    return _indices;
}

auto
CMolecularBasis::max_angular_momentum() const -> int
{
    auto pos = std::ranges::max_element(
        _basis_sets, [&](const auto &lbas, const auto &rbas) { return lbas.max_angular_momentum() < rbas.max_angular_momentum(); });

    return (pos == _basis_sets.end()) ? -1 : pos->max_angular_momentum();
}

auto
CMolecularBasis::max_angular_momentum(const std::vector<int> &atoms) const -> int
{
    auto pos = std::ranges::max_element(atoms, [&](const int i, const int j) {
        return _basis_sets[_indices[i]].max_angular_momentum() < _basis_sets[_indices[j]].max_angular_momentum();
    });

    return (pos == atoms.end()) ? -1 : _basis_sets[_indices[*pos]].max_angular_momentum();
}

auto
CMolecularBasis::basis_functions() const -> std::vector<CBasisFunction>
{
    std::vector<CBasisFunction> bfs;

    std::ranges::for_each(_indices, [&](const int i) {
        auto cbfs = _basis_sets[i].basis_functions();
        bfs.insert(bfs.end(), cbfs.begin(), cbfs.end());
    });

    return bfs;
}

auto
CMolecularBasis::basis_functions(const int angular_momentum) const -> std::vector<CBasisFunction>
{
    std::vector<CBasisFunction> bfs;

    std::ranges::for_each(_indices, [&](const int i) {
        auto cbfs = _basis_sets[i].basis_functions(angular_momentum);
        bfs.insert(bfs.end(), cbfs.begin(), cbfs.end());
    });

    return bfs;
}

auto
CMolecularBasis::basis_functions(const int angular_momentum, const size_t npgtos) const -> std::vector<CBasisFunction>
{
    std::vector<CBasisFunction> bfs;

    std::ranges::for_each(_indices, [&](const int i) {
        auto cbfs = _basis_sets[i].basis_functions(angular_momentum, npgtos);
        bfs.insert(bfs.end(), cbfs.begin(), cbfs.end());
    });

    return bfs;
}

auto
CMolecularBasis::basis_functions(const std::vector<int> &atoms) const -> std::vector<CBasisFunction>
{
    std::vector<CBasisFunction> bfs;

    std::ranges::for_each(atoms, [&](const int i) {
        auto cbfs = _basis_sets[_indices.at(i)].basis_functions();
        bfs.insert(bfs.end(), cbfs.begin(), cbfs.end());
    });

    return bfs;
}

auto
CMolecularBasis::basis_functions(const std::vector<int> &atoms, const int angular_momentum) const -> std::vector<CBasisFunction>
{
    std::vector<CBasisFunction> bfs;

    std::ranges::for_each(atoms, [&](const int i) {
        auto cbfs = _basis_sets[_indices.at(i)].basis_functions(angular_momentum);
        bfs.insert(bfs.end(), cbfs.begin(), cbfs.end());
    });

    return bfs;
}

auto
CMolecularBasis::basis_functions(const std::vector<int> &atoms, const int angular_momentum, const size_t npgtos) const -> std::vector<CBasisFunction>
{
    std::vector<CBasisFunction> bfs;

    std::ranges::for_each(atoms, [&](const int i) {
        auto cbfs = _basis_sets[_indices.at(i)].basis_functions(angular_momentum, npgtos);
        bfs.insert(bfs.end(), cbfs.begin(), cbfs.end());
    });

    return bfs;
}

auto
CMolecularBasis::atomic_indices() const -> std::vector<int>
{
    std::vector<int> atom_indices;

    std::ranges::for_each(std::views::iota(0, static_cast<int>(_indices.size())), [&](const int i) {
        const auto cbfs     = _basis_sets[_indices[i]].basis_functions();
        auto       cindices = std::vector<int>(cbfs.size(), i);
        atom_indices.insert(atom_indices.end(), cindices.begin(), cindices.end());
    });

    return atom_indices;
}

auto
CMolecularBasis::atomic_indices(const int angular_momentum) const -> std::vector<int>
{
    std::vector<int> atom_indices;

    std::ranges::for_each(std::views::iota(0, static_cast<int>(_indices.size())), [&](const int i) {
        const auto cbfs     = _basis_sets[_indices[i]].basis_functions(angular_momentum);
        auto       cindices = std::vector<int>(cbfs.size(), i);
        atom_indices.insert(atom_indices.end(), cindices.begin(), cindices.end());
    });

    return atom_indices;
}

auto
CMolecularBasis::atomic_indices(const int angular_momentum, const size_t npgtos) const -> std::vector<int>
{
    std::vector<int> atom_indices;

    std::ranges::for_each(std::views::iota(0, static_cast<int>(_indices.size())), [&](const int i) {
        const auto cbfs     = _basis_sets[_indices[i]].basis_functions(angular_momentum, npgtos);
        auto       cindices = std::vector<int>(cbfs.size(), i);
        atom_indices.insert(atom_indices.end(), cindices.begin(), cindices.end());
    });

    return atom_indices;
}

auto
CMolecularBasis::atomic_indices(const std::vector<int> &atoms) const -> std::vector<int>
{
    std::vector<int> atom_indices;

    std::ranges::for_each(atoms, [&](const int i) {
        const auto cbfs     = _basis_sets[_indices[i]].basis_functions();
        auto       cindices = std::vector<int>(cbfs.size(), i);
        atom_indices.insert(atom_indices.end(), cindices.begin(), cindices.end());
    });

    return atom_indices;
}

auto
CMolecularBasis::atomic_indices(const std::vector<int> &atoms, const int angular_momentum) const -> std::vector<int>
{
    std::vector<int> atom_indices;

    std::ranges::for_each(atoms, [&](const int i) {
        const auto cbfs     = _basis_sets[_indices[i]].basis_functions(angular_momentum);
        auto       cindices = std::vector<int>(cbfs.size(), i);
        atom_indices.insert(atom_indices.end(), cindices.begin(), cindices.end());
    });

    return atom_indices;
}

auto
CMolecularBasis::atomic_indices(const std::vector<int> &atoms, const int angular_momentum, const size_t npgtos) const -> std::vector<int>
{
    std::vector<int> atom_indices;

    std::ranges::for_each(atoms, [&](const int i) {
        const auto cbfs     = _basis_sets[_indices[i]].basis_functions(angular_momentum, npgtos);
        auto       cindices = std::vector<int>(cbfs.size(), i);
        atom_indices.insert(atom_indices.end(), cindices.begin(), cindices.end());
    });

    return atom_indices;
}

auto
CMolecularBasis::number_of_basis_functions(const int angular_momentum) const -> size_t
{
    return std::accumulate(_indices.begin(), _indices.end(), size_t{0}, [&](const int &sum, const int &i) {
        return sum + _basis_sets[i].number_of_basis_functions(angular_momentum);
    });
}

auto
CMolecularBasis::number_of_basis_functions(const int angular_momentum, const size_t npgtos) const -> size_t
{
    return std::accumulate(_indices.begin(), _indices.end(), size_t{0}, [&](const size_t &sum, const int &i) {
        return sum + _basis_sets[i].number_of_basis_functions(angular_momentum, npgtos);
    });
}

auto
CMolecularBasis::number_of_basis_functions(const std::vector<int> &atoms, const int angular_momentum) const -> size_t
{
    return std::accumulate(atoms.begin(), atoms.end(), size_t{0}, [&](const size_t &sum, const int &i) {
        return sum + _basis_sets[_indices.at(i)].number_of_basis_functions(angular_momentum);
    });
}

auto
CMolecularBasis::number_of_basis_functions(const std::vector<int> &atoms, const int angular_momentum, const size_t npgtos) const -> size_t
{
    return std::accumulate(atoms.begin(), atoms.end(), size_t{0}, [&](const size_t &sum, const int &i) {
        return sum + _basis_sets[_indices.at(i)].number_of_basis_functions(angular_momentum, npgtos);
    });
}

auto
CMolecularBasis::number_of_primitive_functions(const int angular_momentum) const -> size_t
{
    return std::accumulate(_indices.begin(), _indices.end(), size_t{0}, [&](const size_t &sum, const int &i) {
        return sum + _basis_sets[i].number_of_primitive_functions(angular_momentum);
    });
}

auto
CMolecularBasis::number_of_primitive_functions(const std::vector<int> &atoms, const int angular_momentum) const -> size_t
{
    return std::accumulate(atoms.begin(), atoms.end(), size_t{0}, [&](const size_t &sum, const int &i) {
        return sum + _basis_sets[_indices.at(i)].number_of_primitive_functions(angular_momentum);
    });
}

auto
CMolecularBasis::contraction_depths(const int angular_momentum) const -> std::set<size_t>
{
    std::set<size_t> depths;

    std::ranges::for_each(_basis_sets, [&](const auto &abas) {
        auto cdepths = abas.contraction_depths(angular_momentum);
        depths.insert(cdepths.begin(), cdepths.end());
    });

    return depths;
}

auto
CMolecularBasis::contraction_depths(const std::vector<int> &atoms, const int angular_momentum) const -> std::set<size_t>
{
    std::set<size_t> depths;

    std::ranges::for_each(atoms, [&](const int i) {
        auto cdepths = _basis_sets[_indices.at(i)].contraction_depths(angular_momentum);
        depths.insert(cdepths.begin(), cdepths.end());
    });

    return depths;
}

auto
CMolecularBasis::dimensions_of_basis() const -> size_t
{
    return dimensions_of_basis(max_angular_momentum() + 1);
}

auto
CMolecularBasis::dimensions_of_basis(const int angular_momentum) const -> size_t
{
    size_t naos = 0;

    std::ranges::for_each(std::views::iota(0, angular_momentum),
                          [&](const int i) { naos += number_of_basis_functions(i) * tensor::number_of_spherical_components(std::array<int, 1>{i}); });

    return naos;
}

auto
CMolecularBasis::dimensions_of_primitive_basis() const -> size_t
{
    size_t npaos = 0;

    std::ranges::for_each(std::views::iota(0, max_angular_momentum() + 1), [&](const int i) {
        npaos += number_of_primitive_functions(i) * tensor::number_of_spherical_components(std::array<int, 1>{i});
    });

    return npaos;
}

auto
CMolecularBasis::index_map(const int angular_momentum, const size_t npgtos) const -> std::vector<size_t>
{
    std::vector<size_t> ao_indices;

    ao_indices.push_back(number_of_basis_functions(angular_momentum));

    auto offset = dimensions_of_basis(angular_momentum);

    auto bf_index = [&](const auto &bf) {
        if (bf.number_of_primitive_functions() == npgtos) ao_indices.push_back(offset);
        offset++;
    };

    std::ranges::for_each(_indices, [&](const int i) { std::ranges::for_each(_basis_sets[i].basis_functions(angular_momentum), bf_index); });

    return ao_indices;
}

auto
CMolecularBasis::index_map(const std::vector<int> &atoms, const int angular_momentum, const size_t npgtos) const -> std::vector<size_t>
{
    std::vector<size_t> ao_indices;

    ao_indices.push_back(number_of_basis_functions(angular_momentum));

    auto offset = dimensions_of_basis(angular_momentum);

    std::ranges::for_each(std::views::iota(0, static_cast<int>(_indices.size())), [&](const int i) {
        bool not_found = true;
        std::ranges::for_each(atoms, [&](const int j) {
            if (j == i)
            {
                std::ranges::for_each(_basis_sets[_indices[i]].basis_functions(angular_momentum), [&](const auto &bf) {
                    if (bf.number_of_primitive_functions() == npgtos) ao_indices.push_back(offset);
                    offset++;
                });
                not_found = false;
            }
        });
        if (not_found) offset += _basis_sets[_indices[i]].number_of_basis_functions(angular_momentum);
    });

    return ao_indices;
}

auto
CMolecularBasis::main_basis_label() const -> std::string
{
    auto mlabels = _labels_frequency_map();

    auto pos = std::ranges::max_element(mlabels, [&](const auto &lhs, const auto rhs) { return lhs.second < rhs.second; });

    return (pos == mlabels.end()) ? std::string() : pos->first;
}

auto
CMolecularBasis::_labels_frequency_map() const -> std::unordered_map<std::string, int>
{
    std::unordered_map<std::string, int> freqmap;

    std::ranges::for_each(_indices, [&](const int i) { freqmap[_basis_sets[i].get_name()]++; });

    return freqmap;
}

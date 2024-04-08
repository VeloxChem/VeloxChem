#include "MolecularBasis.hpp"

#include <sstream>

#include "AngularMomentum.hpp"
#include "ChemicalElement.hpp"
#include "StringFormat.hpp"

CMolecularBasis::CMolecularBasis(const std::vector<CAtomBasis>& basis_sets, const std::vector<int64_t>& indexes)

    : _basis_sets(basis_sets)

    , _indexes(indexes)
{
}

auto
CMolecularBasis::add(const CAtomBasis& basis) -> void
{
    if (const auto nbases = static_cast<int64_t>(_basis_sets.size()); nbases > 0)
    {
        for (int64_t i = 0; i < nbases; i++)
        {
            if ((_basis_sets[i].getName() == basis.getName()) && (_basis_sets[i].getIdentifier() == basis.getIdentifier()))
            {
                _indexes.push_back(i);

                return;
            }
        }

        _indexes.push_back(nbases);

        _basis_sets.push_back(basis);
    }
    else
    {
        _indexes.push_back(0);

        _basis_sets.push_back(basis);
    }
}

auto
CMolecularBasis::reduceToValenceBasis() const -> CMolecularBasis
{
    if (const auto nbases = _basis_sets.size(); nbases > 0)
    {
        std::vector<CAtomBasis> rbasis_sets;

        for (size_t i = 0; i < nbases; i++)
        {
            rbasis_sets.push_back(_basis_sets[i].reduceToValenceBasis());
        }

        return CMolecularBasis(rbasis_sets, _indexes);
    }
    else
    {
        return CMolecularBasis();
    }
}

auto
CMolecularBasis::slice(const std::vector<int64_t>& atoms) const -> CMolecularBasis
{
    CMolecularBasis mbasis;

    for (const auto atom : atoms)
    {
        mbasis.add(_basis_sets[_indexes[atom]]);
    }

    return mbasis;
}

auto
CMolecularBasis::getBasisSets() const -> std::vector<CAtomBasis>
{
    return _basis_sets;
}

auto
CMolecularBasis::getBasisSetsIndexes() const -> std::vector<int64_t>
{
    return _indexes;
}

auto
CMolecularBasis::getMaxAngularMomentum() const -> int64_t
{
    if (const auto nbases = _basis_sets.size(); nbases > 0)
    {
        auto mang = _basis_sets[0].getMaxAngularMomentum();

        for (size_t i = 1; i < nbases; i++)
        {
            if (const auto angmom = _basis_sets[i].getMaxAngularMomentum(); angmom > mang)
            {
                mang = angmom;
            }
        }

        return mang;
    }
    else
    {
        return -1;
    }
}

auto
CMolecularBasis::getMaxAngularMomentum(const std::vector<int64_t>& atoms) const -> int64_t
{
    if (const auto natoms = atoms.size(); natoms > 0)
    {
        auto mang = _basis_sets[_indexes[atoms[0]]].getMaxAngularMomentum();

        for (size_t i = 1; i < natoms; i++)
        {
            if (const auto angmom = _basis_sets[_indexes[atoms[i]]].getMaxAngularMomentum(); angmom > mang)
            {
                mang = angmom;
            }
        }

        return mang;
    }
    else
    {
        return -1;
    }
}

auto
CMolecularBasis::getBasisFunctions() const -> std::vector<CBasisFunction>
{
    if (const auto nindexes = _indexes.size(); nindexes > 0)
    {
        std::vector<CBasisFunction> gtos;

        for (size_t i = 0; i < nindexes; i++)
        {
            for (const auto& gto : _basis_sets[_indexes[i]].getBasisFunctions())
            {
                gtos.push_back(gto);
            }
        }

        return gtos;
    }
    else
    {
        return std::vector<CBasisFunction>();
    }
}

auto
CMolecularBasis::getBasisFunctions(const int64_t angmom) const -> std::vector<CBasisFunction>
{
    if (const auto nindexes = _indexes.size(); nindexes > 0)
    {
        std::vector<CBasisFunction> gtos;

        for (size_t i = 0; i < nindexes; i++)
        {
            for (const auto& gto : _basis_sets[_indexes[i]].getBasisFunctions(angmom))
            {
                gtos.push_back(gto);
            }
        }

        return gtos;
    }
    else
    {
        return std::vector<CBasisFunction>();
    }
}

auto
CMolecularBasis::getBasisFunctions(const int64_t angmom, const int64_t npgtos) const -> std::vector<CBasisFunction>
{
    if (const auto nindexes = _indexes.size(); nindexes > 0)
    {
        std::vector<CBasisFunction> gtos;

        for (size_t i = 0; i < nindexes; i++)
        {
            for (const auto& gto : _basis_sets[_indexes[i]].getBasisFunctions(angmom, npgtos))
            {
                gtos.push_back(gto);
            }
        }

        return gtos;
    }
    else
    {
        return std::vector<CBasisFunction>();
    }
}

auto
CMolecularBasis::getBasisFunctions(const std::vector<int64_t>& atoms) const -> std::vector<CBasisFunction>
{
    if (const auto natoms = atoms.size(); natoms > 0)
    {
        std::vector<CBasisFunction> gtos;

        for (size_t i = 0; i < natoms; i++)
        {
            for (const auto& gto : _basis_sets[_indexes[atoms[i]]].getBasisFunctions())
            {
                gtos.push_back(gto);
            }
        }

        return gtos;
    }
    else
    {
        return std::vector<CBasisFunction>();
    }
}

auto
CMolecularBasis::getBasisFunctions(const std::vector<int64_t>& atoms, const int64_t angmom) const -> std::vector<CBasisFunction>
{
    if (const auto natoms = atoms.size(); natoms > 0)
    {
        std::vector<CBasisFunction> gtos;

        for (size_t i = 0; i < natoms; i++)
        {
            for (const auto& gto : _basis_sets[_indexes[atoms[i]]].getBasisFunctions(angmom))
            {
                gtos.push_back(gto);
            }
        }

        return gtos;
    }
    else
    {
        return std::vector<CBasisFunction>();
    }
}

auto
CMolecularBasis::getBasisFunctions(const std::vector<int64_t>& atoms, const int64_t angmom, const int64_t npgtos) const -> std::vector<CBasisFunction>
{
    if (const auto natoms = atoms.size(); natoms > 0)
    {
        std::vector<CBasisFunction> gtos;

        for (size_t i = 0; i < natoms; i++)
        {
            for (const auto& gto : _basis_sets[_indexes[atoms[i]]].getBasisFunctions(angmom, npgtos))
            {
                gtos.push_back(gto);
            }
        }

        return gtos;
    }
    else
    {
        return std::vector<CBasisFunction>();
    }
}

auto
CMolecularBasis::getAtomicIndexes() const -> std::vector<int64_t>
{
    if (const auto nindexes = _indexes.size(); nindexes > 0)
    {
        std::vector<int64_t> atom_indexes;

        for (size_t i = 0; i < nindexes; i++)
        {
            if (const auto ngtos = (_basis_sets[_indexes[i]].getBasisFunctions()).size(); ngtos > 0)
            {
                for (size_t j = 0; j < ngtos; j++)
                    atom_indexes.push_back(i);
            }
        }

        return atom_indexes;
    }
    else
    {
        return std::vector<int64_t>();
    }
}

auto
CMolecularBasis::getAtomicIndexes(const int64_t angmom) const -> std::vector<int64_t>
{
    if (const auto nindexes = _indexes.size(); nindexes > 0)
    {
        std::vector<int64_t> atom_indexes;

        for (size_t i = 0; i < nindexes; i++)
        {
            if (const auto ngtos = (_basis_sets[_indexes[i]].getBasisFunctions(angmom)).size(); ngtos > 0)
            {
                for (size_t j = 0; j < ngtos; j++)
                    atom_indexes.push_back(i);
            }
        }

        return atom_indexes;
    }
    else
    {
        return std::vector<int64_t>();
    }
}

auto
CMolecularBasis::getAtomicIndexes(const int64_t angmom, const int64_t npgtos) const -> std::vector<int64_t>
{
    if (const auto nindexes = _indexes.size(); nindexes > 0)
    {
        std::vector<int64_t> atom_indexes;

        for (size_t i = 0; i < nindexes; i++)
        {
            if (const auto ngtos = (_basis_sets[_indexes[i]].getBasisFunctions(angmom, npgtos)).size(); ngtos > 0)
            {
                for (size_t j = 0; j < ngtos; j++)
                    atom_indexes.push_back(i);
            }
        }

        return atom_indexes;
    }
    else
    {
        return std::vector<int64_t>();
    }
}

auto
CMolecularBasis::getAtomicIndexes(const std::vector<int64_t>& atoms) const -> std::vector<int64_t>
{
    if (const auto natoms = atoms.size(); natoms > 0)
    {
        std::vector<int64_t> atom_indexes;

        for (size_t i = 0; i < natoms; i++)
        {
            if (const auto ngtos = (_basis_sets[_indexes[atoms[i]]].getBasisFunctions()).size(); ngtos > 0)
            {
                for (size_t j = 0; j < ngtos; j++)
                    atom_indexes.push_back(atoms[i]);
            }
        }

        return atom_indexes;
    }
    else
    {
        return std::vector<int64_t>();
    }
}

auto
CMolecularBasis::getAtomicIndexes(const std::vector<int64_t>& atoms, const int64_t angmom) const -> std::vector<int64_t>
{
    if (const auto natoms = atoms.size(); natoms > 0)
    {
        std::vector<int64_t> atom_indexes;

        for (size_t i = 0; i < natoms; i++)
        {
            if (const auto ngtos = (_basis_sets[_indexes[atoms[i]]].getBasisFunctions(angmom)).size(); ngtos > 0)
            {
                for (size_t j = 0; j < ngtos; j++)
                    atom_indexes.push_back(atoms[i]);
            }
        }

        return atom_indexes;
    }
    else
    {
        return std::vector<int64_t>();
    }
}

auto
CMolecularBasis::getAtomicIndexes(const std::vector<int64_t>& atoms, const int64_t angmom, const int64_t npgtos) const -> std::vector<int64_t>
{
    if (const auto natoms = atoms.size(); natoms > 0)
    {
        std::vector<int64_t> atom_indexes;

        for (size_t i = 0; i < natoms; i++)
        {
            if (const auto ngtos = (_basis_sets[_indexes[atoms[i]]].getBasisFunctions(angmom, npgtos)).size(); ngtos > 0)
            {
                for (size_t j = 0; j < ngtos; j++)
                    atom_indexes.push_back(atoms[i]);
            }
        }

        return atom_indexes;
    }
    else
    {
        return std::vector<int64_t>();
    }
}

auto
CMolecularBasis::getNumberOfBasisFunctions(const int64_t angmom) const -> int64_t
{
    if (const auto nindexes = _indexes.size(); nindexes > 0)
    {
        int64_t ncgtos = 0;

        for (size_t i = 0; i < nindexes; i++)
        {
            ncgtos += _basis_sets[_indexes[i]].getNumberOfBasisFunctions(angmom);
        }

        return ncgtos;
    }
    else
    {
        return 0;
    }
}

auto
CMolecularBasis::getNumberOfBasisFunctions(const int64_t angmom, const int64_t npgtos) const -> int64_t
{
    if (const auto nindexes = _indexes.size(); nindexes > 0)
    {
        int64_t ncgtos = 0;

        for (size_t i = 0; i < nindexes; i++)
        {
            ncgtos += _basis_sets[_indexes[i]].getNumberOfBasisFunctions(angmom, npgtos);
        }

        return ncgtos;
    }
    else
    {
        return 0;
    }
}

auto
CMolecularBasis::getNumberOfBasisFunctions(const std::vector<int64_t>& atoms, const int64_t angmom) const -> int64_t
{
    if (const auto natoms = atoms.size(); natoms > 0)
    {
        int64_t ncgtos = 0;

        for (size_t i = 0; i < natoms; i++)
        {
            ncgtos += _basis_sets[_indexes[atoms[i]]].getNumberOfBasisFunctions(angmom);
        }

        return ncgtos;
    }
    else
    {
        return 0;
    }
}

auto
CMolecularBasis::getNumberOfBasisFunctions(const std::vector<int64_t>& atoms, const int64_t angmom, const int64_t npgtos) const -> int64_t
{
    if (const auto natoms = atoms.size(); natoms > 0)
    {
        int64_t ncgtos = 0;

        for (size_t i = 0; i < natoms; i++)
        {
            ncgtos += _basis_sets[_indexes[atoms[i]]].getNumberOfBasisFunctions(angmom, npgtos);
        }

        return ncgtos;
    }
    else
    {
        return 0;
    }
}

auto
CMolecularBasis::getNumberOfPrimitiveFunctions(const int64_t angmom) const -> int64_t
{
    if (const auto nindexes = _indexes.size(); nindexes > 0)
    {
        int64_t npgtos = 0;

        for (size_t i = 0; i < nindexes; i++)
        {
            npgtos += _basis_sets[_indexes[i]].getNumberOfPrimitiveFunctions(angmom);
        }

        return npgtos;
    }
    else
    {
        return 0;
    }
}

auto
CMolecularBasis::getNumberOfPrimitiveFunctions(const std::vector<int64_t>& atoms, const int64_t angmom) const -> int64_t
{
    if (const auto natoms = atoms.size(); natoms > 0)
    {
        int64_t npgtos = 0;

        for (size_t i = 0; i < natoms; i++)
        {
            npgtos += _basis_sets[_indexes[atoms[i]]].getNumberOfPrimitiveFunctions(angmom);
        }

        return npgtos;
    }
    else
    {
        return 0;
    }
}

auto
CMolecularBasis::getContractionDepths(const int64_t angmom) const -> std::set<int64_t>
{
    if (const auto nbases = _basis_sets.size(); nbases > 0)
    {
        std::set<int64_t> depths;

        for (size_t i = 0; i < nbases; i++)
        {
            for (const auto npgtos : _basis_sets[i].getContractionDepths(angmom))
            {
                depths.insert(npgtos);
            }
        }

        return depths;
    }
    else
    {
        return std::set<int64_t>();
    }
}

auto
CMolecularBasis::getContractionDepths(const std::vector<int64_t>& atoms, const int64_t angmom) const -> std::set<int64_t>
{
    if (const auto natoms = atoms.size(); natoms > 0)
    {
        std::set<int64_t> depths;

        for (size_t i = 0; i < natoms; i++)
        {
            for (const auto npgtos : _basis_sets[_indexes[atoms[i]]].getContractionDepths(angmom))
            {
                depths.insert(npgtos);
            }
        }

        return depths;
    }
    else
    {
        return std::set<int64_t>();
    }
}

auto
CMolecularBasis::getDimensionsOfBasis() const -> int64_t
{
    return getDimensionsOfBasis(getMaxAngularMomentum() + 1);
}

auto
CMolecularBasis::getDimensionsOfBasis(const int64_t angmom) const -> int64_t
{
    if (angmom > 0)
    {
        int64_t naos = 0;

        for (int64_t i = 0; i < angmom; i++)
        {
            naos += getNumberOfBasisFunctions(i) * angmom::to_SphericalComponents(i);
        }

        return naos;
    }
    else
    {
        return 0;
    }
}

auto
CMolecularBasis::getDimensionsOfPrimitiveBasis() const -> int64_t
{
    if (const auto mang = getMaxAngularMomentum(); mang >= 0)
    {
        int64_t npaos = 0;

        for (int64_t i = 0; i <= mang; i++)
        {
            npaos += getNumberOfPrimitiveFunctions(i) * angmom::to_SphericalComponents(i);
        }

        return npaos;
    }
    else
    {
        return 0;
    }
}

auto
CMolecularBasis::getIndexMap(const int64_t angmom, const int64_t npgtos) const -> std::vector<int64_t>
{
    if (const auto nindexes = _indexes.size(); nindexes > 0)
    {
        std::vector<int64_t> ao_indexes;

        ao_indexes.push_back(getNumberOfBasisFunctions(angmom));

        int64_t offset = getDimensionsOfBasis(angmom);

        for (size_t i = 0; i < nindexes; i++)
        {
            for (const auto& gto : _basis_sets[_indexes[i]].getBasisFunctions(angmom))
            {
                if (gto.getNumberOfPrimitiveFunctions() == npgtos)
                {
                    ao_indexes.push_back(offset);
                }

                offset++;
            }
        }

        return ao_indexes;
    }
    else
    {
        return std::vector<int64_t>();
    }
}

auto
CMolecularBasis::getIndexMap(const std::vector<int64_t>& atoms, const int64_t angmom, const int64_t npgtos) const -> std::vector<int64_t>
{
    if (const auto natoms = atoms.size(); natoms > 0)
    {
        std::vector<int64_t> ao_indexes;

        ao_indexes.push_back(getNumberOfBasisFunctions(angmom));

        int64_t offset = getDimensionsOfBasis(angmom);

        for (size_t i = 0; i < _indexes.size(); i++)
        {
            bool not_found = true;

            for (size_t j = 0; j < natoms; j++)
            {
                if (atoms[j] == i)
                {
                    for (const auto& gto : _basis_sets[_indexes[i]].getBasisFunctions(angmom))
                    {
                        if (gto.getNumberOfPrimitiveFunctions() == npgtos)
                        {
                            ao_indexes.push_back(offset);
                        }

                        offset++;
                    }

                    not_found = false;
                }
            }

            if (not_found)
            {
                offset += _basis_sets[_indexes[i]].getNumberOfBasisFunctions(angmom);
            }
        }

        return ao_indexes;
    }
    else
    {
        return std::vector<int64_t>();
    }
}

void
CMolecularBasis::setLabel(const std::string& label)
{
    _label = label;
}

std::string
CMolecularBasis::getLabel() const
{
    return _label;
}

auto
CMolecularBasis::printBasis(const std::string& title) const -> std::string
{
    std::string str = "Molecular Basis (" + title + ")";

    std::stringstream ss;

    ss << fstr::format(str, 60, fmt_t::center) << "\n";

    str = std::string(str.size() + 2, '=');

    ss << fstr::format(str, 60, fmt_t::center) << "\n\n";

    if (const auto freqmap = _getLabelsFrequencyMap(); !freqmap.empty())
    {
        // determine main basis set in molecular basis

        std::string mlabel;

        int64_t mcount = 0;

        for (const auto& [label, count] : freqmap)
        {
            if (count > mcount)
            {
                mcount = count;

                mlabel = label;
            }
        }

        // print main basis set information

        str.assign("Basis: ");

        str.append(mlabel);

        ss << fstr::format(str, 60, fmt_t::left) << "\n\n";

        ss << "  Atom ";

        ss << fstr::format(std::string("Contracted GTOs"), 26, fmt_t::left);

        ss << fstr::format(std::string("Primitive GTOs"), 30, fmt_t::left);

        ss << "\n\n";

        for (const auto& basis : _basis_sets)
        {
            if (basis.getName() == mlabel)
            {
                std::string lbl("  ");

                if (CChemicalElement elem; elem.setAtomType(basis.getIdentifier()))
                {
                    lbl.append(elem.getName());
                }

                ss << fstr::format(lbl, 6, fmt_t::left);

                ss << fstr::format(basis.getContractionString(), 26, fmt_t::left);

                ss << fstr::format(basis.getPrimitivesString(), 30, fmt_t::left);

                ss << "\n";
            }
        }

        ss << "\n";

        // print additional basis set(s) information

        if (const auto nbases = freqmap.size(); nbases > 1)
        {
            if (nbases > 2)
            {
                ss << fstr::format(std::string("Replacement Basis Sets:"), 84, fmt_t::left) << "\n\n";
            }
            else
            {
                ss << fstr::format(std::string("Replacement Basis Set:"), 84, fmt_t::left) << "\n\n";
            }

            ss << "  Atom No. ";

            ss << fstr::format(std::string("Basis Set"), 20, fmt_t::left);

            ss << fstr::format(std::string("Contracted GTOs"), 26, fmt_t::left);

            ss << fstr::format(std::string("Primitive GTOs"), 30, fmt_t::left);

            ss << "\n\n";

            for (size_t i = 0; i < _indexes.size(); i++)
            {
                if (const auto label = _basis_sets[_indexes[i]].getName(); label != mlabel)
                {
                    ss << "  " << fstr::format(std::to_string(i + 1), 5, fmt_t::left);

                    if (CChemicalElement elem; elem.setAtomType(_basis_sets[_indexes[i]].getIdentifier()))
                    {
                        ss << fstr::format(elem.getName(), 4, fmt_t::left);
                    }

                    ss << fstr::format(label, 20, fmt_t::left);

                    ss << fstr::format(_basis_sets[_indexes[i]].getContractionString(), 26, fmt_t::left);

                    ss << fstr::format(_basis_sets[_indexes[i]].getPrimitivesString(), 30, fmt_t::left);

                    ss << "\n";
                }
            }

            ss << "\n";
        }

        // print molecular basis size

        str.assign("Contracted Basis Functions : ");

        str.append(std::to_string(getDimensionsOfBasis()));

        ss << fstr::format(str, 60, fmt_t::left) << "\n";

        str.assign("Primitive Basis Functions  : ");

        str.append(std::to_string(getDimensionsOfPrimitiveBasis()));

        ss << fstr::format(str, 60, fmt_t::left) << "\n";

        ss << "\n";
    }

    return ss.str();
}

auto
CMolecularBasis::printBasis() const -> std::string
{
    return printBasis("Atomic Basis");
}

auto
CMolecularBasis::_getLabelsFrequencyMap() const -> std::unordered_map<std::string, int64_t>
{
    std::unordered_map<std::string, int64_t> freqmap;

    for (const auto index : _indexes)
    {
        freqmap[_basis_sets[index].getName()]++;
    }

    return freqmap;
}

#ifndef Molecule_hpp
#define Molecule_hpp

#include <set>
#include <string>
#include <vector>

#include "Point.hpp"

/// @brief Class CMolecule stores geometrical data of molecule and provides set
/// of methods for handling of this data.
class CMolecule
{
   public:
    /// @brief The default constructor.
    CMolecule();

    /// @brief The constructor with vector of chemical elements and vector of atoms coordinates.
    /// @param identifiers The vector of chemical element identifiers.
    /// @param coordinates The vector of Cartesian coordinates of atoms.
    /// @param unit The unit used to define coordinates of atoms.
    CMolecule(const std::vector<int> &identifiers, const std::vector<TPoint<double>> &coordinates, const std::string &unit);

    /// @brief The constructor with two molecular fragments.
    /// @param molecule_one The first molecule to merge.
    /// @param molecule_two The second molecule to merge.
    CMolecule(const CMolecule &molecule_one, const CMolecule &molecule_two);

    /// @brief The default copy constructor.
    /// @param other The molecule to be copied.
    CMolecule(const CMolecule &other);

    /// @brief The default move constructor.
    /// @param other The molecule to be moved.
    CMolecule(CMolecule &&other) noexcept;

    /// @brief The default copy assignment operator.
    /// @param other the molecule to be copy assigned.
    /// @return The assigned molecule.
    auto operator=(const CMolecule &other) -> CMolecule &;

    /// @brief The default move assignment operator.
    /// @param other The molecule to be move assigned.
    /// @return The assigned molecule.
    auto operator=(CMolecule &&other) noexcept -> CMolecule &;

    /// @brief The equality operator.
    /// @param other The molecule to be compared.
    /// @return True if molecules are equal, False otherwise.
    auto operator==(const CMolecule &other) const -> bool;

    /// @brief The equality operator.
    /// @param other The molecule to be compared.
    /// @return True if molecules are not equal, False otherwise.
    auto operator!=(const CMolecule &other) const -> bool;

    /// @brief Adds atom to molecule.
    /// @param identifier The chemical element identifier.
    /// @param coordinates The coordinates of atom.
    /// @param unit The unit used to define coordinates of atoms.
    auto add_atom(const int identifier, const TPoint<double> &coordinates, const std::string &unit) -> void;

    /// @brief Slices given set of atoms into new molecule.
    /// @param atoms The vector of atom indices to be slinced.
    /// @return The molcule constructed from sliced atoms.
    auto slice(const std::vector<int> &atoms) const -> CMolecule;

    /// @brief Sets charge of molecule.
    /// @param charge The charge of molecule.
    auto set_charge(const double charge) -> void;

    /// @brief Sets spin multiplicity of molecule.
    /// @param multiplicity The multiplicity (2S+1) of molecule.
    auto set_multiplicity(const int multiplicity) -> void;

    /// @brief Gets charge of molecule.
    /// @return The charge of molecule.
    auto get_charge() const -> double;

    /// @brief Gets spin multiplicity of molecule.
    /// @return The multiplicity of molecule.
    auto get_multiplicity() const -> int;

    /// @brief Gets total number of atoms in molecule.
    /// @return The number of atoms in molecule.
    auto number_of_atoms() const -> int;

    /// @brief Gets number of atoms belonging to specific chemical element in
    /// molecule.
    /// @param identifier The chemical element identifier.
    /// @return The number of atoms in molecule.
    auto number_of_atoms(const int identifier) const -> int;

    /// @brief Gets number of atoms belonging to specific chemical element and in
    /// predefined list of atoms in molecule.
    /// @param iatom The index of first atom in list of atoms.
    /// @param natoms The number of atoms in list of atoms.
    /// @param identifier The chemical element identifier.
    /// @return The number of atoms in molecule.
    auto number_of_atoms(const int iatom, const int natoms, const int identifier) const -> int;

    /// @brief Gets set of unique chemical element identifiers in molecule.
    /// @return The set unique chemical element identifiers.
    auto elemental_composition() const -> std::set<int>;

    /// @brief Gets a number of electrons in molecule.
    /// @return The number of electrons.
    auto number_of_electrons() const -> int;

    /// @brief Gets vector of chemical element identifiers.
    /// @return The vector og chemical element identifiers.
    auto identifiers() const -> std::vector<int>;

    /// @brief Gets vector Cartesian coordinates of atoms in molecule.
    /// @param unit The unit used to define coordinates of atoms.
    /// @return The vector of atom coordinates.
    auto coordinates(const std::string &unit = std::string("au")) const -> std::vector<TPoint<double>>;

    /// @brief Gets charges of all atoms in molecule.
    /// @return The vector of atomic charges of molecule.
    auto charges() const -> std::vector<double>;

    /// @brief Gets masses of all atoms in molecule.
    /// @return The vector of atomic masses of molecule.
    auto masses() const -> std::vector<double>;

    /// @brief Gets labels of all atoms in molecule.
    /// @return The vector of atomic labels of molecule.
    auto labels() const -> std::vector<std::string>;

    /// @brief Gets label of specific atom.
    /// @param iatom The index of atom.
    /// @return The label of atom.
    auto label(const int iatom) const -> std::string;

    /// @brief Gets coordinates of specific atom.
    /// @param iatom The index of atom.
    /// @param unit The unit used to define coordinates of atoms.
    /// @return The coordinates of atom.
    auto atom_coordinates(const int iatom, const std::string &unit = std::string("au")) const -> TPoint<double>;

    /// @brief Gets indices of atoms with given atomic label.
    /// @param label The label of requested atom type.
    /// @return The vector of atom indices.
    auto atom_indices(const std::string &label) const -> std::vector<int>;

    /// @brief Gets nuclear repulsion energy for molecule assuming point charge
    /// model for nucleus.
    /// @return The nuclear repulsion energy.
    auto nuclear_repulsion_energy() const -> double;

    /// @brief Checks if any pair of atoms in molecule is closer than given
    /// minimal distance.
    /// @param distance The minimal distance between two atoms.
    /// @return True if non of atom pairs fail proximity check, False otherwise.
    auto check_proximity(const double distance) const -> bool;

    /// @brief Computes vector of distances to closest neighbouring atom for each
    /// atom.
    /// @return The vector of distances between atoms.
    auto min_distances() const -> std::vector<double>;

   private:
    /// @brief The charge of molecule.
    double _charge{0.0};

    /// @brief The multiplicity of electronic ground state.
    int _multiplicity{1};

    /// @brief The vector of Cartesian coordinates of atoms.
    std::vector<TPoint<double>> _coordinates;

    /// @brief The vector of chemical element identifiers of atoms.
    std::vector<int> _identifiers;

    /// @brief Checks if coordinates of atoms are given in Angstrom.
    /// @param unit The unit used to define coordinates of atoms.
    /// @return True if unit is Angstrom, False otherwise.
    auto _is_angstrom(const std::string &unit) const -> bool;
};

#endif /* Molecule_hpp */

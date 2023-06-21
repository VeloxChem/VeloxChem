#ifndef Molecule_hpp
#define Molecule_hpp

#include <cstdint>

#include <vector>
#include <string>
#include <set>
#include <array>

#include "Point.hpp"

/**
 Class CMolecule stores geometrical data  of molecule and provides set of methods
 for handling of this data.

 @author Z. Rinkevicius
 */
class CMolecule
{
    /**
     The charge of molecule.
     */
    double _charge{0.0};

    /**
     The multiplicity of electronic ground state.
     */
    int64_t _multiplicity{1};

    /**
     The vector of Cartesian coordinates of atoms.
     */
    std::vector<TPoint3D> _coordinates;

    /**
     The vector of chemical element identifiers of atoms.
     */
    std::vector<int64_t> _identifiers;
    
    /**
     Checks if units are given in Angstroms.

     @param units the units of Cartesian coordinates of atoms.
     @return True if units are Angstroms, False otherwise.
     */
    auto
    _isAnngstroms(const std::string& units) const -> bool;

   public:
    
    /**
     Creates an empty molecule.
     */
    CMolecule() = default;

    /**
     Creates a molecule.

     @param identifiers the vector of chemical element identifiers.
     @param coordinates the vector of Cartesian coordinates of atoms.
     @param units the units of Cartesian coordinates of atoms.
     */
    CMolecule(const std::vector<int64_t>&  identifiers,
              const std::vector<TPoint3D>& coordinates,
              const std::string&           units);
    
    /**
     Creates a molecule.

     @param labels the vector of chemical element labels.
     @param coordinates the vector of Cartesian coordinates of atoms.
     @param units the units of Cartesian coordinates of atoms.
     */
    CMolecule(const std::vector<std::string>& labels,
              const std::vector<TPoint3D>&    coordinates,
              const std::string&              units);
    
    /**
     Creates a molecule by merging two molecular fragments.

     @param molfrag_one the first molecular fragment to merge.
     @param molfrag_two the second molecular fragment to merge.
     */
    CMolecule(const CMolecule& molfrag_one,
              const CMolecule& molfrag_two);

    /**
     Adds atom to molecule using given atom label and coordinates.

     @param label the label of atom.
     @param coordinates the coordinates of atom.
     @param units the units of Cartesian coordinates of atom.
    */
    auto
    addAtom(const std::string& label,
            const TPoint3D&    coordinates,
            const std::string& units) -> void;
    
    /**
     Adds atom to molecule using given chemical element identifier and coordinates.

     @param identifier the chemical element identifier of atom.
     @param coordinates the coordinates of atom.
     @param units the units of Cartesian coordinates of atom.
    */
    auto
    addAtom(const int64_t      identifier,
            const TPoint3D&    coordinates,
            const std::string& units) -> void;
    
    /**
     Sets charge of molecule.

     @param charge the charge of molecule.
     */
    auto
    setCharge(const double charge) -> void;

    /**
     Sets spin multiplicity of molecule.

     @param multiplicity the multiplicity (2S+1) of molecule.
     */
    auto
    setMultiplicity(const int64_t multiplicity) -> void;
    
    /**
     Gets charge of molecule.

     @return the charge of molecule.
     */
    auto
    getCharge() const -> double;

    /**
     Gets spin multiplicity of molecule.

     @return the multiplicity of molecules.
     */
    auto
    getMultiplicity() const -> int64_t;

    /**
     Gets total number of atoms in molecule.

     @return the total number of atoms.
     */
    auto
    getNumberOfAtoms() const -> int64_t;

    /**
     Gets number of atoms belonging to specific chemical element in molecule.

     @param identifier the chemical element number.
     @return the number of atoms.
     */
    auto
    getNumberOfAtoms(const int64_t identifier) const -> int64_t;

    /**
    Gets number of atoms belonging to specific chemical element in list of atoms
    in molecule.

    @param iatom the index of first atom in list of atoms.
    @param natoms the number of atoms in list of atoms.
    @param identifier the chemical element number.
    @return the number of atoms.
    */
    auto
    getNumberOfAtoms(const int64_t iatom,
                     const int64_t natoms,
                     const int64_t identifier) const -> int64_t;
    
    /**
     Gets set of unique chemical elements in molecule.

     @return the set of unique chemical elements.
     */
    auto
    getElementalComposition() const -> std::set<int64_t>;

    /**
     Gets a number of electrons in molecule.

     @return the number of electrons.
     */
    auto
    getNumberOfElectrons() const -> int64_t;

    /**
     Gets constant pointer to vector of chemical element identifiers.

     @return the vector of chemical element identifiers.
     */
    auto
    getIdsElemental() const -> std::vector<int64_t>;

    /**
     Gets vector Cartesian  coordinates of atoms in molecule.

     @return the vector of Cartesian coordinates.
     */
    auto
    getCoordinates(const std::string& units) const -> std::vector<TPoint3D>;

    /**
     Gets charges of all atoms in molecule.

     @return the vector of atomic charges.
     */
    auto
    getCharges() const -> std::vector<double>;

    /**
     Gets masses of all atoms in molecule.

     @return the vector of atomic masses.
     */
    auto
    getMasses() const -> std::vector<double>;

    /**
     Gets atomic of all atoms in molecule.

     @return the vector of atomic labels.
     */
    auto
    getLabels() const -> std::vector<std::string>;
    
    /**
     Gets label of specific atom.

     @param iatom the index of atom.
     @return the label of atom.
     */
    auto
    getLabel(const int64_t iatom) const -> std::string;
    
    /**
     Gets vector of atom coordinates of specific atom.

     @param iatom the index of atom.
     @return the vector of atom indexes.
    */
    auto
    getAtomCoordinates(const int64_t      iatom,
                       const std::string& units = std::string("au")) const -> TPoint3D;
    
    /**
     Gets indexes of atoms with given atomic label.

     @param label the label of requested atom type.
     @return the vector of atom indexes.
    */
    auto
    getAtomIndexes(const std::string& label) const -> std::vector<int64_t>;

    /**
     Gets nuclear repulsion energy for molecule assuming point charge model for
     nucleus.

     @return the nuclear repulsion energy.
     */
    auto
    getNuclearRepulsionEnergy() const -> double;
    
    /**
     Checks if any pair of atoms in molecule is closer than given minimal
     distance. Prints error message to output stream for first pair of atoms,
     which are to close.

     @param distance the minimal distance.
     @return true if proximity condition is not violated, false otherwise.
     */
    auto
    checkProximity(const double distance) const -> bool;
    
    /**
     Prints geometry of molecule as table to output stream.

     @return the output string.
     */
     auto
     printGeometry() const -> std::string;
};

#endif /* Molecule_hpp */

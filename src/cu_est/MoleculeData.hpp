#ifndef MoleculeData_hpp
#define MoleculeData_hpp

#include <cstdint>
#include <vector>

#include "Molecule.hpp"

/// @brief Class CMoleculeData stores molecular data in format used in cuEST library.
class CMoleculeData {
    
public:
    /// @brief The default constructor.
    CMoleculeData() {};

    /// @brief The  constructor with molecule.
    /// @param molecule The molecule to create molecule data.
    CMoleculeData(const CMolecule& molecule);
    
private:
    /// @brief The vecotr of number of AO shells per atom.
    std::vector<double> _coords;
    
    /// @brief The vecotr of number of AO shells per atom.
    std::vector<double> _coords;
};

#endif /* MoleculeData_hpp */

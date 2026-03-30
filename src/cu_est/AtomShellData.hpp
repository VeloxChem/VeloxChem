#ifndef AtomShellData_hpp
#define AtomShellData_hpp

#include <cstdint>
#include <vector>

#include "Molecule.hpp"
#include "MolecularBasis.hpp"

#ifdef ENABLE_CUEST
#include "cuest.h"
#endif

/// @brief Class CAtomShellData stores AOs data in format used in cuEST library.
class CAtomShellData {
    
public:
    /// @brief The default constructor.
    CAtomShellData() {};

#ifdef ENABLE_CUEST
    /// @brief The  constructor with cuEST handle, molecule, and molecular basis.
    /// @param cuHandle The cuEST handle required for creating AO shell.
    /// @param molecule The molecule to create atom shell data.
    /// @param basis The molecular basis to create atom shell data.
    CAtomShellData(      cuestHandle_t    cuHandle,
                   const CMolecule&       molecule,
                   const CMolecularBasis& basis);
#endif
    
    /// @brief The custom destructor.
    ~CAtomShellData();
    
private:
    /// @brief The vecotr of number of AO shells per atom.
    std::vector<uint64_t> _aoShellsPerAtom;
    
#ifdef ENABLE_CUEST
    /// @brief The vector of pointers to AO shells.
    std::vector<cuestAOShell_t> _aoShells;
#endif
};

#endif /* AtomShellData_hpp */

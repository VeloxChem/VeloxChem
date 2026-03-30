#include "AtomShellData.hpp"

#ifdef ENABLE_CUEST
CAtomShellData::CAtomShellData(      cuestHandle_t    cuHandle,
                               const CMolecule&       molecule,
                               const CMolecularBasis& basis)
{
    if (const auto natoms = molecule.number_of_atoms(); natoms > 0)
    {
        _aoShellsPerAtom = std::vector<uint64_t>(natoms, 0);
        
        _aoShells = std::vector<cuestAOShell_t>(basis.number_of_basis_functions());
    
        cuestAOShellParameters_t aoshell_parameters;

        cuestParametersCreate(CUEST_AOSHELL_PARAMETERS, &aoshell_parameters);
        
        size_t idx = 0;
        
        for (int i = 0; i < natoms; i++)
        {
            const auto bgtos = basis.basis_functions(std::vector<int>({i, }));
            
            _aoShellsPerAtom[i] = (uint64_t) bgtos.size();
            
            for (const auto& bgto : bgtos)
            {
                const auto pexps = bgto.get_exponents();
                
                const auto pnorms = bgto.get_normalization_factors();
                
                cuestAOShellCreate(cuHandle, 1,
                                   (uint64_t) bgto.get_angular_momentum(),
                                   (uint64_t) bgto.number_of_primitive_functions(),
                                   pexps.data(),
                                   pnorms.data(),
                                   aoshell_parameters,
                                   &_aoShells[idx]);
                
                idx++;
            }
        }
        
        cuestParametersDestroy(CUEST_AOSHELL_PARAMETERS, aoshell_parameters);
    }
}
#endif

CAtomShellData::~CAtomShellData()
{
    //for (size_t i = 0; i < _aoShells.size(); i++)
    //{
    //    delete _aoShells[i];
    //}
}

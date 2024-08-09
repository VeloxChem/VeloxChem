#ifndef Codata_hpp
#define Codata_hpp

namespace units {  // units

// CODATA 2018
// https://physics.nist.gov/cuu/Constants/Table/allascii.txt

/// Gets conversion factor from Bohr to Angstrom.
constexpr auto
bohr_in_angstrom() -> double
{
    // Bohr radius: 0.5291 772 109 03 e-10 [m]
    return 0.529177210903;
}

/// Gets conversion factor from Hartree to eV.
constexpr auto
hartree_in_ev() -> double
{
    // Hartree-eV relationship: 27.211 386 245 988
    return 27.211386245988;
}

}  // namespace units

#endif /* Codata_hpp */

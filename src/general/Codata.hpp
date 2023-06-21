#ifndef Codata_hpp
#define Codata_hpp

namespace units {  // units namespace

// CODATA 2018
// https://physics.nist.gov/cuu/Constants/Table/allascii.txt

/**
 Gets Bohr value in Angstroms.

 @return the conversion factor.
 */
inline auto
getBohrValueInAngstroms() -> double
{
    // Bohr radius: 0.5291 772 109 03 e-10 [m]

    return 0.529177210903;
};

/**
 Gets Hartree value in electronvolts.

 @return the conversion factor.
 */
inline auto
getHartreeValueInElectronVolts() -> double
{
    // hartree-electron volt relationship: 27.211 386 245 988

    return 27.211386245988;
};

}  // namespace units

#endif /* Codata_hpp */

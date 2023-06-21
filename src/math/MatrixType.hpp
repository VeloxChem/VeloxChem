#ifndef MatrixType_hpp
#define MatrixType_hpp

/**
 Enumerate class mat_t:

 Defines supported matrix types:
 mat::symm   - the symmetric square matrix
 mat::antisymm - the antisymmetric square matrix
 mat::gen - the general square or rectangular matrix
 */
enum class mat_t
{
    symm,
    antisymm,
    gen
};

#endif /* MatrixType_hpp */

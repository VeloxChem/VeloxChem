//
//  Tabula — custom-recursion molecular-integral machinery.
//  Basis-function block data, as the fused two-center kernels consume it.
//

#ifndef TabulaKernelBlockData_hpp
#define TabulaKernelBlockData_hpp

namespace tabula {  // tabula namespace

/// @brief A view of one basis-function block's primitive data, as a fused
/// two-center kernel consumes it — raw pointers into arrays the driver
/// fetched once per block. Operator-agnostic: shared by every two-center
/// kernel (overlap, kinetic, …).
struct KernelBlockData
{
    /// @brief The primitive exponents, primitive-major: `[prim·ncgtos + cgto]`.
    const double *exponents;
    /// @brief The primitive normalization factors, same layout.
    const double *norms;
    /// @brief The x coordinate of each contracted GTO.
    const double *x;
    /// @brief The y coordinate of each contracted GTO.
    const double *y;
    /// @brief The z coordinate of each contracted GTO.
    const double *z;
    /// @brief The number of contracted GTOs in the block.
    int ncgtos;
    /// @brief The number of primitives per contracted GTO.
    int nprims;
};

}  // namespace tabula

#endif /* TabulaKernelBlockData_hpp */

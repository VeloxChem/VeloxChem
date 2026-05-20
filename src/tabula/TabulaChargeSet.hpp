//
//  Tabula — custom-recursion molecular-integral machinery.
//  A set of point charges an external-center integral is evaluated against —
//  the molecule's own nuclei, or external charges (a QM/MM embedding).
//

#ifndef TabulaChargeSet_hpp
#define TabulaChargeSet_hpp

namespace tabula {  // tabula namespace

/// @brief A pointer view of the point charges an external-center kernel sums
/// over — raw pointers into arrays the driver owns for the whole compute.
struct ChargeSet
{
    /// @brief The charge magnitudes `Z_N`.
    const double *magnitudes;
    /// @brief The x / y / z coordinate of each charge.
    const double *x, *y, *z;
    /// @brief The number of charges.
    int count;
};

}  // namespace tabula

#endif /* TabulaChargeSet_hpp */

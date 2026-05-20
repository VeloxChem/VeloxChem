//
//  Tabula — custom-recursion molecular-integral machinery.
//  A set of point dipoles a charge-dipole integral is evaluated against — the
//  polarising dipoles of a QM/MM environment.
//

#ifndef TabulaDipoleSet_hpp
#define TabulaDipoleSet_hpp

namespace tabula {  // tabula namespace

/// @brief A pointer view of the point dipoles a charge-dipole kernel sums over
/// — positions and per-site moment components, into arrays the driver owns for
/// the whole compute. The charge-dipole counterpart of `ChargeSet`.
struct DipoleSet
{
    /// @brief The x / y / z coordinate of each dipole site.
    const double *x, *y, *z;
    /// @brief The x / y / z component of each site's moment `d_N`.
    const double *mx, *my, *mz;
    /// @brief The number of dipoles.
    int count;
};

}  // namespace tabula

#endif /* TabulaDipoleSet_hpp */

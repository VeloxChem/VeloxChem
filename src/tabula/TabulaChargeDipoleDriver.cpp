//
//  Tabula — custom-recursion molecular-integral machinery.
//  Two-center charge-dipole driver — supplies the charge-dipole kernel over the
//  shared external-center core, and builds the point-dipole set.
//

#include "TabulaChargeDipoleDriver.hpp"

#include <cmath>
#include <limits>

#include "TabulaChargeDipoleKernel.hpp"
#include "TabulaDipoleSet.hpp"

namespace tabula {  // tabula namespace

namespace {  // unnamed namespace

// per-thread load balance of the most recent charge-dipole compute
ThreadBalance g_balance;

/// @brief Owning storage for a dipole set — the `DipoleSet` view borrows it.
struct DipoleArrays
{
    std::vector<double> x, y, z;     // site positions
    std::vector<double> mx, my, mz;  // site dipole components

    auto view() const -> DipoleSet
    {
        return DipoleSet{x.data(), y.data(), z.data(), mx.data(), my.data(), mz.data(), static_cast<int>(x.size())};
    }

    /// @brief `Σ_N |d_N|` — the factor a screening estimate would scale by.
    auto charge_factor() const -> double
    {
        double s = 0.0;
        for (std::size_t n = 0; n < x.size(); n++)
        {
            s += std::sqrt(mx[n] * mx[n] + my[n] * my[n] + mz[n] * mz[n]);
        }
        return s;
    }
};

/// @brief Point dipoles (moments + au coordinates) as a dipole set.
auto
from_dipoles(const std::vector<std::array<double, 3>>& moments, const std::vector<std::array<double, 3>>& coordinates)
    -> DipoleArrays
{
    DipoleArrays a;
    const auto   n = coordinates.size();
    a.x.reserve(n);
    a.y.reserve(n);
    a.z.reserve(n);
    a.mx.reserve(n);
    a.my.reserve(n);
    a.mz.reserve(n);
    for (std::size_t i = 0; i < n; i++)
    {
        a.x.push_back(coordinates[i][0]);
        a.y.push_back(coordinates[i][1]);
        a.z.push_back(coordinates[i][2]);
        a.mx.push_back(moments[i][0]);
        a.my.push_back(moments[i][1]);
        a.mz.push_back(moments[i][2]);
    }
    return a;
}

/// @brief Screening estimate — keep-all for now (sound: never under-shoots, so
/// the block-sparse path matches dense). A tighter dipole-field bound is a
/// later optimization; until then a positive threshold keeps every block.
auto
dipole_estimate(int, int, const AtomSpan&, const AtomSpan&, double, double) -> double
{
    return std::numeric_limits<double>::max();
}

/// @brief The external-center operator for a dipole set — the charge-dipole
/// kernel with the dipoles bound into the closure. The view `ds` borrows
/// `dipoles`, which must outlive the compute (both are locals of the caller).
auto
make_op(const DipoleArrays& dipoles) -> ExternalCenterOperator
{
    const DipoleSet ds = dipoles.view();
    return ExternalCenterOperator{
        [ds](int l_a, int l_c, const KernelBlockData& bra, int bra_begin, int bra_end, const KernelBlockData& ket,
             double* spherical) { charge_dipole_kernel(l_a, l_c, bra, bra_begin, bra_end, ket, ds, spherical); },
        dipole_estimate, dipoles.charge_factor()};
}

}  // namespace

auto
ChargeDipoleDriver::compute(const CMolecule&                          molecule,
                            const CMolecularBasis&                    basis,
                            const std::vector<std::array<double, 3>>& moments,
                            const std::vector<std::array<double, 3>>& coordinates,
                            const double                              threshold,
                            KernelProfile                            *profile) const -> DenseMatrix
{
    const auto dipoles = from_dipoles(moments, coordinates);
    const auto op      = make_op(dipoles);
    return external_center_compute(molecule, basis, threshold, op, profile, g_balance);
}

auto
ChargeDipoleDriver::computeSparse(const CMolecule&                          molecule,
                                  const CMolecularBasis&                    basis,
                                  const std::vector<std::array<double, 3>>& moments,
                                  const std::vector<std::array<double, 3>>& coordinates,
                                  const double                              threshold,
                                  KernelProfile                            *profile) const -> BlockSparseMatrix
{
    const auto dipoles = from_dipoles(moments, coordinates);
    const auto op      = make_op(dipoles);
    return external_center_compute_sparse(molecule, basis, threshold, op, profile, g_balance);
}

auto
charge_dipole_thread_balance() -> ThreadBalance
{
    return g_balance;
}

}  // namespace tabula

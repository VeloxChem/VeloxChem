//
//  Tabula — custom-recursion molecular-integral machinery.
//  Two-center charge-dipole driver — supplies the charge-dipole kernel over the
//  shared external-center core, and builds the point-dipole set.
//

#include "TabulaChargeDipoleDriver.hpp"

#include <cmath>

#include "MathConst.hpp"
#include "TabulaChargeDipoleFieldKernel.hpp"
#include "TabulaChargeDipoleKernel.hpp"
#include "TabulaDipoleSet.hpp"

namespace tabula {  // tabula namespace

namespace {  // unnamed namespace

// per-thread load balance of the most recent charge-dipole compute
ThreadBalance g_balance;

// Auto-screen defaults for the dense path: a molecule of at least this many
// atoms gets a conservative screen by default, so the dense `compute` does not
// sum every (mostly negligible) shell-pair block over all dipoles on a large,
// extended system.
constexpr int    auto_screen_min_atoms = 100;
constexpr double auto_screen_threshold = 1.0e-12;

/// @brief Resolves the effective screening threshold. A negative `threshold`
/// (the default) selects automatically — a conservative screen for large
/// molecules, exact (0) otherwise; `0` forces exact dense, `> 0` is explicit.
/// `auto_screen_threshold` is sound and effectively lossless.
auto
resolve_threshold(const CMolecule& molecule, const double threshold) -> double
{
    if (threshold >= 0.0) return threshold;
    return (molecule.number_of_atoms() >= auto_screen_min_atoms) ? auto_screen_threshold : 0.0;
}

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

/// @brief A deliberately conservative upper bound on the charge-dipole
/// magnitude of the shell-pair sub-block (bra span `A`, ket span `B`, separated
/// by `r²`), scaled by `charge_factor = Σ_N |d_N|`.
///
/// The field `(a|(r−N)/|r−N|³|c)` is `∇_N` of the nuclear potential, which in
/// the recursion reaches one order higher (the `+1` Q-shift) with the same
/// near-field form. So this is the nuclear bound (`(2π/ζ)·(2ζ_max·r)^m·F_m ≤ 1`,
/// bra–ket overlap decay, spherical prefactors) with the geometric order raised
/// by one — sound (over-counts, never under-shoots), if loose.
auto
dipole_estimate(const int l_a, const int l_c, const AtomSpan& A, const AtomSpan& B, const double r2, const double charge_factor)
    -> double
{
    const int    m        = l_a + l_c + 1;  // one higher than nuclear (the ∇_Q shift)
    const double r        = std::sqrt(r2);
    const double zeta_min = A.exp_min + B.exp_min;
    const double zeta_max = A.exp_max + B.exp_max;
    const double rho_min  = A.exp_min * B.exp_min / (A.exp_min + B.exp_min);  // slowest decay

    double estimate = charge_factor * A.cmax * B.cmax;

    const double geom = 2.0 * zeta_max * r;
    for (int t = 0; t < m; t++) estimate *= geom;

    estimate *= 2.0 * mathconst::pi_value() / zeta_min;  // point-source kernel, F_m ≤ 1
    estimate *= std::exp(-rho_min * r2);                  // bra–ket overlap decay

    const double ia = 0.5 / A.exp_min;
    for (int t = 0; t < l_a; t++) estimate *= ia;

    const double ib = 0.5 / B.exp_min;
    for (int t = 0; t < l_c; t++) estimate *= ib;

    return estimate;
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
    return external_center_compute(molecule, basis, resolve_threshold(molecule, threshold), op, profile, g_balance);
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
ChargeDipoleDriver::computeField(const CMolecule&                          molecule,
                                 const CMolecularBasis&                    basis,
                                 const DenseMatrix&                        density,
                                 const std::vector<std::array<double, 3>>& coordinates,
                                 const double                              threshold) const
    -> std::vector<std::array<double, 3>>
{
    std::vector<double> px, py, pz;
    px.reserve(coordinates.size());
    py.reserve(coordinates.size());
    pz.reserve(coordinates.size());
    for (const auto& c : coordinates)
    {
        px.push_back(c[0]);
        py.push_back(c[1]);
        pz.push_back(c[2]);
    }
    const int np = static_cast<int>(coordinates.size());

    return external_center_field_compute(molecule, basis, density, &charge_dipole_field_kernel, dipole_estimate,
                                         px.data(), py.data(), pz.data(), np, threshold, g_balance);
}

auto
charge_dipole_thread_balance() -> ThreadBalance
{
    return g_balance;
}

}  // namespace tabula

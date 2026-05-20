//
//  Tabula — custom-recursion molecular-integral machinery.
//  Two-center nuclear-attraction driver — supplies the nuclear-attraction
//  kernel and screening estimate over the shared external-center core, and
//  builds the point-charge set.
//

#include "TabulaNuclearAttractionDriver.hpp"

#include <cmath>
#include <string>

#include "MathConst.hpp"
#include "TabulaChargeSet.hpp"
#include "TabulaNuclearAttractionKernel.hpp"

namespace tabula {  // tabula namespace

namespace {  // unnamed namespace

// per-thread load balance of the most recent nuclear-attraction compute
ThreadBalance g_balance;

// Auto-screen defaults for the dense path: a molecule of at least this many
// atoms gets a conservative screen by default, so the dense `compute` does not
// form every (mostly negligible) atom-pair block on a large, extended system.
constexpr int    auto_screen_min_atoms = 100;
constexpr double auto_screen_threshold = 1.0e-12;

/// @brief Resolves the effective screening threshold. A negative `threshold`
/// (the default) selects automatically — a conservative screen for large
/// molecules, exact (0) otherwise; `0` forces exact dense, `> 0` is explicit.
/// `auto_screen_threshold` is sound and ~1e-14 relative, so the auto-screened
/// dense matrix is effectively lossless.
auto
resolve_threshold(const CMolecule& molecule, const double threshold) -> double
{
    if (threshold >= 0.0) return threshold;
    return (molecule.number_of_atoms() >= auto_screen_min_atoms) ? auto_screen_threshold : 0.0;
}

/// @brief A deliberately conservative upper bound on the nuclear-attraction
/// magnitude of the shell-pair sub-block (bra span `A`, ket span `B`, separated
/// by `r²`), scaled by `charge_factor = Σ_N |Z_N|`.
///
/// Dense `compute` is the validated default; this gates the optional
/// block-sparse path. It bounds **every** charge's contribution by the
/// worst-case charge sitting at the bra–ket overlap centre — the near-field
/// nuclear seed `(2π/ζ)·(−2ζ)^m·F_m`, with `F_m ≤ 1`, the geometric factor
/// `(2ζ_max·r)^m ≥ (2ζ)^m·R_AC^m` (`R_AC ≤ r`), and the bra–ket overlap decay
/// `exp(−ρ_min·r²)` (the screenable factor) — then multiplies by `Σ|Z_N|`. It
/// over-counts (most charges are far and contribute less), so it never
/// under-shoots; screening on it is sound, if loose.
auto
nuclear_estimate(const int l_a, const int l_c, const AtomSpan& A, const AtomSpan& B, const double r2, const double charge_factor)
    -> double
{
    const int    m        = l_a + l_c;
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

/// @brief Owning storage for a charge set — the `ChargeSet` view borrows it.
struct ChargeArrays
{
    std::vector<double> mag, x, y, z;

    auto view() const -> ChargeSet
    {
        return ChargeSet{mag.data(), x.data(), y.data(), z.data(), static_cast<int>(mag.size())};
    }

    /// @brief `Σ_N |Z_N|` — the charge factor the screening estimate scales by.
    auto charge_factor() const -> double
    {
        double s = 0.0;
        for (const double z : mag) s += std::abs(z);
        return s;
    }
};

/// @brief The molecule's nuclei as a charge set (charges + au coordinates).
auto
from_nuclei(const CMolecule& molecule) -> ChargeArrays
{
    ChargeArrays a;
    a.mag             = molecule.charges();
    const auto coords = molecule.coordinates(std::string("au"));
    a.x.reserve(coords.size());
    a.y.reserve(coords.size());
    a.z.reserve(coords.size());
    for (const auto& p : coords)
    {
        const auto c = p.coordinates();
        a.x.push_back(c[0]);
        a.y.push_back(c[1]);
        a.z.push_back(c[2]);
    }
    return a;
}

/// @brief Explicit external charges as a charge set.
auto
from_external(const std::vector<double>& magnitudes, const std::vector<std::array<double, 3>>& coordinates)
    -> ChargeArrays
{
    ChargeArrays a;
    a.mag = magnitudes;
    a.x.reserve(coordinates.size());
    a.y.reserve(coordinates.size());
    a.z.reserve(coordinates.size());
    for (const auto& c : coordinates)
    {
        a.x.push_back(c[0]);
        a.y.push_back(c[1]);
        a.z.push_back(c[2]);
    }
    return a;
}

/// @brief The external-center operator for a charge set — the nuclear kernel
/// with the charges bound into the closure, the screening estimate, and the
/// charge factor. The view `cs` borrows `charges`, which must outlive the
/// compute (it does — both are locals of the calling driver method).
auto
make_op(const ChargeArrays& charges) -> ExternalCenterOperator
{
    const ChargeSet cs = charges.view();
    return ExternalCenterOperator{
        [cs](int l_a, int l_c, const KernelBlockData& bra, int bra_begin, int bra_end, const KernelBlockData& ket,
             double* spherical) { nuclear_attraction_kernel(l_a, l_c, bra, bra_begin, bra_end, ket, cs, spherical); },
        nuclear_estimate, charges.charge_factor()};
}

}  // namespace

auto
NuclearAttractionDriver::compute(const CMolecule&       molecule,
                                 const CMolecularBasis& basis,
                                 const double           threshold,
                                 KernelProfile         *profile) const -> DenseMatrix
{
    const auto charges = from_nuclei(molecule);
    const auto op      = make_op(charges);
    return external_center_compute(molecule, basis, resolve_threshold(molecule, threshold), op, profile, g_balance);
}

auto
NuclearAttractionDriver::computeSparse(const CMolecule&       molecule,
                                       const CMolecularBasis& basis,
                                       const double           threshold,
                                       KernelProfile         *profile) const -> BlockSparseMatrix
{
    const auto charges = from_nuclei(molecule);
    const auto op      = make_op(charges);
    return external_center_compute_sparse(molecule, basis, threshold, op, profile, g_balance);
}

auto
NuclearAttractionDriver::compute(const CMolecule&                          molecule,
                                 const CMolecularBasis&                    basis,
                                 const std::vector<double>&                magnitudes,
                                 const std::vector<std::array<double, 3>>& coordinates,
                                 const double                              threshold,
                                 KernelProfile                            *profile) const -> DenseMatrix
{
    const auto charges = from_external(magnitudes, coordinates);
    const auto op      = make_op(charges);
    return external_center_compute(molecule, basis, resolve_threshold(molecule, threshold), op, profile, g_balance);
}

auto
NuclearAttractionDriver::computeSparse(const CMolecule&                          molecule,
                                       const CMolecularBasis&                    basis,
                                       const std::vector<double>&                magnitudes,
                                       const std::vector<std::array<double, 3>>& coordinates,
                                       const double                              threshold,
                                       KernelProfile                            *profile) const -> BlockSparseMatrix
{
    const auto charges = from_external(magnitudes, coordinates);
    const auto op      = make_op(charges);
    return external_center_compute_sparse(molecule, basis, threshold, op, profile, g_balance);
}

auto
nuclear_attraction_thread_balance() -> ThreadBalance
{
    return g_balance;
}

}  // namespace tabula

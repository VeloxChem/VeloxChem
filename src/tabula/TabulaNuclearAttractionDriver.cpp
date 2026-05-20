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

}  // namespace

auto
NuclearAttractionDriver::compute(const CMolecule&       molecule,
                                 const CMolecularBasis& basis,
                                 const double           threshold,
                                 KernelProfile         *profile) const -> DenseMatrix
{
    const auto                   charges = from_nuclei(molecule);
    const ExternalCenterOperator op{nuclear_attraction_kernel, nuclear_estimate, charges.charge_factor()};
    return external_center_compute(molecule, basis, threshold, op, charges.view(), profile, g_balance);
}

auto
NuclearAttractionDriver::computeSparse(const CMolecule&       molecule,
                                       const CMolecularBasis& basis,
                                       const double           threshold,
                                       KernelProfile         *profile) const -> BlockSparseMatrix
{
    const auto                   charges = from_nuclei(molecule);
    const ExternalCenterOperator op{nuclear_attraction_kernel, nuclear_estimate, charges.charge_factor()};
    return external_center_compute_sparse(molecule, basis, threshold, op, charges.view(), profile, g_balance);
}

auto
NuclearAttractionDriver::compute(const CMolecule&                          molecule,
                                 const CMolecularBasis&                    basis,
                                 const std::vector<double>&                magnitudes,
                                 const std::vector<std::array<double, 3>>& coordinates,
                                 const double                              threshold,
                                 KernelProfile                            *profile) const -> DenseMatrix
{
    const auto                   charges = from_external(magnitudes, coordinates);
    const ExternalCenterOperator op{nuclear_attraction_kernel, nuclear_estimate, charges.charge_factor()};
    return external_center_compute(molecule, basis, threshold, op, charges.view(), profile, g_balance);
}

auto
NuclearAttractionDriver::computeSparse(const CMolecule&                          molecule,
                                       const CMolecularBasis&                    basis,
                                       const std::vector<double>&                magnitudes,
                                       const std::vector<std::array<double, 3>>& coordinates,
                                       const double                              threshold,
                                       KernelProfile                            *profile) const -> BlockSparseMatrix
{
    const auto                   charges = from_external(magnitudes, coordinates);
    const ExternalCenterOperator op{nuclear_attraction_kernel, nuclear_estimate, charges.charge_factor()};
    return external_center_compute_sparse(molecule, basis, threshold, op, charges.view(), profile, g_balance);
}

auto
nuclear_attraction_thread_balance() -> ThreadBalance
{
    return g_balance;
}

}  // namespace tabula

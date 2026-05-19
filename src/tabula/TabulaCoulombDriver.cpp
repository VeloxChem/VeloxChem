//
//  Tabula — custom-recursion molecular-integral machinery.
//  Two-center Coulomb driver — the Coulomb kernel and screening estimate
//  over the shared two-center core.
//

#include "TabulaCoulombDriver.hpp"

#include <algorithm>
#include <cmath>

#include "MathConst.hpp"
#include "TabulaCoulombKernel.hpp"

namespace tabula {  // tabula namespace

namespace {  // unnamed namespace

// per-thread load balance of the most recent Coulomb compute
ThreadBalance g_balance;

/// @brief A conservative upper bound on the Coulomb magnitude of the
/// shell-pair sub-block formed by bra span `A` and ket span `B`, separated by
/// `r² = r2` — see `docs/screening-predicates.md`.
auto
coulomb_estimate(const int l_a, const int l_c, const AtomSpan &A, const AtomSpan &B, const double r2) -> double
{
    const int    m       = l_a + l_c;
    const double r       = std::sqrt(r2);
    const double rho_min = A.exp_min * B.exp_min / (A.exp_min + B.exp_min);

    double estimate = A.cmax * B.cmax;

    // the geometric r^m factor of the MD recursion
    for (int t = 0; t < m; t++) estimate *= r;

    // the Coulomb kernel bound (2π/ρ_min)·min(1/(2m+1), C_m·T^{−m−1/2}),
    // T = ρ_min·r², C_m = ½·Γ(m+½) — see screening-predicates.md
    const double T          = rho_min * r2;
    const double bounded    = 1.0 / static_cast<double>(2 * m + 1);
    double       asymptotic = 0.5 * std::tgamma(static_cast<double>(m) + 0.5) / std::sqrt(T);
    for (int t = 0; t < m; t++) asymptotic /= T;
    estimate *= (2.0 * mathconst::pi_value() / rho_min) * std::min(bounded, asymptotic);

    const double pe = mathconst::pi_value() / (A.exp_min + B.exp_min);
    estimate *= pe * std::sqrt(pe);

    const double ia = 0.5 / A.exp_min;
    for (int t = 0; t < l_a; t++) estimate *= ia;

    const double ib = 0.5 / B.exp_min;
    for (int t = 0; t < l_c; t++) estimate *= ib;

    return estimate;
}

}  // namespace

auto
CoulombDriver::compute(const CMolecule       &molecule,
                       const CMolecularBasis &basis,
                       const double           threshold,
                       KernelProfile         *profile) const -> DenseMatrix
{
    const TwoCenterOperator op{coulomb_kernel, coulomb_estimate};
    return two_center_compute(molecule, basis, threshold, op, profile, g_balance);
}

auto
CoulombDriver::computeSparse(const CMolecule       &molecule,
                             const CMolecularBasis &basis,
                             const double           threshold,
                             KernelProfile         *profile) const -> BlockSparseMatrix
{
    const TwoCenterOperator op{coulomb_kernel, coulomb_estimate};
    return two_center_compute_sparse(molecule, basis, threshold, op, profile, g_balance);
}

auto
coulomb_thread_balance() -> ThreadBalance
{
    return g_balance;
}

}  // namespace tabula

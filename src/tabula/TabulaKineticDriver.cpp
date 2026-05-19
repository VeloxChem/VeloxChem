//
//  Tabula — custom-recursion molecular-integral machinery.
//  Two-center kinetic-energy driver — the kinetic kernel and screening
//  estimate over the shared two-center core.
//

#include "TabulaKineticDriver.hpp"

#include <cmath>

#include "MathConst.hpp"
#include "TabulaKineticKernel.hpp"

namespace tabula {  // tabula namespace

namespace {  // unnamed namespace

// per-thread load balance of the most recent kinetic compute
ThreadBalance g_balance;

/// @brief A conservative upper bound on the kinetic-energy magnitude of the
/// shell-pair sub-block formed by bra span `A` and ket span `B`, separated by
/// `r² = r2` — see `docs/screening-predicates.md`.
auto
kinetic_estimate(const int l_a, const int l_c, const AtomSpan &A, const AtomSpan &B, const double r2) -> double
{
    const int    m       = l_a + l_c;
    const double r       = std::sqrt(r2);
    const double rho_min = A.exp_min * B.exp_min / (A.exp_min + B.exp_min);
    const double rho_max = A.exp_max * B.exp_max / (A.exp_max + B.exp_max);

    double estimate = A.cmax * B.cmax;

    // the geometric r^m factor of the MD recursion, then the kinetic kernel
    // bound ρ_max·(2m+3 + 2·ρ_max·r²)·e^{−ρ_min·r²} — see screening-predicates.md
    for (int t = 0; t < m; t++) estimate *= r;

    estimate *= rho_max * (static_cast<double>(2 * m + 3) + 2.0 * rho_max * r2);

    estimate *= std::exp(-rho_min * r2);

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
KineticDriver::compute(const CMolecule       &molecule,
                       const CMolecularBasis &basis,
                       const double           threshold,
                       KernelProfile         *profile) const -> DenseMatrix
{
    const TwoCenterOperator op{kinetic_kernel, kinetic_estimate};
    return two_center_compute(molecule, basis, threshold, op, profile, g_balance);
}

auto
KineticDriver::computeSparse(const CMolecule       &molecule,
                             const CMolecularBasis &basis,
                             const double           threshold,
                             KernelProfile         *profile) const -> BlockSparseMatrix
{
    const TwoCenterOperator op{kinetic_kernel, kinetic_estimate};
    return two_center_compute_sparse(molecule, basis, threshold, op, profile, g_balance);
}

auto
kinetic_thread_balance() -> ThreadBalance
{
    return g_balance;
}

}  // namespace tabula

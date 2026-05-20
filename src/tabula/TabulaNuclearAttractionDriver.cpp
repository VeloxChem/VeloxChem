//
//  Tabula — custom-recursion molecular-integral machinery.
//  Two-center nuclear-attraction driver — supplies the nuclear-attraction
//  kernel and screening estimate over the shared external-center core, and
//  builds the point-charge set.
//

#include "TabulaNuclearAttractionDriver.hpp"

#include <limits>
#include <string>

#include "TabulaChargeSet.hpp"
#include "TabulaNuclearAttractionKernel.hpp"

namespace tabula {  // tabula namespace

namespace {  // unnamed namespace

// per-thread load balance of the most recent nuclear-attraction compute
ThreadBalance g_balance;

/// @brief Phase-1 placeholder screening estimate — returns a value above any
/// threshold, so no atom span is dropped. Sound (it never under-estimates); a
/// real point-source estimate replaces it in the screening pass.
auto
nuclear_estimate(int, int, const AtomSpan&, const AtomSpan&, double) -> double
{
    return std::numeric_limits<double>::max();
}

/// @brief Owning storage for a charge set — the `ChargeSet` view borrows it.
struct ChargeArrays
{
    std::vector<double> mag, x, y, z;

    auto view() const -> ChargeSet
    {
        return ChargeSet{mag.data(), x.data(), y.data(), z.data(), static_cast<int>(mag.size())};
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
    const ExternalCenterOperator op{nuclear_attraction_kernel, nuclear_estimate};
    const auto                   charges = from_nuclei(molecule);
    return external_center_compute(molecule, basis, threshold, op, charges.view(), profile, g_balance);
}

auto
NuclearAttractionDriver::computeSparse(const CMolecule&       molecule,
                                       const CMolecularBasis& basis,
                                       const double           threshold,
                                       KernelProfile         *profile) const -> BlockSparseMatrix
{
    const ExternalCenterOperator op{nuclear_attraction_kernel, nuclear_estimate};
    const auto                   charges = from_nuclei(molecule);
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
    const ExternalCenterOperator op{nuclear_attraction_kernel, nuclear_estimate};
    const auto                   charges = from_external(magnitudes, coordinates);
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
    const ExternalCenterOperator op{nuclear_attraction_kernel, nuclear_estimate};
    const auto                   charges = from_external(magnitudes, coordinates);
    return external_center_compute_sparse(molecule, basis, threshold, op, charges.view(), profile, g_balance);
}

auto
nuclear_attraction_thread_balance() -> ThreadBalance
{
    return g_balance;
}

}  // namespace tabula

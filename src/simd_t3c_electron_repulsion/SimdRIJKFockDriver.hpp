//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//
//  Redistribution and use in source and binary forms, with or without modification,
//  are permitted provided that the following conditions are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this
//     list of conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//  3. Neither the name of the copyright holder nor the names of its contributors
//     may be used to endorse or promote products derived from this software without
//     specific prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
//  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
//  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
//  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
//  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
//  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
//  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#ifndef SimdRIJKFockDriver_hpp
#define SimdRIJKFockDriver_hpp

#include <cstddef>
#include <utility>
#include <vector>

#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "PackedMatrix.hpp"
#include "SimdRIFockCommon.hpp"
#include "SimdRIFockDriver.hpp"
#include "SparseTensor.hpp"
#include "TripleSparsityPattern.hpp"

/// @brief The times of the phases of one Fock build of the direct mode.
/// @note The direct mode repeats every phase on every iteration, so a phase which
/// does not widen with the threads bounds the whole calculation however many cores
/// are given to it. Run a calculation at one thread and again at many, with
/// VLX_RIJK_PROFILE set, and the phase whose time does not fall between the two is
/// the one worth working on. The whole is timed as well as the parts, so that what
/// the parts do not account for is visible rather than assumed.
/// @note A build is made of the exchange pass, the fitting and the Coulomb pass,
/// which are three calls so that a communicator can gather the fitting between the
/// first and the last. The times of the three are gathered here and reported once,
/// at the end of the last of them, so that a build still reads as one build. The
/// total is the time of the calls and not of the gap between them, which belongs to
/// whoever is dividing the work.
struct CDirectTimes
{
    double allocate    = 0.0;
    double integrals_a = 0.0;
    double transform   = 0.0;
    double closure     = 0.0;
    double copies      = 0.0;
    double solve       = 0.0;
    double exchange    = 0.0;
    double integrals_b = 0.0;
    double coulomb     = 0.0;
    double total       = 0.0;
};

/// @brief Class CSimdRIJKFockDriver builds the Fock matrices of the resolution of
/// the identity for one molecule and one pair of bases, one matrix per call.
///
/// @note The driver holds what does not change between the calls: the inverted
/// Cholesky factor of the metric of the fitting basis and the B vectors. Both are
/// formed once by prepare. A call then contracts a density into the Y vector and
/// the Coulomb matrix, and transforms the B vectors with the orbitals into the W
/// matrices and adds their exchange.
///
/// @note The W matrices are formed for a range of the auxiliary basis at a time,
/// their exchange is added, and the same storage is reused for the next range.
/// They are rebuilt on every call, as the orbitals change, so holding all of them
/// would cost the memory of the whole auxiliary basis and save nothing. Only the
/// B vectors have to be resident, which is what the memory check is about.
///
/// @note The exchange is added with a factor the caller passes, so that a hybrid
/// functional scales it by its fraction of exact exchange. A pure functional is
/// served by a driver which never forms the B vectors and is not this one.
class CSimdRIJKFockDriver
{
   public:
    /// @brief The default constructor.
    CSimdRIJKFockDriver() = default;

    /// @brief Gets the memory the driver holds for a molecule and its bases.
    /// @param molecule The molecule to compute the Fock matrices of.
    /// @param basis The molecular basis.
    /// @param aux_basis The auxiliary molecular basis.
    /// @param threshold The screening threshold.
    /// @param aux_atoms The atoms of the auxiliary basis to answer for, or none of
    /// them to answer for all of them. A rank of a communicator asks for its own
    /// share, as answering the memory of the whole molecule would put every rank on
    /// the direct way for a calculation each of them holds a fitting share of.
    /// @return The memory of the B vectors in bytes.
    /// @note The sparsity pattern of the B vectors is described to answer this,
    /// which is what the driver would do anyway and is a small part of forming
    /// them, so the answer is the memory they will take rather than an estimate of
    /// it. The W matrices are not counted, as one range of them is held at a time
    /// and is small beside the B vectors.
    /// @param range_separated True for a hybrid range separated functional, whose
    /// build holds a second set of B vectors of the attenuated operator. The two
    /// sets are on one sparsity pattern, so the memory is twice the one set.
    auto required_memory(const CMolecule        &molecule,
                         const CMolecularBasis  &basis,
                         const CMolecularBasis  &aux_basis,
                         const double            threshold,
                         const std::vector<int> &aux_atoms       = {},
                         const bool              range_separated = false) const -> size_t;

    /// @brief Forms the metric a way of building asks for, and the way it is for.
    /// @param molecule The molecule to compute the metric of.
    /// @param aux_basis The auxiliary molecular basis.
    /// @param metric_threshold The threshold below which a direction of the metric
    /// carries nothing and is dropped.
    /// @param use_inverse_square_root Whether to invert the square root of the
    /// metric rather than take the factor or its inverse.
    /// @param mode The way of building the metric is for, which must be named: the
    /// two ways want different metrics, so there is nothing to form for a way which
    /// has not been chosen.
    /// @return The metric, and the way it is for, which is the way asked for. The
    /// direct way takes either metric: it solves the Cholesky factor where it has
    /// one and multiplies by the inverted square root where it was asked for one or
    /// where the fitting basis has no factor to be had.
    /// @note prepare forms the metric with this, and a caller which prepares the
    /// ranks of a communicator forms it once with this and hands it to them, so
    /// that the fallbacks are decided in one place rather than raced for on every
    /// rank.
    auto make_metric(const CMolecule       &molecule,
                     const CMolecularBasis &aux_basis,
                     const double           metric_threshold,
                     const bool             use_inverse_square_root,
                     const rimode           mode) const -> std::pair<CPackedMatrix, rimode>;

    /// @brief Forms the metric of the Coulomb operator and that of the attenuated
    /// one, both inverted, for a hybrid range separated functional.
    /// @param molecule The molecule to compute the metrics of.
    /// @param aux_basis The auxiliary molecular basis.
    /// @param metric_threshold The eigenvalues of a metric at or below which a
    /// direction is dropped, when it is inverted through its square root.
    /// @param use_inverse_square_root True to invert the square roots rather than
    /// the Cholesky factors.
    /// @param mode The way of building the metrics are for, which must be the one
    /// which holds the B vectors.
    /// @param omega The range separation parameter, which must be positive.
    /// @return The inverted metric of the Coulomb operator and that of the
    /// attenuated one.
    /// @note The two operators are formed in one call of the two-center range
    /// separated driver rather than in two sweeps of the fitting basis, and each
    /// matrix is inverted by the route asked for with the same fallback. The way
    /// which answers is not returned, as unlike make_metric above there is only one
    /// it can be: the range separated way holds its B vectors.
    auto make_metric_rs(const CMolecule       &molecule,
                        const CMolecularBasis &aux_basis,
                        const double           metric_threshold,
                        const bool             use_inverse_square_root,
                        const rimode           mode,
                        const double           omega) const -> std::pair<CPackedMatrix, CPackedMatrix>;

    /// @brief Forms the inverted factor of the metric and the B vectors.
    /// @param molecule The molecule to compute the Fock matrices of.
    /// @param basis The molecular basis.
    /// @param aux_basis The auxiliary molecular basis.
    /// @param threshold The screening threshold.
    /// @param memory_budget The memory the driver may hold, in bytes.
    /// @note The memory is checked against the budget before the integrals are
    /// computed, from the sparsity pattern alone, so a calculation which cannot fit
    /// is told so rather than dying in the allocator after minutes of work.
    /// @param metric_threshold The eigenvalues of the metric at or below which a
    /// direction is dropped, when the metric is inverted through its square root.
    /// @param use_inverse_square_root True to invert the square root of the metric
    /// rather than its Cholesky factor.
    /// @note The Cholesky factor is the cheaper of the two by an order of
    /// magnitude and is tried first. A fitting basis which is close to linearly
    /// dependent has no Cholesky factor to invert, and the square root is inverted
    /// instead, with a warning. Setting the flag takes that way from the start.
    /// @param aux_atoms The atoms of the auxiliary basis this driver forms the B
    /// vectors of, or none of them for all of them. This is how the work is divided
    /// over a communicator: each rank is given a share of the atoms and answers a
    /// share of the Fock matrix, which the ranks sum. The direct way refuses a
    /// division, as its triangular solve reaches across the whole auxiliary basis.
    /// @param min_parts The fewest parts the direct way is to sweep the auxiliary
    /// basis in. The parts are cut to fit the memory of a build, and a machine with
    /// memory to spare gives one of them -- which one rank then sweeps for the
    /// Coulomb matrix while every other rank waits. A caller dividing the work over a
    /// communicator asks here for at least as many parts as it has ranks. It costs
    /// nothing to ask: the parts are a division of the same atoms either way, so the
    /// integrals of a sweep are the same integrals however they are grouped.
    /// @param metric The metric to build with, or an empty matrix to form it here.
    /// A metric given must be the one make_metric answers for the mode given, and
    /// the mode must then be named rather than automatic, as the fallbacks which
    /// change the mode have already been taken where the metric was formed.
    /// @param mode Which way the Fock matrices are formed, or automatic to hold the
    /// B vectors when they fit in the budget and to form the integrals again on
    /// every call when they do not.
    /// @note The B vectors of a large molecule do not fit in the memory of any one
    /// machine, and the direct mode is what makes such a molecule reachable. It
    /// forms the integrals once for every batch of occupied orbitals and once more
    /// for the Coulomb matrix, so it is several times slower for each Fock matrix
    /// and asks for a hundredth of the memory. The automatic choice takes the held
    /// form wherever it fits.
    auto prepare(const CMolecule       &molecule,
                 const CMolecularBasis &basis,
                 const CMolecularBasis &aux_basis,
                 const double           threshold,
                 const size_t           memory_budget,
                 const double           metric_threshold        = 1.0e-12,
                 const bool             use_inverse_square_root = false,
                 const rimode           mode                    = rimode::automatic,
                 const std::vector<int> &aux_atoms              = {},
                 const CPackedMatrix    &metric                 = CPackedMatrix(),
                 const size_t            min_parts              = 1,
                 const double            omega                  = 0.0,
                 const CPackedMatrix    &metric_erf             = CPackedMatrix()) -> void;

    /// @brief Computes the Fock matrix of a density and a set of orbitals.
    /// @param density The density matrix, in the packed format, symmetric for a
    /// self consistent field calculation and general for a response one.
    /// @param coefficients The molecular orbital coefficients of the occupied
    /// orbitals, as a general matrix of one row per basis function and one column
    /// per orbital.
    /// @param exchange_scaling_factor The factor the exchange is scaled by, which
    /// is one for the exchange of a field calculation and the fraction of exact
    /// exchange of a hybrid functional. The exchange is not formed at all when it
    /// is zero.
    /// @return The Fock matrix, in the packed format as a symmetric matrix.
    /// @note The matrix is twice the Coulomb matrix less the scaled exchange,
    /// which is the convention of a closed shell calculation, whose density is
    /// that of one spin.
    /// @param erf_exchange_scaling_factor The factor the exchange of the attenuated
    /// operator is scaled by, which is the erf coefficient of a hybrid range
    /// separated functional and zero for every other calculation. It is subtracted
    /// as the plain exchange is, so the caller passes what the four-center way is
    /// passed and the sign is the same.
    /// @note The two exchanges are added inside one pass over the ranges of the
    /// auxiliary basis, and the storage of the W matrices is filled again from the
    /// attenuated B vectors rather than doubled: a range's plain exchange is added
    /// before its attenuated W matrices are formed, so nothing of the first is
    /// needed when the second is being built.
    auto compute(const CPackedMatrix &density,
                 const CPackedMatrix &coefficients,
                 const double         exchange_scaling_factor,
                 const double         erf_exchange_scaling_factor = 0.0) -> CPackedMatrix;

    /// @brief Computes the Fock matrices of the two spins of an open shell.
    /// @param density The total density, which is that of both spins added, in the
    /// packed format as a symmetric matrix.
    /// @param coefficients_alpha The coefficients of the occupied orbitals of the
    /// alpha spin, as a general matrix of one row per basis function and one
    /// column per orbital.
    /// @param coefficients_beta The same for the beta spin, which has a number of
    /// columns of its own: the two spins of an open shell do not occupy the same
    /// number of orbitals.
    /// @param exchange_scaling_factor The factor the exchange is scaled by, which
    /// is one for a field calculation and the fraction of exact exchange of a
    /// hybrid functional. Neither exchange is formed at all when it is zero.
    /// @return The Fock matrix of the alpha spin and that of the beta spin, in the
    /// packed format as symmetric matrices.
    /// @note The Coulomb matrix is formed once from the total density and is **not**
    /// doubled, where the closed shell call above doubles the Coulomb of one spin's
    /// density. Each spin's exchange is then subtracted from its own matrix.
    /// @note The two spins are done inside one pass over the ranges of the
    /// auxiliary basis rather than in two passes, so a range's B vectors are
    /// touched once and serve both. The storage of the W matrices is therefore
    /// two, one per spin, as the two have different numbers of columns, and a
    /// range holds a matrix of the basis by the occupied orbitals of each spin.
    /// @note The way which forms the integrals again on every call is not served
    /// here. Its fitting is accumulated from the integrals during the sweep which
    /// builds the exchange, and two spins there means two exchanges and one
    /// fitting summed over both inside that sweep, which is a different piece of
    /// work from this one and is refused rather than approximated.
    /// @param erf_exchange_scaling_factor The factor the exchange of the attenuated
    /// operator is scaled by, as above. Each spin's attenuated exchange goes into
    /// its own matrix, as its plain exchange does.
    auto compute(const CPackedMatrix &density,
                 const CPackedMatrix &coefficients_alpha,
                 const CPackedMatrix &coefficients_beta,
                 const double         exchange_scaling_factor,
                 const double         erf_exchange_scaling_factor = 0.0) -> std::pair<CPackedMatrix, CPackedMatrix>;

    /// @brief Computes the exchange of a range of the orbitals, and the right hand
    /// side of the fitting it closes on the way.
    /// @param coefficients The molecular orbital coefficients of all the occupied
    /// orbitals, as a general matrix of one row per basis function and one column
    /// per orbital.
    /// @param exchange_scaling_factor The factor the exchange is scaled by.
    /// @param ofirst The first orbital of the range.
    /// @param olast One past the last orbital of the range.
    /// @return The scaled exchange, negated as it enters the Fock matrix, and the
    /// right hand side of the fitting of this range of the orbitals.
    /// @note Every term of both is a sum over the orbitals, so the ranks of a
    /// communicator take a range each and their exchange matrices and right hand
    /// sides add. This is the index the direct way divides over: the auxiliary
    /// basis is not one, as the triangular solve of this pass reaches across all
    /// of it.
    /// @note Only the direct way builds this way, as the way which holds the B
    /// vectors has no pass to divide.
    auto compute_exchange(const CPackedMatrix &coefficients,
                          const double         exchange_scaling_factor,
                          const size_t         ofirst,
                          const size_t         olast) -> std::pair<CPackedMatrix, std::vector<double>>;

    /// @brief Solves the metric against the right hand side of the fitting.
    /// @param gamma The right hand side, of one value per auxiliary basis function,
    /// summed over every orbital there is.
    /// @return The coefficients of the fitting.
    /// @note This is what couples the two passes, and is why they are two calls: the
    /// solve reaches across the whole auxiliary basis and needs a right hand side
    /// which is complete, so the ranks of a communicator have to add theirs
    /// together before it. It costs the square of the auxiliary basis and every rank
    /// holds the factor, so each of them solves it rather than one solving and
    /// sending.
    auto solve_fitting(std::vector<double> gamma) -> std::vector<double>;

    /// @brief Adds the Coulomb matrix of the given parts of the auxiliary basis to
    /// a matrix.
    /// @param gamma The coefficients of the fitting, of one value per auxiliary
    /// basis function.
    /// @param parts The parts to sweep, as their indices, which number_of_parts
    /// bounds. The ranks of a communicator take some each and their matrices add.
    /// @param matrix The matrix to add to, which is the Fock matrix being built.
    /// @note The Coulomb matrix enters twice, as the density is that of one spin,
    /// and this adds it that way.
    auto compute_coulomb(const std::vector<double> &gamma,
                         const std::vector<int>    &parts,
                         CPackedMatrix             &matrix) -> void;

    /// @brief Gets the work each atom of the auxiliary basis carries.
    /// @param molecule The molecule to compute the Fock matrices of.
    /// @param basis The molecular basis.
    /// @param aux_basis The auxiliary molecular basis.
    /// @param threshold The screening threshold.
    /// @return The memory of the B vectors of each atom, in bytes, of one value per
    /// atom of the molecule and zero for an atom the screening leaves nothing of.
    /// @note This is what a communicator should divide by, and not the count of the
    /// atoms. An auxiliary function on an atom in the middle of a molecule survives
    /// screening against far more atom pairs than one on an atom at the edge, so
    /// ranks given equal counts of atoms are given unequal counts of values, and it
    /// is the values which cost. The weights sum to what required_memory answers for
    /// the whole molecule, so a rank may take the memory of its share from them
    /// rather than describing a second pattern.
    auto aux_atom_weights(const CMolecule       &molecule,
                          const CMolecularBasis &basis,
                          const CMolecularBasis &aux_basis,
                          const double           threshold) const -> std::vector<double>;

    /// @brief Gets the number of auxiliary basis functions this driver holds the B
    /// vectors of, which is the whole auxiliary basis unless the atoms were divided.
    /// @return The number of functions a build sweeps.
    /// @note This is what says whether a division divided anything. A rank given a
    /// share of the atoms answers a share of the Fock matrix whether or not the work
    /// was divided -- a function it holds nothing of contributes nothing either way
    /// -- so an energy cannot tell the two apart, and this can.
    auto number_of_aux_functions() const -> size_t;

    /// @brief Gets the number of parts the Coulomb pass of the direct mode is
    /// divided over.
    /// @return The number of parts, which is zero in the mode which holds the B
    /// vectors.
    /// @note This is what compute_coulomb indexes into, and what a caller dividing
    /// the Coulomb pass over the ranks of a communicator deals out.
    auto number_of_parts() const -> size_t;

    /// @brief Gets the number of parts the exchange pass of the direct mode sweeps
    /// the auxiliary basis in.
    /// @return The number of parts, which is zero in the mode which holds the B
    /// vectors.
    /// @note Every rank sweeps every one of these, so they are cut by memory alone
    /// and are not divided over anything. Reported so that a caller can see that
    /// asking for more parts of the Coulomb pass has not added any here.
    auto number_of_sweep_parts() const -> size_t;

    /// @brief Checks that the driver has been prepared.
    /// @return True if the driver is ready to form a Fock matrix.
    auto is_prepared() const -> bool;

    /// @brief Gets the way the driver forms the Fock matrices.
    /// @return The mode, which is never automatic once the driver is prepared.
    auto get_mode() const -> rimode;

    /// @brief Gets the B vectors the driver holds.
    /// @return The B vectors.
    auto get_bq_vectors() const -> const CSparseTensor &;

    /// @brief Gets the B vectors of the attenuated operator.
    /// @return The B vectors, which are empty unless the driver was prepared for a
    /// hybrid range separated functional.
    auto get_bq_vectors_erf() const -> const CSparseTensor &;

    /// @brief Gets the range separation parameter the driver was prepared at.
    /// @return The parameter, or zero where there is no attenuated set. A caller
    /// which holds a driver someone else prepared asks this to know whether it can
    /// serve the functional it has, rather than finding out in the build.
    auto get_omega() const -> double;

    /// @brief Gets the inverted Cholesky factor of the metric the driver holds.
    /// @return The inverted factor.
    auto get_metric() const -> const CPackedMatrix &;

    /// @brief Sets the density of the B vectors at which the exchange half
    /// transformation expands them into a square.
    /// @param threshold The density. A threshold of zero or less expands always,
    /// which is the default and what every calculation has done; one above one
    /// walks the values always; between the two the driver counts the values of
    /// each range and chooses.
    /// @note Forwarded to the driver of the B vectors, which is where the choice
    /// is made. It is exposed here because that driver is held privately and a
    /// calculation reaches only this one.
    auto set_dense_threshold(const double threshold) -> void;

    /// @brief Gets the density at which the exchange half transformation expands
    /// the B vectors.
    /// @return The density.
    auto get_dense_threshold() const -> double;

    /// @brief The fraction of the B vectors this driver holds which is actually
    /// filled, against what a dense tensor of the same dimensions would hold.
    /// @return The density, or zero where the driver holds no B vectors.
    /// @note This is the quantity the exchange half transformation compares against
    /// the threshold. Reported rather than only decided on, so a calculation can say
    /// how sparse its tensor was instead of leaving it to be inferred from a timing.
    auto bq_density() const -> double;

    /// @brief Gets the inverted metric of the attenuated operator.
    /// @return The inverted metric, which is empty unless the driver was prepared
    /// for a hybrid range separated functional.
    /// @note It is a different matrix from the plain one and the two are not
    /// interchangeable: a quantity fitted in one operator's metric and contracted
    /// in the other's is a fitting of neither.
    auto get_metric_erf() const -> const CPackedMatrix &;

   private:
    /// @brief Computes the Fock matrix by forming the integrals again on every
    /// call, holding no B vectors.
    /// @param coefficients The molecular orbital coefficients.
    /// @param exchange_scaling_factor The factor the exchange is scaled by.
    /// @return The Fock matrix, twice the Coulomb less the scaled exchange.
    /// @note The density is not needed and is not taken: the Coulomb matrix of this
    /// way comes from the orbitals, through the fitting they close.
    auto _compute_direct(const CPackedMatrix &coefficients,
                         const double         exchange_scaling_factor) -> CPackedMatrix;

    /// @brief Applies the metric the direct mode holds to a set of right hand
    /// sides, in place.
    /// @param values The right hand sides, as a row major array of one row per
    /// auxiliary basis function and ncols columns, overwritten by the result.
    /// @param nrows The number of auxiliary basis functions.
    /// @param ncols The number of right hand sides.
    /// @param transposed True to apply the transpose, which the Cholesky factor has
    /// and the inverted square root, being symmetric, does not.
    /// @note What is applied depends on which metric the driver was given. The
    /// Cholesky factor is solved against, which is cheaper than multiplying by an
    /// inverse and is why the direct way keeps the factor where it can; the inverted
    /// square root is multiplied by. Both close the same sum: solving the factor
    /// gives B with B^T B equal to A^T V^-1 A, and so does multiplying by the root,
    /// since the root is its own transpose.
    auto _apply_metric(double *values, const size_t nrows, const size_t ncols, const bool transposed) const -> void;

    /// @brief Multiplies a set of right hand sides by the metric, in place.
    /// @param metric The metric, expanded into a row major square.
    /// @param values The right hand sides, of one row per auxiliary basis function
    /// and ncols columns, overwritten by the product.
    /// @param nrows The number of auxiliary basis functions.
    /// @param ncols The number of right hand sides.
    auto _multiply_metric(const double *metric, double *values, const size_t nrows, const size_t ncols) const -> void;

    /// @brief The most values the buffer of a multiplication by the metric may hold.
    /// @note The product cannot be taken in place, and the right hand sides of a
    /// batch are the auxiliary basis by the basis by the orbitals -- gigabytes -- so
    /// a buffer of that size would be added to the peak of a build. The columns are
    /// taken in chunks against a buffer of this size instead, and what the chunking
    /// costs is a copy for each, which is nothing beside the product.
    static constexpr size_t _metric_buffer = size_t{32} * 1024 * 1024;

    /// @brief Divides the auxiliary basis into the parts the direct mode sweeps.
    /// @param molecule The molecule to compute the integrals of.
    /// @param basis The molecular basis.
    /// @param aux_basis The auxiliary molecular basis.
    /// @param threshold The screening threshold.
    /// @param pattern The sparsity pattern of the whole, to measure the parts from.
    /// @param min_parts The fewest parts to cut, whatever the memory allows.
    /// @return The pattern of each part.
    /// @brief Gets the parts the Coulomb pass is divided over.
    /// @return The parts of the Coulomb pass, which are the parts of the sweep unless
    /// a finer division was asked for.
    auto _coulomb_patterns() const -> const std::vector<CTripleSparsityPattern> &;

    auto _make_parts(const CMolecule              &molecule,
                     const CMolecularBasis        &basis,
                     const CMolecularBasis        &aux_basis,
                     const double                  threshold,
                     const CTripleSparsityPattern &pattern,
                     const size_t                  min_parts) const -> std::vector<CTripleSparsityPattern>;

    /// @brief The fewest auxiliary functions a range of the W matrices holds.
    /// @note The range is what the transformation divides over the threads, and a
    /// range fixed at sixty four -- four times the sixteen cores it was chosen on --
    /// left half of a hundred and twenty eight with nothing. Measured on a node of
    /// that width, the way which holds the B vectors was then between two and three
    /// times slower than the direct way at every basis set of two molecules, by a
    /// ratio which did not move with the size of the calculation. This is only the
    /// floor now: the range is taken from the memory below, as what a call costs
    /// beside its tasks is paid once for the call however long the range is, and
    /// seven calls of five hundred and twelve carried 0.3 seconds a build of that
    /// over one call of the whole basis.
    static constexpr size_t _w_batch = 64;

    /// @brief The memory the W matrices of a range may take together.
    /// @note They are the basis by the occupied orbitals, one for each function of
    /// the range, so this and not the threads is what sets the range of a large
    /// calculation -- and the bound is on the calculation rather than on the
    /// machine. Two gigabytes gave a cluster of three hundred and twenty atoms a
    /// hundred and thirty one functions a call, which is a hundred and twenty two
    /// calls of the transformation for every build, and every one of them ends with
    /// an exchange which adds a triangle of the basis for each thread. Sixteen
    /// makes it sixteen calls. The matrices are held once, where the direct mode
    /// holds its half transformed integrals twice, so this is the smaller of the
    /// two claims on the memory of such a calculation.
    /// @note Sixteen gigabytes is the bound on a driver, and a node given to several
    /// ranks holds several drivers, so this alone let fourteen ranks of a laptop ask
    /// for thirty two gigabytes of W matrices between them on a machine with thirty
    /// six. The triangles of the exchange and the copies of the Kohn-Sham matrix are
    /// bounded by the threads as well, which shrink as the ranks of a node grow, and
    /// do not multiply this way; this one is bounded by the auxiliary basis, which
    /// does not. The share of the memory budget below bounds it for that reason.
    static constexpr size_t _w_batch_memory = size_t{16} * 1024 * 1024 * 1024;

    /// @brief The share of the memory budget a range of the W matrices may take.
    /// @note The budget is what one rank may hold, already divided by the ranks
    /// sharing a host, so a share of it is a bound which knows about the machine
    /// where the constant above does not. A quarter leaves the B vectors, which the
    /// budget was measured against, the room they were given.
    static constexpr size_t _w_batch_divisor = 4;

    /// @brief The memory the direct mode is allowed to reach, taken from the budget
    /// the driver was prepared with.
    /// @note The direct mode holds two things of its own: the half transformed
    /// integrals of a batch of orbitals, held twice over, and the integrals of the
    /// part of the auxiliary basis being swept. Both are live at once, so each is
    /// given half of this and the whole stays within it. A larger batch of orbitals
    /// is fewer passes over the integrals and more memory, which is the trade this
    /// sets.
    size_t _budget = size_t{4} * 1024 * 1024 * 1024;

    /// @brief The way the driver forms the Fock matrices.
    rimode _mode = rimode::automatic;

    /// @brief The molecule, which the direct mode forms the integrals of again on
    /// every call.
    CMolecule _molecule;

    /// @brief The parts the Coulomb pass of the direct mode is divided over, which
    /// is empty when that pass divides over the parts of the sweep itself.
    /// @note The two passes want different numbers of parts. The Coulomb pass is
    /// divided over them, so it wants at least one for every rank; the exchange pass
    /// sweeps every one of them on every rank, and each is another call of the
    /// transformation -- which costs a square of the basis allocated and zeroed for
    /// every thread, whatever the part holds. Cutting the sweep finer to balance the
    /// Coulomb pass cost more than the balance was worth at a thousand functions on
    /// two hundred and fifty six threads, so the two are cut apart.
    std::vector<CTripleSparsityPattern> _coulomb_parts;

    /// @brief The sparsity patterns the direct mode sweeps, one for each part of
    /// the auxiliary basis, described once.
    /// @note The blocks of a pattern are described pair block by pair block and
    /// every auxiliary group within one, so a run of them spans every auxiliary
    /// function. Sweeping runs like that would take the half transform of one
    /// function once for every run. The parts hold disjoint atoms of the auxiliary
    /// basis instead, so each function belongs to exactly one of them and is
    /// transformed once.
    std::vector<CTripleSparsityPattern> _parts;

    /// @brief The metric the direct mode builds with, which is the lower triangular
    /// Cholesky factor of it or the inverted square root of it.
    /// @note Which of the two it is, is read from the matrix rather than remembered
    /// beside it: a Cholesky factor is lower triangular and an inverted square root
    /// is symmetric, so the type of the matrix says how it is to be applied and
    /// cannot disagree with it.
    CPackedMatrix _direct_metric;

    /// @brief The molecular basis.
    CMolecularBasis _basis;

    /// @brief The auxiliary molecular basis.
    CMolecularBasis _aux_basis;

    /// @brief The inverted Cholesky factor of the metric of the auxiliary basis.
    CPackedMatrix _metric;

    /// @brief The B vectors.
    CSparseTensor _bq_vectors;

    /// @brief The B vectors of the attenuated operator, held only for a hybrid range
    /// separated functional and empty otherwise.
    /// @note They are on the same sparsity pattern as the plain ones, so an element
    /// of either is at the same place in the other.
    CSparseTensor _bq_vectors_erf;

    /// @brief The inverted metric of the attenuated operator, held with its B
    /// vectors and empty otherwise.
    CPackedMatrix _metric_erf;

    /// @brief The range separation parameter the attenuated set was formed at, or
    /// zero when there is no attenuated set. This is what a build asks to know
    /// whether it may add an attenuated exchange.
    double _omega = 0.0;

    /// @brief The W matrices of one range of the auxiliary basis.
    std::vector<CPackedMatrix> _w_vectors;

    /// @brief The W matrices of the second spin of an open shell.
    /// @note A second storage rather than the one above reused, because the two
    /// spins of an open shell occupy different numbers of orbitals and the matrices
    /// of a range differ in their columns between them. Reusing one would form the
    /// storage again at every range of every build.
    std::vector<CPackedMatrix> _w_vectors_beta;

    /// @brief The driver of the B vectors and of the matrices formed from them.
    CSimdRIFockDriver _drv;

    /// @brief The dense indices of the auxiliary basis functions this driver holds
    /// the B vectors of, in ascending order.
    /// @note The functions of an atom are not consecutive -- the dense index is
    /// ordered by angular momentum across the whole molecule -- so the share of a
    /// rank is a set and not a range, and a build sweeps the set. Sweeping the range
    /// which covers it instead would have every rank form, zero and multiply a
    /// matrix for every function of every other rank, which divides the memory of
    /// the B vectors and nothing else.
    std::vector<size_t> _aux_functions;

    /// @brief The times of the phases of the build being made, gathered over the
    /// calls it is made of and reported at the end of the last of them.
    CDirectTimes _direct_times;

    /// @brief Whether the B vectors have been formed.
    bool _prepared = false;

    /// @brief Whether this rank was dealt no atoms of the auxiliary basis.
    /// @note Such a rank holds no B vectors and its share of every sum over the
    /// auxiliary basis is zero, so it answers matrices of zeros. It is not the same
    /// as a driver which was asked for the whole molecule: that one is handed an
    /// empty list of atoms too, and the two are told apart by how many ranks are
    /// dividing the work.
    bool _holds_nothing = false;
};

#endif /* SimdRIJKFockDriver_hpp */

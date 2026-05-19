//
//  Tabula — custom-recursion molecular-integral machinery.
//  Matrix storage, independent of VeloxChem's own CMatrix.
//

#ifndef TabulaDenseMatrix_hpp
#define TabulaDenseMatrix_hpp

#include <cstddef>
#include <memory>
#include <new>
#include <utility>
#include <vector>

namespace tabula {  // tabula namespace

/// @brief An allocator that default-initializes — for a trivial type, leaves
/// the storage uninitialized — on a no-argument construct, so a `resize` can
/// grow the buffer without the zero-fill. Value-initializing constructs
/// (`vector(n, value)`) are unaffected.
template <typename T>
struct DefaultInitAllocator : std::allocator<T>
{
    DefaultInitAllocator() = default;

    template <typename U>
    DefaultInitAllocator(const DefaultInitAllocator<U>&) noexcept
    {
    }

    template <typename U>
    struct rebind
    {
        using other = DefaultInitAllocator<U>;
    };

    template <typename U>
    void construct(U* ptr)
    {
        ::new (static_cast<void*>(ptr)) U;
    }

    template <typename U, typename... Args>
    void construct(U* ptr, Args&&... args)
    {
        ::new (static_cast<void*>(ptr)) U(std::forward<Args>(args)...);
    }
};

/// @brief The symmetry of a dense matrix.
enum class Symmetry
{
    /// @brief No symmetry — every element stored and independent.
    general,
    /// @brief Symmetric — M(j,i) == M(i,j).
    symmetric,
    /// @brief Antisymmetric — M(j,i) == -M(i,j). The diagonal is independent
    /// storage and is not constrained to zero.
    antisymmetric
};

/// @brief A dense, row-major matrix of doubles — Tabula's own integral-result
/// storage, independent of VeloxChem's block-structured CMatrix.
///
/// A symmetric or antisymmetric matrix still holds the full `rows x columns`
/// buffer (no packing); the symmetry is a tag the evaluator and `symmetrize`
/// honour. Packed and sparse representations are separate, later types.
class DenseMatrix
{
   public:
    /// @brief Creates an empty (0 x 0) general matrix.
    DenseMatrix();

    /// @brief Creates a matrix.
    /// @param rows The number of rows.
    /// @param columns The number of columns.
    /// @param symmetry The matrix symmetry.
    /// @param initialize When true (the default) the buffer is zeroed; when
    /// false it is left uninitialized — the caller must write every element
    /// before reading, which a full unscreened evaluation does.
    DenseMatrix(const std::size_t rows,
                const std::size_t columns,
                const Symmetry    symmetry   = Symmetry::general,
                const bool        initialize = true);

    /// @brief Element access for assignment.
    /// @param row The row index.
    /// @param column The column index.
    /// @return A reference to the addressed element.
    auto operator()(const std::size_t row, const std::size_t column) -> double&;

    /// @brief Element access.
    /// @param row The row index.
    /// @param column The column index.
    /// @return The value of the addressed element.
    auto operator()(const std::size_t row, const std::size_t column) const -> double;

    /// @brief Gets the number of rows.
    /// @return The number of rows.
    auto rows() const -> std::size_t;

    /// @brief Gets the number of columns.
    /// @return The number of columns.
    auto columns() const -> std::size_t;

    /// @brief Gets the matrix symmetry.
    /// @return The matrix symmetry.
    auto symmetry() const -> Symmetry;

    /// @brief Gets a pointer to the row-major value buffer.
    /// @return The value buffer pointer.
    auto values() -> double*;

    /// @brief Gets a constant pointer to the row-major value buffer.
    /// @return The constant value buffer pointer.
    auto values() const -> const double*;

    /// @brief Sets all elements to zero.
    auto zero() -> void;

    /// @brief Scales all elements by a factor.
    /// @param factor The scaling factor.
    auto scale(const double factor) -> void;

    /// @brief Mirrors the upper triangle into the lower one according to the
    /// matrix symmetry — `M(j,i) = M(i,j)` (symmetric) or `-M(i,j)`
    /// (antisymmetric). The diagonal is left untouched. A no-op for a general
    /// matrix. Requires a square matrix.
    auto symmetrize() -> void;

   private:
    /// @brief The number of rows.
    std::size_t _rows;

    /// @brief The number of columns.
    std::size_t _columns;

    /// @brief The matrix symmetry.
    Symmetry _symmetry;

    /// @brief The row-major value buffer, length `_rows * _columns`.
    std::vector<double, DefaultInitAllocator<double>> _values;
};

}  // namespace tabula

#endif /* TabulaDenseMatrix_hpp */

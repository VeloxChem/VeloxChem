#ifndef Matrix_hpp
#define Matrix_hpp

#include <map>
#include <set>
#include <utility>
#include <vector>

#include "SubMatrix.hpp"

/// @brief Enumerate class mat_t: defines supported matrix types.
/// mat::symmetric   - the symmetric square matrix
/// mat::antisymmetric - the antisymmetric square matrix
/// mat::general - the general square or rectangular matrix
enum class mat_t
{
    symmetric,
    antisymmetric,
    general
};

/// @brief Class CMatrix stores matrix data and provides set of methods
/// for handling of matrix data.
class CMatrix
{
   public:
    /// @brief The default constructor.
    CMatrix();

    /// @brief The constructor with dictionary of submatrices and matrix type.
    /// @param sub_matrices The map of submatrices.
    /// @param mat_type The matrix type.
    CMatrix(const std::map<std::pair<int, int>, CSubMatrix> &sub_matrices, const mat_t mat_type);

    /// @brief The constructor with vector of angular pairs, vector of submatrices, and matrix type.
    /// @param ang_pairs The vector of angular pairs.
    /// @param sub_matrices The vector of submatrices.
    /// @param mat_type The matrix type.
    CMatrix(const std::vector<std::pair<int, int>> &ang_pairs, const std::vector<CSubMatrix> &sub_matrices, const mat_t mat_type);

    /// @brief The default copy constructor.
    /// @param other The matrix to be copied.
    CMatrix(const CMatrix &other);

    /// @brief The default move constructor.
    /// @param other The matrix to be moved.
    CMatrix(CMatrix &&other) noexcept;

    /// @brief The default destructor.
    ~CMatrix();

    /// @brief The default copy assignment operator.
    /// @param other The matrix to be copy assigned.
    /// @return The assigned matrix.
    auto operator=(const CMatrix &other) -> CMatrix &;

    /// @brief The default move assignment operator.
    /// @param other The matrix to be move assigned.
    /// @return The assigned matrix.
    auto operator=(CMatrix &&other) noexcept -> CMatrix &;

    /// @brief The equality operator.
    /// @param other The matrix to be compared.
    /// @return True if matrices are equal, False otherwise.
    auto operator==(const CMatrix &other) const -> bool;

    /// @brief The equality operator.
    /// @param other The matrix to be compared.
    /// @return True if matrices are not equal, False otherwise.
    auto operator!=(const CMatrix &other) const -> bool;

    /// @brief The addition operator.
    /// @param other The matrix to be added.
    /// @return The sum of two matrices.
    auto operator+(const CMatrix &other) const -> CMatrix;

    /// @brief Adds submatrix to matrix.
    /// @param sub_matrix The submatrix to be added.
    /// @param angpair The angular pair associated with submatrix to be added.
    auto add(const CSubMatrix &sub_matrix, const std::pair<int, int> &angpair) -> void;

    /// @brief Adds submatrix to matrix.
    /// @param dimensions The dimensions of submatrix to be added.
    /// @param angpair The angular pair associated with submatrix to be added.
    auto add(const std::array<size_t, 4> &dimensions, const std::pair<int, int> &angpair) -> void;

    /// @brief Set type of matrix.
    /// @param mat_type The matrix type.
    auto set_type(const mat_t mat_type) -> void;

    /// @brief Set matrix values to zero.
    auto zero() -> void;
    
    /// @brief Assigns flat vector of values to matrix.
    /// @param values The flat vector of values.
    auto assign_flat_values(const std::vector<double>& values) -> void;

    /// @brief Scales matrix values by factor.
    /// @param factor The factor to scale matrix values.
    auto scale(const double factor) -> void;

    /// @brief Symmetrizes submatrices  on diagonal of matrix.
    auto symmetrize() -> void;

    /// @brief Gets vector of angular pairs for submatrices.
    /// @return The vector of angular pairs.
    auto angular_pairs() const -> std::vector<std::pair<int, int>>;

    /// @brief Gets type of matrix.
    /// @return The type of matrix.
    auto get_type() const -> mat_t;

    /// @brief Gets vector of submatrices.
    /// @return The vector of submatrices.
    auto sub_matrices() const -> std::vector<CSubMatrix>;

    /// @brief Gets pointer to requested submatrix.
    /// @param angpair The angular pair of requested submatrix.
    /// @return The pointer to submatrix.
    auto sub_matrix(const std::pair<int, int> &angpair) -> CSubMatrix *;

    /// @brief Gets constant pointer to requested submatrix.
    /// @param angpair The angular pair of requested submatrix.
    /// @return The constant pointer to submatrix.
    auto sub_matrix(const std::pair<int, int> &angpair) const -> const CSubMatrix *;

    /// @brief Checks if submatrix with requested angular pair is stored in matrix.
    /// @param angpair The angular pair of requested submatrix.
    /// @return True if submatrix is stored in matrix, False otherwise.
    auto is_angular_order(const std::pair<int, int> &angpair) const -> bool;

    /// @brief Gets number of rows in matrix.
    /// @return The number of rows in matrix.
    auto number_of_rows() const -> size_t;

    /// @brief Gets number of columns in matrix.
    /// @return The number of columns in matrix.
    auto number_of_columns() const -> size_t;

    /// @brief Reconstructs full matrix as submatrix with dimensions (0, 0, naos, naos).
    /// @return The reconstructed matrix into submatrix.
    auto full_matrix() const -> CSubMatrix;

    /// @brief Gets pointer to matrix.
    /// @return The  pointer to matrix.
    auto pointer() -> CMatrix *;

    /// @brief Gets constant pointer to matrix.
    /// @return The constant pointer to matrix.
    auto pointer() const -> const CMatrix *;
    
    /// @brief Gets flatened values of matrix elements into vector.
    /// @return The constant pointer to matrix.
    auto flat_values() const -> std::vector<double>;

   private:
    /// @brief The map of submatrices.
    std::map<std::pair<int, int>, CSubMatrix *> _sub_matrices;

    /// @brief The type of matrix.
    mat_t _mat_type;

    /// @brief Gets set of row angular keys.
    /// @return The set of row angular keys.
    auto _row_angular_keys() const -> std::set<int>;

    /// @brief Gets set of columns angular keys.
    /// @return The set of columns angular keys.
    auto _column_angular_keys() const -> std::set<int>;

    /// @brief Deallocates submatrices in map of submatrices.
    auto _deallocate() -> void;
};

#endif /* Matrix_hpp */

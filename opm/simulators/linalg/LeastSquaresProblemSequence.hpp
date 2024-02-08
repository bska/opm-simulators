/*
  Copyright 2024 Equinor ASA

  This file is part of the Open Porous Media Project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef OPM_LINALG_LEAST_SQUARES_PROBLEM_SEQUENCE_HPP
#define OPM_LINALG_LEAST_SQUARES_PROBLEM_SEQUENCE_HPP

#include <cstddef>
#include <vector>

/// \file Facility for defining and solving a sequnce of full-rank linear
/// least-squares problems.  The intended use case is that in which the
/// coefficient matrix rarely changes, if at all, but in which the system
/// right-hand side will change frequently.  We also expect that the size of
/// each individual LLS problem will be fairly small--e.g., the number of
/// rows corresponds to the number of per-cell connections/intersections and
/// the number of columns corresponds to the number of spatial dimensions.
/// To this end, we use the QR method and represent the orthogonal matrix
/// 'Q' in factored form by calculating and storing Householder reflectors.
///
/// This implementation is very much inspired by the LAPACK routines
///
///   -# xORQRF -- Generates a QR factorisation
///   -# xORMQR -- Apply Q or Q^T (== Q^{-1}) to a vector
///   -# xGETRS -- Solve triangular system of linear equations
///
/// and how these manage their internal data.  We refer to
/// https://www.netlib.org/lapack/lug/ for additional details on those
/// routines.

namespace Opm {

    /// Sequence of linear least-squares problems
    template <typename Scalar>
    class LeastSquaresProblemSequence
    {
    public:
        /// LLS coefficient matrix.  Supports per-row construction.
        class CoefficientMatrix
        {
        public:
            /// Constructor
            ///
            /// \param[in] ncol Number of columns in the coefficient matrix.
            explicit CoefficientMatrix(const int ncol)
                : ncol_{ ncol }
            {}

            /// Add a row to the coefficient matrix.
            ///
            /// Values will be copied into internal storage.
            ///
            /// \param[in] row Pointer to start of row data.  Assumed to
            ///   point to the start of at least \c ncol elements.
            void addRow(const Scalar* row)
            {
                this->a_.insert(this->a_.end(), row, row + this->ncol_);
            }

            friend class LeastSquaresProblemSequence;

        private:
            /// Number of columns in the coefficient matrix.
            int ncol_;

            /// Coefficient matrix values.  Row major (i.e., 'C++') order,
            /// meaning A_{ij} is stored in a_[i*ncol_ + j] for zero-based
            /// indices i and j.
            std::vector<Scalar> a_{};
        };

        /// Single LLS problem.
        class LLSProblem
        {
        public:
            friend class LeastSquaresProblemSequence;

            /// Number of rows in this LLS problem.
            int nrow() const { return this->nrow_; }

            /// Number of columns in this LLS problem
            int ncol() const { return this->ncol_; }

            /// Solve LLS problem with specific system right-hand side
            ///
            /// \param[in] b System right-hand side.  Assumed to point to
            ///   the start of at least nrows() elements, the first nrows()
            ///   elements of which constitute the actual RHS values.
            ///
            /// \return LLS solution vector.  Pointer to the start of the
            ///   solution vector, the first ncols() elements of which
            ///   constitute the LLS solution data.  The elements of this
            ///   sequence may be overwritten, and should be copied to
            ///   client storage for subsequent processing.
            const Scalar* solve(const Scalar* b) const;

        private:
            /// Constructor
            ///
            /// For use by LeastSquaresProblemSequence only
            ///
            /// \param[in] nrow Number of rows in LLS coefficient matrix
            ///
            /// \param[in] ncol Number of columns in LLS coefficient matrix
            ///
            /// \param[in] A QR factorisation of LLS coefficient matrix in
            ///    factored form.  Pointer to start of \code nrow * ncol
            ///    \endcode data values in column major order.  Essential
            ///    part of Householder reflectors below diagonal, upper
            ///    triangular matrix 'R' on and above diagonal.
            ///
            /// \param[in] beta Householder reflector scaling data.  Pointer
            ///    to the start of at least \p ncol data elements.
            ///
            /// \param[out] house Work array for reconstituting full
            ///    Householder reflectors.  Pointer to the start of at least
            ///    \p nrow data elements.
            ///
            /// \param[out] work Internal scratch array for holding results
            ///   of intermediate stages.  Starts as a copy of the system
            ///   RHS, continues as the result of applying the Householder
            ///   reflectors to that RHS, and ends as the final solution to
            ///   the upper triangular system of linear equations.  Pointer
            ///   to the start of at least \p nrow data elements.
            LLSProblem(const int     nrow,
                       const int     ncol,
                       const Scalar* A,
                       const Scalar* beta,
                       Scalar*       house,
                       Scalar*       work)
                : nrow_ { nrow }
                , ncol_ { ncol }
                , A_    { A }
                , beta_ { beta }
                , house_{ house }
                , work_ { work }
            {}

            /// Number of LLS coefficient matrix rows.
            int nrow_{};

            /// Number of LLS coefficient matrix columns.
            int ncol_{};

            /// QR factorisation of LLS coefficient matrix in factored form.
            const Scalar* A_{};

            /// Scaling elements for Householder reflectors.
            const Scalar* beta_{};

            /// Work array for reconstituting individual full Householder
            /// reflectors.  Pointer to start of at least nrow_ elements.
            Scalar* house_{};

            /// Scratch array for intermediate results.  Pointer to the
            /// start of at least nrow_ elements.
            Scalar* work_{};
        };

        /// Insert new coefficient matrix at the end of the current
        /// sequence.
        ///
        /// \param[in] A LLS coefficient matrix.  Matrix will be copied into
        ///   internal storage, translated to column major order, and put
        ///   into a factored form QR factorisation.
        void push_back(const CoefficientMatrix& A);

        /// Retieve solver for individual LLS problem.
        ///
        /// \param[in] i Sequence index for individual LLS problem.
        ///
        /// \return Solver for individual LLS problem, i.e., with a
        ///   predfined coefficient matrix entered in a prior call to
        ///   push_back().  Will solve one or more specific LLS problems by
        ///   varying the system RHS.  Client code is expected to capture
        ///   solution values in own, separate storage if those values will
        ///   be needed for additional processing.
        LLSProblem operator[](const std::size_t i) const;

    private:
        /// Type alias for representing offset into value arrays.
        using Offset = typename std::vector<Scalar>::size_type;

        /// CSR-like start pointers for matrix data.
        std::vector<Offset> startMatrix_{0};

        /// CSR-like start pointers for Householder reflector scaling data
        /// ("beta").
        std::vector<Offset> startHHScale_{0};

        /// Factored form QR factorisations of sequence coefficient matrices.
        std::vector<Scalar> matrixData_{};

        /// Scaling data for Householder reflectors for every problem in the
        /// sequence.
        std::vector<Scalar> householderScale_{};

        /// Work array for calculating individual Householder reflectors.
        /// Also passed as the 'house' constructor argument to the \c
        /// LLSProblem formed in \code operator[]() \endcode.
        mutable std::vector<Scalar> house_{};

        /// Work array for applying a Householder reflector to a vector.
        /// Also passed as the 'work' constructor argument to the \c
        /// LLSProblem formed in \code operator[]() \endcode.
        mutable std::vector<Scalar> work_{};

        /// Copy LLS coefficient matrix data into internal storage
        ///
        /// Translates from row major format to column major format.
        ///
        /// \param[in] A LLS coefficient matrix.
        void internaliseCoeffMatrixData(const CoefficientMatrix& A);

        /// Factorise LLS coefficient matrix into QR form.
        ///
        /// Store 'Q' matrix result in factored form as sequence of
        /// essential part of Housholder reflectors below the matrix
        /// diagonal, the result upper triangular matrix 'R' on and above
        /// the matrix diagonal, and the Housholder scaling vectors in
        /// separate array.
        ///
        /// Mutates matrixData_ and householderScale_.  Uses house_ and
        /// work_ for intermediate results.
        ///
        /// \param[in] i LLS sequence index.
        void factorQR(const std::size_t i);

        /// Compute Householder reflector for a single column of a single
        /// LLS coefficient matrix.
        ///
        /// Reflector value available in house_ upon completion.
        ///
        /// \param[in] A Coefficient matrix.  Pointer to the start of
        ///   unfactored matrix data.
        ///
        /// \param[in] nrow Number of rows in LLS coefficient matrix.
        ///
        /// \param[in] col Zero-based column index specifying the column of
        ///   \p A for which to compute the Householder reflector.
        ///
        /// \return Scaling ("beta") of Householder reflector.
        Scalar computeHHReflector(const Scalar*     A,
                                  const std::size_t nrow,
                                  const std::size_t col);

        /// Apply Householder reflector to remainder of coefficient matrix
        ///
        /// Performs a rank-1 (outer product) update of columns to the right
        /// of the column for which we're currently computing the
        /// factorisation.  When computing the Householder reflection of
        /// column 'j', this function updates the values in A(j:end,j:end).
        /// Following this update, elements A(0:j,j) hold column 'j' of the
        /// upper triangular matrix 'R' in the QR factorisation.
        ///
        /// Mutates work_ to hold the intermediate result
        ///
        ///     w = beta * v' * A(j:end, j:end)
        ///
        /// \param[in,out] A In-progress factorised LLS coefficient matrix.
        ///   Pointer to the start of the matrix data.
        ///
        /// \param[in] beta Scaling for single Householder reflector.
        ///
        /// \param[in] nrow Number of coefficient matrix rows.
        ///
        /// \param[in] ncol Number of coefficient matrix columns.
        ///
        /// \param[in] col Column for which we're currently computing the QR
        ///   factorisation.
        void applyHHReflector(Scalar*           A,
                              const Scalar      beta,
                              const std::size_t nrow,
                              const std::size_t ncol,
                              const std::size_t col);

        /// Store essential part of single Householder reflector back into
        /// the matrix data, below the diagonal.
        ///
        /// \param[out] A In-progress factorised LLS coefficient matrix.
        ///   Pointer to the start of the matrix data.
        ///
        /// \param[in] nrow Number of coefficient matrix rows.
        ///
        /// \param[in] col Zero-based column index specifying the column of
        ///   \p A for which we're currently calculating the QR factorisation.
        void preserveEssentialHHReflector(Scalar*           A,
                                          const std::size_t nrow,
                                          const std::size_t col);
    };

} // namespace Opm

#endif // OPM_LINALG_LEAST_SQUARES_PROBLEM_SEQUENCE_HPP

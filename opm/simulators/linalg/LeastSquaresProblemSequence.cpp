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

#include <config.h>

#include <opm/simulators/linalg/LeastSquaresProblemSequence.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numeric>
#include <vector>

template <typename Scalar>
const Scalar*
Opm::LeastSquaresProblemSequence<Scalar>::LLSProblem::solve(const Scalar* b) const
{
    // w <- b
    std::copy_n(b, this->nrow_, this->work_);

    // w <- Q' * w (== Q' * b == Q \ b)
    for (auto col = 0*this->ncol_; col < this->ncol_; ++col) {
        std::fill_n(this->house_, col, 0.0);
        this->house_[col] = 1.0;
        std::copy(this->A_ + col*this->nrow_ + col + 1,
                  this->A_ + (col + 1)*this->nrow_,
                  this->house_ + col + 1);

        const auto scale = this->beta_[col] *
            std::inner_product(this->work_, this->work_ + this->nrow_,
                               this->house_, 0.0);

        for (auto row = 0*this->nrow_; row < this->nrow_; ++row) {
            this->work_[row] -= scale * this->house_[row];
        }
    }

    // w <- R \ w (== R \ (Q \ b) == least-squares approximation)
    for (auto col = this->ncol_; col > 0; --col) {
        this->work_[col - 1] /= this->A_[(col - 1) * (this->nrow_ + 1)];

        const auto* Acol = &this->A_[0 + (col - 1)*this->nrow_];
        const auto  xc   = this->work_[col - 1];

        for (auto row = col - 1; row > 0; --row) {
            this->work_[row - 1] -= xc * Acol[row - 1];
        }
    }

    return this->work_;
}

// ---------------------------------------------------------------------------

template <typename Scalar>
void Opm::LeastSquaresProblemSequence<Scalar>::push_back(const CoefficientMatrix& A)
{
    this->matrixData_.resize(this->matrixData_.size() + A.a_.size());
    this->householderScale_.resize(this->householderScale_.size() + A.ncol_);

    this->internaliseCoeffMatrixData(A);

    if (const auto nrow = A.a_.size() / A.ncol_; nrow > this->work_.size()) {
        this->work_.resize(nrow);
    }

    this->startMatrix_.push_back(this->matrixData_.size());
    this->startHHScale_.push_back(this->householderScale_.size());

    this->factorQR(this->startMatrix_.size() - 1 - 1);
}

template <typename Scalar>
typename Opm::LeastSquaresProblemSequence<Scalar>::LLSProblem
Opm::LeastSquaresProblemSequence<Scalar>::operator[](const std::size_t i) const
{
    const auto ncol = this->startHHScale_[i + 1] - this->startHHScale_[i];
    const auto nrow = (this->startMatrix_[i + 1] - this->startMatrix_[i]) / ncol;

    return {
        static_cast<int>(nrow), static_cast<int>(ncol),
        this->matrixData_.data() + this->startMatrix_[i],
        this->householderScale_.data() + this->startHHScale_[i],
        const_cast<std::vector<Scalar>&>(this->house_).data(),
        const_cast<std::vector<Scalar>&>(this->work_).data()
    };
}

template <typename Scalar>
void Opm::LeastSquaresProblemSequence<Scalar>::
internaliseCoeffMatrixData(const CoefficientMatrix& A)
{
    auto* columnMajorA = this->matrixData_.data() + this->startMatrix_.back();

    const auto nrow = A.a_.size() / A.ncol_;
    for (auto col = 0*A.ncol_; col < A.ncol_; ++col) {
        for (auto row = 0*nrow; row < nrow; ++row) {
            *columnMajorA++ = A.a_[row*A.ncol_ + col];
        }
    }
}

template <typename Scalar>
void Opm::LeastSquaresProblemSequence<Scalar>::factorQR(const std::size_t i)
{
    const auto ncol = this->startHHScale_[i + 1] - this->startHHScale_[i];
    const auto nrow = (this->startMatrix_[i + 1] - this->startMatrix_[i]) / ncol;

    auto* A    = this->matrixData_.data() + this->startMatrix_[i];
    auto* beta = this->householderScale_.data() + this->startHHScale_[i];

    for (auto col = 0*ncol; col < ncol; ++col) {
        beta[col] = this->computeHHReflector(A, nrow, col);
        this->applyHHReflector(A, beta[col], nrow, ncol, col);
        this->preserveEssentialHHReflector(A, nrow, col);
    }
}

template <typename Scalar>
Scalar Opm::LeastSquaresProblemSequence<Scalar>::
computeHHReflector(const Scalar*     A,
                   const std::size_t nrow,
                   const std::size_t col)
{
    const auto nelem = nrow - col;

    if (nelem > this->house_.size()) {
        // When col == 0 and nrow exceeds the high watermark seen thus far
        // for number of coefficient matrix rows.
        this->house_.resize(nelem);
    }

    std::copy_n(A + col*(nrow + 1), nelem, this->house_.begin());
    const auto nrm2 = std::inner_product(this->house_.begin(),
                                         this->house_.begin() + nelem,
                                         this->house_.begin(), 0.0);

    const auto nrm2_rest = nrm2 - this->house_[0]*this->house_[0];

    if (! (nrm2_rest > 0.0)) {
        std::fill(this->house_.begin(), this->house_.end(), 0.0);
        this->house_[0] = 1.0;
        return 0.0;
    }

    // Normalise input vector => "mu" = 1.
    std::transform(this->house_.begin(), this->house_.begin() + nelem,
                   this->house_.begin(), [nrm = std::sqrt(nrm2)](const auto x)
                   { return x / nrm; });

    const auto sigma = nrm2_rest / nrm2;
    const auto mu = 1.0;

    auto& v1 = this->house_[0];

    if (v1 > 0.0) {
        v1 = - sigma / (v1 + mu);
    }
    else {
        v1 -= mu;
    }

    const auto beta = 2 * v1*v1 / (sigma + v1*v1);

    std::transform(this->house_.begin(), this->house_.begin() + nelem,
                   this->house_.begin(), [scale = v1](const auto x)
                   {
                       return x / scale;
                   });

    return beta;
}

template <typename Scalar>
void Opm::LeastSquaresProblemSequence<Scalar>::
applyHHReflector(Scalar*           A,
                 const Scalar      beta,
                 const std::size_t nrow,
                 const std::size_t ncol,
                 const std::size_t col)
{
    // work <- beta * v' * A(col:end, col:end) ("w")
    for (auto c = 0*ncol; c < ncol - col; ++c) {
        this->work_[c] = beta *
            std::inner_product(this->house_.begin(),
                               this->house_.begin() + nrow - col,
                               A + nrow*(col + c) + col, 0.0);
    }

    // A(col:end, col:end) -= v * w'
    for (auto c = col; c < ncol; ++c) {
        const auto w = this->work_[c - col];
        for (auto r = col; r < nrow; ++r) {
            A[r + c*nrow] -= this->house_[r - col] * w;
        }
    }
}

template <typename Scalar>
void Opm::LeastSquaresProblemSequence<Scalar>::
preserveEssentialHHReflector(Scalar*           A,
                             const std::size_t nrow,
                             const std::size_t col)
{
    const auto nelem = nrow - (col + 1);

    std::copy_n(this->house_.begin() + 1, nelem, A + col*(nrow + 1) + 1);
}

// ===========================================================================
// Explicit Specialisations of LeastSquaresProblemSequence Template
//
// No other code below this separator
// ===========================================================================

template class Opm::LeastSquaresProblemSequence<float>;
template class Opm::LeastSquaresProblemSequence<double>;

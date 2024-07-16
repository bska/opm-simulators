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

#include <dune/common/fvector.hh>

#include <dune/grid/common/partitionset.hh>
#include <dune/grid/common/rangegenerators.hh>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <numeric>

template <typename Scalar>
template <typename GridView, typename ElementMapper, typename TransContainer>
void Opm::NormalisedTotalConcVariation<Scalar>::
initialise(const GridView&       gv,
           const ElementMapper&  emap,
           const TransContainer& transContainer)
{
    using CoeffMat = typename LeastSquaresProblemSequence<Scalar>::CoefficientMatrix;

    // Interior, boundary, and ghost inclusive.
    this->massFrac_.resize(gv.size(0));

    auto addTrans = [this](const Scalar t1, const Scalar t2) {
        this->trans_.push_back(Scalar{1} / (Scalar{1}/t1 + Scalar{1}/t2));
    };

    constexpr auto dimWorld = GridView::dimensionworld;

    auto c = Dune::FieldVector<Scalar, dimWorld>{};
    auto numInternal = typename std::vector<Scalar>::size_type{0};

    for (const auto& cell : elements(gv, Dune::Partitions::interior)) {
        ++numInternal;
        const auto i1 = emap.index(cell);

        auto CC = CoeffMat { dimWorld };
        for (const auto& is : intersections(gv, cell)) {
            c = is.geometry().center();
            transContainer.distanceVector(i1, c);

            CC.addRow(c.data());

            if (! is.neighbor()) {
                continue;
            }

            const auto i2 = emap.index(is.outside());

            addTrans(transContainer.thermalHalfTrans(i1, i2),
                     transContainer.thermalHalfTrans(i2, i1));
        }

        this->llsSeq_.push_back(CC);
    }

    this->concVariation_.resize(numInternal);
}

template <typename Scalar>
void Opm::NormalisedTotalConcVariation<Scalar>::prepareAccumulation()
{
    this->maxMassFrac_ = std::numeric_limits<Scalar>::lowest();
    std::fill(this->massFrac_.begin(), this->massFrac_.end(), Scalar{0});
}

template <typename Scalar>
void Opm::NormalisedTotalConcVariation<Scalar>::
accumulate(const typename std::vector<Scalar>::size_type i,
           const Scalar                                  massFrac)
{
    this->massFrac_[i] = massFrac;
    this->maxMassFrac_ = std::max(this->maxMassFrac_, massFrac);
}

template <typename Scalar>
template <typename GridView, typename ElementMapper>
void Opm::NormalisedTotalConcVariation<Scalar>::
endAccumulation(const GridView& gv, const ElementMapper& emap)
{
    this->maxMassFrac_ = gv.comm().max(this->maxMassFrac_);

    if (! (this->maxMassFrac_ > Scalar{0})) {
        // Nothing to do.
        return;
    }

    auto norm2 = [](const Scalar* x)
    {
        return std::sqrt(std::inner_product(x, x + GridView::dimensionworld, x, Scalar{0}));
    };

    auto i = typename std::vector<Scalar>::size_type{0};
    auto t = this->trans_.begin();
    auto f = std::vector<Scalar>{};

    for (const auto& cell : elements(gv, Dune::Partitions::interior)) {
        f.clear();

        {
            const auto f1 = this->massFrac_[emap.index(cell)];

            for (const auto& is : intersections(gv, cell)) {
                if (! is.neighbor()) {
                    f.push_back(Scalar{0});
                    continue;
                }

                const auto f2 = this->massFrac_[emap.index(is.outside())];

                f.push_back(*t++ * (f1 - f2) / this->maxMassFrac_);
            }
        }

        const auto lls = this->llsSeq_[i];

        assert (lls.nrow() == static_cast<int>(f.size()));
        assert (lls.ncol() == GridView::dimensionworld);

        this->concVariation_[i++] = norm2(lls.solve(f.data()));
    }
}

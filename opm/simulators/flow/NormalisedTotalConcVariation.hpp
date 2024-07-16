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

#ifndef OPM_TOTAL_CONCENTRATION_VARIATION_INCLUDED_HPP
#define OPM_TOTAL_CONCENTRATION_VARIATION_INCLUDED_HPP

#include <opm/simulators/linalg/LeastSquaresProblemSequence.hpp>

#include <vector>

namespace Opm {

    template <typename Scalar>
    class NormalisedTotalConcVariation
    {
    public:
        template <typename GridView, typename ElementMapper, typename TransContainer>
        void initialise(const GridView&       gv,
                        const ElementMapper&  emap,
                        const TransContainer& transContainer);

        void prepareAccumulation();
        void accumulate(const typename std::vector<Scalar>::size_type i,
                        const Scalar                                  massFrac);

        template <typename GridView, typename ElementMapper>
        void endAccumulation(const GridView&      gv,
                             const ElementMapper& emap);

        const std::vector<Scalar>& concVariation() const
        {
            return this->concVariation_;
        }

    private:
        LeastSquaresProblemSequence<Scalar> llsSeq_{};
        std::vector<Scalar> trans_{};
        Scalar maxMassFrac_{};
        std::vector<Scalar> massFrac_{};
        std::vector<Scalar> concVariation_;
    };

} // namespace Opm

#include "NormalisedTotalConcVariation_impl.hpp"

#endif // OPM_TOTAL_CONCENTRATION_VARIATION_INCLUDED_HPP

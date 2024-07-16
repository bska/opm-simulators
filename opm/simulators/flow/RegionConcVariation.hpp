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

#ifndef OPM_REGION_TOTAL_CONCENTRATION_VARIATION_INCLUDED_HPP
#define OPM_REGION_TOTAL_CONCENTRATION_VARIATION_INCLUDED_HPP

#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <cstddef>
#include <functional>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace Opm {

    template <typename Scalar>
    class RegionNormalisedTotalConcVariation
    {
    public:
        using RegDef = std::function<const std::vector<int>&(const std::string&)>;
        using ConcVarIter = typename std::vector<Scalar>::const_iterator;
        using ConcVarSpan = std::pair<ConcVarIter, ConcVarIter>;

        RegionNormalisedTotalConcVariation(const std::vector<std::string>& regions,
                                           const RegDef&                   regDef,
                                           const std::size_t               declaredMaxRegID,
                                           const Parallel::Communication&  comm);

        void prepareAccumulation();

        void addConcVariation(const std::size_t cell,
                              const Scalar      concVar);

        void accumulateParallel(const Parallel::Communication& comm);

        const std::vector<std::string>& regions() const
        {
            return this->regions_;
        }

        ConcVarSpan getConcVariation(const std::size_t region) const;

    private:
        std::vector<std::string> regions_{};

        std::size_t numCells_{};
        std::size_t maxRegID_{};

        std::vector<int> regDef_{};
        std::vector<Scalar> concVar_{};
    };

} // namespace Opm

#endif // OPM_REGION_TOTAL_CONCENTRATION_VARIATION_INCLUDED_HPP

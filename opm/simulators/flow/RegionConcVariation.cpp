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

#include <opm/simulators/flow/RegionConcVariation.hpp>

#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <algorithm>
#include <cstddef>
#include <functional>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

    std::string normalisedRegsetName(const std::string& regSet)
    {
        const auto maxChar = std::string::size_type{3};
        const auto prefix = std::string_view { "FIP" };
        const auto start = prefix.size() * (regSet.find(prefix) == 0);

        return regSet.substr(start, maxChar);
    }

    std::vector<std::string> normalisedRegsetNames(std::vector<std::string> regSets)
    {
        if (std::find(regSets.begin(), regSets.end(), "FIPNUM") == regSets.end()) {
            // Standard region set FIPNUM is always present, even if not
            // explicitly mentioned in the input.
            regSets.push_back("FIPNUM");
        }

        std::transform(regSets.begin(), regSets.end(), regSets.begin(),
                        &normalisedRegsetName);

        return regSets;
    }

    std::vector<std::string> sorted(std::vector<std::string> strings)
    {
        std::sort(strings.begin(), strings.end());
        return strings;
    }

} // Anonymous namespace

// ===========================================================================

template <typename Scalar>
Opm::RegionNormalisedTotalConcVariation<Scalar>::
RegionNormalisedTotalConcVariation(const std::vector<std::string>& regions,
                                   const RegDef&                   regDef,
                                   const std::size_t               declaredMaxRegID,
                                   const Parallel::Communication&  comm)
    : regions_ { sorted(normalisedRegsetNames(regions)) }
{
    for (const auto& region : this->regions_) {
        const auto& reg = regDef("FIP" + region);

        if (this->numCells_ == 0) {
            this->numCells_ = reg.size();
        }
        else if (this->numCells_ != reg.size()) {
            throw std::invalid_argument {
                "Size mismatch in region " + region
            };
        }

        this->regDef_.insert(this->regDef_.end(), reg.begin(), reg.end());

        const auto maxRegID =
            std::accumulate(reg.begin(), reg.end(), declaredMaxRegID,
                            [](const std::size_t i1, const std::size_t i2)
                            { return std::max(i1, i2); });

        this->maxRegID_ = std::max(this->maxRegID_, maxRegID);
    }

    this->maxRegID_ = comm.max(this->maxRegID_);
    this->concVar_.resize(this->regions_.size() * (this->maxRegID_ + 1));
}

template <typename Scalar>
void Opm::RegionNormalisedTotalConcVariation<Scalar>::prepareAccumulation()
{
    std::fill(this->concVar_.begin(), this->concVar_.end(), Scalar{0});
}

template <typename Scalar>
void Opm::RegionNormalisedTotalConcVariation<Scalar>::
addConcVariation(const std::size_t cell, const Scalar concVar)
{
    const auto nReg = this->regions_.size();

    for (auto reg = 0*nReg; reg < nReg; ++reg) {
        const auto regID = this->regDef_[reg*this->numCells_ + cell];
        this->concVar_[reg*(this->maxRegID_ + 1) + regID] += concVar;
    }
}

template <typename Scalar>
void Opm::RegionNormalisedTotalConcVariation<Scalar>::
accumulateParallel(const Parallel::Communication& comm)
{
    comm.sum(this->concVar_.data(), this->concVar_.size());
}

template <typename Scalar>
typename Opm::RegionNormalisedTotalConcVariation<Scalar>::ConcVarSpan
Opm::RegionNormalisedTotalConcVariation<Scalar>::
getConcVariation(const std::size_t region) const
{
    return {
        this->concVar_.begin() + (region + 0)*(this->maxRegID_ + 1),
        this->concVar_.begin() + (region + 1)*(this->maxRegID_ + 1)
    };
}

// ===========================================================================

template class Opm::RegionNormalisedTotalConcVariation<float>;
template class Opm::RegionNormalisedTotalConcVariation<double>;

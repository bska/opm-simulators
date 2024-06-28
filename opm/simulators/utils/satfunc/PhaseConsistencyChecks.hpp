/*
  Copyright 2024 Equinor AS

  This file is part of the Open Porous Media project (OPM).

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

#ifndef OIL_PHASE_CONSISTENCY_CHECKS_HPP_INCLUDED
#define OIL_PHASE_CONSISTENCY_CHECKS_HPP_INCLUDED

#include <opm/simulators/utils/satfunc/SatfuncConsistencyChecks.hpp>

#include <memory>
#include <vector>

namespace Opm::Satfunc::PhaseChecks {

    template <typename Scalar>
    using CheckPtr = std::unique_ptr<typename SatfuncConsistencyChecks<Scalar>::Check>;

    template <typename Scalar>
    std::vector<CheckPtr<Scalar>> createGasChecks();

    template <typename Scalar>
    std::vector<CheckPtr<Scalar>> createOilChecks();

    template <typename Scalar>
    std::vector<CheckPtr<Scalar>> createWaterChecks();

} // namespace Opm::Satfunc::PhaseChecks

#endif // OIL_PHASE_CONSISTENCY_CHECKS_HPP_INCLUDED

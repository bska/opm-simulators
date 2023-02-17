/*
  Copyright 2023 SINTEF Digital, Mathematics and Cybernetics.

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

#ifndef OPM_NONLINEARCONVERGENCERATE_HEADER_INCLUDED
#define OPM_NONLINEARCONVERGENCERATE_HEADER_INCLUDED

#include <array>
#include <cstddef>
#include <optional>
#include <vector>

namespace Opm {
    class ConvergenceReport;
} // namespace Opm

namespace Opm {

    class NonlinearConvergenceRate
    {
    public:
        class Target
        {
        public:
            enum class ValueCategory { Raw, Log10, };
            struct Point { double mb{}; double cnv{}; };

            double distance(const Point&        point,
                            const ValueCategory cat = ValueCategory::Raw) const;

            Target& mb(const double        value,
                       const ValueCategory cat = ValueCategory::Raw)
            {
                this->conv_.mb = Target::metric(value, cat);
                return *this;
            }

            Target& cnv(const double        value,
                        const ValueCategory cat = ValueCategory::Raw)
            {
                this->conv_.cnv = Target::metric(value, cat);
                return *this;
            }

        private:
            Point conv_{ -7.0, -2.0 };

            static double metric(const double x, const ValueCategory cat);
        };

        NonlinearConvergenceRate() = default;
        explicit NonlinearConvergenceRate(const Target& target)
            : target_ { target }
        {}

        Target& defineTarget()
        {
            return this->target_;
        }

        void measureRate(const ConvergenceReport& convRep);
        void clear();

        std::size_t numConvergenceViolations() const
        {
            return this->numViolations_;
        }

        double currentDistance() const
        {
            return this->distance_[this->current_];
        }

        const std::vector<double>& improvements() const
        {
            return this->improvements_;
        }

    private:
        Target target_{};
        std::array<double, 2> distance_{};
        std::size_t current_{0};
        std::size_t numViolations_{0};
        std::optional<std::size_t> optimalNumIter_{};
        std::vector<double> improvements_{};

        std::size_t nextIndex() const;

        void commitDistanceMetric(const double d);
        void estimateOptimalNumberOfIterations(const double d);
        void inferConvergenceBehaviour(const double d);
    };

} // namespace Opm

#endif // OPM_NONLINEARCONVERGENCERATE_HEADER_INCLUDED

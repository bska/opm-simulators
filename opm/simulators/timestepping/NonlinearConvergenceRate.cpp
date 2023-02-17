#include <config.h>

#include <opm/simulators/timestepping/NonlinearConvergenceRate.hpp>
#include <opm/simulators/timestepping/ConvergenceReport.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace {
    double lineDistance(const Opm::NonlinearConvergenceRate::Target::Point&          pt,
                        const Opm::NonlinearConvergenceRate::Target::Point&          conv,
                        const double Opm::NonlinearConvergenceRate::Target::Point::* comp,
                        const bool                                                   outside)
    {
        return outside
            ? std::hypot(pt.mb - conv.mb, pt.cnv - conv.cnv) // L2
            : std::abs(pt.*comp - conv.*comp);               // L_\infty
    }

    Opm::NonlinearConvergenceRate::Target::Point
    getConvergencePoint(const Opm::ConvergenceReport& convRep)
    {
        using Metric = Opm::ConvergenceReport::ReservoirFailure::Type;

        auto pt = Opm::NonlinearConvergenceRate::Target::Point{};

        for (const auto& metric : convRep.reservoirConvergence()) {
            switch (metric.type()) {
            case Metric::MassBalance:
                pt.mb += metric.value();
                break;

            case Metric::Cnv:
                pt.cnv += metric.value();
                break;

            default:
                break;
            }
        }

        return pt;
    }

    std::size_t optimalNumberOfIterations(const double distance,
                                          const double offset,
                                          const double reduction)
    {
        const auto numIt = std::ceil((distance - offset) / reduction);

        return std::max(static_cast<std::size_t>(numIt), std::size_t{1});
    }
} // Anonymous namespace

double
Opm::NonlinearConvergenceRate::Target::
distance(const Point&        pt0,
         const ValueCategory cat) const
{
    const auto pt = Point {
        Target::metric(pt0.mb , cat),
        Target::metric(pt0.cnv, cat)
    };

    const auto outside_1 = pt.mb  > this->conv_.mb;
    const auto outside_2 = pt.cnv > this->conv_.cnv;

    if (! (outside_1 || outside_2)) {
        // We're inside convergence target region.
        return 0.0;
    }

    const auto d1 = lineDistance(pt, this->conv_, &Point::mb , outside_2);
    const auto d2 = lineDistance(pt, this->conv_, &Point::cnv, outside_1);

    return std::min(d1, d2);
}

double
Opm::NonlinearConvergenceRate::Target::metric(const double        x,
                                              const ValueCategory cat)
{
    switch (cat) {
    case ValueCategory::Raw:   return std::log10(x);
    case ValueCategory::Log10: return x;
    }

    throw std::logic_error {
        "Unsupported value category '" +
        std::to_string(static_cast<int>(cat)) + '\''
    };
}

// ------------------------------------------------------------------------

void
Opm::NonlinearConvergenceRate::measureRate(const ConvergenceReport& convRep)
{
    const auto currDist = this->target_
        .distance(getConvergencePoint(convRep));

    if (! this->optimalNumIter_.has_value()) {
        this->estimateOptimalNumberOfIterations(currDist);
    }
    else if (*this->optimalNumIter_ > std::size_t{0}) {
        this->inferConvergenceBehaviour(currDist);
    }

    this->commitDistanceMetric(currDist);
}

std::size_t
Opm::NonlinearConvergenceRate::nextIndex() const
{
    return (this->current_ + 1) % this->distance_.size();
}

void
Opm::NonlinearConvergenceRate::commitDistanceMetric(const double currDist)
{
    this->current_ = this->nextIndex();
    this->distance_[this->current_] = currDist;
}

void
Opm::NonlinearConvergenceRate::estimateOptimalNumberOfIterations(const double currDist)
{
    this->optimalNumIter_ = (currDist > 0.0)
        ? optimalNumberOfIterations(currDist, 0.05, 0.75)
        : 0;
}

void
Opm::NonlinearConvergenceRate::inferConvergenceBehaviour(const double currDist)
{
    const auto improve = this->distance_[this->current_] - currDist;
    const auto estIter = std::ceil(currDist / improve);

    this->numViolations_ += (estIter < 0.0) ||
        (estIter > 4.0 * *this->optimalNumIter_);

    this->improvements_.push_back(improve);
}

void
Opm::NonlinearConvergenceRate::clear()
{
    this->distance_.fill(std::numeric_limits<double>::max());
    this->current_ = 0;
    this->numViolations_ = 0;
    this->optimalNumIter_.reset();
    this->improvements_.clear();
}

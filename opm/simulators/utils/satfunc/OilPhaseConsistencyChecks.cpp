#include <config.h>

#include <opm/simulators/utils/satfunc/OilPhaseConsistencyChecks.hpp>

#include <opm/material/fluidmatrixinteractions/EclEpsScalingPoints.hpp>

#include <cmath>

namespace {
    constexpr auto one_bit = static_cast<unsigned char>(1);
    constexpr auto failed_bit = one_bit << 0u;
    constexpr auto critical_bit = one_bit << 1u;
}

// ===========================================================================

template <typename Scalar>
void Opm::Satfunc::PhaseChecks::Oil::SOCheckBase<Scalar>::
test(const EclEpsScalingPoints<Scalar>& endPoints)
{
    this->flags_ = static_cast<unsigned char>(0);

    this->testImpl(endPoints);
}

template <typename Scalar>
bool Opm::Satfunc::PhaseChecks::Oil::SOCheckBase<Scalar>::isViolated() const
{
    return (this->flags_ & failed_bit) != 0;
}

template <typename Scalar>
bool Opm::Satfunc::PhaseChecks::Oil::SOCheckBase<Scalar>::isCritical() const
{
    return (this->flags_ & critical_bit) != 0;
}

template <typename Scalar>
void Opm::Satfunc::PhaseChecks::Oil::SOCheckBase<Scalar>::setViolated()
{
    this->flags_ |= failed_bit;
}

template <typename Scalar>
void Opm::Satfunc::PhaseChecks::Oil::SOCheckBase<Scalar>::setCritical()
{
    this->flags_ |= critical_bit;
}

// ---------------------------------------------------------------------------

template <typename Scalar>
void Opm::Satfunc::PhaseChecks::Oil::SOcr_GO<Scalar>::
testImpl(const EclEpsScalingPoints<Scalar>& endPoints)
{
    this->sogcr_ = endPoints.Sogcr;

    if (! std::isfinite(this->sogcr_)) {
        this->setViolated();
        this->setCritical();
    }

    const auto low = this->sogcr_ < Scalar{0};
    const auto high = ! (this->sogcr_ < Scalar{1});

    if (low || high) { this->setViolated(); }
    if (low)         { this->setCritical(); }
}

// ---------------------------------------------------------------------------

template <typename Scalar>
void Opm::Satfunc::PhaseChecks::Oil::SOmin_GO<Scalar>::
testImpl(const EclEpsScalingPoints<Scalar>& endPoints)
{
    this->swl_ = endPoints.Swl;
    this->sgu_ = endPoints.Sgu;

    if (! std::isfinite(this->swl_) ||
        ! std::isfinite(this->sgu_))
    {
        this->setViolated();
        this->setCritical();
    }

    if (this->swl_ + this->sgu_ > Scalar{1}) {
        this->setViolated();
    }
}

// ---------------------------------------------------------------------------

template <typename Scalar>
void Opm::Satfunc::PhaseChecks::Oil::SOcr_OW<Scalar>::
testImpl(const EclEpsScalingPoints<Scalar>& endPoints)
{
    this->sowcr_ = endPoints.Sowcr;

    if (! std::isfinite(this->sowcr_)) {
        this->setViolated();
        this->setCritical();
    }

    const auto low = this->sowcr_ < Scalar{0};
    const auto high = ! (this->sowcr_ < Scalar{1});

    if (low || high) { this->setViolated(); }
    if (low)         { this->setCritical(); }
}

// ---------------------------------------------------------------------------

template <typename Scalar>
void Opm::Satfunc::PhaseChecks::Oil::SOmin_OW<Scalar>::
testImpl(const EclEpsScalingPoints<Scalar>& endPoints)
{
    this->sgl_ = endPoints.Sgl;
    this->swu_ = endPoints.Swu;

    if (! std::isfinite(this->sgl_) ||
        ! std::isfinite(this->swu_))
    {
        this->setViolated();
        this->setCritical();
    }

    if (this->sgl_ + this->swu_ > Scalar{1}) {
        this->setViolated();
    }
}

// ---------------------------------------------------------------------------

template <typename Scalar>
void Opm::Satfunc::PhaseChecks::Oil::MobileOil_GO_SGmin<Scalar>::
testImpl(const EclEpsScalingPoints<Scalar>& endPoints)
{
    this->swl_   = endPoints.Swl;
    this->sgl_   = endPoints.Sgl;
    this->sogcr_ = endPoints.Sogcr;

    if (! std::isfinite(this->swl_) ||
        ! std::isfinite(this->sgl_) ||
        ! std::isfinite(this->sogcr_))
    {
        this->setViolated();
        this->setCritical();
    }

    if (! (this->sogcr_ < Scalar{1} - this->swl_ - this->sgl_)) {
        this->setViolated();
    }
}

// ---------------------------------------------------------------------------

template <typename Scalar>
void Opm::Satfunc::PhaseChecks::Oil::MobileOil_GO_SGcr<Scalar>::
testImpl(const EclEpsScalingPoints<Scalar>& endPoints)
{
    this->swl_   = endPoints.Swl;
    this->sgcr_  = endPoints.Sgcr;
    this->sogcr_ = endPoints.Sgcr;

    if (! std::isfinite(this->swl_) ||
        ! std::isfinite(this->sgcr_) ||
        ! std::isfinite(this->sogcr_))
    {
        this->setViolated();
        this->setCritical();
    }

    if (! (this->sogcr_ < Scalar{1} - this->swl_ - this->sgcr_)) {
        this->setViolated();
    }
}

// ---------------------------------------------------------------------------

template <typename Scalar>
void Opm::Satfunc::PhaseChecks::Oil::MobileOil_OW_SWmin<Scalar>::
testImpl(const EclEpsScalingPoints<Scalar>& endPoints)
{
    this->swl_   = endPoints.Swl;
    this->sgl_   = endPoints.Sgl;
    this->sowcr_ = endPoints.Sowcr;

    if (! std::isfinite(this->swl_) ||
        ! std::isfinite(this->sgl_) ||
        ! std::isfinite(this->sogcr_))
    {
        this->setViolated();
        this->setCritical();
    }

    if (! (this->sowcr_ < Scalar{1} - this->swl_ - this->sgl_)) {
        this->setViolated();
    }
}

// ---------------------------------------------------------------------------

template <typename Scalar>
void Opm::Satfunc::PhaseChecks::Oil::MobileOil_GO_SGcr<Scalar>::
testImpl(const EclEpsScalingPoints<Scalar>& endPoints)
{
    this->sgl_   = endPoints.Sgl;
    this->swcr_  = endPoints.Swcr;
    this->sowcr_ = endPoints.Sowcr;

    if (! std::isfinite(this->sgl_) ||
        ! std::isfinite(this->swcr_) ||
        ! std::isfinite(this->sowcr_))
    {
        this->setViolated();
        this->setCritical();
    }

    if (! (this->sowcr_ < Scalar{1} - this->swcr_ - this->sgl_)) {
        this->setViolated();
    }
}

// ===========================================================================
// Explicit Specialisations of Individual Check Templates
//
// No other code below this separator
// ===========================================================================

template class Opm::Satfunc::PhaseChecks::Oil::SOCheckBase<float>;
template class Opm::Satfunc::PhaseChecks::Oil::SOCheckBase<double>;

// ---------------------------------------------------------------------------

template class Opm::Satfunc::PhaseChecks::Oil::SOcr_GO<float>;
template class Opm::Satfunc::PhaseChecks::Oil::SOcr_GO<double>;

// ---------------------------------------------------------------------------

template class Opm::Satfunc::PhaseChecks::Oil::SOmin_GO<float>;
template class Opm::Satfunc::PhaseChecks::Oil::SOmin_GO<double>;

// ---------------------------------------------------------------------------

template class Opm::Satfunc::PhaseChecks::Oil::MobileOil_GO_SGmin<float>;
template class Opm::Satfunc::PhaseChecks::Oil::MobileOil_GO_SGmin<double>;

// ---------------------------------------------------------------------------

template class Opm::Satfunc::PhaseChecks::Oil::MobileOil_GO_SGcr<float>;
template class Opm::Satfunc::PhaseChecks::Oil::MobileOil_GO_SGcr<double>;

// ---------------------------------------------------------------------------

template class Opm::Satfunc::PhaseChecks::Oil::SOcr_OW<float>;
template class Opm::Satfunc::PhaseChecks::Oil::SOcr_OW<double>;

// ---------------------------------------------------------------------------

template class Opm::Satfunc::PhaseChecks::Oil::SOmin_OW<float>;
template class Opm::Satfunc::PhaseChecks::Oil::SOmin_OW<double>;

// ---------------------------------------------------------------------------

template class Opm::Satfunc::PhaseChecks::Oil::MobileOil_OW_SWmin<float>;
template class Opm::Satfunc::PhaseChecks::Oil::MobileOil_OW_SWmin<double>;

// ---------------------------------------------------------------------------

template class Opm::Satfunc::PhaseChecks::Oil::MobileOil_OW_SWcr<float>;
template class Opm::Satfunc::PhaseChecks::Oil::MobileOil_OW_SWcr<double>;

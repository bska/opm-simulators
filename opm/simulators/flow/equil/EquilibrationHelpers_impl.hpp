// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
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

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/

#include <opm/common/TimingMacros.hpp>

#include <opm/common/utility/numeric/RootFinders.hpp>

#include <opm/material/fluidmatrixinteractions/EclMaterialLawManager.hpp>
#include <opm/material/fluidstates/SimpleModularFluidState.hpp>
#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>

#include <opm/simulators/flow/equil/EquilibrationHelpers.hpp>

#include <fmt/format.h>

namespace Opm { namespace EQUIL {

using FluidSystemSimple = BlackOilFluidSystem<double>;

// Adjust oil pressure according to gas saturation and cap pressure
using SatOnlyFluidState = SimpleModularFluidState<double,
                                                  /*numPhases=*/3,
                                                  /*numComponents=*/3,
                                                  FluidSystemSimple,
                                                  /*storePressure=*/false,
                                                  /*storeTemperature=*/false,
                                                  /*storeComposition=*/false,
                                                  /*storeFugacity=*/false,
                                                  /*storeSaturation=*/true,
                                                  /*storeDensity=*/false,
                                                  /*storeViscosity=*/false,
                                                  /*storeEnthalpy=*/false>;

namespace Miscibility {

template<class FluidSystem>
RsVD<FluidSystem>::RsVD(const int pvtRegionIdx,
                        const std::vector<Scalar>& depth,
                        const std::vector<Scalar>& rs)
    : pvtRegionIdx_(pvtRegionIdx)
    , rsVsDepth_(depth, rs)
{}

template<class FluidSystem>
typename RsVD<FluidSystem>::Scalar
RsVD<FluidSystem>::
operator()(const Scalar depth,
           const Scalar press,
           const Scalar temp,
           const Scalar satGas) const
{
    const auto sat_rs = satRs(press, temp);

    if (satGas > std::sqrt(std::numeric_limits<Scalar>::epsilon())) {
        return sat_rs;
    }

    if (this->rsVsDepth_.xMin() > depth) {
        return std::min(sat_rs, this->rsVsDepth_.valueAt(0));
    }

    if (this->rsVsDepth_.xMax() < depth) {
        return std::min(sat_rs, this->rsVsDepth_.valueAt(rsVsDepth_.numSamples() - 1));
    }

    return std::min(sat_rs, this->rsVsDepth_.eval(depth, /*extrapolate=*/false));
}

template<class FluidSystem>
typename RsVD<FluidSystem>::Scalar
RsVD<FluidSystem>::
satRs(const Scalar press, const Scalar temp) const
{
    return FluidSystem::oilPvt().saturatedGasDissolutionFactor(pvtRegionIdx_, temp, press);
}

template<class FluidSystem>
PBVD<FluidSystem>::PBVD(const int pvtRegionIdx,
                        const std::vector<Scalar>& depth,
                        const std::vector<Scalar>& pbub)
    : pvtRegionIdx_(pvtRegionIdx)
    , pbubVsDepth_ (depth, pbub)
{}

template<class FluidSystem>
typename PBVD<FluidSystem>::Scalar
PBVD<FluidSystem>::
operator()(const Scalar depth,
           const Scalar cellPress,
           const Scalar temp,
           const Scalar satGas) const
{
    auto press = cellPress;

    if (satGas <= 0.0) {
        if (pbubVsDepth_.xMin() > depth) {
            press = pbubVsDepth_.valueAt(0);
        }
        else if (pbubVsDepth_.xMax() < depth) {
            press = pbubVsDepth_.valueAt(pbubVsDepth_.numSamples() - 1);
        }
        else {
            press = pbubVsDepth_.eval(depth, /*extrapolate=*/false);
        }
    }

    return satRs(std::min(press, cellPress), temp);
}

template<class FluidSystem>
typename PBVD<FluidSystem>::Scalar
PBVD<FluidSystem>::
satRs(const Scalar press, const Scalar temp) const
{
    return FluidSystem::oilPvt().saturatedGasDissolutionFactor(pvtRegionIdx_, temp, press);
}

template<class FluidSystem>
PDVD<FluidSystem>::PDVD(const int pvtRegionIdx,
                        const std::vector<Scalar>& depth,
                        const std::vector<Scalar>& pdew)
    : pvtRegionIdx_(pvtRegionIdx)
    , pdewVsDepth_(depth, pdew)
{}

template<class FluidSystem>
typename PDVD<FluidSystem>::Scalar
PDVD<FluidSystem>::
operator()(const Scalar depth,
           const Scalar cellPress,
           const Scalar temp,
           const Scalar satOil) const
{
    auto press = cellPress;

    if (satOil <= 0.0) {
        if (pdewVsDepth_.xMin() > depth) {
            press = pdewVsDepth_.valueAt(0);
        }
        else if (pdewVsDepth_.xMax() < depth) {
            press = pdewVsDepth_.valueAt(pdewVsDepth_.numSamples() - 1);
        }
        else {
            press = pdewVsDepth_.eval(depth, /*extrapolate=*/false);
        }
    }

    return satRv(std::min(press, cellPress), temp);
}

template<class FluidSystem>
typename PDVD<FluidSystem>::Scalar
PDVD<FluidSystem>::
satRv(const Scalar press, const Scalar temp) const
{
    return FluidSystem::gasPvt().saturatedOilVaporizationFactor(pvtRegionIdx_, temp, press);
}

template<class FluidSystem>
RvVD<FluidSystem>::RvVD(const int pvtRegionIdx,
                        const std::vector<Scalar>& depth,
                        const std::vector<Scalar>& rv)
    : pvtRegionIdx_(pvtRegionIdx)
    , rvVsDepth_   (depth, rv)
{}

template<class FluidSystem>
typename RvVD<FluidSystem>::Scalar
RvVD<FluidSystem>::
operator()(const Scalar depth,
           const Scalar press,
           const Scalar temp,
           const Scalar satOil) const
{
    if (satOil < - std::sqrt(std::numeric_limits<Scalar>::epsilon())) {
        throw std::logic_error {
            "Must not pass negative oil saturation"
        };
    }

    const auto sat_rv = satRv(press, temp);

    if (satOil > std::sqrt(std::numeric_limits<Scalar>::epsilon())) {
        return sat_rv;
    }

    if (rvVsDepth_.xMin() > depth) {
        return std::min(sat_rv, rvVsDepth_.valueAt(0));
    }

    if (rvVsDepth_.xMax() < depth) {
        return std::min(sat_rv, rvVsDepth_.valueAt(rvVsDepth_.numSamples() - 1));
    }

    return std::min(sat_rv, rvVsDepth_.eval(depth, /*extrapolate=*/false));
}

template<class FluidSystem>
typename RvVD<FluidSystem>::Scalar
RvVD<FluidSystem>::
satRv(const Scalar press, const Scalar temp) const
{
    return FluidSystem::gasPvt().saturatedOilVaporizationFactor(pvtRegionIdx_, temp, press);
}

template<class FluidSystem>
RvwVD<FluidSystem>::RvwVD(const int pvtRegionIdx,
                          const std::vector<Scalar>& depth,
                          const std::vector<Scalar>& rvw)
    : pvtRegionIdx_(pvtRegionIdx)
    , rvwVsDepth_(depth, rvw)
{}

template<class FluidSystem>
typename RvwVD<FluidSystem>::Scalar
RvwVD<FluidSystem>::
operator()(const Scalar depth,
           const Scalar press,
           const Scalar temp,
           const Scalar satWat) const
{
    if (satWat < - std::sqrt(std::numeric_limits<double>::epsilon())) {
        throw std::logic_error {
            "Must not pass negative water saturation"
        };
    }

    const auto sat_rvw = satRvw(press, temp);

    if (satWat > std::sqrt(std::numeric_limits<double>::epsilon())) {
        return sat_rvw; // saturated Rvw
    }

    if (rvwVsDepth_.xMin() > depth) {
        return std::min(sat_rvw, this->rvwVsDepth_.valueAt(0));
    }

    if (rvwVsDepth_.xMax() < depth) {
        return std::min(sat_rvw, this->rvwVsDepth_.valueAt(rvwVsDepth_.numSamples() - 1));
    }

    return std::min(sat_rvw, this->rvwVsDepth_.eval(depth, /*extrapolate=*/false));
}

template<class FluidSystem>
typename RvwVD<FluidSystem>::Scalar
RvwVD<FluidSystem>::
satRvw(const Scalar press, const Scalar temp) const
{
    return FluidSystem::gasPvt().saturatedWaterVaporizationFactor(pvtRegionIdx_, temp, press);
}

template<class FluidSystem>
RsSatAtContact<FluidSystem>::
RsSatAtContact(const int pvtRegionIdx, const Scalar pContact,  const Scalar T_contact)
    : pvtRegionIdx_(pvtRegionIdx)
{
    rsSatContact_ = satRs(pContact, T_contact);
}

template<class FluidSystem>
typename RsSatAtContact<FluidSystem>::Scalar
RsSatAtContact<FluidSystem>::
operator()(const Scalar /* depth */,
           const Scalar press,
           const Scalar temp,
           const Scalar satGas) const
{
    if (satGas > std::sqrt(std::numeric_limits<Scalar>::epsilon())) {
        return satRs(press, temp);
    }

    return std::min(satRs(press, temp), this->rsSatContact_);
}

template<class FluidSystem>
typename RsSatAtContact<FluidSystem>::Scalar
RsSatAtContact<FluidSystem>::
satRs(const Scalar press, const Scalar temp) const
{
    return FluidSystem::oilPvt()
        .saturatedGasDissolutionFactor(pvtRegionIdx_, temp, press);
}

template<class FluidSystem>
RvSatAtContact<FluidSystem>::
RvSatAtContact(const int pvtRegionIdx, const Scalar pContact, const Scalar T_contact)
    : pvtRegionIdx_(pvtRegionIdx)
{
    rvSatContact_ = satRv(pContact, T_contact);
}

template<class FluidSystem>
typename RvSatAtContact<FluidSystem>::Scalar
RvSatAtContact<FluidSystem>::
operator()(const Scalar /*depth*/,
           const Scalar press,
           const Scalar temp,
           const Scalar satOil) const
{
    if (satOil > std::sqrt(std::numeric_limits<Scalar>::epsilon())) {
        return satRv(press, temp);
    }

    return std::min(satRv(press, temp), rvSatContact_);
}

template<class FluidSystem>
typename RvSatAtContact<FluidSystem>::Scalar
RvSatAtContact<FluidSystem>::
satRv(const Scalar press, const Scalar temp) const
{
    return FluidSystem::gasPvt()
        .saturatedOilVaporizationFactor(this->pvtRegionIdx_, temp, press);;
}

template<class FluidSystem>
RvwSatAtContact<FluidSystem>::
RvwSatAtContact(const int pvtRegionIdx, const Scalar pContact, const Scalar T_contact)
    : pvtRegionIdx_(pvtRegionIdx)
{
    rvwSatContact_ = satRvw(pContact, T_contact);
}

template<class FluidSystem>
typename RvwSatAtContact<FluidSystem>::Scalar
RvwSatAtContact<FluidSystem>::
operator()(const Scalar /*depth*/,
           const Scalar press,
           const Scalar temp,
           const Scalar satWat) const
{
    if (satWat > std::sqrt(std::numeric_limits<Scalar>::epsilon())) {
        return satRvw(press, temp);
    }

    return std::min(satRvw(press, temp), rvwSatContact_);
}

template<class FluidSystem>
typename RvwSatAtContact<FluidSystem>::Scalar
RvwSatAtContact<FluidSystem>::
satRvw(const Scalar press, const Scalar temp) const
{
    return FluidSystem::gasPvt()
        .saturatedWaterVaporizationFactor(pvtRegionIdx_, temp, press);;
}

} // namespace Miscibility

template<class Scalar>
EquilReg<Scalar>::EquilReg(const EquilRecord& rec,
                           std::shared_ptr<Miscibility::RsFunction<Scalar>> rs,
                           std::shared_ptr<Miscibility::RsFunction<Scalar>> rv,
                           std::shared_ptr<Miscibility::RsFunction<Scalar>> rvw,
                           const TabulatedFunction& tempVdTable,
                           const TabulatedFunction& saltVdTable,
                           const int pvtIdx)
    : rec_    (rec)
    , rs_     (rs)
    , rv_     (rv)
    , rvw_     (rvw)
    , tempVdTable_ (tempVdTable)
    , saltVdTable_ (saltVdTable)
    , pvtIdx_ (pvtIdx)
{}

template<class Scalar>
Scalar EquilReg<Scalar>::datum() const
{
    return this->rec_.datumDepth();
}

template<class Scalar>
Scalar EquilReg<Scalar>::pressure() const
{
    return this->rec_.datumDepthPressure();
}

template<class Scalar>
Scalar EquilReg<Scalar>::zwoc() const
{
    return this->rec_.waterOilContactDepth();
}

template<class Scalar>
Scalar EquilReg<Scalar>::pcowWoc() const
{
    return this->rec_.waterOilContactCapillaryPressure();
}

template<class Scalar>
Scalar EquilReg<Scalar>::zgoc() const
{
    return this->rec_.gasOilContactDepth();
}

template<class Scalar>
Scalar EquilReg<Scalar>::pcgoGoc() const
{
    return this->rec_.gasOilContactCapillaryPressure();
}

template<class Scalar>
int EquilReg<Scalar>::equilibrationAccuracy() const
{
    return this->rec_.initializationTargetAccuracy();
}

template<class Scalar>
const typename EquilReg<Scalar>::CalcDissolution&
EquilReg<Scalar>::dissolutionCalculator() const
{
    return *this->rs_;
}

template<class Scalar>
const typename EquilReg<Scalar>::CalcEvaporation&
EquilReg<Scalar>::evaporationCalculator() const
{
    return *this->rv_;
}

template<class Scalar>
const typename EquilReg<Scalar>::CalcWaterEvaporation&
EquilReg<Scalar>::waterEvaporationCalculator() const
{
    return *this->rvw_;
}

template<class Scalar>
const typename EquilReg<Scalar>::TabulatedFunction&
EquilReg<Scalar>::saltVdTable() const
{
    return saltVdTable_;
}

template<class Scalar>
const typename EquilReg<Scalar>::TabulatedFunction&
EquilReg<Scalar>::tempVdTable() const
{
    return tempVdTable_;
}

template<class Scalar>
int EquilReg<Scalar>::pvtIdx() const
{
    return this->pvtIdx_;
}

template<class FluidSystem, class MaterialLawManager>
PcEq<FluidSystem,MaterialLawManager>::
PcEq(const MaterialLawManager& materialLawManager,
     const int phase,
     const int cell,
     const Scalar targetPc)
    : materialLawManager_(materialLawManager)
    , phase_(phase)
    , cell_(cell)
    , targetPc_(targetPc)
{
    // This code should only be called for water or gas phase
    assert(phase_ != FluidSystem::oilPhaseIdx);
}

template<class FluidSystem, class MaterialLawManager>
typename PcEq<FluidSystem,MaterialLawManager>::Scalar
PcEq<FluidSystem,MaterialLawManager>::operator()(const Scalar s) const
{
    SatOnlyFluidState fluidState;
    fluidState.setSaturation(FluidSystem::waterPhaseIdx, 0.0);
    fluidState.setSaturation(FluidSystem::oilPhaseIdx, 0.0);
    fluidState.setSaturation(FluidSystem::gasPhaseIdx, 0.0);
    fluidState.setSaturation(phase_, s);

    // The capillaryPressure code only uses the wetting phase saturation
    if (phase_ == FluidSystem::gasPhaseIdx) {
        fluidState.setSaturation(FluidSystem::waterPhaseIdx, 1.0 - s);
        fluidState.setSaturation(FluidSystem::oilPhaseIdx, 1.0 - s);
    }

    std::array<Scalar, FluidSystem::numPhases> pc{0.0};
    {
        const auto& matParams = materialLawManager_.materialLawParams(cell_);

        using MaterialLaw = typename MaterialLawManager::MaterialLaw;
        MaterialLaw::capillaryPressures(pc, matParams, fluidState);
    }

    const Scalar sign = (phase_ == FluidSystem::waterPhaseIdx)? -1.0 : 1.0;
    const Scalar pcPhase = pc[FluidSystem::oilPhaseIdx] + sign*pc[phase_];

    return pcPhase - targetPc_;
}

template<class FluidSystem, class MaterialLawManager>
PcEqSum<FluidSystem,MaterialLawManager>::
PcEqSum(const MaterialLawManager& materialLawManager,
        const int phase1,
        const int phase2,
        const int cell,
        const Scalar targetPc)
    : materialLawManager_(materialLawManager)
    , phase1_(phase1)
    , phase2_(phase2)
    , cell_(cell)
    , targetPc_(targetPc)
{}

template<class FluidSystem, class MaterialLawManager>
typename PcEqSum<FluidSystem,MaterialLawManager>::Scalar
PcEqSum<FluidSystem,MaterialLawManager>::operator()(const Scalar s) const
{
    SatOnlyFluidState fluidState{};
    fluidState.setSaturation(FluidSystem::waterPhaseIdx, 0.0);
    fluidState.setSaturation(FluidSystem::oilPhaseIdx, 0.0);
    fluidState.setSaturation(FluidSystem::gasPhaseIdx, 0.0);
    fluidState.setSaturation(phase1_, s);
    fluidState.setSaturation(phase2_, 1.0 - s);

    std::array<Scalar, FluidSystem::numPhases> pc {0.0};
    {
        const auto& matParams = materialLawManager_.materialLawParams(cell_);

        using MaterialLaw = typename MaterialLawManager::MaterialLaw;
        MaterialLaw::capillaryPressures(pc, matParams, fluidState);
    }

    const Scalar sign1 = (phase1_ == FluidSystem::waterPhaseIdx)? -1.0 : 1.0;
    const Scalar pc1 = pc[FluidSystem::oilPhaseIdx] + sign1*pc[phase1_];

    const Scalar sign2 = (phase2_ == FluidSystem::waterPhaseIdx)? -1.0 : 1.0;
    const Scalar pc2 = pc[FluidSystem::oilPhaseIdx] + sign2*pc[phase2_];

    return pc1 + pc2 - targetPc_;
}

template <class FluidSystem, class MaterialLawManager>
typename FluidSystem::Scalar
minSaturations(const MaterialLawManager& materialLawManager,
               const int phase, const int cell)
{
    if (phase == FluidSystem::oilPhaseIdx) {
        throw std::runtime_error {
            "Minimum saturation not implemented for oil phase."
        };
    }

    if ((phase != FluidSystem::gasPhaseIdx) &&
        (phase != FluidSystem::waterPhaseIdx))
    {
        throw std::runtime_error("Unknown phase");
    }

    const auto& scaledDrainageInfo =
        materialLawManager.oilWaterScaledEpsInfoDrainage(cell);

    if (phase == FluidSystem::gasPhaseIdx) {
        return scaledDrainageInfo.Sgl;
    }

    return scaledDrainageInfo.Swl;
}

template <class FluidSystem, class MaterialLawManager>
typename FluidSystem::Scalar
maxSaturations(const MaterialLawManager& materialLawManager,
               const int phase, const int cell)
{
    if (phase == FluidSystem::oilPhaseIdx) {
        throw std::runtime_error {
            "Maximum saturation not implemented for oil phase."
        };
    }

    if ((phase != FluidSystem::gasPhaseIdx) &&
        (phase != FluidSystem::waterPhaseIdx))
    {
        throw std::runtime_error("Unknown phase");
    }

    const auto& scaledDrainageInfo =
        materialLawManager.oilWaterScaledEpsInfoDrainage(cell);

    if (phase == FluidSystem::gasPhaseIdx) {
        return scaledDrainageInfo.Sgu;
    }

    return scaledDrainageInfo.Swu;
}

template <class FluidSystem, class MaterialLawManager>
typename FluidSystem::Scalar
satFromPc(const MaterialLawManager& materialLawManager,
          const int phase,
          const int cell,
          const typename FluidSystem::Scalar targetPc,
          const bool increasing)
{
    using Scalar = typename FluidSystem::Scalar;

    // Find minimum and maximum saturations.
    const Scalar s0 = increasing
        ? maxSaturations<FluidSystem>(materialLawManager, phase, cell)
        : minSaturations<FluidSystem>(materialLawManager, phase, cell);

    const Scalar s1 = increasing
        ? minSaturations<FluidSystem>(materialLawManager, phase, cell)
        : maxSaturations<FluidSystem>(materialLawManager, phase, cell);

    // Create the equation f(s) = pc(s) - targetPc
    const auto f = PcEq<FluidSystem, MaterialLawManager> {
        materialLawManager, phase, cell, targetPc
    };

    const auto f0 = f(s0);
    const auto f1 = f(s1);
    if (! std::isfinite(f0 + f1)) {
        throw std::logic_error {
            fmt::format("The capillary pressure values "
                        "{} and {} are not both finite", f0, f1)
        };
    }

    if (f0 <= Scalar{0}) { return s0; }
    if (f1 >= Scalar{0}) { return s1; }

    const auto tol = static_cast<Scalar>(1e-10);

    // Should at least converge in 2 times bisection but some safety here:
    const int maxIter = -2*static_cast<int>(std::log2(tol)) + 10;
    int usedIterations = -1;

    return RegulaFalsiBisection<ThrowOnError>::
        solve(f, s0, s1, maxIter, tol, usedIterations);
}

template<class FluidSystem, class MaterialLawManager>
typename FluidSystem::Scalar
satFromSumOfPcs(const MaterialLawManager& materialLawManager,
                const int phase1,
                const int phase2,
                const int cell,
                const typename FluidSystem::Scalar targetPc)
{
    using Scalar = typename FluidSystem::Scalar;

    // Create the equation f(s) = pc1(s) + pc2(1-s) - targetPc
    const auto f = PcEqSum<FluidSystem, MaterialLawManager> {
        materialLawManager, phase1, phase2, cell, targetPc
    };

    const auto s0 = minSaturations<FluidSystem>(materialLawManager, phase1, cell);
    const auto f0 = f(s0);
    if (f0 <= Scalar{0}) { return s0; }

    const auto s1 = maxSaturations<FluidSystem>(materialLawManager, phase1, cell);
    const auto f1 = f(s1);
    if (f1 >= Scalar{0}) { return s1; }

    assert ((f0 > Scalar{0}) && (f1 < Scalar{0}));
    const auto tol = static_cast<Scalar>(1e-10);

    // should at least converge in 2 times bisection but some safety here:
    const int maxIter = -2*static_cast<int>(std::log2(tol)) + 10;

    int usedIterations = -1;

    return RegulaFalsiBisection<ThrowOnError>::
        solve(f, s0, s1, maxIter, tol, usedIterations);
}

template<class FluidSystem, class MaterialLawManager>
typename FluidSystem::Scalar
satFromDepth(const MaterialLawManager& materialLawManager,
             const typename FluidSystem::Scalar cellDepth,
             const typename FluidSystem::Scalar contactDepth,
             const int phase,
             const int cell,
             const bool increasing)
{
    if (cellDepth <= contactDepth) {
        return increasing
            ? maxSaturations<FluidSystem>(materialLawManager, phase, cell)
            : minSaturations<FluidSystem>(materialLawManager, phase, cell);
    }

    return increasing
        ? minSaturations<FluidSystem>(materialLawManager, phase, cell)
        : maxSaturations<FluidSystem>(materialLawManager, phase, cell);
}

template<class FluidSystem, class MaterialLawManager>
bool isConstPc(const MaterialLawManager& materialLawManager,
               const int phase,
               const int cell)
{
    using Scalar = typename FluidSystem::Scalar;

    // Create the equation f(s) = pc(s);
    const auto f = PcEq<FluidSystem, MaterialLawManager> {
        materialLawManager, phase, cell, 0
    };

    const auto f0 = f(minSaturations<FluidSystem>(materialLawManager, phase, cell));
    const auto f1 = f(maxSaturations<FluidSystem>(materialLawManager, phase, cell));

    return std::abs(f0 - f1) < std::numeric_limits<Scalar>::epsilon();
}

}} // namespace Opm::EQUIL

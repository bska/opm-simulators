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

#include <cstddef>
#include <string>

namespace Opm::Satfunc::PhaseChecks::Oil {

    template <typename Scalar>
    class SOCheckBase : public typename SatfuncConsistencyChecks<Scalar>::Check
    {
    public:
        void test(const EclEpsScalingPointsInfo<Scalar>&) override;
        bool isViolated() const override;
        bool isCritical() const override;

    protected:
        void setViolated();
        void setCritical();

    private:
        unsigned char flags_{0};

        virtual void testImpl(const EclEpsScalingPointsInfo<Scalar>&) = 0;
    };

    // -----------------------------------------------------------------------

    template <typename Scalar>
    class SOcr_GO : public SOCheckBase<Scalar>
    {
    public:
        std::size_t numExportedCheckValues() const override { return 1 };

        void exportCheckValues(Scalar* exportedCheckValues) const override
        {
            exportedCheckValues[0] = this->sogcr_;
        }

        std::string description() const override
        {
            return { "Non-negative critical oil saturation in G/O system" };
        }

        std::string condition() const override
        {
            return { "0 <= SOGCR < 1" };
        }

        void columnNames(std::string* headers) const override
        {
            headers[0] = "SOGCR";
        }

    private:
        Scalar sogcr_;

        void testImpl(const EclEpsScalingPointsInfo<Scalar>&) override;
    };

    template <typename Scalar>
    class SOmin_GO : public SOCheckBase<Scalar>
    {
    public:
        std::size_t numExportedCheckValues() const override { return 3 };

        void exportCheckValues(double* exportedCheckValues) const override
        {
            exportedCheckValues[0] = this->swl_;
            exportedCheckValues[1] = this->sgu_;
            exportedCheckValues[2] = this->swl_ + this->sgu_;
        }

        std::string description() const override
        {
            return { "Non-negative minimum oil saturation in G/O system" };
        }

        std::string condition() const override
        {
            return { "SWL + SGU <= 1" };
        }

        void columnNames(std::string* headers) const override
        {
            headers[0] = "SWL";
            headers[1] = "SGU";
            headers[2] = "SWL + SGU";
        }

    private:
        Scalar swl_;
        Scalar sgu_;

        void testImpl(const EclEpsScalingPointsInfo<Scalar>&) override;
    };

    template <typename Scalar>
    class MobileOil_GO_SGmin : public SOCheckBase<Scalar>
    {
    public:
        std::size_t numExportedCheckValues() const override { return 4 };

        void exportCheckValues(double* exportedCheckValues) const override
        {
            exportedCheckValues[0] = this->swl_;
            exportedCheckValues[1] = this->sgl_;
            exportedCheckValues[2] = this->sogcr_;
            exportedCheckValues[3] = Scalar{1} - (this->swl_ + this->sgl_);
        }

        std::string description() const override
        {
            return { "Mobile oil saturation in G/O system at minimum gas saturation" };
        }

        std::string condition() const override
        {
            return { "SOGCR < 1 - SWL - SGL" };
        }

        void columnNames(std::string* headers) const override
        {
            headers[0] = "SWL";
            headers[1] = "SGL";
            headers[2] = "SOGCR";
            headers[3] = "1 - SWL - SGL";
        }

    private:
        Scalar swl_;
        Scalar sgl_;
        Scalar sogcr_;

        void testImpl(const EclEpsScalingPointsInfo<Scalar>&) override;
    };

    template <typename Scalar>
    class MobileOil_GO_SGcr : public SOCheckBase<Scalar>
    {
    public:
        std::size_t numExportedCheckValues() const override { return 4 };

        void exportCheckValues(double* exportedCheckValues) const override
        {
            exportedCheckValues[0] = this->swl_;
            exportedCheckValues[1] = this->sgcr_;
            exportedCheckValues[2] = this->sogcr_;
            exportedCheckValues[3] = Scalar{1} - (this->swl_ + this->sgcr_);
        }

        std::string description() const override
        {
            return { "Mobile oil saturation in G/O system at critical gas saturation" };
        }

        std::string condition() const override
        {
            return { "SOGCR < 1 - SWL - SGCR" };
        }

        void columnNames(std::string* headers) const override
        {
            headers[0] = "SWL";
            headers[1] = "SGCR";
            headers[2] = "SOGCR";
            headers[3] = "1 - SWL - SGCR";
        }

    private:
        Scalar swl_;
        Scalar sgcr_;
        Scalar sogcr_;

        void testImpl(const EclEpsScalingPointsInfo<Scalar>&) override;
    };

    // -----------------------------------------------------------------------

    template <typename Scalar>
    class SOcr_OW : public SOCheckBase<Scalar>
    {
    public:
        std::size_t numExportedCheckValues() const override { return 1 };

        void exportCheckValues(Scalar* exportedCheckValues) const override
        {
            exportedCheckValues[0] = this->sowcr_;
        }

        std::string description() const override
        {
            return { "Non-negative critical oil saturation in O/W system" };
        }

        std::string condition() const override
        {
            return { "0 <= SOWCR < 1" };
        }

        void columnNames(std::string* headers) const override
        {
            headers[0] = "SOWCR";
        }

    private:
        Scalar sowcr_;

        void testImpl(const EclEpsScalingPointsInfo<Scalar>&) override;
    };

    template <typename Scalar>
    class SOmin_OW : public SOCheckBase<Scalar>
    {
    public:
        std::size_t numExportedCheckValues() const override { return 3 };

        void exportCheckValues(double* exportedCheckValues) const override
        {
            exportedCheckValues[0] = this->sgl_;
            exportedCheckValues[1] = this->swu_;
            exportedCheckValues[2] = this->sgl_ + this->swu_;
        }

        std::string description() const override
        {
            return { "Non-negative minimum oil saturation in G/O system" };
        }

        std::string condition() const override
        {
            return { "SGL + SWU <= 1" };
        }

        void columnNames(std::string* headers) const override
        {
            headers[0] = "SGL";
            headers[1] = "SWU";
            headers[2] = "SGL + SUU";
        }

    private:
        Scalar sgl_;
        Scalar swu_;

        void testImpl(const EclEpsScalingPointsInfo<Scalar>&) override;
    };

    template <typename Scalar>
    class MobileOil_OW_SWmin : public SOCheckBase<Scalar>
    {
    public:
        std::size_t numExportedCheckValues() const override { return 4 };

        void exportCheckValues(double* exportedCheckValues) const override
        {
            exportedCheckValues[0] = this->swl_;
            exportedCheckValues[1] = this->sgl_;
            exportedCheckValues[2] = this->sowcr_;
            exportedCheckValues[3] = Scalar{1} - (this->swl_ + this->sgl_);
        }

        std::string description() const override
        {
            return { "Mobile oil saturation in O/W system at minimum water saturation" };
        }

        std::string condition() const override
        {
            return { "SOWCR < 1 - SWL - SGL" };
        }

        void columnNames(std::string* headers) const override
        {
            headers[0] = "SWL";
            headers[1] = "SGL";
            headers[2] = "SOWCR";
            headers[3] = "1 - SWL - SGL";
        }

    private:
        Scalar swl_;
        Scalar sgl_;
        Scalar sowcr_;

        void testImpl(const EclEpsScalingPointsInfo<Scalar>&) override;
    };

    template <typename Scalar>
    class MobileOil_OW_SWcr : public SOCheckBase<Scalar>
    {
    public:
        std::size_t numExportedCheckValues() const override { return 4 };

        void exportCheckValues(double* exportedCheckValues) const override
        {
            exportedCheckValues[0] = this->sgl_;
            exportedCheckValues[1] = this->swcr_;
            exportedCheckValues[2] = this->sogcr_;
            exportedCheckValues[3] = Scalar{1} - this->swcr_ - this->sgl_;
        }

        std::string description() const override
        {
            return { "Mobile oil saturation in O/W system at critical water saturation" };
        }

        std::string condition() const override
        {
            return { "SOWCR < 1 - SWCR - SGL" };
        }

        void columnNames(std::string* headers) const override
        {
            headers[0] = "SGL";
            headers[1] = "SWCR";
            headers[2] = "SOWCR";
            headers[3] = "1 - SWCR - SGL";
        }

    private:
        Scalar sgl_;
        Scalar swcr_;
        Scalar sowcr_;

        void testImpl(const EclEpsScalingPointsInfo<Scalar>&) override;
    };

} // namespace Opm::Satfunc::PhaseChecks::Oil

#endif // OIL_PHASE_CONSISTENCY_CHECKS_HPP_INCLUDED

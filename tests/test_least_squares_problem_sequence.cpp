/*
  Copyright 2024 Equinor.

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

#include <config.h>

#define BOOST_TEST_MODULE Linear_Least_Squares_Problem_Sequence

#include <boost/test/unit_test.hpp>

#include <opm/simulators/linalg/LeastSquaresProblemSequence.hpp>

#include <array>
#include <cstddef>

using LLSSeq = Opm::LeastSquaresProblemSequence<double>;

BOOST_AUTO_TEST_SUITE(Cartesian_Geometry)

BOOST_AUTO_TEST_SUITE(Two_Space_Dimensions)

namespace {
    constexpr int ncol()
    {
        return 2;               // Number of space dimensions
    }

    template <typename CentroidVectors>
    LLSSeq makeLLSequence(const int n, CentroidVectors&& c)
    {
        auto llsseq = LLSSeq{};

        for (auto i = 0*n; i < n; ++i) {
            auto a = LLSSeq::CoefficientMatrix { ncol() };

            for (const auto& row : c) {
                a.addRow(row.data());
            }

            llsseq.push_back(a);
        }

        return llsseq;
    }

    LLSSeq unitSquare(const int n = 1)
    {
        return makeLLSequence(n, std::array {
            std::array { -0.5,  0.0 },
            std::array {  0.5,  0.0 },
            std::array {  0.0, -0.5 },
            std::array {  0.0,  0.5 },
            });
    }

#if 0
    // Isoceles triangle with base 1 and height 4.  Centroid in (0, 4/3),
    // face centroids in (0,0), (-1/2, 2), (1/2, 2).
    LLSSeq isoscelesTriangle(const int n = 1)
    {
        const auto fourThirds = 4.0 / 3.0;
        const auto quarter    = 1.0 / 4.0;

        return makeLLSequence(n, std::array {
            std::array {      0.0, 0.0 - fourThirds },
            std::array {  quarter, 2.0 - fourThirds },
            std::array { -quarter, 2.0 - fourThirds },
            });
    }
#endif
}

BOOST_AUTO_TEST_CASE(Diagonal_Gradient)
{
    const auto llsProblems = unitSquare();

    const auto llsProblem = llsProblems[0];

    BOOST_CHECK_EQUAL(llsProblem.nrow(), 4);
    BOOST_CHECK_EQUAL(llsProblem.ncol(), 2);

    const auto flux = std::array { -0.5, 0.5, -0.5, 0.5 };
    const auto* const gradient = llsProblem.solve(flux.data());

    BOOST_CHECK_CLOSE(gradient[0], 1.0, 1.0e-8);
    BOOST_CHECK_CLOSE(gradient[1], 1.0, 1.0e-8);
}

BOOST_AUTO_TEST_CASE(Multiple_Diagonal_Gradients)
{
    const auto n = 5;

    const auto llsProblems = unitSquare(n);
    const auto flux = std::array { -0.5, 0.5, -0.5, 0.5 };

    for (auto p = 0*n; p < n; ++p) {
        const auto llsProblem = llsProblems[p];
        const auto* const gradient = llsProblem.solve(flux.data());

        BOOST_TEST_MESSAGE("LLS " << (p + 1) << '/' << n);

        BOOST_CHECK_CLOSE(gradient[0], 1.0, 1.0e-8);
        BOOST_CHECK_CLOSE(gradient[1], 1.0, 1.0e-8);
    }
}

BOOST_AUTO_TEST_SUITE_END()     // Two_Space_Dimensions

BOOST_AUTO_TEST_SUITE_END()     // Cartesian_Geometry

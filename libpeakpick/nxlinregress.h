/*
 * <n times linear regression. - Use as many linear functions as possible to fit a data set.>
 * Copyright (C) 2018  Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#pragma once

#include <map>
#include <vector>

#include "mathhelper.h"
#include "peakpick.h"

namespace PeakPick {

struct MultiRegression {
    std::vector<LinearRegression> regressions;
    double sum_err = 0;
    std::vector<int> start;
};

class JacobVector {

public:
    inline JacobVector(std::vector<int> start, std::vector<int> end)
        : m_value(start)
        , m_end(end)
    {
    }

    inline bool Up()
    {
        for (unsigned int i = m_end.size() - 1; i != 0; --i) {
            if (m_value[i] < m_end[i]) {
                m_value[i]++;
                return true;
            } else if (m_value[i] == m_end[i] && i) {
                std::vector<int> initial = m_value;
                if (initial[i - 1] < m_end[i - 1]) {
                    initial[i - 1]++;
                    for (unsigned int j = i; j < m_end.size(); ++j) {
                        /* FIXME: this assigns [i] on every pass instead of [j], so for 4+ functions
                         * the entries above i keep their old (maxed) values and the vector goes
                         * non-monotonic, which is what makes segments come out empty. Left as-is
                         * deliberately: correcting it makes whole families of splits valid that are
                         * currently discarded, which changes results for 4+ functions (reachable
                         * from the regression analysis dialog). That needs to be a separate,
                         * measured change - it is not a no-op cleanup. */
                        initial[i] = initial[i - 1] + 2;
                    }
                    m_value = initial;
                    return true;
                }

            } else if (!i) {
                return false;
            }
        }
        return false;
    }

    inline std::vector<int> Value() const { return m_value; }

private:
    std::vector<int> m_end;
    std::vector<int> m_value;
};

inline std::map<double, MultiRegression> LeastSquares(const Vector& x, const Vector& y, unsigned int functions)
{
    std::map<double, MultiRegression> regressions;

    if (x.size() != y.size())
        return regressions;

    if (functions == 1) {
        LinearRegression regression = LeastSquares(x, y);
        MultiRegression reg;
        reg.regressions.push_back(regression);
        reg.sum_err = regression.sum_err;
        reg.start.push_back(0);
        reg.start.push_back(x.size() - 1);
        regressions[reg.sum_err] = reg;
    } else {
        /* The segment starts are staggered 0, 2, 4, ... 2*(functions-1) and the ends mirror them
         * from the far side. With fewer points than that span, the starts run past the data and the
         * ends go negative, so the segment loop below indexes x[j] out of bounds. There is no way to
         * place this many lines on this few points anyway - report no fit. */
        if (x.size() < 2 * static_cast<int>(functions - 1))
            return regressions;

        std::vector<int> starts, ends;
        for (unsigned int i = 0; i < functions; i++) {
            int start = 2 * i;
            starts.push_back(start);
            int end = x.size() - 2 * i;
            ends.insert(ends.begin(), end);
        }
        JacobVector vector(starts, ends);

        while (true) {
            if (vector.Value()[0] != 0)
                break;
            MultiRegression reg;
            double sum = 0;
            std::vector<int> work = vector.Value();
            work.push_back(x.size());
            bool valid = true;

            for (unsigned int i = 0; i < work.size() - 1; ++i) {

                std::vector<double> x_i, y_i;
                for (int j = work[i]; j < work[i + 1]; ++j) {

                    x_i.push_back(x[j]);
                    y_i.push_back(y[j]);
                }
                if (x_i.empty()) {
                    /* A split can be degenerate (work[i] >= work[i+1]): JacobVector walks the last
                     * index up to x.size(), which leaves the final segment with no points at all.
                     * Such a configuration is rejected below, so stop here instead of fitting an
                     * empty segment - that fit is meaningless, and reading &x_i[0] out of an empty
                     * vector is undefined behaviour. */
                    valid = false;
                    break;
                }
                LinearRegression regression = LeastSquares(x_i, y_i);
                reg.regressions.push_back(regression);

                reg.start.push_back(work[i]);
                reg.start.push_back(work[i + 1] - 1);

                sum += regression.sum_err;
            }
            reg.sum_err = sum;
            if (!std::isnan(sum) && valid)
                regressions[reg.sum_err] = reg;
            if (regressions.size() > 10)
                regressions.erase(--(regressions.end()));

            if (!vector.Up())
                break;
        }
    }
    return regressions;
}
}

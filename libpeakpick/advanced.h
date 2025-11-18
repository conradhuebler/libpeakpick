/*
 * <Advanced algorithms for libpeakpick>
 * Copyright (C) 2024  Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "analyse.h"
#include "mathhelper.h"
#include "spectrum.h"

typedef Eigen::VectorXd Vector;

namespace PeakPick {

/*!
 * \brief Calculate prominence of a peak
 * Prominence measures how much a peak stands out from its surroundings
 */
inline double CalculatePeakProminence(const spectrum* spec, unsigned int peak_idx)
{
    if (peak_idx >= spec->size())
        return 0.0;

    double peak_height = spec->Y(peak_idx);

    // Find lowest contour line on left side
    double left_min = peak_height;
    for (int i = peak_idx - 1; i >= 0; --i) {
        if (spec->Y(i) > peak_height)
            break; // Found higher peak
        left_min = std::min(left_min, spec->Y(i));
    }

    // Find lowest contour line on right side
    double right_min = peak_height;
    for (unsigned int i = peak_idx + 1; i < spec->size(); ++i) {
        if (spec->Y(i) > peak_height)
            break; // Found higher peak
        right_min = std::min(right_min, spec->Y(i));
    }

    // Prominence is height above the highest valley
    double base = std::max(left_min, right_min);
    return peak_height - base;
}

/*!
 * \brief Improved peak picking with prominence threshold
 * \param spec Input spectrum
 * \param threshold Minimum height threshold
 * \param min_prominence Minimum prominence (default 0.0 = disabled)
 * \param min_distance Minimum distance between peaks (default 1)
 * \return Vector of detected peaks
 */
inline std::vector<Peak> PickPeaksAdvanced(const spectrum* spec,
    double threshold,
    double min_prominence = 0.0,
    unsigned int min_distance = 1,
    unsigned int start = 0,
    unsigned int end = 0)
{
    std::vector<Peak> peaks;

    if (spec->size() < 3)
        return peaks;

    if (end == 0 || end > spec->size())
        end = spec->size();

    // Find local maxima
    std::vector<unsigned int> candidates;
    for (unsigned int i = start + 1; i < end - 1; ++i) {
        double y_curr = spec->Y(i);

        // Check if above threshold
        if (y_curr < threshold)
            continue;

        // Check if local maximum
        if (y_curr > spec->Y(i - 1) && y_curr > spec->Y(i + 1)) {
            candidates.push_back(i);
        }
    }

    // Filter by prominence
    std::vector<unsigned int> prominent_peaks;
    for (unsigned int idx : candidates) {
        double prominence = CalculatePeakProminence(spec, idx);
        if (prominence >= min_prominence) {
            prominent_peaks.push_back(idx);
        }
    }

    // Filter by minimum distance
    std::vector<bool> keep(prominent_peaks.size(), true);
    for (size_t i = 0; i < prominent_peaks.size(); ++i) {
        if (!keep[i])
            continue;

        for (size_t j = i + 1; j < prominent_peaks.size(); ++j) {
            if (!keep[j])
                continue;

            unsigned int dist = (prominent_peaks[j] > prominent_peaks[i]) ? (prominent_peaks[j] - prominent_peaks[i]) : (prominent_peaks[i] - prominent_peaks[j]);

            if (dist < min_distance) {
                // Keep the higher peak
                if (spec->Y(prominent_peaks[j]) > spec->Y(prominent_peaks[i])) {
                    keep[i] = false;
                    break;
                } else {
                    keep[j] = false;
                }
            }
        }
    }

    // Create Peak structures
    for (size_t i = 0; i < prominent_peaks.size(); ++i) {
        if (!keep[i])
            continue;

        Peak peak;
        peak.max = prominent_peaks[i];

        // Find peak start (left valley)
        peak.start = peak.max;
        for (int j = peak.max - 1; j >= static_cast<int>(start); --j) {
            if (spec->Y(j) < spec->Y(peak.start)) {
                peak.start = j;
            }
            // Stop at significant rise
            if (j > 0 && spec->Y(j) < spec->Y(j - 1))
                break;
        }

        // Find peak end (right valley)
        peak.end = peak.max;
        for (unsigned int j = peak.max + 1; j < end; ++j) {
            if (spec->Y(j) < spec->Y(peak.end)) {
                peak.end = j;
            }
            // Stop at significant rise
            if (j < spec->size() - 1 && spec->Y(j) < spec->Y(j + 1))
                break;
        }

        peak.int_start = peak.start;
        peak.int_end = peak.end;

        peaks.push_back(peak);
    }

    std::cout << "PickPeaksAdvanced: Found " << peaks.size() << " peaks (prominence >= " << min_prominence << ")" << std::endl;

    return peaks;
}

/*!
 * \brief Estimate noise level using Median Absolute Deviation (MAD)
 * \param spec Input spectrum
 * \return Estimated noise standard deviation
 */
inline double EstimateNoise(const spectrum& spec)
{
    if (spec.size() == 0)
        return 0.0;

    // Calculate first differences
    std::vector<double> differences;
    differences.reserve(spec.size() - 1);

    for (unsigned int i = 1; i < spec.size(); ++i) {
        differences.push_back(std::abs(spec.Y(i) - spec.Y(i - 1)));
    }

    // Calculate median
    std::sort(differences.begin(), differences.end());
    double median = differences[differences.size() / 2];

    // MAD-based noise estimate (assuming Gaussian noise)
    // Factor 1.4826 makes MAD consistent estimator for standard deviation
    return 1.4826 * median;
}

/*!
 * \brief Asymmetric Least Squares (AsLS) baseline correction
 * \param spec Input spectrum
 * \param lambda Smoothness parameter (typical: 1e2 to 1e9, larger = smoother)
 * \param p Asymmetry parameter (typical: 0.001 to 0.1, smaller = more asymmetric)
 * \param max_iter Maximum iterations
 * \return Baseline-corrected spectrum
 */
inline spectrum BaselineAsLS(const spectrum& spec,
    double lambda = 1e6,
    double p = 0.01,
    unsigned int max_iter = 10)
{
    unsigned int n = spec.size();
    if (n < 3) {
        std::cerr << "Warning: BaselineAsLS - spectrum too small" << std::endl;
        return spec;
    }

    std::cout << "BaselineAsLS: Computing baseline (lambda=" << lambda << ", p=" << p << ")" << std::endl;

    // Build difference matrix D (second derivative)
    typedef Eigen::Triplet<double> T;
    std::vector<T> triplets;
    triplets.reserve(3 * (n - 2));

    for (unsigned int i = 0; i < n - 2; ++i) {
        triplets.push_back(T(i, i, 1.0));
        triplets.push_back(T(i, i + 1, -2.0));
        triplets.push_back(T(i, i + 2, 1.0));
    }

    Eigen::SparseMatrix<double> D(n - 2, n);
    D.setFromTriplets(triplets.begin(), triplets.end());

    // D^T * D
    Eigen::SparseMatrix<double> DTD = D.transpose() * D;

    // Initial weights
    Vector w = Vector::Ones(n);
    Vector y = spec.y();
    Vector z = y;

    // Iterative reweighting
    for (unsigned int iter = 0; iter < max_iter; ++iter) {
        // Build weighted matrix W
        Eigen::SparseMatrix<double> W(n, n);
        W.reserve(Eigen::VectorXi::Constant(n, 1));
        for (unsigned int i = 0; i < n; ++i) {
            W.insert(i, i) = w(i);
        }

        // Solve: (W + lambda * D^T * D) * z = W * y
        Eigen::SparseMatrix<double> A = W + lambda * DTD;
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
        solver.compute(A);

        if (solver.info() != Eigen::Success) {
            std::cerr << "Warning: BaselineAsLS - solver failed at iteration " << iter << std::endl;
            break;
        }

        z = solver.solve(W * y);

        // Update weights: penalize points above baseline
        Vector d = y - z;
        double max_change = 0.0;
        for (unsigned int i = 0; i < n; ++i) {
            double w_old = w(i);
            w(i) = (d(i) < 0) ? p : (1.0 - p);
            max_change = std::max(max_change, std::abs(w(i) - w_old));
        }

        // Check convergence
        if (max_change < 1e-6) {
            std::cout << "BaselineAsLS: Converged after " << iter + 1 << " iterations" << std::endl;
            break;
        }
    }

    // Return baseline-corrected spectrum
    Vector corrected = y - z;
    return spectrum(spec.x(), corrected);
}

/*!
 * \brief Cubic spline interpolation coefficients
 */
struct CubicSpline {
    Vector x, y, a, b, c, d;
    int n;

    inline CubicSpline(const Vector& x_data, const Vector& y_data)
        : x(x_data)
        , y(y_data)
        , n(x_data.size())
    {
        if (n < 2) {
            std::cerr << "Warning: CubicSpline requires at least 2 points" << std::endl;
            return;
        }

        // Natural cubic spline (second derivatives = 0 at endpoints)
        a = y;
        b = Vector::Zero(n);
        c = Vector::Zero(n);
        d = Vector::Zero(n);

        if (n == 2) {
            // Linear interpolation
            b(0) = (y(1) - y(0)) / (x(1) - x(0));
            return;
        }

        // Build tridiagonal system
        Vector h(n - 1);
        for (int i = 0; i < n - 1; ++i) {
            h(i) = x(i + 1) - x(i);
        }

        Vector alpha(n - 1);
        for (int i = 1; i < n - 1; ++i) {
            alpha(i) = 3.0 / h(i) * (y(i + 1) - y(i)) - 3.0 / h(i - 1) * (y(i) - y(i - 1));
        }

        // Solve tridiagonal system
        Vector l(n), mu(n), z(n);
        l(0) = 1.0;
        mu(0) = 0.0;
        z(0) = 0.0;

        for (int i = 1; i < n - 1; ++i) {
            l(i) = 2.0 * (x(i + 1) - x(i - 1)) - h(i - 1) * mu(i - 1);
            mu(i) = h(i) / l(i);
            z(i) = (alpha(i) - h(i - 1) * z(i - 1)) / l(i);
        }

        l(n - 1) = 1.0;
        z(n - 1) = 0.0;
        c(n - 1) = 0.0;

        for (int j = n - 2; j >= 0; --j) {
            c(j) = z(j) - mu(j) * c(j + 1);
            b(j) = (y(j + 1) - y(j)) / h(j) - h(j) * (c(j + 1) + 2.0 * c(j)) / 3.0;
            d(j) = (c(j + 1) - c(j)) / (3.0 * h(j));
        }
    }

    inline double evaluate(double x_val) const
    {
        if (n < 2)
            return 0.0;

        // Find interval
        int i = 0;
        for (int j = 0; j < n - 1; ++j) {
            if (x_val >= x(j) && x_val <= x(j + 1)) {
                i = j;
                break;
            }
        }

        // Clamp to boundaries
        if (x_val < x(0))
            i = 0;
        if (x_val > x(n - 1))
            i = n - 2;

        // Evaluate cubic polynomial
        double dx = x_val - x(i);
        return a(i) + b(i) * dx + c(i) * dx * dx + d(i) * dx * dx * dx;
    }
};

/*!
 * \brief Resample spectrum using cubic spline interpolation
 * \param spec Input spectrum
 * \param new_x New X grid
 * \return Resampled spectrum
 */
inline spectrum resampleCubic(const spectrum& spec, const Vector& new_x)
{
    CubicSpline spline(spec.x(), spec.y());

    Vector new_y(new_x.size());
    for (int i = 0; i < new_x.size(); ++i) {
        new_y(i) = spline.evaluate(new_x(i));
    }

    return spectrum(new_x, new_y);
}

} // namespace PeakPick

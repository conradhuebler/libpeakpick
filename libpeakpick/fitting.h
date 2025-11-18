/*
 * <Multi-peak fitting algorithms>
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
#include <unsupported/Eigen/NonLinearOptimization>

#include <cmath>
#include <iostream>
#include <vector>

#include "analyse.h"
#include "mathhelper.h"
#include "spectrum.h"

typedef Eigen::VectorXd Vector;

namespace PeakPick {

enum class FitType {
    Gaussian,
    Lorentzian,
    PseudoVoigt
};

struct FitResult {
    Vector parameters;          // Fitted parameters
    Eigen::MatrixXd covariance; // Covariance matrix
    Vector residuals;           // Residuals
    double chi_squared;         // Chi-squared value
    double reduced_chi_squared; // Reduced chi-squared
    int iterations;             // Number of iterations
    bool converged;             // Convergence flag

    std::vector<Peak> fitted_peaks; // Peaks with fitted parameters
};

/*!
 * \brief Multi-peak fitting functor for Levenberg-Marquardt
 */
template <typename _Scalar, int NX = Eigen::Dynamic, int NY = Eigen::Dynamic>
struct MultiPeakFitFunctor {
    typedef _Scalar Scalar;
    enum {
        InputsAtCompileTime = NX,
        ValuesAtCompileTime = NY
    };
    typedef Eigen::Matrix<Scalar, InputsAtCompileTime, 1> InputType;
    typedef Eigen::Matrix<Scalar, ValuesAtCompileTime, 1> ValueType;
    typedef Eigen::Matrix<Scalar, ValuesAtCompileTime, InputsAtCompileTime> JacobianType;

    const Vector& m_x;
    const Vector& m_y;
    int m_num_peaks;
    FitType m_fit_type;
    int m_params_per_peak;

    MultiPeakFitFunctor(const Vector& x, const Vector& y, int num_peaks, FitType type)
        : m_x(x)
        , m_y(y)
        , m_num_peaks(num_peaks)
        , m_fit_type(type)
    {
        // Gaussian: position, height, width (3)
        // Lorentzian: position, height, width (3)
        // PseudoVoigt: position, height, width, eta (4)
        m_params_per_peak = (type == FitType::PseudoVoigt) ? 4 : 3;
    }

    int inputs() const { return m_num_peaks * m_params_per_peak; }
    int values() const { return m_x.size(); }

    inline double evaluate_peak(double x, const Vector& params, int peak_idx) const
    {
        int offset = peak_idx * m_params_per_peak;
        double pos = params(offset);
        double height = params(offset + 1);
        double width = params(offset + 2);

        if (m_fit_type == FitType::Gaussian) {
            return height * exp(-0.5 * pow((x - pos) / width, 2));
        } else if (m_fit_type == FitType::Lorentzian) {
            return height / (1.0 + pow((x - pos) / width, 2));
        } else { // PseudoVoigt
            double eta = params(offset + 3);
            double gauss = exp(-0.5 * pow((x - pos) / width, 2));
            double lorentz = 1.0 / (1.0 + pow((x - pos) / width, 2));
            return height * (eta * lorentz + (1.0 - eta) * gauss);
        }
    }

    int operator()(const InputType& params, ValueType& fvec) const
    {
        for (int i = 0; i < m_x.size(); ++i) {
            double model = 0.0;
            for (int p = 0; p < m_num_peaks; ++p) {
                model += evaluate_peak(m_x(i), params, p);
            }
            fvec(i) = m_y(i) - model;
        }
        return 0;
    }
};

/*!
 * \brief Fit multiple peaks simultaneously
 * \param spec Input spectrum
 * \param peaks Initial peak guesses
 * \param fit_type Type of peak shape to fit
 * \param max_iter Maximum iterations
 * \return Fit result with parameters and statistics
 */
inline FitResult FitMultiplePeaks(const spectrum& spec,
    const std::vector<Peak>& peaks,
    FitType fit_type = FitType::Gaussian,
    int max_iter = 100)
{
    FitResult result;
    result.converged = false;

    if (peaks.empty()) {
        std::cerr << "Warning: FitMultiplePeaks - no peaks provided" << std::endl;
        return result;
    }

    int num_peaks = peaks.size();
    int params_per_peak = (fit_type == FitType::PseudoVoigt) ? 4 : 3;
    int num_params = num_peaks * params_per_peak;

    // Initial parameter guess
    Vector params(num_params);
    for (size_t i = 0; i < peaks.size(); ++i) {
        int offset = i * params_per_peak;
        params(offset) = spec.X(peaks[i].max);     // Position
        params(offset + 1) = spec.Y(peaks[i].max); // Height
        params(offset + 2) = (spec.X(peaks[i].end) - spec.X(peaks[i].start)) / 4.0; // Width estimate
        if (fit_type == FitType::PseudoVoigt) {
            params(offset + 3) = 0.5; // eta (mixing parameter)
        }
    }

    // Extract region for fitting
    unsigned int fit_start = peaks[0].start;
    unsigned int fit_end = peaks[peaks.size() - 1].end;
    for (const auto& p : peaks) {
        fit_start = std::min(fit_start, p.start);
        fit_end = std::max(fit_end, p.end);
    }

    Vector x_fit(fit_end - fit_start);
    Vector y_fit(fit_end - fit_start);
    for (unsigned int i = fit_start; i < fit_end; ++i) {
        x_fit(i - fit_start) = spec.X(i);
        y_fit(i - fit_start) = spec.Y(i);
    }

    // Setup optimizer
    MultiPeakFitFunctor<double> functor(x_fit, y_fit, num_peaks, fit_type);
    Eigen::NumericalDiff<MultiPeakFitFunctor<double>> numDiff(functor);
    Eigen::LevenbergMarquardt<Eigen::NumericalDiff<MultiPeakFitFunctor<double>>> lm(numDiff);
    lm.parameters.maxfev = max_iter;

    std::cout << "FitMultiplePeaks: Fitting " << num_peaks << " peaks with " << num_params << " parameters" << std::endl;

    Eigen::LevenbergMarquardtSpace::Status status = lm.minimize(params);

    result.converged = (status == Eigen::LevenbergMarquardtSpace::RelativeErrorTooSmall ||
                        status == Eigen::LevenbergMarquardtSpace::RelativeReductionTooSmall);
    result.iterations = 0; // Iteration count not available in this Eigen version
    result.parameters = params;

    // Calculate residuals and chi-squared
    result.residuals = Vector(x_fit.size());
    for (int i = 0; i < x_fit.size(); ++i) {
        double model = 0.0;
        for (int p = 0; p < num_peaks; ++p) {
            model += functor.evaluate_peak(x_fit(i), params, p);
        }
        result.residuals(i) = y_fit(i) - model;
    }

    result.chi_squared = result.residuals.squaredNorm();
    int dof = x_fit.size() - num_params;
    result.reduced_chi_squared = (dof > 0) ? result.chi_squared / dof : 0.0;

    // Approximate covariance matrix from Jacobian
    Eigen::MatrixXd jacobian(x_fit.size(), num_params);
    numDiff.df(params, jacobian);
    result.covariance = (jacobian.transpose() * jacobian).inverse();

    // Update peaks with fitted parameters
    result.fitted_peaks = peaks;
    for (size_t i = 0; i < peaks.size(); ++i) {
        int offset = i * params_per_peak;
        result.fitted_peaks[i].deconv_x = params(offset);
        result.fitted_peaks[i].deconv_y = params(offset + 1);
    }

    std::cout << "FitMultiplePeaks: " << (result.converged ? "Converged" : "Did not converge")
              << " after " << result.iterations << " iterations"
              << " (χ²=" << result.chi_squared << ", reduced χ²=" << result.reduced_chi_squared << ")" << std::endl;

    return result;
}

/*!
 * \brief Extract parameter errors from covariance matrix
 * \param result Fit result
 * \return Vector of parameter standard deviations
 */
inline Vector ExtractParameterErrors(const FitResult& result)
{
    Vector errors(result.parameters.size());
    for (int i = 0; i < result.parameters.size(); ++i) {
        errors(i) = sqrt(std::max(0.0, result.covariance(i, i)));
    }
    return errors;
}

} // namespace PeakPick

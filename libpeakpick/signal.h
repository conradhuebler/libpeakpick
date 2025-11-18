/*
 * <Signal processing algorithms>
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
#include <unsupported/Eigen/FFT>

#include <cmath>
#include <complex>
#include <iostream>
#include <vector>

#include "analyse.h"
#include "spectrum.h"

typedef Eigen::VectorXd Vector;

namespace PeakPick {

/*!
 * \brief Calculate first derivative of spectrum
 * \param spec Input spectrum
 * \return First derivative spectrum
 */
inline spectrum FirstDerivative(const spectrum& spec)
{
    if (spec.size() < 2) {
        std::cerr << "Warning: FirstDerivative - spectrum too small" << std::endl;
        return spec;
    }

    Vector deriv(spec.size());

    // Forward difference at start
    deriv(0) = (spec.Y(1) - spec.Y(0)) / (spec.X(1) - spec.X(0));

    // Central difference in middle
    for (unsigned int i = 1; i < spec.size() - 1; ++i) {
        double dx = spec.X(i + 1) - spec.X(i - 1);
        double dy = spec.Y(i + 1) - spec.Y(i - 1);
        deriv(i) = dy / dx;
    }

    // Backward difference at end
    int n = spec.size() - 1;
    deriv(n) = (spec.Y(n) - spec.Y(n - 1)) / (spec.X(n) - spec.X(n - 1));

    return spectrum(spec.x(), deriv);
}

/*!
 * \brief Calculate second derivative of spectrum
 * \param spec Input spectrum
 * \return Second derivative spectrum
 */
inline spectrum SecondDerivative(const spectrum& spec)
{
    if (spec.size() < 3) {
        std::cerr << "Warning: SecondDerivative - spectrum too small" << std::endl;
        return spec;
    }

    Vector deriv2(spec.size());

    // Three-point formula for interior points
    for (unsigned int i = 1; i < spec.size() - 1; ++i) {
        double dx = spec.X(i) - spec.X(i - 1); // Assume uniform spacing
        double d2y = spec.Y(i - 1) - 2.0 * spec.Y(i) + spec.Y(i + 1);
        deriv2(i) = d2y / (dx * dx);
    }

    // Boundary points (forward/backward differences)
    deriv2(0) = deriv2(1);
    deriv2(spec.size() - 1) = deriv2(spec.size() - 2);

    return spectrum(spec.x(), deriv2);
}

/*!
 * \brief Pick peaks using second derivative (peaks are minima in 2nd derivative)
 * \param spec Input spectrum
 * \param threshold Minimum height in original spectrum
 * \param deriv_threshold Maximum value in 2nd derivative (should be negative)
 * \return List of detected peaks
 */
inline std::vector<Peak> PickPeaksSecondDerivative(const spectrum& spec,
    double threshold = 0.0,
    double deriv_threshold = -0.1)
{
    spectrum deriv2 = SecondDerivative(spec);
    std::vector<Peak> peaks;

    std::cout << "PickPeaksSecondDerivative: Searching for minima in 2nd derivative" << std::endl;

    for (unsigned int i = 2; i < spec.size() - 2; ++i) {
        // Check if local minimum in 2nd derivative
        if (deriv2.Y(i) < deriv2.Y(i - 1) &&
            deriv2.Y(i) < deriv2.Y(i + 1) &&
            deriv2.Y(i) < deriv_threshold) {

            // Check if above threshold in original spectrum
            if (spec.Y(i) > threshold) {
                Peak peak;
                peak.max = i;

                // Find peak boundaries (where 2nd deriv crosses zero)
                peak.start = i;
                for (int j = i - 1; j >= 0; --j) {
                    if (deriv2.Y(j) > 0) {
                        peak.start = j;
                        break;
                    }
                }

                peak.end = i;
                for (unsigned int j = i + 1; j < spec.size(); ++j) {
                    if (deriv2.Y(j) > 0) {
                        peak.end = j;
                        break;
                    }
                }

                peak.int_start = peak.start;
                peak.int_end = peak.end;

                peaks.push_back(peak);
            }
        }
    }

    std::cout << "PickPeaksSecondDerivative: Found " << peaks.size() << " peaks" << std::endl;
    return peaks;
}

/*!
 * \brief Apply low-pass filter to spectrum
 * \param spec Input spectrum
 * \param cutoff_freq Cutoff frequency (normalized to Nyquist)
 * \return Filtered spectrum
 */
inline spectrum LowPassFilter(const spectrum& spec, double cutoff_freq = 0.1)
{
    if (cutoff_freq <= 0.0 || cutoff_freq >= 1.0) {
        std::cerr << "Warning: LowPassFilter - cutoff_freq must be in (0, 1), using 0.1" << std::endl;
        cutoff_freq = 0.1;
    }

    Eigen::FFT<double> fft;
    std::vector<std::complex<double>> freq_domain;

    // Forward FFT
    std::vector<double> time_domain(spec.size());
    for (unsigned int i = 0; i < spec.size(); ++i) {
        time_domain[i] = spec.Y(i);
    }

    fft.fwd(freq_domain, time_domain);

    // Apply filter in frequency domain
    int cutoff_index = static_cast<int>(freq_domain.size() * cutoff_freq);
    for (size_t i = cutoff_index; i < freq_domain.size() - cutoff_index; ++i) {
        freq_domain[i] = std::complex<double>(0.0, 0.0);
    }

    // Inverse FFT
    std::vector<double> filtered;
    fft.inv(filtered, freq_domain);

    Vector filtered_y(filtered.size());
    for (size_t i = 0; i < filtered.size(); ++i) {
        filtered_y(i) = filtered[i];
    }

    std::cout << "LowPassFilter: Applied with cutoff = " << cutoff_freq << std::endl;

    return spectrum(spec.x(), filtered_y);
}

/*!
 * \brief Compute power spectrum (FFT magnitude)
 * \param spec Input spectrum
 * \return Power spectrum
 */
inline spectrum PowerSpectrum(const spectrum& spec)
{
    Eigen::FFT<double> fft;
    std::vector<std::complex<double>> freq_domain;

    std::vector<double> time_domain(spec.size());
    for (unsigned int i = 0; i < spec.size(); ++i) {
        time_domain[i] = spec.Y(i);
    }

    fft.fwd(freq_domain, time_domain);

    // Calculate magnitude
    Vector power(freq_domain.size() / 2); // Only positive frequencies
    Vector freq(freq_domain.size() / 2);

    double step = spec.Step();
    double sampling_rate = 1.0 / step;

    for (size_t i = 0; i < freq_domain.size() / 2; ++i) {
        power(i) = std::abs(freq_domain[i]);
        freq(i) = i * sampling_rate / freq_domain.size();
    }

    return spectrum(freq, power);
}

/*!
 * \brief Convolve two spectra
 * \param spec1 First spectrum
 * \param spec2 Second spectrum (kernel)
 * \return Convolved spectrum
 */
inline spectrum Convolve(const spectrum& spec1, const spectrum& spec2)
{
    // Simple convolution (not FFT-based for simplicity)
    int n1 = spec1.size();
    int n2 = spec2.size();
    int n_out = n1 + n2 - 1;

    Vector result = Vector::Zero(n_out);

    for (int i = 0; i < n1; ++i) {
        for (int j = 0; j < n2; ++j) {
            result(i + j) += spec1.Y(i) * spec2.Y(j);
        }
    }

    // Create X-axis for result
    Vector x_out(n_out);
    double step = spec1.Step();
    for (int i = 0; i < n_out; ++i) {
        x_out(i) = spec1.XMin() + i * step;
    }

    return spectrum(x_out, result);
}

} // namespace PeakPick

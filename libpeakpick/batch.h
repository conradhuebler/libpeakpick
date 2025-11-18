/*
 * <Batch processing framework>
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

#include <functional>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "advanced.h"
#include "analyse.h"
#include "baseline.h"
#include "logger.h"
#include "peakpick.h"
#include "spectrum.h"

namespace PeakPick {

/*!
 * \brief Results from batch spectrum processing
 */
struct BatchResults {
    std::vector<spectrum> spectra;
    std::vector<std::vector<Peak>> all_peaks;
    std::vector<std::string> filenames;
    std::vector<bool> success;
    int total_processed = 0;
    int total_failed = 0;
};

/*!
 * \brief Batch spectrum processor with fluent interface
 */
class BatchProcessor {
private:
    std::vector<spectrum> spectra_;
    std::vector<std::string> filenames_;
    std::vector<std::function<void(spectrum&)>> operations_;
    bool verbose_ = true;

public:
    BatchProcessor(bool verbose = true)
        : verbose_(verbose) {}

    /*!
     * \brief Add spectrum from file
     */
    BatchProcessor& addFile(const std::string& filename)
    {
        filenames_.push_back(filename);
        if (verbose_)
            Logger::info("Added file: " + filename, "BatchProcessor");
        return *this;
    }

    /*!
     * \brief Add multiple files
     */
    BatchProcessor& addFiles(const std::vector<std::string>& filenames)
    {
        for (const auto& fn : filenames) {
            addFile(fn);
        }
        return *this;
    }

    /*!
     * \brief Add spectrum directly
     */
    BatchProcessor& addSpectrum(const spectrum& spec, const std::string& name = "")
    {
        spectra_.push_back(spec);
        filenames_.push_back(name.empty() ? "spectrum_" + std::to_string(spectra_.size()) : name);
        return *this;
    }

    /*!
     * \brief Apply smoothing to all spectra
     */
    BatchProcessor& smooth(unsigned int points)
    {
        operations_.push_back([points](spectrum& spec) {
            SmoothFunction(&spec, points);
        });
        if (verbose_)
            Logger::info("Added smoothing operation (points=" + std::to_string(points) + ")", "BatchProcessor");
        return *this;
    }

    /*!
     * \brief Apply normalization
     */
    BatchProcessor& normalize(double min = 0.0, double max = 1.0)
    {
        operations_.push_back([min, max](spectrum& spec) {
            Normalise(&spec, min, max);
        });
        if (verbose_)
            Logger::info("Added normalization operation", "BatchProcessor");
        return *this;
    }

    /*!
     * \brief Apply AsLS baseline correction
     */
    BatchProcessor& baselineCorrect(double lambda = 1e6, double p = 0.01)
    {
        operations_.push_back([lambda, p](spectrum& spec) {
            spec = BaselineAsLS(spec, lambda, p);
        });
        if (verbose_)
            Logger::info("Added AsLS baseline correction", "BatchProcessor");
        return *this;
    }

    /*!
     * \brief Apply custom operation
     */
    BatchProcessor& apply(std::function<void(spectrum&)> operation)
    {
        operations_.push_back(operation);
        if (verbose_)
            Logger::info("Added custom operation", "BatchProcessor");
        return *this;
    }

    /*!
     * \brief Process all spectra and return results
     */
    BatchResults process()
    {
        BatchResults results;

        if (verbose_)
            Logger::info("Starting batch processing of " + std::to_string(filenames_.size()) + " files", "BatchProcessor");

        // Load files if not already loaded
        if (spectra_.empty() && !filenames_.empty()) {
            for (const auto& fn : filenames_) {
                spectrum spec = loadFromFile(fn);
                spectra_.push_back(spec);
            }
        }

        results.filenames = filenames_;
        results.spectra.reserve(spectra_.size());
        results.success.resize(spectra_.size(), true);

        // Process each spectrum
        for (size_t i = 0; i < spectra_.size(); ++i) {
            if (verbose_) {
                int progress = static_cast<int>((i + 1) * 100.0 / spectra_.size());
                Logger::info("Processing " + filenames_[i] + " (" + std::to_string(progress) + "%)", "BatchProcessor");
            }

            try {
                spectrum processed = spectra_[i];

                // Apply all operations
                for (const auto& op : operations_) {
                    op(processed);
                }

                results.spectra.push_back(processed);
                results.total_processed++;
            } catch (const std::exception& e) {
                Logger::error("Failed to process " + filenames_[i] + ": " + e.what(), "BatchProcessor");
                results.success[i] = false;
                results.total_failed++;
                results.spectra.push_back(spectrum()); // Empty placeholder
            } catch (...) {
                Logger::error("Failed to process " + filenames_[i], "BatchProcessor");
                results.success[i] = false;
                results.total_failed++;
                results.spectra.push_back(spectrum());
            }
        }

        if (verbose_) {
            Logger::info("Batch processing complete: " +
                             std::to_string(results.total_processed) + " succeeded, " +
                             std::to_string(results.total_failed) + " failed",
                "BatchProcessor");
        }

        return results;
    }

    /*!
     * \brief Process and find peaks in all spectra
     */
    BatchResults processAndFindPeaks(double threshold, double min_prominence = 0.0, unsigned int min_distance = 1)
    {
        BatchResults results = process();

        if (verbose_)
            Logger::info("Finding peaks in processed spectra", "BatchProcessor");

        results.all_peaks.reserve(results.spectra.size());

        for (size_t i = 0; i < results.spectra.size(); ++i) {
            if (!results.success[i]) {
                results.all_peaks.push_back(std::vector<Peak>());
                continue;
            }

            try {
                auto peaks = PickPeaksAdvanced(&results.spectra[i], threshold, min_prominence, min_distance);
                results.all_peaks.push_back(peaks);

                if (verbose_)
                    Logger::info("Found " + std::to_string(peaks.size()) + " peaks in " + results.filenames[i], "BatchProcessor");
            } catch (...) {
                Logger::error("Failed to find peaks in " + results.filenames[i], "BatchProcessor");
                results.all_peaks.push_back(std::vector<Peak>());
            }
        }

        return results;
    }

    /*!
     * \brief Reset processor
     */
    void reset()
    {
        spectra_.clear();
        filenames_.clear();
        operations_.clear();
    }
};

} // namespace PeakPick

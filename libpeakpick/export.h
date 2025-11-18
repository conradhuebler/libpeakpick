/*
 * <Export formats for spectra and peaks>
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

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "analyse.h"
#include "spectrum.h"

namespace PeakPick {

struct CSVOptions {
    bool include_header = true;
    char delimiter = ',';
    int precision = 6;
};

/*!
 * \brief Export spectrum to CSV format
 * \param spec Spectrum to export
 * \param filename Output filename
 * \param options CSV formatting options
 * \return true if successful
 */
inline bool ExportCSV(const spectrum& spec, const std::string& filename, const CSVOptions& options = CSVOptions())
{
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file for CSV export: " << filename << std::endl;
        return false;
    }

    file << std::fixed << std::setprecision(options.precision);

    // Header
    if (options.include_header) {
        file << "X" << options.delimiter << "Y" << std::endl;
    }

    // Data
    for (unsigned int i = 0; i < spec.size(); ++i) {
        file << spec.X(i) << options.delimiter << spec.Y(i) << std::endl;
    }

    file.close();
    std::cout << "Exported spectrum to CSV: " << filename << " (" << spec.size() << " points)" << std::endl;
    return true;
}

/*!
 * \brief Export peaks to CSV format
 * \param peaks Peak list to export
 * \param filename Output filename
 * \param options CSV formatting options
 * \return true if successful
 */
inline bool ExportPeaksCSV(const std::vector<Peak>& peaks, const std::string& filename, const CSVOptions& options = CSVOptions())
{
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file for peaks CSV export: " << filename << std::endl;
        return false;
    }

    file << std::fixed << std::setprecision(options.precision);

    // Header
    if (options.include_header) {
        file << "Index" << options.delimiter
             << "Start" << options.delimiter
             << "Max" << options.delimiter
             << "End" << options.delimiter
             << "Position" << options.delimiter
             << "Height" << options.delimiter
             << "Area" << std::endl;
    }

    // Data
    for (size_t i = 0; i < peaks.size(); ++i) {
        file << i << options.delimiter
             << peaks[i].start << options.delimiter
             << peaks[i].max << options.delimiter
             << peaks[i].end << options.delimiter
             << peaks[i].deconv_x << options.delimiter
             << peaks[i].deconv_y << options.delimiter
             << peaks[i].integ_num << std::endl;
    }

    file.close();
    std::cout << "Exported " << peaks.size() << " peaks to CSV: " << filename << std::endl;
    return true;
}

/*!
 * \brief Export spectrum to JSON format
 * \param spec Spectrum to export
 * \param filename Output filename
 * \param include_metadata Include spectrum statistics
 * \return true if successful
 */
inline bool ExportJSON(const spectrum& spec, const std::string& filename, bool include_metadata = true)
{
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file for JSON export: " << filename << std::endl;
        return false;
    }

    file << std::fixed << std::setprecision(6);

    file << "{" << std::endl;

    // Metadata
    if (include_metadata) {
        file << "  \"metadata\": {" << std::endl;
        file << "    \"size\": " << spec.size() << "," << std::endl;
        file << "    \"mean\": " << spec.Mean() << "," << std::endl;
        file << "    \"max\": " << spec.Max() << "," << std::endl;
        file << "    \"min\": " << spec.Min() << "," << std::endl;
        file << "    \"stddev\": " << spec.StdDev() << "," << std::endl;
        file << "    \"x_min\": " << spec.XMin() << "," << std::endl;
        file << "    \"x_max\": " << spec.XMax() << std::endl;
        file << "  }," << std::endl;
    }

    // Data
    file << "  \"data\": {" << std::endl;
    file << "    \"x\": [";
    for (unsigned int i = 0; i < spec.size(); ++i) {
        if (i > 0)
            file << ", ";
        file << spec.X(i);
    }
    file << "]," << std::endl;

    file << "    \"y\": [";
    for (unsigned int i = 0; i < spec.size(); ++i) {
        if (i > 0)
            file << ", ";
        file << spec.Y(i);
    }
    file << "]" << std::endl;
    file << "  }" << std::endl;

    file << "}" << std::endl;

    file.close();
    std::cout << "Exported spectrum to JSON: " << filename << " (" << spec.size() << " points)" << std::endl;
    return true;
}

/*!
 * \brief Export peaks to JSON format
 * \param peaks Peak list to export
 * \param filename Output filename
 * \return true if successful
 */
inline bool ExportPeaksJSON(const std::vector<Peak>& peaks, const std::string& filename)
{
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file for peaks JSON export: " << filename << std::endl;
        return false;
    }

    file << std::fixed << std::setprecision(6);

    file << "{" << std::endl;
    file << "  \"peaks\": [" << std::endl;

    for (size_t i = 0; i < peaks.size(); ++i) {
        file << "    {" << std::endl;
        file << "      \"index\": " << i << "," << std::endl;
        file << "      \"start\": " << peaks[i].start << "," << std::endl;
        file << "      \"max\": " << peaks[i].max << "," << std::endl;
        file << "      \"end\": " << peaks[i].end << "," << std::endl;
        file << "      \"position\": " << peaks[i].deconv_x << "," << std::endl;
        file << "      \"height\": " << peaks[i].deconv_y << "," << std::endl;
        file << "      \"area_numerical\": " << peaks[i].integ_num << "," << std::endl;
        file << "      \"area_analytical\": " << peaks[i].integ_analyt << std::endl;
        file << "    }";
        if (i < peaks.size() - 1)
            file << ",";
        file << std::endl;
    }

    file << "  ]," << std::endl;
    file << "  \"num_peaks\": " << peaks.size() << std::endl;
    file << "}" << std::endl;

    file.close();
    std::cout << "Exported " << peaks.size() << " peaks to JSON: " << filename << std::endl;
    return true;
}

/*!
 * \brief Import spectrum from CSV
 * \param filename Input filename
 * \param has_header Whether file has header row
 * \param delimiter CSV delimiter
 * \return Loaded spectrum
 */
inline spectrum ImportCSV(const std::string& filename, bool has_header = true, char delimiter = ',')
{
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open CSV file: " << filename << std::endl;
        return spectrum();
    }

    std::vector<double> x_data, y_data;
    std::string line;

    // Skip header if present
    if (has_header) {
        std::getline(file, line);
    }

    // Read data
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#')
            continue;

        std::stringstream ss(line);
        std::string item;
        std::vector<double> row;

        while (std::getline(ss, item, delimiter)) {
            try {
                row.push_back(std::stod(item));
            } catch (...) {
                continue;
            }
        }

        if (row.size() >= 2) {
            x_data.push_back(row[0]);
            y_data.push_back(row[1]);
        }
    }

    file.close();

    if (x_data.empty()) {
        std::cerr << "Warning: No data read from CSV file" << std::endl;
        return spectrum();
    }

    Vector x = Vector::Map(&x_data[0], x_data.size());
    Vector y = Vector::Map(&y_data[0], y_data.size());

    std::cout << "Imported spectrum from CSV: " << filename << " (" << x_data.size() << " points)" << std::endl;

    return spectrum(x, y);
}

} // namespace PeakPick

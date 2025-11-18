/*
 * <Command-line interface for libpeakpick>
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

#include "libpeakpick/advanced.h"
#include "libpeakpick/analyse.h"
#include "libpeakpick/batch.h"
#include "libpeakpick/comparison.h"
#include "libpeakpick/export.h"
#include "libpeakpick/peakpick.h"
#include "libpeakpick/signal.h"

#include <iostream>
#include <string>
#include <vector>

void printUsage(const char* program_name)
{
    std::cout << "libpeakpick CLI - Spectral Analysis Tool" << std::endl;
    std::cout << "==========================================" << std::endl << std::endl;
    std::cout << "Usage: " << program_name << " <command> [options]" << std::endl << std::endl;
    std::cout << "Commands:" << std::endl;
    std::cout << "  analyze <file>           Analyze single spectrum" << std::endl;
    std::cout << "  batch <files...>         Batch process multiple files" << std::endl;
    std::cout << "  compare <file1> <file2>  Compare two spectra" << std::endl;
    std::cout << "  convert <file> <format>  Convert spectrum to format (csv/json)" << std::endl;
    std::cout << "  peaks <file>             Find and export peaks" << std::endl;
    std::cout << "  help                     Show this help message" << std::endl << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "  --smooth <N>             Apply Savitzky-Golay smoothing (N points)" << std::endl;
    std::cout << "  --baseline <lambda>      Apply AsLS baseline correction" << std::endl;
    std::cout << "  --normalize              Normalize to [0,1]" << std::endl;
    std::cout << "  --threshold <val>        Peak detection threshold" << std::endl;
    std::cout << "  --prominence <val>       Minimum peak prominence" << std::endl;
    std::cout << "  --output <file>          Output filename" << std::endl;
    std::cout << "  --verbose                Verbose output" << std::endl << std::endl;
    std::cout << "Examples:" << std::endl;
    std::cout << "  " << program_name << " analyze data.txt --smooth 5 --normalize" << std::endl;
    std::cout << "  " << program_name << " batch *.txt --baseline 1e6 --output results/" << std::endl;
    std::cout << "  " << program_name << " peaks data.txt --threshold 0.5 --prominence 1.0" << std::endl;
    std::cout << "  " << program_name << " compare spec1.txt spec2.txt" << std::endl;
}

int cmdAnalyze(int argc, char* argv[])
{
    if (argc < 3) {
        std::cerr << "Error: analyze command requires input file" << std::endl;
        return 1;
    }

    std::string filename = argv[2];
    bool do_smooth = false;
    bool do_baseline = false;
    bool do_normalize = false;
    int smooth_points = 5;
    double baseline_lambda = 1e6;
    std::string output_file;

    // Parse options
    for (int i = 3; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "--smooth" && i + 1 < argc) {
            do_smooth = true;
            smooth_points = std::stoi(argv[++i]);
        } else if (arg == "--baseline" && i + 1 < argc) {
            do_baseline = true;
            baseline_lambda = std::stod(argv[++i]);
        } else if (arg == "--normalize") {
            do_normalize = true;
        } else if (arg == "--output" && i + 1 < argc) {
            output_file = argv[++i];
        }
    }

    std::cout << "Analyzing: " << filename << std::endl;

    // Load spectrum
    PeakPick::spectrum spec = loadFromFile(filename);

    if (spec.size() == 0) {
        std::cerr << "Error: Could not load spectrum from " << filename << std::endl;
        return 1;
    }

    std::cout << "Loaded spectrum: " << spec.size() << " points" << std::endl;
    std::cout << "Statistics:" << std::endl;
    std::cout << "  Mean:   " << spec.Mean() << std::endl;
    std::cout << "  Max:    " << spec.Max() << " at " << spec.PosMax() << std::endl;
    std::cout << "  Min:    " << spec.Min() << " at " << spec.PosMin() << std::endl;
    std::cout << "  StdDev: " << spec.StdDev() << std::endl;
    std::cout << "  Noise:  " << PeakPick::EstimateNoise(spec) << std::endl;

    // Apply processing
    if (do_smooth) {
        std::cout << "\nApplying smoothing (" << smooth_points << " points)..." << std::endl;
        PeakPick::SmoothFunction(&spec, smooth_points);
    }

    if (do_baseline) {
        std::cout << "Applying AsLS baseline correction (lambda=" << baseline_lambda << ")..." << std::endl;
        spec = PeakPick::BaselineAsLS(spec, baseline_lambda);
    }

    if (do_normalize) {
        std::cout << "Normalizing..." << std::endl;
        PeakPick::Normalise(&spec);
    }

    // Save if output specified
    if (!output_file.empty()) {
        PeakPick::ExportJSON(spec, output_file);
    }

    std::cout << "\nAnalysis complete." << std::endl;
    return 0;
}

int cmdBatch(int argc, char* argv[])
{
    if (argc < 3) {
        std::cerr << "Error: batch command requires input files" << std::endl;
        return 1;
    }

    std::vector<std::string> files;
    bool do_smooth = false;
    bool do_baseline = false;
    bool do_normalize = false;
    int smooth_points = 5;
    double baseline_lambda = 1e6;
    std::string output_dir = ".";

    // Parse files and options
    for (int i = 2; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "--smooth" && i + 1 < argc) {
            do_smooth = true;
            smooth_points = std::stoi(argv[++i]);
        } else if (arg == "--baseline" && i + 1 < argc) {
            do_baseline = true;
            baseline_lambda = std::stod(argv[++i]);
        } else if (arg == "--normalize") {
            do_normalize = true;
        } else if (arg == "--output" && i + 1 < argc) {
            output_dir = argv[++i];
        } else if (arg[0] != '-') {
            files.push_back(arg);
        }
    }

    std::cout << "Batch processing " << files.size() << " files..." << std::endl;

    PeakPick::BatchProcessor processor;
    processor.addFiles(files);

    if (do_smooth)
        processor.smooth(smooth_points);
    if (do_baseline)
        processor.baselineCorrect(baseline_lambda);
    if (do_normalize)
        processor.normalize();

    auto results = processor.process();

    std::cout << "\nBatch processing complete:" << std::endl;
    std::cout << "  Succeeded: " << results.total_processed << std::endl;
    std::cout << "  Failed:    " << results.total_failed << std::endl;

    return 0;
}

int cmdCompare(int argc, char* argv[])
{
    if (argc < 4) {
        std::cerr << "Error: compare command requires two input files" << std::endl;
        return 1;
    }

    std::string file1 = argv[2];
    std::string file2 = argv[3];

    std::cout << "Comparing: " << file1 << " vs " << file2 << std::endl;

    PeakPick::spectrum spec1 = loadFromFile(file1);
    PeakPick::spectrum spec2 = loadFromFile(file2);

    if (spec1.size() == 0 || spec2.size() == 0) {
        std::cerr << "Error: Could not load one or both spectra" << std::endl;
        return 1;
    }

    double rmsd = PeakPick::RMSD(spec1, spec2);
    double corr = PeakPick::CorrelationCoefficient(spec1, spec2);
    double cos_sim = PeakPick::CosineSimilarity(spec1, spec2);
    double angle = PeakPick::SpectralAngle(spec1, spec2);
    double mae = PeakPick::MAE(spec1, spec2);

    std::cout << "\nComparison Metrics:" << std::endl;
    std::cout << "  RMSD:                 " << rmsd << std::endl;
    std::cout << "  Correlation:          " << corr << std::endl;
    std::cout << "  Cosine Similarity:    " << cos_sim << std::endl;
    std::cout << "  Spectral Angle:       " << angle << "°" << std::endl;
    std::cout << "  Mean Absolute Error:  " << mae << std::endl;

    return 0;
}

int cmdPeaks(int argc, char* argv[])
{
    if (argc < 3) {
        std::cerr << "Error: peaks command requires input file" << std::endl;
        return 1;
    }

    std::string filename = argv[2];
    double threshold = 0.5;
    double prominence = 0.0;
    std::string output_file;

    for (int i = 3; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "--threshold" && i + 1 < argc) {
            threshold = std::stod(argv[++i]);
        } else if (arg == "--prominence" && i + 1 < argc) {
            prominence = std::stod(argv[++i]);
        } else if (arg == "--output" && i + 1 < argc) {
            output_file = argv[++i];
        }
    }

    std::cout << "Finding peaks in: " << filename << std::endl;

    PeakPick::spectrum spec = loadFromFile(filename);

    if (spec.size() == 0) {
        std::cerr << "Error: Could not load spectrum" << std::endl;
        return 1;
    }

    auto peaks = PeakPick::PickPeaksAdvanced(&spec, threshold, prominence, 5);

    std::cout << "Found " << peaks.size() << " peaks" << std::endl;

    for (size_t i = 0; i < peaks.size(); ++i) {
        std::cout << "  Peak " << i + 1 << ": position=" << spec.X(peaks[i].max)
                  << ", height=" << spec.Y(peaks[i].max) << std::endl;
    }

    if (!output_file.empty()) {
        PeakPick::ExportPeaksJSON(peaks, output_file);
        std::cout << "\nPeaks exported to: " << output_file << std::endl;
    }

    return 0;
}

int main(int argc, char* argv[])
{
    if (argc < 2) {
        printUsage(argv[0]);
        return 0;
    }

    std::string command = argv[1];

    if (command == "analyze") {
        return cmdAnalyze(argc, argv);
    } else if (command == "batch") {
        return cmdBatch(argc, argv);
    } else if (command == "compare") {
        return cmdCompare(argc, argv);
    } else if (command == "peaks") {
        return cmdPeaks(argc, argv);
    } else if (command == "help" || command == "--help" || command == "-h") {
        printUsage(argv[0]);
        return 0;
    } else {
        std::cerr << "Unknown command: " << command << std::endl;
        printUsage(argv[0]);
        return 1;
    }
}

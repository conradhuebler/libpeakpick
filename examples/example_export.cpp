// Export and import example
#include "libpeakpick/advanced.h"
#include "libpeakpick/export.h"
#include "libpeakpick/peakpick.h"
#include <iostream>

int main() {
    // Load and process spectrum
    PeakPick::spectrum spec = loadFromFile("samples/lorentzian_1");

    // Normalize
    PeakPick::Normalise(&spec);

    // Export to different formats
    std::cout << "Exporting spectrum..." << std::endl;

    // Export as CSV
    PeakPick::CSVOptions csv_opts;
    csv_opts.precision = 8;
    PeakPick::ExportCSV(spec, "spectrum_output.csv", csv_opts);

    // Export as JSON with metadata
    PeakPick::ExportJSON(spec, "spectrum_output.json", true);

    // Find peaks
    auto peaks = PeakPick::PickPeaksAdvanced(&spec, 0.1, 0.5, 10);

    // Export peaks
    PeakPick::ExportPeaksCSV(peaks, "peaks_output.csv", csv_opts);
    PeakPick::ExportPeaksJSON(peaks, "peaks_output.json");

    std::cout << "\nExports complete!" << std::endl;

    // Import back
    std::cout << "\nImporting CSV..." << std::endl;
    PeakPick::spectrum imported = PeakPick::ImportCSV("spectrum_output.csv");
    std::cout << "Imported spectrum: " << imported.size() << " points" << std::endl;

    return 0;
}

// Batch processing example
#include "libpeakpick/batch.h"
#include "libpeakpick/export.h"
#include "libpeakpick/logger.h"

int main() {
    // Set logging level
    PeakPick::Logger::setLevel(PeakPick::LogLevel::Info);

    // Create batch processor
    PeakPick::BatchProcessor processor;

    // Configure processing pipeline
    auto results = processor
                       .addFile("samples/lorentzian_1")
                       .addFile("samples/lorentzian_2")
                       .smooth(3)
                       .normalize()
                       .processAndFindPeaks(0.1, 0.5, 5);

    // Export results
    for (size_t i = 0; i < results.spectra.size(); ++i) {
        if (results.success[i]) {
            std::string output = "output_" + std::to_string(i) + ".json";
            PeakPick::ExportJSON(results.spectra[i], output);
            PeakPick::ExportPeaksJSON(results.all_peaks[i], "peaks_" + std::to_string(i) + ".json");
        }
    }

    return 0;
}

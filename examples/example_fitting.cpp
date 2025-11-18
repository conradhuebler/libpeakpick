// Multi-peak fitting example
#include "libpeakpick/advanced.h"
#include "libpeakpick/fitting.h"
#include "libpeakpick/peakpick.h"
#include <iostream>

int main() {
    // Load spectrum
    PeakPick::spectrum spec = loadFromFile("samples/lorentzian_1");

    // Find peaks
    auto peaks = PeakPick::PickPeaksAdvanced(&spec, 0.1, 0.5, 10);

    std::cout << "Found " << peaks.size() << " peaks" << std::endl;

    // Fit multiple peaks simultaneously
    if (peaks.size() > 0) {
        auto result = PeakPick::FitMultiplePeaks(spec, peaks, PeakPick::FitType::Gaussian);

        std::cout << "\nFit Results:" << std::endl;
        std::cout << "  Converged: " << (result.converged ? "Yes" : "No") << std::endl;
        std::cout << "  Iterations: " << result.iterations << std::endl;
        std::cout << "  Chi-squared: " << result.chi_squared << std::endl;
        std::cout << "  Reduced chi-squared: " << result.reduced_chi_squared << std::endl;

        // Print fitted parameters
        auto errors = PeakPick::ExtractParameterErrors(result);
        std::cout << "\nFitted Parameters:" << std::endl;
        for (size_t i = 0; i < peaks.size(); ++i) {
            std::cout << "  Peak " << i + 1 << ":" << std::endl;
            std::cout << "    Position: " << result.parameters(i * 3) << " ± " << errors(i * 3) << std::endl;
            std::cout << "    Height:   " << result.parameters(i * 3 + 1) << " ± " << errors(i * 3 + 1) << std::endl;
            std::cout << "    Width:    " << result.parameters(i * 3 + 2) << " ± " << errors(i * 3 + 2) << std::endl;
        }
    }

    return 0;
}

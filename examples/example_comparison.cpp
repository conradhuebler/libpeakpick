// Spectrum comparison example
#include "libpeakpick/comparison.h"
#include "libpeakpick/peakpick.h"
#include <iostream>

int main() {
    // Load two spectra
    PeakPick::spectrum spec1 = loadFromFile("samples/lorentzian_1");
    PeakPick::spectrum spec2 = loadFromFile("samples/lorentzian_2");

    std::cout << "Comparing two spectra..." << std::endl;

    // Calculate similarity metrics
    double rmsd = PeakPick::RMSD(spec1, spec2);
    double corr = PeakPick::CorrelationCoefficient(spec1, spec2);
    double cos_sim = PeakPick::CosineSimilarity(spec1, spec2);
    double angle = PeakPick::SpectralAngle(spec1, spec2);
    double mae = PeakPick::MAE(spec1, spec2);

    std::cout << "\nSimilarity Metrics:" << std::endl;
    std::cout << "  RMSD:              " << rmsd << std::endl;
    std::cout << "  Correlation:       " << corr << std::endl;
    std::cout << "  Cosine Similarity: " << cos_sim << std::endl;
    std::cout << "  Spectral Angle:    " << angle << "°" << std::endl;
    std::cout << "  MAE:               " << mae << std::endl;

    // Align spectra
    std::cout << "\nAligning spectra..." << std::endl;
    PeakPick::spectrum aligned = PeakPick::AlignSpectra(spec1, spec2, 5.0, 0.5);

    // Recalculate after alignment
    double rmsd_aligned = PeakPick::RMSD(spec1, aligned);
    std::cout << "RMSD after alignment: " << rmsd_aligned << std::endl;

    return 0;
}

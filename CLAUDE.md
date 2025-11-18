# CLAUDE.md - AI Assistant Guide for libpeakpick

## Project Overview

**libpeakpick** is a header-only C++ library for peak picking and spectral analysis. It provides tools for analyzing spectroscopic data including baseline correction, peak detection, Gaussian/Lorentzian fitting, and various signal processing operations.

- **Language**: C++11
- **License**: GPLv3
- **Build System**: CMake
- **Architecture**: Header-only library (no compilation required for library itself)
- **Dependencies**: Eigen (included as git submodule)
- **Parallelization**: OpenMP support for performance-critical operations
- **Status**: Pre-alpha (active development, API may change)

### Key Features
- Spectral analysis and statistics (mean, max, min, stddev)
- Peak picking and detection (basic and advanced with prominence)
- Baseline correction algorithms (including AsLS)
- Gaussian/Lorentzian/PseudoVoigt curve fitting (single and multi-peak)
- Savitzky-Golay smoothing and signal processing
- FFT-based operations and derivatives
- Spectrum comparison metrics (RMSD, correlation, cosine similarity)
- CSV and JSON import/export
- Batch processing framework
- Command-line interface tool
- Logging system with color support
- Linear regression and deconvolution
- Python bindings guide

## Repository Structure

```
libpeakpick/
├── libpeakpick/           # Header-only library files (main implementation)
│   ├── analyse.h          # Analysis algorithms (IMPROVED)
│   ├── advanced.h         # Advanced algorithms (NEW: AsLS, prominence, splines)
│   ├── baseline.h         # Baseline correction methods (IMPROVED)
│   ├── batch.h            # Batch processing framework (NEW)
│   ├── comparison.h       # Spectrum comparison metrics (NEW)
│   ├── deconvulate.h      # Deconvolution algorithms
│   ├── export.h           # CSV/JSON import/export (NEW)
│   ├── fitting.h          # Multi-peak fitting (NEW)
│   ├── glfit.h            # Gaussian/Lorentzian fitting
│   ├── logger.h           # Logging system (NEW)
│   ├── mathhelper.h       # Mathematical utilities
│   ├── nxlinregress.h     # N-dimensional linear regression
│   ├── peakpick.h         # Main header with file loading utilities
│   ├── savitzky.h         # Savitzky-Golay filter implementation
│   ├── signal.h           # Signal processing (FFT, derivatives, filters) (NEW)
│   ├── spectrum.h         # Core spectrum class and operations (IMPROVED)
│   └── utilities.h        # Utility functions (save, SNR, FWHM, etc.)
├── tests/                 # Unit tests (CMake/CTest based)
│   ├── CMakeLists.txt     # Test build configuration
│   ├── test_spectrum.cpp  # Spectrum class tests
│   ├── test_mathhelper.cpp # Math helper function tests
│   ├── test_analyse.cpp   # Analysis function tests
│   ├── test_baseline.cpp  # Baseline correction tests
│   ├── test_savitzky.cpp  # Savitzky-Golay filter tests
│   ├── test_peakpick.cpp  # File I/O tests
│   ├── test_integration.cpp # Integration tests
│   └── test_advanced.cpp  # Tests for advanced algorithms (NEW)
├── examples/              # Example programs (NEW)
│   ├── CMakeLists.txt     # Examples build configuration
│   ├── example_batch.cpp  # Batch processing example
│   ├── example_comparison.cpp # Spectrum comparison example
│   ├── example_export.cpp # Import/export formats example
│   └── example_fitting.cpp # Multi-peak fitting example
├── src/                   # Test/example source files
│   ├── testLorentzian.h
│   └── testLorentzian.cpp
├── samples/               # Sample spectral data files
│   ├── lorentzian_1
│   └── lorentzian_2
├── eigen/                 # Eigen library (git submodule)
├── main.cpp               # Example usage and demos
├── peakpick_cli.cpp       # Command-line interface tool (NEW)
├── CMakeLists.txt         # CMake build configuration
├── .clang-format          # Code formatting rules
├── PYTHON_BINDINGS.md     # Guide for Python bindings (NEW)
└── README.md              # Basic usage instructions
```

## Development Workflow

### Building the Project

```bash
# Clone with submodules
git clone --recursive https://github.com/conradhuebler/libpeakpick.git

# Or initialize submodules after cloning
git submodule update --init --recursive

# Build the example executable
mkdir build
cd build
cmake ..
make
```

### Using as a Library

Since this is a header-only library, integration is simple:

1. Include the library directory in your CMakeLists.txt:
   ```cmake
   include_directories(${CMAKE_CURRENT_BINARY_DIR} libpeakpick/)
   ```

2. Include the required headers in your C++ code:
   ```cpp
   #include "libpeakpick/peakpick.h"
   #include "libpeakpick/spectrum.h"
   // Include other headers as needed
   ```

### Running Tests

The project includes comprehensive unit tests using CMake's CTest framework:

```bash
# Build with tests enabled (default)
mkdir build
cd build
cmake ..
make

# Run all tests
ctest

# Run tests with verbose output
ctest --verbose

# Run specific test
./tests/test_spectrum

# Run tests in parallel
ctest -j4
```

To disable building tests:
```bash
cmake -DBUILD_TESTS=OFF ..
```

### Test Coverage

The test suite includes:
- **test_spectrum**: Tests for the spectrum class (constructors, accessors, operations)
- **test_mathhelper**: Tests for mathematical functions (mean, stddev, Gaussian, Lorentzian, etc.)
- **test_analyse**: Tests for analysis functions (peak finding, integration, normalization)
- **test_baseline**: Tests for baseline correction structures and algorithms
- **test_savitzky**: Tests for Savitzky-Golay filter coefficients
- **test_peakpick**: Tests for file I/O operations
- **test_integration**: End-to-end integration tests for complete workflows
- **test_advanced**: Tests for advanced algorithms (prominence, AsLS, splines, noise estimation)

## Recent Improvements (2024)

### Critical Bug Fixes

**⚠️ IMPORTANT**: The following critical bugs have been fixed. These are **drop-in replacements** - your existing code should work without changes, but verify behavior with your datasets:

1. **XtoIndex() - spectrum.h** (FIXED)
   - **Problem**: Buggy while-loop and incorrect index modification caused wrong results
   - **Fix**: Complete rewrite with proper bounds checking and binary search fallback
   - **Runtime Notice**: No warnings - fully backward compatible
   - **Action**: **Test your code** that uses `XtoIndex()` or `Y(double x)` to ensure results match expectations

2. **SmoothFunction() - analyse.h** (FIXED)
   - **Problem**: Lost boundary points, started at i=1, stored at i-1
   - **Fix**: Proper boundary handling, now processes all points including edges
   - **Runtime Notice**: Warnings for invalid parameters (size 0, unsupported points)
   - **Action**: **Verify smoothed spectra** - results near boundaries will be different (more accurate)

3. **Normalise() - analyse.h** (FIXED)
   - **Problem**: `min` parameter was ignored, only scaled to max
   - **Fix**: Now properly normalizes to [min, max] range
   - **Runtime Notice**: Warning if spectrum has no range (max == min)
   - **Action**: **Check normalization** - if you relied on old behavior, update `min` parameter

4. **FitBaseLine() - baseline.h** (FIXED)
   - **Problem**: Used Qt type `qreal`, incomplete convergence check
   - **Fix**: Proper convergence criteria and divergence detection
   - **Runtime Notice**: Warnings for optimization failures
   - **Action**: **Monitor baseline fits** - may converge differently

### New Advanced Algorithms (advanced.h)

Include `#include "libpeakpick/advanced.h"` to access these:

1. **PickPeaksAdvanced()** - Improved peak detection with prominence filtering
   ```cpp
   std::vector<Peak> PickPeaksAdvanced(const spectrum* spec,
       double threshold,
       double min_prominence = 0.0,  // NEW: Filter by prominence
       unsigned int min_distance = 1, // NEW: Minimum distance between peaks
       unsigned int start = 0,
       unsigned int end = 0);
   ```
   - Filters noise better than original `PickPeaks()`
   - Runtime message shows number of peaks found with prominence threshold

2. **BaselineAsLS()** - Asymmetric Least Squares baseline correction
   ```cpp
   spectrum BaselineAsLS(const spectrum& spec,
       double lambda = 1e6,    // Smoothness (1e2 to 1e9)
       double p = 0.01,        // Asymmetry (0.001 to 0.1)
       unsigned int max_iter = 10);
   ```
   - Modern baseline correction method
   - Runtime messages show convergence progress
   - Returns baseline-corrected spectrum

3. **EstimateNoise()** - Robust noise estimation using MAD
   ```cpp
   double EstimateNoise(const spectrum& spec);
   ```
   - Uses Median Absolute Deviation for robust noise estimation
   - Useful for automatic threshold selection in peak picking

4. **resampleCubic()** - Cubic spline interpolation
   ```cpp
   spectrum resampleCubic(const spectrum& spec, const Vector& new_x);
   ```
   - Higher quality than linear interpolation (utilities.h)
   - Smoother interpolation for plotting and analysis

5. **CalculatePeakProminence()** - Peak prominence calculation
   ```cpp
   double CalculatePeakProminence(const spectrum* spec, unsigned int peak_idx);
   ```
   - Measures how much a peak stands out from surroundings
   - Used internally by PickPeaksAdvanced()

### Migration Guide for Existing Code

**No changes required** for basic functionality - all fixes are drop-in replacements. However:

1. **Verify Results**: Run your analysis pipelines and compare results
2. **Update Tests**: If you have unit tests, some values may differ slightly
3. **Check Warnings**: Runtime warnings indicate potential issues
4. **Consider New Methods**: `PickPeaksAdvanced()` and `BaselineAsLS()` may improve results

**Example: Migrating to advanced peak detection**
```cpp
// Old code (still works)
std::vector<Peak> peaks = PickPeaks(&spec, threshold);

// New code (better noise rejection)
std::vector<Peak> peaks = PickPeaksAdvanced(&spec, threshold,
    0.5,  // min_prominence - adjust for your data
    5);   // min_distance - minimum 5 points between peaks
```

### Multi-Peak Fitting (fitting.h)

Advanced curve fitting for overlapping peaks with error estimation and multiple peak types.

**Include**: `#include "libpeakpick/fitting.h"`

**Features**:
- Fit multiple peaks simultaneously
- Support for Gaussian, Lorentzian, and Pseudo-Voigt profiles
- Error estimation via covariance matrix
- Chi-squared goodness-of-fit statistics
- Convergence diagnostics

**Usage Example**:
```cpp
// Load spectrum and find peaks
PeakPick::spectrum spec = loadFromFile("data.txt");
std::vector<Peak> peaks = PeakPick::PickPeaksAdvanced(&spec, 0.5, 1.0);

// Fit multiple Gaussian peaks
auto result = PeakPick::FitMultiplePeaks(spec, peaks,
    PeakPick::FitType::Gaussian,  // or Lorentzian, PseudoVoigt
    100);  // max iterations

// Check convergence
if (result.converged) {
    std::cout << "Chi-squared: " << result.chi_squared << std::endl;
    std::cout << "Reduced chi-squared: " << result.reduced_chi_squared << std::endl;

    // Extract parameter errors
    Vector errors = PeakPick::ExtractParameterErrors(result);

    // Access fitted peaks
    for (const auto& peak : result.fitted_peaks) {
        std::cout << "Peak at " << peak.max << " with area " << peak.integ_num << std::endl;
    }
}
```

**API Reference**:
```cpp
enum class FitType {
    Gaussian,      // Gaussian profile
    Lorentzian,    // Lorentzian profile
    PseudoVoigt    // Pseudo-Voigt (mixed Gaussian/Lorentzian)
};

struct FitResult {
    Vector parameters;           // Fitted parameters (height, center, width for each peak)
    Eigen::MatrixXd covariance; // Covariance matrix for error estimation
    Vector residuals;           // Fit residuals
    double chi_squared;         // Chi-squared statistic
    double reduced_chi_squared; // Reduced chi-squared (chi²/dof)
    int iterations;             // Number of iterations performed
    bool converged;             // Convergence status
    std::vector<Peak> fitted_peaks; // Peak structures with fitted parameters
};

FitResult FitMultiplePeaks(const spectrum& spec,
    const std::vector<Peak>& peaks,
    FitType fit_type = FitType::Gaussian,
    int max_iter = 100);

Vector ExtractParameterErrors(const FitResult& result);
```

### Spectrum Comparison (comparison.h)

Quantitative metrics for comparing spectra and peak lists.

**Include**: `#include "libpeakpick/comparison.h"`

**Features**:
- Multiple similarity metrics (RMSD, correlation, cosine similarity)
- Spectral alignment for drift correction
- Peak matching between spectra
- Spectral angle calculation

**Usage Example**:
```cpp
// Load two spectra
PeakPick::spectrum spec1 = loadFromFile("sample1.txt");
PeakPick::spectrum spec2 = loadFromFile("sample2.txt");

// Calculate similarity metrics
double rmsd = PeakPick::RMSD(spec1, spec2);
double corr = PeakPick::CorrelationCoefficient(spec1, spec2);
double cos_sim = PeakPick::CosineSimilarity(spec1, spec2);
double angle = PeakPick::SpectralAngle(spec1, spec2);
double mae = PeakPick::MAE(spec1, spec2);

std::cout << "RMSD: " << rmsd << std::endl;
std::cout << "Correlation: " << corr << std::endl;
std::cout << "Cosine Similarity: " << cos_sim << std::endl;
std::cout << "Spectral Angle: " << angle << "°" << std::endl;

// Align spectra (correct for X-axis drift)
PeakPick::spectrum aligned = PeakPick::AlignSpectra(spec1, spec2,
    10.0,  // search range
    0.1);  // step size

// Compare peak lists
std::vector<Peak> peaks1 = PeakPick::PickPeaksAdvanced(&spec1, 0.5);
std::vector<Peak> peaks2 = PeakPick::PickPeaksAdvanced(&spec2, 0.5);
int matching_peaks = PeakPick::ComparePeaks(peaks1, peaks2, 1.0);  // tolerance
```

**API Reference**:
```cpp
double RMSD(const spectrum& spec1, const spectrum& spec2);
double CorrelationCoefficient(const spectrum& spec1, const spectrum& spec2);
double CosineSimilarity(const spectrum& spec1, const spectrum& spec2);
double MAE(const spectrum& spec1, const spectrum& spec2);
double SpectralAngle(const spectrum& spec1, const spectrum& spec2);
spectrum AlignSpectra(const spectrum& spec1, const spectrum& spec2,
    double search_range = 10.0, double step = 0.1);
int ComparePeaks(const std::vector<Peak>& peaks1,
    const std::vector<Peak>& peaks2, double tolerance = 1.0);
```

### Export Formats (export.h)

Import and export spectra and peaks in standard formats.

**Include**: `#include "libpeakpick/export.h"`

**Features**:
- CSV import/export with customizable delimiters
- JSON import/export with metadata
- Peak list export
- Header and precision control

**Usage Example**:
```cpp
// Export spectrum to CSV
PeakPick::CSVOptions csv_opts;
csv_opts.include_header = true;
csv_opts.delimiter = ',';
csv_opts.precision = 6;
PeakPick::ExportCSV(spec, "output.csv", csv_opts);

// Export to JSON with metadata
PeakPick::ExportJSON(spec, "output.json", true);

// Export peaks to CSV
std::vector<Peak> peaks = PeakPick::PickPeaksAdvanced(&spec, 0.5);
PeakPick::ExportPeaksCSV(peaks, "peaks.csv", csv_opts);
PeakPick::ExportPeaksJSON(peaks, "peaks.json");

// Import from CSV
PeakPick::spectrum imported = PeakPick::ImportCSV("data.csv", true, ',');
```

**API Reference**:
```cpp
struct CSVOptions {
    bool include_header = true;
    char delimiter = ',';
    int precision = 6;
};

bool ExportCSV(const spectrum& spec, const std::string& filename,
    const CSVOptions& options = CSVOptions());
bool ExportPeaksCSV(const std::vector<Peak>& peaks, const std::string& filename,
    const CSVOptions& options = CSVOptions());
bool ExportJSON(const spectrum& spec, const std::string& filename,
    bool include_metadata = true);
bool ExportPeaksJSON(const std::vector<Peak>& peaks, const std::string& filename);
spectrum ImportCSV(const std::string& filename, bool has_header = true,
    char delimiter = ',');
```

### Signal Processing (signal.h)

FFT-based operations, derivatives, and filtering.

**Include**: `#include "libpeakpick/signal.h"`

**Features**:
- First and second derivatives
- Peak detection via second derivative
- Low-pass filtering
- FFT power spectrum
- Convolution

**Usage Example**:
```cpp
// Calculate derivatives
PeakPick::spectrum first_deriv = PeakPick::FirstDerivative(spec);
PeakPick::spectrum second_deriv = PeakPick::SecondDerivative(spec);

// Find peaks using second derivative
std::vector<Peak> peaks = PeakPick::PickPeaksSecondDerivative(spec,
    0.0,    // intensity threshold
    -0.1);  // derivative threshold (negative for maxima)

// Apply low-pass filter
PeakPick::spectrum filtered = PeakPick::LowPassFilter(spec, 0.1);

// Calculate power spectrum (FFT magnitude)
PeakPick::spectrum power = PeakPick::PowerSpectrum(spec);

// Convolve two spectra
PeakPick::spectrum convolved = PeakPick::Convolve(spec1, spec2);
```

**API Reference**:
```cpp
spectrum FirstDerivative(const spectrum& spec);
spectrum SecondDerivative(const spectrum& spec);
std::vector<Peak> PickPeaksSecondDerivative(const spectrum& spec,
    double threshold = 0.0, double deriv_threshold = -0.1);
spectrum LowPassFilter(const spectrum& spec, double cutoff_freq = 0.1);
spectrum PowerSpectrum(const spectrum& spec);
spectrum Convolve(const spectrum& spec1, const spectrum& spec2);
```

### Logging System (logger.h)

Centralized logging with configurable levels and color output.

**Include**: `#include "libpeakpick/logger.h"`

**Features**:
- Multiple log levels (Debug, Info, Warning, Error)
- Color-coded console output
- Context information (function names)
- Global level control
- Convenient macros

**Usage Example**:
```cpp
// Set logging level
PeakPick::Logger::setLevel(PeakPick::LogLevel::Info);

// Enable/disable colors
PeakPick::Logger::setColors(true);

// Log messages
PeakPick::Logger::info("Processing started", "main");
PeakPick::Logger::warning("Low SNR detected", "analyzeSpectrum");
PeakPick::Logger::error("File not found", "loadData");

// Use convenient macros (automatically includes function name)
LOG_INFO("Peak detection complete");
LOG_WARNING("Threshold may be too low");
LOG_ERROR("Invalid parameter");
```

**API Reference**:
```cpp
enum class LogLevel {
    Debug = 0,
    Info = 1,
    Warning = 2,
    Error = 3,
    None = 4
};

class Logger {
    static void setLevel(LogLevel level);
    static LogLevel getLevel();
    static void setColors(bool enable);
    static void log(LogLevel level, const std::string& message,
        const std::string& context = "");
    static void debug(const std::string& message, const std::string& context = "");
    static void info(const std::string& message, const std::string& context = "");
    static void warning(const std::string& message, const std::string& context = "");
    static void error(const std::string& message, const std::string& context = "");
};

// Macros
#define LOG_DEBUG(msg) PeakPick::Logger::debug(msg, __FUNCTION__)
#define LOG_INFO(msg) PeakPick::Logger::info(msg, __FUNCTION__)
#define LOG_WARNING(msg) PeakPick::Logger::warning(msg, __FUNCTION__)
#define LOG_ERROR(msg) PeakPick::Logger::error(msg, __FUNCTION__)
```

### Batch Processing (batch.h)

Fluent interface for processing multiple spectra with consistent parameters.

**Include**: `#include "libpeakpick/batch.h"`

**Features**:
- Fluent/chainable API
- Automatic error handling
- Success/failure tracking
- Custom operation support
- Peak finding across batches

**Usage Example**:
```cpp
// Create batch processor
PeakPick::BatchProcessor processor;

// Add files
processor.addFile("sample1.txt")
         .addFile("sample2.txt")
         .addFile("sample3.txt");

// Or add multiple files
std::vector<std::string> files = {"data1.txt", "data2.txt", "data3.txt"};
processor.addFiles(files);

// Chain operations (fluent interface)
auto results = processor
    .smooth(5)              // Savitzky-Golay smoothing
    .normalize(0.0, 1.0)    // Normalize to [0,1]
    .baselineCorrect(1e6, 0.01)  // AsLS baseline correction
    .process();

// Check results
std::cout << "Processed: " << results.total_processed << std::endl;
std::cout << "Failed: " << results.total_failed << std::endl;

// Process and find peaks
auto results_with_peaks = processor
    .smooth(5)
    .normalize()
    .processAndFindPeaks(0.5,  // threshold
                         1.0,  // prominence
                         5);   // min_distance

// Access results
for (size_t i = 0; i < results_with_peaks.spectra.size(); ++i) {
    if (results_with_peaks.success[i]) {
        std::cout << "File: " << results_with_peaks.filenames[i] << std::endl;
        std::cout << "Peaks found: " << results_with_peaks.all_peaks[i].size() << std::endl;
    }
}

// Apply custom operation
processor.apply([](PeakPick::spectrum& spec) {
    // Custom processing
    spec.center();
    PeakPick::SmoothFunction(&spec, 3);
});
```

**API Reference**:
```cpp
struct BatchResults {
    std::vector<spectrum> spectra;
    std::vector<std::vector<Peak>> all_peaks;
    std::vector<std::string> filenames;
    std::vector<bool> success;
    int total_processed = 0;
    int total_failed = 0;
};

class BatchProcessor {
    BatchProcessor& addFile(const std::string& filename);
    BatchProcessor& addFiles(const std::vector<std::string>& filenames);
    BatchProcessor& smooth(unsigned int points);
    BatchProcessor& normalize(double min = 0.0, double max = 1.0);
    BatchProcessor& baselineCorrect(double lambda = 1e6, double p = 0.01);
    BatchProcessor& apply(std::function<void(spectrum&)> operation);
    BatchResults process();
    BatchResults processAndFindPeaks(double threshold,
        double min_prominence = 0.0, unsigned int min_distance = 1);
};
```

### Command-Line Interface (peakpick_cli)

Standalone CLI tool for common spectral analysis tasks.

**Building**:
```bash
cd build
cmake ..
make
./peakpick_cli help
```

**Commands**:
- `analyze <file>` - Analyze single spectrum with statistics
- `batch <files...>` - Batch process multiple files
- `compare <file1> <file2>` - Compare two spectra
- `peaks <file>` - Find and export peaks
- `help` - Show usage information

**Options**:
- `--smooth <N>` - Apply Savitzky-Golay smoothing (N points)
- `--baseline <lambda>` - Apply AsLS baseline correction
- `--normalize` - Normalize to [0,1]
- `--threshold <val>` - Peak detection threshold
- `--prominence <val>` - Minimum peak prominence
- `--output <file>` - Output filename or directory
- `--verbose` - Verbose output

**Usage Examples**:
```bash
# Analyze spectrum with smoothing and normalization
./peakpick_cli analyze data.txt --smooth 5 --normalize

# Batch process with baseline correction
./peakpick_cli batch *.txt --baseline 1e6 --output results/

# Find peaks with prominence filtering
./peakpick_cli peaks data.txt --threshold 0.5 --prominence 1.0 --output peaks.json

# Compare two spectra
./peakpick_cli compare sample1.txt sample2.txt
```

### Examples Directory

Four complete example programs demonstrating library features.

**Building Examples**:
```bash
cd build
cmake -DBUILD_EXAMPLES=ON ..
make
```

**Available Examples**:

1. **example_batch.cpp** - Batch processing demonstration
   ```bash
   ./examples/example_batch
   ```
   Demonstrates: BatchProcessor, fluent API, error handling

2. **example_fitting.cpp** - Multi-peak fitting demonstration
   ```bash
   ./examples/example_fitting
   ```
   Demonstrates: FitMultiplePeaks, FitType options, error estimation

3. **example_comparison.cpp** - Spectrum comparison
   ```bash
   ./examples/example_comparison
   ```
   Demonstrates: RMSD, correlation, alignment, peak matching

4. **example_export.cpp** - Import/export formats
   ```bash
   ./examples/example_export
   ```
   Demonstrates: CSV/JSON export, import, peak export

### Python Bindings

See `PYTHON_BINDINGS.md` for comprehensive guide on creating Python bindings using pybind11.

**Quick Start**:
```python
# After building Python bindings (see PYTHON_BINDINGS.md)
import libpeakpick as pp
import numpy as np

# Create spectrum
x = np.linspace(0, 100, 1000)
y = np.exp(-((x - 50)**2) / 100) + 0.1 * np.random.randn(1000)
spec = pp.Spectrum(x, y)

# Analyze
print(f"Mean: {spec.mean()}")
print(f"Max: {spec.max()}")

# Find peaks
peaks = pp.pick_peaks_advanced(spec, threshold=0.5, min_prominence=1.0)
print(f"Found {len(peaks)} peaks")

# Export
pp.export_json(spec, "output.json")
```

## Code Conventions

### Code Style

- **Formatter**: clang-format (WebKit-based style)
- **Indentation**: 4 spaces
- **Braces**: Opening brace on same line for functions
- **Alignment**: Don't align after open brackets
- **Functions**: Allow short functions on single line
- **Comments**: Trailing comments not aligned

To format code:
```bash
clang-format -i <file>
```

### Naming Conventions

- **Namespace**: `PeakPick`
- **Classes**: lowercase (e.g., `spectrum`)
- **Member variables**: `m_` prefix (e.g., `m_x`, `m_y`)
- **Methods**: PascalCase (e.g., `Mean()`, `StdDev()`, `PosMax()`)
- **Functions**: camelCase for utilities (e.g., `loadFromFile`)
- **Types**: PascalCase for aliases (e.g., `Vector` for `Eigen::VectorXd`)

### File Organization

- All library headers use `#pragma once` for include guards
- Each header includes necessary Eigen headers
- Headers are self-contained (include their dependencies)
- All implementation is inline (header-only)
- Copyright header in every file with GPL license notice

### Error Handling

- Use exceptions for critical errors (e.g., `throw 2` for size mismatches)
- Print to `std::cout` for non-critical issues (e.g., "Unable to open file")
- No extensive error recovery mechanisms (library expects valid inputs)

## Key Components

### Core Classes and Types

#### `Vector` Type
```cpp
typedef Eigen::VectorXd Vector;
```
Used throughout the library for storing spectral data.

#### `PeakPick::spectrum` Class
The central class for spectral data manipulation:

```cpp
// Construction
spectrum(const Vector& x, const Vector& y);  // From x,y data
spectrum(const Vector& y, double start, double end);  // Auto-generate x-axis
spectrum(const spectrum& other);  // Copy constructor

// Analysis methods
double Mean();      // Calculate mean value
double StdDev();    // Standard deviation
double Max();       // Maximum value
double Min();       // Minimum value
int PosMax();       // Position of maximum
int PosMin();       // Position of minimum
void print();       // Print spectrum data
```

### Important Enumerations

```cpp
namespace PeakPick {
    enum {
        Innovative = 1,
        Liberal = 2,
        Conservative = 3
    };
}
```
These likely control peak picking sensitivity/threshold levels.

### Utility Functions

#### File I/O (peakpick.h)
```cpp
// Load spectrum from file
PeakPick::spectrum loadFromFile(const std::string& filename, double min = 0, double max = 0);
```

#### Enhanced Utilities (utilities.h)
```cpp
// Save spectrum to file
bool saveToFile(const spectrum& spec, const std::string& filename, bool save_x = true);

// Calculate signal-to-noise ratio
double calculateSNR(const spectrum& spec,
    unsigned int peak_region_start, unsigned int peak_region_end,
    unsigned int noise_region_start, unsigned int noise_region_end);

// Find baseline points automatically
std::vector<unsigned int> findBaselinePoints(const spectrum& spec,
    unsigned int num_points, double percentile = 0.1);

// Resample spectrum to new X grid
spectrum resample(const spectrum& spec, const Vector& new_x);

// Calculate full width at half maximum (FWHM)
double calculateFWHM(const spectrum& spec, const Peak& peak);

// Spectrum arithmetic
spectrum subtract(const spectrum& spec1, const spectrum& spec2);
spectrum add(const spectrum& spec1, const spectrum& spec2);
spectrum scale(const spectrum& spec, double factor);
```

## Dependencies

### Eigen Library
- **Version**: Included as git submodule (from eigen-git-mirror)
- **Used for**: Matrix/vector operations, non-linear optimization
- **Location**: `eigen/` directory
- **Key features used**:
  - `Eigen::VectorXd` for data storage
  - Dense and Sparse matrix operations
  - `Eigen::NonLinearOptimization` for curve fitting
  - Template-based numerical algorithms

### OpenMP (Optional)
- Automatically detected by CMake
- Enables parallelization for performance improvements
- Not required but recommended for large datasets

## Compiler Settings

### GCC-Specific Flags (from CMakeLists.txt)
The project uses extensive warning flags for code quality:
- `-Wall -Wextra -pedantic` - Standard warnings
- `-Wcast-align -Wcast-qual` - Casting checks
- `-Wformat -Wformat-security` - Format string security
- `-Winit-self -Winvalid-pch` - Initialization checks
- `-Wunused-*` - Unused code detection
- `-Wno-deprecated-declarations` - Allow deprecated features
- `-Wno-error=enum-compare` - Enum comparison as warnings

### C++ Standard
- **Minimum**: C++11
- Set via: `set_property(TARGET libpeakpick PROPERTY CXX_STANDARD 11)`

## Git Workflow

### Branching Strategy
- **Main branch**: `master` (default)
- **Feature branches**: Use `claude/` prefix for AI-generated changes
- **Current branch**: `claude/claude-md-mi4cqde0bp3i783x-01WmEES7VUX38bStjP8t9fsG`

### Commit Message Style
Based on recent commits, the project uses:
- Short, informal messages
- Lowercase, no periods
- Focus on what changed, not why
- Examples:
  - "some speed up"
  - "fix crash"
  - "change ..."
  - "fix compilation on windows"
  - "update submodule"

### Submodule Management
```bash
# Update Eigen submodule
git submodule update --remote eigen

# After pulling changes
git submodule update --init --recursive
```

## Common Development Tasks

### Adding New Analysis Methods

1. Determine appropriate header file (e.g., `analyse.h` for general analysis)
2. Implement as inline function or method
3. Add to `spectrum` class if it's a core operation
4. Follow existing patterns (use `Vector` type, inline implementation)
5. Add example usage to `main.cpp` or create test in `src/`

### Adding New Fitting Algorithms

1. Look at `glfit.h` and `baseline.h` for patterns
2. Use Eigen's NonLinearOptimization module
3. Create functor struct with `InputType`, `ValueType`, `operator()`
4. Implement fit function that uses Eigen's LevenbergMarquardt solver

### Writing Unit Tests

When adding new features, always add corresponding tests:

1. **Create test file**: Add `test_<feature>.cpp` in `tests/` directory
2. **Follow test structure**:
   ```cpp
   #include "libpeakpick/<header>.h"
   #include <iostream>
   #include <cmath>

   #define TEST_ASSERT(condition, message) \
       if (!(condition)) { \
           std::cerr << "FAILED: " << message << std::endl; \
           return 1; \
       }

   #define TEST_ASSERT_NEAR(val1, val2, epsilon, message) \
       if (std::abs((val1) - (val2)) > (epsilon)) { \
           std::cerr << "FAILED: " << message << std::endl; \
           return 1; \
       }

   int main() {
       // Test code here
       std::cout << "All tests passed!" << std::endl;
       return 0;
   }
   ```

3. **Add to CMakeLists.txt**: Update `tests/CMakeLists.txt` to include new test
4. **Build and run**:
   ```bash
   cd build
   make
   ctest --verbose
   ```

5. **Test naming**: Use descriptive test names and verify edge cases
6. **Use sample data**: Leverage files in `samples/` directory for realistic tests

### Debugging Tips

- Use `-DCMAKE_BUILD_TYPE=Debug` for debug builds
- Eigen provides `.size()`, `.rows()`, `.cols()` for dimension checks
- Use `.print()` method on spectrum objects for quick inspection
- Check for dimension mismatches (common source of `throw 2` errors)

## Platform Compatibility

### Supported Platforms
- **Linux**: Primary development platform
- **Windows**: Supported (including x32 with MinGW)
- **macOS**: Should work (uses standard C++11 and CMake)

### Known Platform Issues
- Windows compilation fixes were needed (see commits)
- Ensure OpenMP is available or disabled appropriately

## External Usage Examples

The library is used in production by:
- [SupraFit](https://github.com/conradhuebler/SupraFit) - Supramolecular titration data fitting
- [QBit](https://github.com/conradhuebler/QBit) - Quantum chemistry analysis tools

Reference these for real-world usage patterns.

## Important Notes for AI Assistants

### When Making Changes

1. **Preserve header-only architecture**: All implementations must be inline
2. **Maintain backward compatibility**: Library is used by other projects
3. **Test thoroughly**: Changes can affect curve fitting accuracy
4. **Update main.cpp**: Add examples for new features
5. **Check Eigen documentation**: Ensure correct usage of Eigen operations
6. **Consider performance**: This library may process large spectral datasets
7. **Follow GPL v3**: All contributions must be GPL-compatible

### Code Quality Guidelines

- Run clang-format before committing
- Ensure all warnings are addressed (extensive warning flags enabled)
- Use const correctness (`const Vector&` for read-only parameters)
- Prefer inline implementations for header-only library
- Use Eigen's optimized operations (avoid manual loops where possible)
- Add comments for complex mathematical operations
- Include mathematical formulas in comments where applicable

### Common Pitfalls to Avoid

1. **Dimension mismatches**: Always check x and y vectors have same size
2. **Missing includes**: Each header should include its Eigen dependencies
3. **Breaking API**: Users include these headers directly, changes affect them
4. **Memory inefficiency**: Use Eigen's references and maps to avoid copies
5. **Thread safety**: If adding state, consider thread safety for OpenMP
6. **Platform-specific code**: Keep implementations cross-platform

### Testing Strategy

- **Unit testing**: Add tests in `src/` directory
- **Integration testing**: Add examples in `main.cpp`
- **Sample data**: Use files in `samples/` for realistic tests
- **Compile testing**: Verify with strict warning flags
- **Cross-platform**: Test on Windows and Linux if possible

## File Format Specifications

### Input Data Format
The `loadFromFile` function expects:
- Plain text files
- One value per line (y-values)
- Lines starting with `#` are comments/metadata
- Optional metadata: `#start = <value>` and `#end = <value>`
- X-values are auto-generated if not provided

## Performance Considerations

- **Eigen optimizations**: Eigen uses expression templates and vectorization
- **OpenMP**: Enable for parallel processing of large datasets
- **Memory layout**: Eigen::VectorXd uses column-major storage
- **Copy avoidance**: Use references (`const Vector&`) to avoid copies
- **Compiler optimizations**: Use `-O2` or `-O3` for release builds

## Mathematical Background

The library implements common spectroscopy operations:
- **Peak picking**: Identifies local maxima in spectral data
- **Baseline correction**: Removes background signals
- **Curve fitting**: Fits Gaussian/Lorentzian functions to peaks
- **Savitzky-Golay**: Polynomial smoothing filter preserving peak shapes
- **Deconvolution**: Separates overlapping spectral features

Understanding these concepts helps when modifying algorithms.

## Contact and Contribution

- **Author**: Conrad Hübler <Conrad.Huebler@gmx.net>
- **Repository**: https://github.com/conradhuebler/libpeakpick
- **License**: GNU GPL v3 - all contributions must maintain this license
- **Development status**: Pre-alpha - expect API changes

When contributing, follow the existing patterns and maintain compatibility with SupraFit and QBit projects.

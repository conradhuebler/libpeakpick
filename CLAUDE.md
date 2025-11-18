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
- Peak picking and detection
- Baseline correction algorithms
- Gaussian and Lorentzian curve fitting
- Savitzky-Golay smoothing
- Linear regression
- Deconvolution
- File I/O for spectral data

## Repository Structure

```
libpeakpick/
├── libpeakpick/           # Header-only library files (main implementation)
│   ├── analyse.h          # Analysis algorithms
│   ├── baseline.h         # Baseline correction methods
│   ├── deconvulate.h      # Deconvolution algorithms
│   ├── glfit.h            # Gaussian/Lorentzian fitting
│   ├── mathhelper.h       # Mathematical utilities
│   ├── nxlinregress.h     # N-dimensional linear regression
│   ├── peakpick.h         # Main header with file loading utilities
│   ├── savitzky.h         # Savitzky-Golay filter implementation
│   └── spectrum.h         # Core spectrum class and operations
├── src/                   # Test/example source files
│   ├── testLorentzian.h
│   └── testLorentzian.cpp
├── samples/               # Sample spectral data files
│   ├── lorentzian_1
│   └── lorentzian_2
├── eigen/                 # Eigen library (git submodule)
├── main.cpp               # Example usage and demos
├── CMakeLists.txt         # CMake build configuration
├── .clang-format          # Code formatting rules
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

```bash
# Build and run the example executable
./libpeakpick
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

```cpp
// Load spectrum from file (in peakpick.h)
PeakPick::spectrum loadFromFile(const std::string& filename, double min = 0, double max = 0);
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

### Testing Changes

1. Add test case to `main.cpp` or create new file in `src/`
2. Use sample data from `samples/` directory
3. Build and run: `./libpeakpick`
4. Verify output matches expectations

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

# Python Bindings for libpeakpick

## Overview

While libpeakpick is a C++ library, you can create Python bindings using **pybind11**.

## Quick Setup

### 1. Install pybind11

```bash
pip install pybind11
```

### 2. Create bindings file (`python/bindings.cpp`)

```cpp
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include "libpeakpick/peakpick.h"
#include "libpeakpick/advanced.h"
#include "libpeakpick/fitting.h"
#include "libpeakpick/comparison.h"
#include "libpeakpick/export.h"

namespace py = pybind11;

PYBIND11_MODULE(libpeakpick_py, m) {
    m.doc() = "Python bindings for libpeakpick";

    // Spectrum class
    py::class_<PeakPick::spectrum>(m, "Spectrum")
        .def(py::init<const Vector&, const Vector&>())
        .def("Mean", &PeakPick::spectrum::Mean)
        .def("Max", &PeakPick::spectrum::Max)
        .def("Min", &PeakPick::spectrum::Min)
        .def("StdDev", &PeakPick::spectrum::StdDev)
        .def("size", &PeakPick::spectrum::size)
        .def("X", &PeakPick::spectrum::X)
        .def("Y", py::overload_cast<unsigned int>(&PeakPick::spectrum::Y, py::const_));

    // Peak structure
    py::class_<PeakPick::Peak>(m, "Peak")
        .def(py::init<>())
        .def_readwrite("start", &PeakPick::Peak::start)
        .def_readwrite("max", &PeakPick::Peak::max)
        .def_readwrite("end", &PeakPick::Peak::end)
        .def_readwrite("deconv_x", &PeakPick::Peak::deconv_x)
        .def_readwrite("deconv_y", &PeakPick::Peak::deconv_y);

    // Functions
    m.def("loadFromFile", &loadFromFile, "Load spectrum from file");
    m.def("PickPeaksAdvanced", &PeakPick::PickPeaksAdvanced,
        py::arg("spec"), py::arg("threshold"),
        py::arg("min_prominence") = 0.0,
        py::arg("min_distance") = 1,
        py::arg("start") = 0,
        py::arg("end") = 0);
    m.def("EstimateNoise", &PeakPick::EstimateNoise);
    m.def("BaselineAsLS", &PeakPick::BaselineAsLS,
        py::arg("spec"), py::arg("lambda") = 1e6,
        py::arg("p") = 0.01, py::arg("max_iter") = 10);
    m.def("RMSD", &PeakPick::RMSD);
    m.def("CorrelationCoefficient", &PeakPick::CorrelationCoefficient);
    m.def("ExportCSV", &PeakPick::ExportCSV);
    m.def("ExportJSON", &PeakPick::ExportJSON);
}
```

### 3. Create `setup.py`

```python
from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup

ext_modules = [
    Pybind11Extension(
        "libpeakpick_py",
        ["python/bindings.cpp"],
        include_dirs=[".", "eigen/"],
        cxx_std=11,
    ),
]

setup(
    name="libpeakpick",
    version="1.0.0",
    author="Conrad Hübler",
    description="Python bindings for libpeakpick",
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
)
```

### 4. Build and Install

```bash
pip install .
```

## Usage Example

```python
import libpeakpick_py as lpp
import numpy as np

# Load spectrum
spec = lpp.loadFromFile("data.txt")

# Get statistics
print(f"Mean: {spec.Mean()}")
print(f"Max: {spec.Max()}")
print(f"Noise: {lpp.EstimateNoise(spec)}")

# Baseline correction
spec_corrected = lpp.BaselineAsLS(spec, lambda_=1e6, p=0.01)

# Find peaks
peaks = lpp.PickPeaksAdvanced(spec_corrected,
                               threshold=0.5,
                               min_prominence=1.0,
                               min_distance=5)

print(f"Found {len(peaks)} peaks")

for i, peak in enumerate(peaks):
    print(f"Peak {i+1}: position={peak.deconv_x}, height={peak.deconv_y}")

# Export
lpp.ExportJSON(spec_corrected, "output.json", True)
```

## Integration with NumPy

```python
# Convert to NumPy arrays
x = np.array([spec.X(i) for i in range(spec.size())])
y = np.array([spec.Y(i) for i in range(spec.size())])

# Plot with matplotlib
import matplotlib.pyplot as plt
plt.plot(x, y)
plt.show()
```

## Integration with Pandas

```python
import pandas as pd

# Create DataFrame from spectrum
df = pd.DataFrame({
    'X': [spec.X(i) for i in range(spec.size())],
    'Y': [spec.Y(i) for i in range(spec.size())]
})

# Export to CSV
df.to_csv("spectrum.csv", index=False)
```

## Further Resources

- [pybind11 documentation](https://pybind11.readthedocs.io/)
- [Eigen-pybind11 integration](https://pybind11.readthedocs.io/en/stable/advanced/cast/eigen.html)

## Notes

- Requires pybind11 >= 2.6.0
- Eigen must be in include path
- C++11 compiler required
- For better NumPy integration, use `pybind11/eigen.h`

# libpeakpick
Small header-only library to performe some peak picking and spectra analysis

## Download and requirements
git clones automatically Eigen. Eigen is used as non-linear optimimization tool.

Include directory like this in your CMakeLists.txt
```
include_directories(${CMAKE_CURRENT_BINARY_DIR} libpeakpick/)
```

and use the header files directly in your project. Example usage can be found in
- [SupraFit](https://github.com/conradhuebler/SupraFit)
- [QBit](https://github.com/conradhuebler/QBit)


## Note
This is pure pre-pre alpha.

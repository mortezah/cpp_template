# a template starting project for cpp

* this template uses `CMake` and `CMakePresets` to set up a cross platform builds
* dependencies are handle via `vcpkg`s in manifest mode (i.e. the packages listed in `vcpkg.json` will be automatically downloaded, installed, and made available to the project during the configuration stage)
* this includes a testing template that used a combination of `gtest` and `ctest`


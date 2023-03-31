# palabos

Based on release v2.0r0.

The modifications implemented in this fork allow to gather lists of values in a vector using MPI.
This feature can be used to evaluate the strain rate at the solid/fluid interface to determine where to apply erosion or sedimentation for reshaping with [SEROS](https://github.com/AIT-LKR/SEROS).
More information on the palabos library is available at [https://palabos.unige.ch/](https://palabos.unige.ch/).

## Linked software
The test framework [Catch2](https://github.com/catchorg/Catch2) was used during development of this fork. This is no crucial part to use palabos.
For running [SEROS](https://github.com/AIT-LKR/SEROS), the tool [IsoSurfaceExtraction](https://github.com/AIT-LKR/IsoSurfaceExtraction) was created as fork of [https://github.com/mkazhdan/IsoSurfaceExtraction](https://github.com/mkazhdan/IsoSurfaceExtraction) and linked with this version of palabos as external library. This is also no crucial part to use palabos.

## Acknowledgement
The work presented here was conducted in the context of the Clean Sky 2 project COMBO3D. This project has received funding from the Clean Sky 2 Joint Undertaking under the European Union’s Horizon 2020 research and innovation programme under grant agree- ment No 831989. This publication reflects only the author’s views and the European Union is not liable for any use that may be made of the information contained therein.

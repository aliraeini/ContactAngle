# Changelog

This repository started as a fork of the
[original code by AlRatrout et al](https://github.com/AhmedAlratrout/ContactAngle-Curvature-Roughness)
and has been **heavily refactored** since (build system, image I/O, scripts and C++ sources).
Regression tracking started only recently (October 2026): earlier changes were not checked against
reference results, so results of this fork may differ from the original code.
The reference results are in [ref/](ref/) (`tutorial_subvolume_stats.txt`, made by `ref/post_process.py`
from the tutorial) and in `test.py` (`dropInCylinder`); they are the results of the code before the refactoring
and are checked after each change.

## Unreleased

### Changed: unit of the curvature `Kc`

`Kc` is now in **1/voxel**, i.e. the voxel size in the image header (`ElementSize`) no longer
affects the output: `AllRunContAngle` sets the voxel size to 1 before the surface is extracted.

Before, `AllRunContAngle` multiplied the voxel size of the header by 1e6 (`reset s 1e6`), so the curvature was in
1/(1e6 x header length unit), e.g. 1/micrometre for headers in metres. To get physical units from the new output:

    Kc [1/m] = Kc [1/voxel] / voxel size [m]


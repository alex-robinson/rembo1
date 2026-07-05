# rembo1

## Dependencies (fesm-utils compatibility shift)

rembo1 now relies on **fesm-utils** for all its shared modules — both the utility
modules (`ncio`, `nml`, `interp1D`, `gaussian_filter`, `ncio_transpose`) and the
coordinate/grid modules (formerly the standalone `coordinates` library, now folded
into `libfesmutils`). It no longer vendors its own copies, nor links `coordinates`
separately. This requires **fesm-utils dev at `4a8511e` (2026-07-05) or later**, which
folds the coords module into `libfesmutils` and hoists `utils/` to the repo root — the
layout this build's paths assume (`fesm-utils/include-serial`, `fesm-utils/lis/lis-serial`).
Set up a `fesm-utils` symlink/checkout as for the other fesm models (the build links
`libfesmutils`).

## Build

make clean
make rembo-static
make rembo



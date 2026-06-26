# rembo1

## Dependencies (fesm-utils compatibility shift)

rembo1 now relies on **fesm-utils** for its shared utility modules (`ncio`, `nml`,
`interp1D`, `gaussian_filter`, `ncio_transpose`) — it no longer vendors its own copies.
This requires **fesm-utils dev at `3f415cc` (2026-06-26) or later**, which adds
`ncio_transpose` and a double-precision-capable `interp1D`. Set up a `fesm-utils`
symlink/checkout as for the other fesm models (the build links `libfesmutils`, listed
ahead of the `coordinates` library so its modules take precedence).

## Build

make clean
make rembo-static
make rembo



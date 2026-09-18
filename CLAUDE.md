# fmmlib3d_em_matlab

Matlab/Octave toolbox for 3D electromagnetic (Maxwell) dipole interactions: direct and FMM-accelerated field evaluation plus FMM-accelerated boundary integral equation demos (AEFIE, MFIE, Muller) on triangulated surfaces.

## Conventions that must not drift

- Green's function convention: the dyadic kernels use `exp(ikr)/r` with no `1/(4 pi)` scaling. Every routine header repeats this; keep new code consistent and do not "fix" it toward the 4 pi-normalized form.
- FMM precision flag `iprec` (values -2 .. 5) maps to tolerances 0.5, 0.5e-1, 0.5e-2, 0.5e-3, 0.5e-6, 0.5e-9, 0.5e-12, 0.5e-15 (see the header of `emfmm3dpart_matlab.m`).

## Dependencies

- fmmlib3d (the scalar Helmholtz FMM library) is external; `startup.m` adds its `matlab/` directory to the path (`addpath('../fmmlib3d/matlab')`). The FMM-accelerated routines (`emfmm3dpart_matlab.m`, `emfmm3dtria_*.m`) call its `hfmm3dpart`.
- Building the MEX file needs mwrap and gfortran (> 4.4.0 strongly recommended), plus Matlab or Octave.

## Wrapper duality

Direct evaluation exists in two flavors that must stay numerically equivalent:

- MEX-backed: `em3dpartdirect.m`, `em3dccpartdirect.m` call the compiled Fortran in `emtools.mex` (generated from `emtools.mw` by mwrap; sources in `src/`).
- Pure Matlab: `em3dpartdirect_matlab.m`, `emfmm3dpart_matlab.m` need no MEX build and use `em3dipole3et.m` / `em3dipole3mt.m` for the dipole fields.

When changing a kernel, change both and rerun `test_emfmm3dpart` / `test_emfmm3dparttarg` in Matlab or Octave.

## Building

`./makefile` has `linux`, `macosx` (Intel and arm64 Matlab), `windows` (64-bit cross-build only) and `clean` targets; each delegates to `makefile.mwrap` with a `TARGET=` from `mwrap.inc`, which is the authoritative list of platform configurations.

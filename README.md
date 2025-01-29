# NSInterpolate

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://bamburgh.github.io/NSInterpolate.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://bamburgh.github.io/NSInterpolate.jl/dev/)
[![Build Status](https://github.com/bamburgh/NSInterpolate.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/bamburgh/NSInterpolate.jl/actions/workflows/CI.yml?query=branch%3Amain)

 NSInterpolate is a geophysical gridding program based on NSInterp (Naprstek & Smith, 2019).

 This is a Julia version based on the original C# version 1.8 and available on github at:

     https://github.com/TomasNaprstek/Naprstek-Smith-Interpolation

and published in:

	T. Naprstek & R. S. Smith. (2019). A new method for interpolating linear
	features in aeromagnetic data. Geophysics, 84(3), JM15-JM24.
	https://library.seg.org/doi/10.1190/geo2018-0156.1

Thank you to Tomas Naprstek and Richard Smith for all their excellent work!

## Changes

The only changes, other than those resulting from language translation, are at input
and output.

Currently NSInterpolate.jl accepts a number of parameters to control the
code. These may be passed directly to NSinterp, or may be stored in a JSON file
in which case, the name of the JSON file should be passed to NSinterp.

The output grid is written to a netCDF4 file.

## Warnings

This is the first release of the Julia code and will likely contain bugs. In particular,
the code makes some complex use of array indexing and the change from C#, with 0-indexing,
to Julia, with 1-indexing, may have introduced some errors.

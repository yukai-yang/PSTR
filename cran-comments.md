## Summary

This is an update to the CRAN package "PSTR".

Version 2.1.1 improves numerical robustness, refines the internal implementation, fixes multi-switch transition-function evaluation, and modernises several utilities while keeping the statistical methodology and the user-facing workflow unchanged.

## Main changes

* Fixed transition-function evaluation for models with multiple switches (`im > 1`), which could previously lead to recycling warnings and dimension mismatch errors in `EstPSTR()`.
* Improved numerical robustness in nonlinear least squares estimation and variance computation.
* Updated bootstrap routines and internal parallel helpers.
* Improved printing and diagnostic plotting utilities.
* Updated documentation, examples, and vignette.

Existing user code following documented usage should continue to work as before.

## Reverse dependencies

There are currently no reverse dependencies on CRAN.

## R CMD check results

Platform: aarch64-apple-darwin20  
R version: 4.5.2  

0 errors | 0 warnings | 0 notes
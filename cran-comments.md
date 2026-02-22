This is a major update of the CRAN package "PSTR".

Version 2.0.0 improves the internal structure and numerical robustness of the package.
The core implementation has been reorganised using an R6 class structure,
while exported user-level functions remain available for standard workflows.

## Main changes in version 2.0.0

* Internal refactoring to an R6-based architecture.
* Clearer separation of specification, estimation, and evaluation routines.
* Improved numerical stability in matrix inversion (SVD-based pseudoinverse when necessary).
* Rewritten bootstrap routines to improve internal consistency.
* Modernised console output and printing.
* Standardised internal storage of evaluation results.

The statistical methodology and default user interface remain unchanged.
Existing user code following documented usage should continue to work as before.

## Reverse dependencies

There are currently no reverse dependencies on CRAN.

## R CMD check results

Platform: aarch64-apple-darwin20  
R version: 4.5.2  

0 ERRORs | 0 WARNINGs | 0 NOTEs.

## Additional notes

The corresponding methodological reference is a working paper
without a DOI. The full reference and link are provided in the package
documentation and can be accessed via `?PSTR`.
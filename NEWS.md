<!-- README.md is generated from README.Rmd. Please edit that file -->

# PSTR 2.0.0 ‘Green Panel’

## New features

- The `PSTR` class now follows an R6 design with reference semantics.
- The main modelling steps are implemented as class methods.
- Exported wrapper functions remain available for backward
  compatibility.
- The command line output is redesigned using the `cli` package.
- The print layout is modernised using `knitr::kable`.
- Separate print views are available for summary, tests, estimation, and
  evaluation.
- Internal getters and setters are added to support a clean R6 workflow.

## Improvements

- Evaluation and bootstrap routines are rewritten under the R6
  structure.
- Result objects are stored with consistent naming (`tv`, `ht`,
  `wcb_tv`, `wcb_ht`).
- Numerical robustness is improved in near-singular matrix inversion
  cases.
- An SVD-based pseudoinverse is used when necessary.
- Bootstrap routines use deep cloning to avoid overwriting the original
  object.
- Internal print helpers are reorganised into dedicated private methods.

## Notes

- This release includes a substantial internal refactor.
- Direct access to model components is discouraged.
- Users should rely on documented print modes and exported interfaces.

# New Features in PSTR 1.3.0 ‘Yellow Panel’

- fix the errors in R version 4.
- change documentation.

# New Features in PSTR 1.2.5 ‘Orange Panel’

- fix the bug in the function “plot_transition”.
- fix the bug in the function “plot_coefficients”.

# New Features in PSTR 1.2.4 ‘Orange Panel’

- improve the plotting function “plot_transition” and “plot_response”.
- add new plotting function “plot_coefficients”.
- change documentation.

# New Features in PSTR 1.2.3 ‘Orange Panel’

- fix some bugs found in the previous version.

# New Features in PSTR 1.2.2 ‘Orange Panel’

- fix some bugs found in the previous version.

# New Features in PSTR 1.2.1 ‘Orange Panel’

- Link to the working paper.

# New Features in PSTR 1.2.0 ‘Orange Panel’

- New function “plot_response” which use the ggplot2 and plotly packages
  to plot the curve or surface of the expected response (contribution)
  against some variable and the transition variable.
- The estimation function “EstPSTR” will return the estimated
  time-invariante individule effects.
- The width of the separating lines in printing the results can adjust
  itself automatically to fit the width of the console.
- Rewrite the “version” function.

# New Features in PSTR 1.1.0 ‘Red Panel’

- Better documentation.
- In the function “EstPSTR”, now the initial value for “gamma” is
  allowed. The user can freely choose either “gamma” or “delta” to input
  for the estimation.
- Estimation of a linear panel regression model with fixed effects.
- The initial values can now be chosen automatically and will be
  returned after estimation.
- New function “plot_transition” which use the ggplot2 package to plot
  the estimated transition function.
- The package ggplot2 imported.
- A new data set “sunspot” included.
- Some bugs fixed.

# New Featurs in PSTR 1.0.1

- All the functions
- Documetation including README, LICENSE and etc.

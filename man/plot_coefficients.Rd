% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{plot_coefficients}
\alias{plot_coefficients}
\title{Plot the coefficients, the standard errors and the p-values against the transition variable.}
\usage{
plot_coefficients(obj, vars, length.out = 100, color = "blue", size = 1.5)
}
\arguments{
\item{obj}{an object of the class PSTR returned from some functions in the package. Note that the corresponding PSTR model must be estimated first.}

\item{vars}{a vector of column numbers or names (character strings) specifying which variables in the nonlinear part to use.}

\item{length.out}{a scalar of desired length (number of points) for building the x-axis. 100 by default.}

\item{color}{the color of the lines.}

\item{size}{the size of the lines.}
}
\value{
A list of plottable objects from the \code{ggplot2} package.
}
\description{
This function plots the curves of the coefficients, the standard errors and the p-values against the transition variable.
}
\details{
The curves of the coefficients, the standard errors and the p-values against the transition variable are functions
\deqn{f_1(x) = \beta_{0j} + \beta_{1j}g(x ; \gamma, c)}
\deqn{f_2(x) = se(f_1(x))}
\deqn{f_3(x) = 1 - Prob\{ X < [f_1(x)/f_2(x)]^2 \} }
where \eqn{x} is a variable taking the position of the transition variable,
\eqn{se} stands for the cluster-robust and heteroskedasticity-consistent standard error of the estimate \eqn{f_1(x)} at \eqn{x},
\eqn{X} is a random variable following chi-square distribution with degrees of freedom one.

More than one variable can be put in \code{vars}.

The return value is a list of the same length as \code{vars}, whose elements are plottable objects.
}
\examples{
\donttest{
pstr = NewPSTR(Hansen99, dep='inva', indep=4:20, indep_k=c('vala','debta','cfa','sales'),
    tvars=c('vala','debta','cfa','sales'), iT=14) # create a new PSTR object

# estimate the PSTR model first
pstr = EstPSTR(use=pstr, im=1, iq=1, useDelta=TRUE, par=c(.63,0), method='CG')

# plot the curve and surfaces
ret = plot_coefficients(pstr, vars=1:4, length.out=100, color="dodgerblue4", size=2)
ret[[1]]
ret[[1]] + ggplot2::scale_x_log10()
}

}
\seealso{
Functions which return an object of the class PSTR can be input into this function

\code{\link{EstPSTR}}
}
\author{
Yukai Yang, \email{yukai.yang@statistik.uu.se}
}
\keyword{utils}

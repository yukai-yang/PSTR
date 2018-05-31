#' A balanced panel of 565 US firms observed for the years 1973–1987
#'
#' A dataset containing a balanced panel data of annual observations
#' over the period 1973-1987 (15 years) for 560 US firms for the variables described below.
#'
#' The structure of the dataset is such that the time index runs “fast”, while the firm index runs “slow”;
#' that is, first all 14 observations for the first firm are given, then the 14 observations for the second firm, etc.
#'
#' Since we used one year lagged variables of "vala", "debta", "cfa" and "cfa" as regressors,
#' the records in 1973 are skipped.
#'
#' All values are nominal and millions of dollars except where otherwise noted. Stocks are end of year.
#'
#' @format A tibble with 7840 rows and 20 variables:
#' \describe{
#'   \item{cusip}{Committee on Uniform Security Identication Procedures firm code number, the first 6 digits (CNUM)}
#'   \item{year}{2-digit year of the data}
#'   \item{inva}{investment to assets ratio}
#'   \item{dt_75}{dummy variable for 1975}
#'   \item{dt_76}{dummy variable for 1976}
#'   \item{dt_77}{dummy variable for 1977}
#'   \item{dt_78}{dummy variable for 1978}
#'   \item{dt_79}{dummy variable for 1979}
#'   \item{dt_80}{dummy variable for 1980}
#'   \item{dt_81}{dummy variable for 1981}
#'   \item{dt_82}{dummy variable for 1982}
#'   \item{dt_83}{dummy variable for 1983}
#'   \item{dt_84}{dummy variable for 1984}
#'   \item{dt_85}{dummy variable for 1985}
#'   \item{dt_86}{dummy variable for 1986}
#'   \item{dt_87}{dummy variable for 1987}
#'   \item{vala}{lagged total market value to assets ratio ("Tobin's Q")}
#'   \item{debta}{lagged long term debt to assets ratio}
#'   \item{cfa}{lagged cash flow to assets ratio}
#'   \item{sales}{lagged sales during the year (million USD)}
#' }
#' @source \url{http://www.ssc.wisc.edu/~bhansen/progs/joe_99.html}
"Hansen99"


#' Transformed Wolf annual sunspot numbers for the years 1710-1979
#'
#' A dataset containing the transformed Wolf annual sunspot numbers for the years 1710-1979.
#'
#' Each column of the data matrix is a lagged transformed sunspot observations from lag order 0 to 10.
#'
#' The data were transformed by using the formula
#' \deqn{y_t = 2 \left\{ (1 + x_t)^{1/2} -1 \right\}}
#' see Ghaddar and Tong (1981)
#'
#' @format A tibble with 270 rows and 11 variables:
#' \describe{
#'   \item{spot_0}{transformed sunspot}
#'   \item{spot_1}{transformed sunspot, lag one}
#'   \item{spot_2}{transformed sunspot, lag two}
#'   \item{spot_3}{transformed sunspot, lag three}
#'   \item{spot_4}{transformed sunspot, lag four}
#'   \item{spot_5}{transformed sunspot, lag five}
#'   \item{spot_6}{transformed sunspot, lag six}
#'   \item{spot_7}{transformed sunspot, lag seven}
#'   \item{spot_8}{transformed sunspot, lag eight}
#'   \item{spot_9}{transformed sunspot, lag nine}
#'   \item{spot_10}{transformed sunspot, lag ten}
#' }
#' 
#' @section References:
#' Ghaddar, D. K. and Tong, H. (1981) Data transformation and self‐exciting threshold autoregression, Applied Statistics, 30, 238–48.
#' 
#' @source \url{http://sidc.oma.be/html/sunspot.html}
"sunspot"

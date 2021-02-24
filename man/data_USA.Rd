\name{data_USA}
\alias{data_USA}
\docType{data}
\title{
Historical data for the United States
}
\description{
data_USA contains:
- NBER based Recession Indicators for the United States from 1953:04 to 2020:01
- 10Years TB for the United States from 1953:04 to 2020:01
- 3Months TB for the United States from 1953:04 to 2020:01
- Yield Spread (10Years TB minus 3Months TB) for the United States from 1975:03 to 2019:05
}
\usage{data("data_USA")}

\format{
  A data frame with 268 observations on the following 5 variables.
  \describe{
    \item{\code{Date}}{Vector of dates.}
    \item{\code{X10Y}}{Historical 10 years Treasury bond.}
    \item{\code{X3M}}{Historical 3 months Treasury bond.}
    \item{\code{Spread}}{Historical yield spread.}
    \item{\code{NBER}}{Historical binary variable related to historical recessions.}
  }
}

\source{
https://fred.stlouisfed.org/
}
\examples{
data("data_USA")
head("data_USA")
}
\keyword{datasets}

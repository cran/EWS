\name{data_USA}
\alias{data_USA}
\docType{data}
\title{
Historical data for the United States
}
\description{
data_USA contains:
- NBER based Recession Indicators for the United States from 1975:03 to 2019:05
- Yield Spread (10Years TB minus 3Months TB) for the United States from 1975:03 to 2019:05
}
\usage{data("data_USA")}

\format{
  A data frame with 531 observations on the following 4 variables.
  \describe{
    \item{\code{country}}{USA.}
    \item{\code{Date}}{Vector of dates.}
    \item{\code{YIESPR}}{Historical yield spread.}
    \item{\code{NBER}}{Historical binary variable related to historical recessions.}
  }
}

\source{
https://fred.stlouisfed.org/
}

\keyword{datasets}

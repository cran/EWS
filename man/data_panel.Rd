\name{data_panel}
\alias{data_panel}
\docType{data}
\title{
Historical data for 13 OECD countries
}
\description{
data_USA contains:
- OECD based Recession Indicators for 13 OECD countries from the Peak through the Trough from 1975:03 to 2019:05
- Yield Spread (10Years TB minus 3Months TB) for 13 OECD countries from 1975:03 to 2019:05

List of countries: Australia, Belgium, Canada, France, Germany, Italy, Japan, the Netherlands, New Zealand, Sweden, Switzerland, the United Kinkdom, the United States.
}
\usage{data("data_panel")}
\format{
  A data frame with 6903 observations on the following 4 variables.
  \describe{
    \item{\code{country}}{List of countries.}
    \item{\code{Date}}{Vector of dates.}
    \item{\code{YIESPR}}{historical yield spread for the 13 OECD countries.}
    \item{\code{OECD}}{Historical binary variable related to historical recessions for the 13 OECD countries.}
  }
}

\source{
https://fred.stlouisfed.org/
}
\examples{
data("data_panel")
head("data_panel")
}
\keyword{datasets}

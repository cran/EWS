\name{Vector_Error}
\alias{Vector_Error}
\title{
  Vector of Errors
}
\description{
The function measures the estimation errors from the logistic estimation, and stores them in a vector. This function is used to initialize a shock in impulse response analysis as in Koop, Pesaran and Potter (1996).
}
\usage{
Vector_Error(Dicho_Y, Exp_X, Intercept, Nb_Id, Lag, type_model)
}
\arguments{
  \item{Dicho_Y}{
    Vector of the binary time series.
  }
  \item{Exp_X}{
    Vector or Matrix of explanatory time series.
  }
  \item{Intercept}{
  Boolean value: TRUE for an estimation with intercept, and FALSE otherwise.
  }
  \item{Nb_Id}{
    Number of individuals studied for a panel approach. Nb_Id=1 in the univariate case.
  }
  \item{Lag}{
    Number of lags used for the estimation.
  }
  \item{type_model}{
    Model number: 1, 2, 3 or 4.
  }
}

\value{
  A numeric vector containing estimation errors.
}

\references{

Kauppi, Heikki, and Pentti Saikkonen. "Predicting US recessions with dynamic binary response models." The Review of Economics and Statistics 90.4 (2008): 777-791.

Koop, Gary, M. Hashem Pesaran, and Simon M. Potter. "Impulse response analysis in nonlinear multivariate models." Journal of econometrics 74.1 (1996): 119-147.

Lajaunie, Quentin. Generalized Impulse Response Function for Dichotomous Models. No. 2852. Orleans Economics Laboratory/Laboratoire d'Economie d'Orleans (LEO), University of Orleans, 2021.

}
\author{
Jean-Baptiste Hasse and Quentin Lajaunie
}

\examples{

# NOT RUN {

# Import data
data("data_USA")

# Data process
Var_Y <- as.vector(data_USA$NBER)
Var_X <- as.vector(data_USA$Spread)

# Estimate the estimation errors
results <- Vector_Error(Dicho_Y = Var_Y, Exp_X = Var_X, Intercept = TRUE,
                      Nb_Id = 1, Lag = 1, type_model = 4)

# print results
results

#}
}


\keyword{Econometrics}
\keyword{Dichotomous}
\keyword{Shock}

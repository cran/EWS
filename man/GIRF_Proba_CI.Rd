\name{GIRF_Proba_CI}
\alias{GIRF_Proba_CI}
\title{
Confidence Intervals for the Probability - GIRF Analysis
}
\description{
From the results of the Simulation_GIRF function, this function calculates the values of the upper and lower bounds of the confidence intervals, as well as the average of the different response functions for the probability.
}
\usage{
GIRF_Proba_CI(results_simul_GIRF, CI_bounds, n_simul, horizon_forecast)
}
\arguments{
  \item{results_simul_GIRF}{
    Matrix containing results of the Simulation_GIRF function.
  }
  \item{CI_bounds}{
    Numeric variable between 0 and 1 for the size of the confidence intervals.
  }
  \item{n_simul}{
    Numeric variable equal to the total number of replications.
  }
  \item{horizon_forecast}{
    Numeric variable corresponding to the horizon target for the GIRF analysis.
  }
}



\value{
A list with:

  \item{horizon}{Numeric vector containing the horizon targer.}
  \item{Simulation_CI_proba_shock}{a matrix containing the set of simulations of probabilities with shock belonging to the confidence interval.}
  \item{Simulation_CI_proba}{a matrix containing the set of simulations of probabilities belonging to the confidence interval.}
  \item{CI_proba_shock}{a matrix containing three columns: the lower bound of the CI, the average of the IRFs, and the upper bound of the CI for the probabilities with shock.}
  \item{CI_proba}{a matrix containing three columns: the lower bound of the CI, the average of the IRFs, and the upper bound of the CI for the probabilities.}
}

\references{
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

# Simulation for the GIRF analysis
results_simulation <- Simul_GIRF(Var_Y, Var_X, TRUE, 1, 1, 2, 0.95, 3, "AM")

# Confidence intervals for the index
results <- GIRF_Proba_CI(results_simulation, 0.95, 2, 3)

# print results
results

#}
}

\keyword{Econometrics}
\keyword{Dichotomous}
\keyword{IRF}
\keyword{Confidence-Intervals}

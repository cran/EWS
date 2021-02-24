\name{Simul_GIRF}
\alias{Simul_GIRF}
\title{
GIRF Simulations
}
\description{
This function calls the BlockBootstrap function of the EWS package and then calculates response functions for each simulation. It then measures the confidence intervals as in Lajaunie (2021). The response functions are based on the 4 specifications proposed by Kauppi & Saikkonen (2008).
}
\usage{
Simul_GIRF(Dicho_Y, Exp_X, Int, Lag, t_mod, n_simul, centile_shock, horizon, OC)
}
\arguments{
  \item{Dicho_Y}{
Vector of the binary time series.
}
  \item{Exp_X}{
Vector or Matrix of explanatory time series.
}
  \item{Int}{
Boolean value: TRUE for an estimation with intercept, and FALSE otherwise.
}
  \item{Lag}{
Number of lags used for the estimation.
}
  \item{t_mod}{
Model number: 1, 2, 3 or 4.

-> 1 for the static model:

\deqn{P_{t-1}(Y_{t}) = F(\pi_{t})=F(\alpha + \beta'X_{t})}


-> 2 for the dynamic model with lag binary variable:

\deqn{P_{t-1}(Y_{t}) = F(\pi_{t})=F(\alpha + \beta'X_{t} + \gamma Y_{t-l})}


-> 3 for the dynamic model with lag index variable:

\deqn{P_{t-1}(Y_{t}) = F(\pi_{t})=F(\alpha + \beta'X_{t} + \eta \pi_{t-l})}


-> 4 for the dynamic model with both lag binary variable and lag index variable:

\deqn{P_{t-1}(Y_{t}) = F(\pi_{t})=F(\alpha + \beta'X_{t} + \eta \pi_{t-l} + \gamma Y_{t-l})}
}
  \item{n_simul}{
    Numeric variable equal to the total number of replications.
  }
  \item{centile_shock}{
    Numeric variable corresponding to the centile of the shock following Koop, Pesaran and Potter (1996).
  }
  \item{horizon}{
    Numeric variable corresponding to the horizon target for the GIRF analysis.
  }
  \item{OC}{
    Either a numeric variable equal to the optimal cut-off (threshold) or a character variable of the method chosen to calculate the optimal cut-off ("NSR", "CSA", "AM").
  }
}

\value{
A matrix containing the GIRF analysis for each replication. For each replication, the function returns 7 colomns with:

  \item{column 1}{horizon}
  \item{column 2}{index}
  \item{column 3}{index with shock}
  \item{column 4}{probability associated to the index}
  \item{column 5}{probability associated to the index with shock}
  \item{column 6}{binary variable associated to the index}
  \item{column 7}{binary variable associated to the index with shock}

The matrix contains \eqn{7 \times S} colomns, where \eqn{S} denotes the number of replications.
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

# Simulations
results <- Simul_GIRF(Dicho_Y = Var_Y, Exp_X = Var_X, Int = TRUE, Lag = 1, t_mod = 1 ,
                      n_simul = 2 , centile_shock = 0.95, horizon = 3, OC = "AM")

# print results
results

#}
}

\keyword{Econometrics}
\keyword{Dichotomous}
\keyword{IRF}
\keyword{Confidence-Intervals}


\name{GIRF_Dicho}
\alias{GIRF_Dicho}
\title{
GIRF for Dichotomous models
}
\description{
This function estimates the response functions of dichotomous models in a univariate framework using the method proposed by Lajaunie (2021). The response functions are based on the 4 specifications proposed by Kauppi & Saikkonen (2008).
}
\usage{
GIRF_Dicho(Dicho_Y, Exp_X, Lag, Int, t_mod, horizon, shock_size, OC)
}
\arguments{

  \item{Dicho_Y}{
Vector of the binary time series.
}
  \item{Exp_X}{
Vector or Matrix of explanatory time series.
}
  \item{Lag}{
Number of lags used for the estimation.
}
  \item{Int}{
Boolean value: TRUE for an estimation with intercept, and FALSE otherwise.
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
  \item{horizon}{
Numeric variable corresponding to the horizon target for the GIRF analysis.
}
  \item{shock_size}{
Numeric variable equal to the size of the shock. It can be estimated with the Vector_Error function.
}
  \item{OC}{
Numeric variable equal to the Optimal Cut-off (threshold). This threshold can be considered arbitrarily, with a value between 0 and 1, or it can be estimated with one of the functions EWS_AM_Criterion, EWS_CSA_Criterion, or EWS_NSR_Criterion.
}
}

\value{

Matrix with 7 columns:

  \item{column 1}{horizon}
  \item{column 2}{index}
  \item{column 3}{index with shock}
  \item{column 4}{probability associated to the index}
  \item{column 5}{probability associated to the index with shock}
  \item{column 6}{binary variable associated to the index}
  \item{column 7}{binary variable associated to the index with shock}
}

\references{
Kauppi, Heikki, and Pentti Saikkonen. "Predicting US recessions with dynamic binary response models." The Review of Economics and Statistics 90.4 (2008): 777-791.

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

# Estimate the logit regression
Logistic_results <- Logistic_Estimation(Dicho_Y = Var_Y, Exp_X = Var_X, Intercept = TRUE,
                      Nb_Id = 1, Lag = 1, type_model = 1)

# Vector of probabilities
vector_proba <- as.vector(rep(0,length(Var_Y)-1))
vector_proba <- Logistic_results$prob

# Vector of binary variables
Lag <- 1
vector_binary <- as.vector(rep(0,length(Var_Y)-1))
vector_binary <- Var_Y[(1+Lag):length(Var_Y)]

# optimal cut-off that maximizes the AM criterion
Threshold_AM <- EWS_AM_Criterion(Var_Proba = vector_proba, Dicho_Y = vector_binary,
                      cutoff_interval = 0.0001)

# Estimate the estimation errors
Residuals <- Vector_Error(Dicho_Y = Var_Y, Exp_X = Var_X, Intercept = TRUE,
                      Nb_Id = 1, Lag = 1, type_model = 1)

# Initialize the shock
size_shock <- quantile(Residuals,0.95)

# GIRF Analysis
results <- GIRF_Dicho(Dicho_Y = Var_Y, Exp_X = Var_X, Lag = 1, Int = TRUE, t_mod = 1,
                      horizon = 3, shock_size = size_shock, OC = Threshold_AM)

# print results
results

#}
}


\keyword{Econometrics}
\keyword{Dichotomous}
\keyword{IRF}

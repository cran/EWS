\name{EWS_NSR_Criterion}
\alias{EWS_NSR_Criterion}
\title{
NSR Threshold - optimal cut-off
}
\description{
This function provides a method to compute the optimal NSR (Noise to Signal Ratio) criterion proposed by Kaminsky, Lizondo and Reinhart (1998). As defined in Candelon, Dumitrescu and Hurlin (2012), the NSR represents the ratio of the false alarms (type II error) to the number of crises correctly identified by the EWS for a given cut-off. The optimal cut-off minimizes the NSR criterion.
}
\usage{
EWS_NSR_Criterion(Var_Proba, Dicho_Y, cutoff_interval)
}
\arguments{
  \item{Var_Proba}{
    Vector containing the estimated probabilities obtained with the Logistic Estimation function.
    }
  \item{Dicho_Y}{
    Vector of the binary time series.
  }
  \item{cutoff_interval}{
    Numeric variable between 0 and 1.
  }
}


\value{
A numeric variable containing the optimal cut-off that minimizes the NSR criterion.
}

\references{

Candelon, Bertrand, Elena-Ivona Dumitrescu, and Christophe Hurlin. "How to evaluate an early-warning system: Toward a unified statistical framework for assessing financial crises forecasting methods." IMF Economic Review 60.1 (2012): 75-113.

Kaminsky, Graciela, Saul Lizondo, and Carmen M. Reinhart. "Leading indicators of currency crises." IMF Staff Papers 45.1 (1998): 1-48.
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
                      Nb_Id = 1, Lag = 1, type_model = 4)

# Vector of probabilities
vector_proba <- as.vector(rep(0,length(Var_Y)-1))
vector_proba <- Logistic_results$prob

# Vector of binary variables
Lag <- 1
vector_binary <- as.vector(rep(0,length(Var_Y)-1))
vector_binary <- Var_Y[(1+Lag):length(Var_Y)]

# optimal cut-off that minimizes the NSR criterion
results <- EWS_NSR_Criterion(Var_Proba = vector_proba, Dicho_Y = vector_binary,
                      cutoff_interval = 0.0001)

# print results
results

#}
}

\keyword{Econometrics}
\keyword{Dichotomous}
\keyword{Threshold}

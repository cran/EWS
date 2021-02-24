\name{EWS_AM_Criterion}
\alias{EWS_AM_Criterion}
\title{
AM Threshold - optimal cut-off
}
\description{
This function provides a method to compute the optimal AM (Accuracy Measure) criterion. As defined in Candelon, Dumitrescu and Hurlin (2012), this approach consists in aggregating the number of crisis and calm periods correctly identified by the EWS. The optimal cut-off maximizes the number of correctly identified periods.
}
\usage{
EWS_AM_Criterion(Var_Proba, Dicho_Y, cutoff_interval)
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
A numeric variable containing the optimal cut-off that maximizes the higher proportion of calm and crisis periods correctly identified.
}

\references{

Candelon, Bertrand, Elena-Ivona Dumitrescu, and Christophe Hurlin. "How to evaluate an early-warning system: Toward a unified statistical framework for assessing financial crises forecasting methods." IMF Economic Review 60.1 (2012): 75-113.

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
                      Nb_Id = 1, Lag = 1, type_model = 4)

# Vector of probabilities
vector_proba <- as.vector(rep(0,length(Var_Y)-1))
vector_proba <- Logistic_results$prob

# Vector of binary variables
Lag <- 1
vector_binary <- as.vector(rep(0,length(Var_Y)-1))
vector_binary <- Var_Y[(1+Lag):length(Var_Y)]

# optimal cut-off that maximizes the AM criterion
results <- EWS_AM_Criterion(Var_Proba = vector_proba, Dicho_Y = vector_binary,
                      cutoff_interval = 0.0001)

# print results
results

#}
}

\keyword{Econometrics}
\keyword{Dichotomous}
\keyword{Threshold}

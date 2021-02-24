\name{BlockBootstrapp}
\alias{BlockBootstrapp}
\title{
Block Bootstrapp
}
\description{
This function enables the estimation of the block size for resampling. The size of the blocks is computed as in Hall, Horowitz and Jing (1995). Then, the function returns in a matrix the new resampled input variables. These variables are then used to determine the confidence intervals of the response functions proposed by Lajaunie (2021).
}
\usage{
BlockBootstrapp(Dicho_Y, Exp_X, Intercept, n_simul)
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
  \item{n_simul}{
    Numeric variable equal to the total number of replications.
  }
}


\value{
A matrix containing the replications of the new resampled input variables. The matrix contains \eqn{n \times S} colomns, where \eqn{n} denotes the number of input variables, and \eqn{S} denotes the number of replications.
}

\references{

Hall, Peter, Joel L. Horowitz, and Bing-Yi Jing. "On blocking rules for the bootstrap with dependent data." Biometrika 82.3 (1995): 561-574.


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

# Resample
results <- BlockBootstrapp(Dicho_Y = Var_Y, Exp_X = Var_X, Intercept = TRUE, n_simul = 100)

# print results
results

#}
}

\keyword{Econometrics}
\keyword{Bootstrapp}
\keyword{Confidence-Intervals}

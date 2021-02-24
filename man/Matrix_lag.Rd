\name{Matrix_lag}
\alias{Matrix_lag}
\title{
Matrix Lag - data processing
}
\description{
Compute a lagged version of a time series, shifting the time base back by a given number of observations defined by the user. The user must enter three parameters for this function: the matrix, the number of lags, and of boolean variable calls 'beginning'. If 'beginning'=TRUE, then the lag will be applied at the beginning of the matrix whereas if 'beginning'=FALSE, then the lag will be applied at the end of the matrix.
}
\usage{
Matrix_lag(Matrix_target, Nb_lag, beginning)
}
\arguments{
  \item{Matrix_target}{
Initial Matrix
}
  \item{Nb_lag}{
Number of lag
}
  \item{beginning}{
Boolean variable. If 'place'=TRUE, the lag is applied at the beginning of the matrix. If 'place'=FALSE, the lag is applied at the end of the matrix.
}

}


\value{
A numeric Matrix.
}

\examples{

# NOT RUN {

# Initialize the following matrix
Matrix_example <- matrix(data=(1:10), nrow=5, ncol=2)

# Use Matrix_lag
new_matrix <- Matrix_lag(Matrix_target = Matrix_example, Nb_lag = 1, beginning = TRUE)

new_matrix

# Results:
#> new_matrix
#     [,1] [,2]
#[1,]    2    7
#[2,]    3    8
#[3,]    4    9
#[4,]    5   10

#}
}


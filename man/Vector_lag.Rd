\name{Vector_lag}
\alias{Vector_lag}
\title{
Vector lag - data processing
}
\description{
Compute a lagged version of a time series, shifting the time base back by a given number of observations defined by the user. The user must enter three parameters for this function: the vector, the number of lags, and a boolean variable named 'beginning'. If 'beginning'=TRUE, then the lag will be applied at the beginning of the vector whereas if 'beginning'=FALSE, then the lag will be applied at the end of the vector.
}
\usage{
Vector_lag(Vector_target, Nb_lag, beginning)
}
\arguments{
  \item{Vector_target}{
Initial vector
}
  \item{Nb_lag}{
Number of lag
}
  \item{beginning}{
Boolean variable. If 'beginning'=TRUE, the lag is applied at the beginning of the vector. If 'beginning'=FALSE, the lag is applied at the end of the vector.
}
}


\value{
A numeric Vector.
}

\examples{

# NOT RUN {

# Initialize the following vector
vector_example <- as.vector(1:10)

# Use Vector_lag
new_vector <- Vector_lag(Vector_target = vector_example, Nb_lag = 2, beginning = TRUE)

new_vector
# Results:
#> new_vector
#[1]  3  4  5  6  7  8  9 10

#}
}

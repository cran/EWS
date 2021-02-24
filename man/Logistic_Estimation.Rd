\name{Logistic_Estimation}
\alias{Logistic_Estimation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Logistic Estimation for Dichotomous Analysis
}
\description{
This function provides methods for estimating the four dichotomous models as in Kauppi & Saikkonen (2008). Based on a logit approach, models are estimated in a univariate or a balanced panel framework as in Candelon, Dumitrescu and Hurlin (2014). This estimation has been used in recent papers such in Ben Naceur, Candelon and Lajaunie (2019) and Hasse and Lajaunie (2020).
}
\usage{
Logistic_Estimation(Dicho_Y, Exp_X, Intercept, Nb_Id, Lag, type_model)
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

-> 1 for the static model:

\deqn{P_{t-1}(Y_{t}) = F(\pi_{t})=F(\alpha + \beta'X_{t})}


-> 2 for the dynamic model with lag binary variable:

\deqn{P_{t-1}(Y_{t}) = F(\pi_{t})=F(\alpha + \beta'X_{t} + \gamma Y_{t-l})}


-> 3 for the dynamic model with lag index variable:

\deqn{P_{t-1}(Y_{t}) = F(\pi_{t})=F(\alpha + \beta'X_{t} + \eta \pi_{t-l})}


-> 4 for the dynamic model with both lag binary variable and lag index variable:

\deqn{P_{t-1}(Y_{t}) = F(\pi_{t})=F(\alpha + \beta'X_{t} + \eta \pi_{t-l} + \gamma Y_{t-l})}

}

}

\value{

A list with:

  \item{Estimation}{a dataframe containing the coefficients of the logitic estimation, the Standard Error for each coefficient, the Z-score and the associated critical probability}
  \item{AIC}{a numeric vector containing the Akaike information criterion}
  \item{BIC}{a numeric vector containing the Bayesian information criterion}
  \item{R2}{a numeric vector containing the Pseudo R Square}
  \item{index}{a numeric vector containing the estimated index}
  \item{prob}{a numeric vector containing the estimated probabilities}
  \item{LogLik}{a numeric vector containing the Log likelihood value of the estimation}
  \item{VCM}{a numeric matrix of the Variance Covariance of the estimation}
}
\references{

Candelon, Bertrand, Elena-Ivona Dumitrescu, and Christophe Hurlin. "Currency crisis early warning systems: Why they should be dynamic." International Journal of Forecasting 30.4 (2014): 1016-1029.

Hasse, Jean-Baptiste, Lajaunie Quentin. "Does the Yield Curve Signal Recessions? New Evidence from an International Panel Data Analysis." (2020)

Kauppi, Heikki, and Pentti Saikkonen. "Predicting US recessions with dynamic binary response models." The Review of Economics and Statistics 90.4 (2008): 777-791.

Naceur, Sami Ben, Bertrand Candelon, and Quentin Lajaunie. "Taming financial development to reduce crises." Emerging Markets Review 40 (2019): 100618.
}
\author{
Jean-Baptiste Hasse and Quentin Lajaunie
}
\note{
For the panel estimation, data must be stacked one after the other for each country or for each individual.
}


\examples{

# NOT RUN {

# Import data
data("data_USA")

# Data process
Var_Y <- as.vector(data_USA$NBER)
Var_X <- as.vector(data_USA$Spread)

# Estimate the logit regression
results <- Logistic_Estimation(Dicho_Y = Var_Y, Exp_X = Var_X, Intercept = TRUE,
                      Nb_Id = 1, Lag = 1, type_model = 1)

# print results
results

# }

}

\keyword{Econometrics}
\keyword{Panel}
\keyword{Dichotomous}






# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%% Function vector lag  %%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Input:
# 	- Vector_target: a vector of size [T]
# 	- Nb_lag: number of lags
# 	- beginning: TRUE: lag at the beginning
# 		           FALSE: lag at the end
# Output: new vector of size [T-nb_lag]


Vector_lag <- function(Vector_target, Nb_lag, beginning){
  if (beginning==TRUE){
    Vector_target[1:Nb_lag] <- NA
    Vector_target <- Vector_target[!is.na(Vector_target)]
    results <- as.vector(Vector_target)
    return(results)
  }else{
    size_vector <- length(Vector_target)
    Vector_target[(size_vector+1-Nb_lag):size_vector] <- NA
    Vector_target <- Vector_target[!is.na(Vector_target)]
    results <- as.vector(Vector_target)
    return(results)
  }
}



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%% Function matrix lag  %%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Input:
# 	- Matrix_target: a matrix of size [T,n]
# 	- Nb_lag: number of lags
# 	- beginning: TRUE: lag at the beginning
# 		           FALSE: lag at the end
# Output: Matrix of size [T-nb_lag,n]


Matrix_lag <- function(Matrix_target, Nb_lag, beginning){

  ncol_matrix <- ncol(Matrix_target)
  nrow_matrix <- nrow(Matrix_target) - Nb_lag

  Var_transition<- matrix(0, nrow_matrix , ncol_matrix )

  if (beginning==TRUE){
    Matrix_target[1:Nb_lag,] <- NA
    Matrix_target <- Matrix_target[!is.na(Matrix_target)]
    for (compteur_col in 1:ncol_matrix){
      Var_transition[,compteur_col] <- Matrix_target[ ( (compteur_col-1) * nrow_matrix + 1 ):( (compteur_col) * nrow_matrix ) ]
    }
    results <- as.matrix(Var_transition)
    return(results)
  }else{
    Matrix_target[(nrow_matrix +1):(nrow_matrix + Nb_lag),] <- NA
    Matrix_target <- Matrix_target[!is.na(Matrix_target)]
    for (compteur_col in 1:ncol_matrix){
      Var_transition[,compteur_col] <- Matrix_target[ ( (compteur_col-1) * nrow_matrix + 1 ):( (compteur_col) * nrow_matrix ) ]
    }
    results <- as.matrix(Var_transition)
    return(results)
  }
}

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%% Logistic Estimation Function - 4 models %%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Input:
# 	- Dicho_Y: Binary variable
# 	- Exp_X: Matrix of explanatory variables
#	- Intercept: Boolean variable that is equal to TRUE (with intercept) or FALSE (without)
# 	- Nb_Id: number of individuals
# 	- Lag: number of lag for the logistic estimation
# 	- type_model:
#			1: static
#			2: dynamic with the lag binary variable
#			3: dynamiC with the lag index variable
#			4: dynamiC with both lag binary and lag index variable
# Output:
#		- "coeff": Coefficients of the logit regression
#		- "AIC" : AIC criteria value
#		- "BIC" : BIC criteria value
#		- "R2" : R2 value
#		- "LogLik" : LogLikelihood value
#		- "VarCov" : matrix VarCov


Logistic_Estimation <- function(Dicho_Y, Exp_X, Intercept, Nb_Id, Lag, type_model)
{


  # --------------------------------------------------------------------------------
  # ------------------------------- Data Processing --------------------------------
  # --------------------------------------------------------------------------------

  T_size <- length(Dicho_Y) 					# number of periods for all sample
  Nb_periode <- T_size/Nb_Id 				# number of periods for each country


  # Add the intercept in the matrix of explanatory variables if Intercept=TRUE
  if (Intercept==TRUE){
    Cste <- as.vector(rep(1,T_size))
    Exp_X <- cbind(Cste, Exp_X)
  }else{
  }

  # For the model 2 and the model 4, we add the binary crisis in the matrix of explanatory variables
  if (type_model==2){
    Exp_X <- cbind(Exp_X, Dicho_Y)
  }else if (type_model==4){
    Exp_X <- cbind(Exp_X, Dicho_Y)
  }else{
  }

  # Creation of Var_X_Reg and Var_Y_Reg
  # 	- Var_X_Reg: vector or matrix of explanatory variables for the logit regression
  # 	- Var_Y_Reg: vector of the binary variable for the logit regression


  if (Nb_Id == 1){ # If there is only one Id (univariate case)

    if ( length(Exp_X) == length(Dicho_Y) ) {

      # Initialization
      Var_X_Reg <- as.vector( rep ( 0 , (Nb_periode- Lag) ) )
      Var_Y_Reg <- as.vector( rep ( 0 , (Nb_periode- Lag) ) )

      # data processing with lag
      Var_X_Reg[ 1 : (Nb_periode - Lag) ] <- Vector_lag( (Exp_X[1:Nb_periode]) , Lag , FALSE )
      Var_Y_Reg[ 1 : (Nb_periode - Lag) ] <- Vector_lag( (Dicho_Y[1:Nb_periode]) , Lag , TRUE )

    } else {

      # Initialization
      Var_X_Reg <- matrix(data = 0, nrow = (Nb_periode- Lag), ncol = (ncol(Exp_X)) )
      Var_Y_Reg <- as.vector(rep(0, (Nb_periode- Lag)))

      # data processing with lag
      Var_X_Reg[1:(Nb_periode- Lag),] <- Matrix_lag( (Exp_X[1:Nb_periode,]) , Lag , FALSE )
      Var_Y_Reg[1:(Nb_periode- Lag)] <- Vector_lag( (Dicho_Y[1:Nb_periode]) , Lag , TRUE )
    }

  }else{ # balanced panel case

    # Creation of the fixed effects matrix
    Fixed_Effects <- matrix( data=0, ncol = Nb_Id-1 , nrow=T_size )

    for ( compteur in 1:(Nb_Id-1) )	{
      Fixed_Effects[((compteur-1)*Nb_periode+1):(compteur*Nb_periode),compteur] <- 1
    }

    # New matrix of explanatory variables
    Exp_X <- cbind( Exp_X , Fixed_Effects )

    # Initialization
    Var_X_Reg <- matrix( data = 0 , nrow = (Nb_periode- Lag) * Nb_Id , ncol = ( ncol(Exp_X) ) )
    Var_Y_Reg <- as.vector( rep ( 0 , Nb_Id * (Nb_periode- Lag) ) )

    # data processing with lag
    for (compteur_Id in 1:Nb_Id){
      Var_X_Reg[ ( 1 + ( (compteur_Id-1) * (Nb_periode- Lag) ) ) : ( (compteur_Id) * (Nb_periode- Lag) ),] <- Matrix_lag( (Exp_X[ ( 1 + ( (compteur_Id-1) * Nb_periode ) ) :( (compteur_Id) * Nb_periode ),]) , Lag , FALSE )
      Var_Y_Reg[ ( 1 + ( (compteur_Id-1) * (Nb_periode- Lag) ) ) : ( compteur_Id * (Nb_periode- Lag) ) ] <- Vector_lag( (Dicho_Y[( 1 + ( (compteur_Id-1) * Nb_periode ) ) :( (compteur_Id) * Nb_periode )]) , Lag , TRUE )
    }
  }


  # --------------------------------------------------------------------------------
  # ----------------------------- Logistic Estimation ------------------------------
  # --------------------------------------------------------------------------------


  # ----------------------------------- Model 1 -----------------------------------

  if (type_model==1) {


    if (length(Var_X_Reg)==length(Var_Y_Reg)){
      # Initialization - simple OLS
      coeff_initialize <- solve(t(Var_X_Reg)%*%(Var_X_Reg))%*%t(Var_X_Reg)%*%(Var_Y_Reg)
      coeff_initialize <- coeff_initialize/0.25
    }else{
      # Initialization - simple OLS
      coeff_initialize <- as.vector(rep(0,ncol(Var_X_Reg)))
      coeff_initialize <- solve(t(Var_X_Reg)%*%(Var_X_Reg))%*%t(Var_X_Reg)%*%(Var_Y_Reg)

      if (Intercept==TRUE) {
        coeff_initialize[1] <- (coeff_initialize[1]-0.4)/0.25

        for (compteur in 2:ncol(Var_X_Reg)) {
          coeff_initialize[compteur] <- coeff_initialize[compteur]/0.25
        }
      }else {
        for (compteur in 1:ncol(Var_X_Reg)) {
          coeff_initialize[compteur] <- coeff_initialize[compteur]/0.25
        }
      }
    }

    F_mod1 <- function(par) {

      # Sample Size
      T_size_function <- length(Var_Y_Reg)

      # Number of Coeff
      last <- length(par)

      # Vector of parameters except the lagged index
      beta <- as.vector(rep(0,last))
      beta <- par[1:(last)]

      # =================================
      # === Construction of the index ===
      # =================================

      # Initialization of the vector of probability
      expf <- as.vector(rep(0,T_size_function))

      for (compteur in 1:T_size_function) {
        if (length(Var_X_Reg)==T_size_function) {
          expf[compteur] <- exp(Var_X_Reg[compteur] %*% beta)/(1+exp(Var_X_Reg[compteur] %*% beta))
        }else{
          expf[compteur] <- exp(Var_X_Reg[compteur,] %*% beta)/(1+exp(Var_X_Reg[compteur,] %*% beta))
        }
      }

      prob <- as.vector(rep(0,T_size_function))

      for (compteur in 1:T_size_function) {
        prob[compteur] <- min(max(expf[compteur],0.0000001),0.999999)
      }

      # Vector of individual log-likelihood
      lik <- Var_Y_Reg*log(prob)+(1-Var_Y_Reg)*log(1-prob)

      # log-likelihood maximization function
      return(-sum(lik))
    }

    if (length(Var_X_Reg)==length(Var_Y_Reg)) {
      results <- optim(coeff_initialize, F_mod1, gr = NULL,
                       lower = -1, upper = 1, method = "Brent", control = list(maxit = 50000, factr = TRUE, abstol=0.00001), hessian = FALSE)
    }else{
      results <- optim(coeff_initialize, F_mod1, gr = NULL, method = "Nelder-Mead", control = list(maxit = 50000, factr = TRUE, abstol=0.00001), hessian = FALSE)
    }

    # Estimated parameters
    if (length(Var_X_Reg)==length(Var_Y_Reg)){
      Beta <- results$par[1]
    }else{
      Beta <- as.vector(rep(0,ncol(Var_X_Reg)))
      Beta <- results$par[1:ncol(Var_X_Reg)]
    }

    # Approximation Hessian Matrix
    hessc <- hessian(func=F_mod1, x=Beta , "Richardson")

    # Initialisation of the vector of index
    T_size <- length(Var_Y_Reg)

    # Initialization of the vector of probability
    expf <- as.vector(rep(0,T_size))

    # Initialization of the vector of density
    pdf <- as.vector(rep(0,T_size))

    for (compteur in 1:T_size) {
      if (length(Var_X_Reg)==length(Var_Y_Reg)) {
        expf[compteur] <- exp(Var_X_Reg[compteur] %*% Beta)/(1+exp(Var_X_Reg[compteur] %*% Beta))
        pdf[compteur] <- exp(Var_X_Reg[compteur] %*% Beta)/((1+exp(Var_X_Reg[compteur] %*% Beta))^2)
      }else{
        expf[compteur] <- exp(Var_X_Reg[compteur,] %*% Beta)/(1+exp(Var_X_Reg[compteur,] %*% Beta))
        pdf[compteur] <- exp(Var_X_Reg[compteur,] %*% Beta)/((1+exp(Var_X_Reg[compteur,] %*% Beta))^2)
      }
    }

    # The vector of estimated probabilities
    prob <- as.vector(rep(0,T_size))
    for (compteur in 1:T_size) {
      prob[compteur] <- min(max(expf[compteur],0.0000001),0.999999)
    }

    # Matrix of explanatory variables
    if (length(Var_X_Reg)==length(Var_Y_Reg)) {
      X <- as.vector(rep(0,T_size-1))
      X <- Var_X_Reg[2:T_size]
    }else{
      X <- matrix(0,nrow=T_size-1,ncol=(ncol(Var_X_Reg)))
      X <- Var_X_Reg[2:T_size,]
    }

    # Matrix of gradient
    vect_gradient <- as.vector(rep(0,T_size-1))
    vect_gradient <- (Var_Y_Reg[2:T_size]-prob[2:T_size])/(prob[2:T_size]*(1-prob[2:T_size]))*pdf[2:T_size]

    if (length(Var_X_Reg)==length(Var_Y_Reg)) {
      gradient <- as.vector(rep(0,(T_size-1)))
      gradient <- vect_gradient
    }else{
      gradient <- matrix(0,nrow=(T_size-1),ncol=(ncol(Var_X_Reg)))

      for (compteur in 1:(ncol(Var_X_Reg))) {
        gradient[,compteur] <- vect_gradient
      }
    }

    gradient <- gradient * X

    # Matrix of covariance of gradient
    I=t(gradient)%*%gradient/(T_size-1)

    # Bandwith parameter
    bdw <- floor(4*(T_size/100)^(2/9))

    for (compteur in 1:bdw) {
      u=abs(compteur/bdw)
      if ((0.5<u)&&(u<=1)) {
        w=2*(1-abs(u))^3
      }else if((0<u)&&(u<=0.5)) {
        w=1-6*abs(u)^2+6*abs(u)^3
      }

      # Matrix of estimated autovariance of gradient

      if (length(Var_X_Reg)==length(Var_Y_Reg)){
        Covariance <- t(gradient[(1+compteur):(T_size-1)]) %*% gradient[1:(T_size-1-compteur)]/(T_size-compteur)
      }else{
        Covariance <- t(gradient[(1+compteur):(T_size-1),]) %*% gradient[1:(T_size-1-compteur),]/(T_size-compteur)
      }

      # Matrix of asymptotic correction
      I <- I + w * (Covariance + t(Covariance))
    }


    # ===============
    # === Results ===
    # ===============

    # Estimated parameters
    if (length(Var_X_Reg)==length(Var_Y_Reg)) {
      result_param <- Beta
      ind <- as.vector(rep(0,T_size))
      ind <- Var_X_Reg * Beta
    }else{
      result_param <- as.vector(rep(0,ncol(Var_X_Reg)))
      result_param <- Beta
      ind <- as.vector(rep(0,T_size))
      ind <- Var_X_Reg %*% Beta
    }

    # Asymptotic Matrix of Var-Cov
    V <- solve(hessc)

    # Asymptotic standard errors (non robust)
    Std <- t(diag(sqrt(abs(V))))

    # Robust var-cov
    VCM <- V %*% I %*% V * T_size

    # Robust standard errors
    Rob_Std <- t(diag(sqrt(abs(VCM))))

    # Log-Likelihood
    loglikelihood <- -results$value

    if (length(Var_X_Reg)==length(Var_Y_Reg)){

      # AIC information criteria
      AIC <- -2*loglikelihood + 1
      # BIC information criteria
      BIC <- -2*loglikelihood + log(T_size)

    }else{

      # AIC information criteria
      AIC <- -2*loglikelihood + ncol(Var_X_Reg)

      # BIC information criteria
      BIC <- -2*loglikelihood + (ncol(Var_X_Reg))*log(T_size)
    }

    # R2
    Lc <- sum(Var_Y_Reg*log(mean(Var_Y_Reg))+(1-Var_Y_Reg)*log(1-mean(Var_Y_Reg)))
    R2 <- 1 - ((loglikelihood -ncol(Var_X_Reg))/Lc)

    # initialize the coefficients matrix results
    if (Intercept==TRUE){
      nb_Var_X <- ncol(Var_X_Reg)
      nameVarX <- as.vector(rep(0,nb_Var_X-1))
      nameVarX <- c(1:(nb_Var_X-1))
    }else{
      if (length(Var_X_Reg)==length(Var_Y_Reg)){
        nb_Var_X <- 1
        nameVarX <- 1
      }else{
        nb_Var_X <- ncol(Var_X_Reg)
        nameVarX <- as.vector(rep(0,nb_Var_X))
        nameVarX <- c(1:(nb_Var_X))
      }
    }

    name <- as.vector(rep(0,nb_Var_X))
    Estimate <- as.vector(rep(0,nb_Var_X))
    Std.Error <- as.vector(rep(0,nb_Var_X))
    zvalue <- as.vector(rep(0,nb_Var_X))
    Pr <- as.vector(rep(0,nb_Var_X))

    if (Intercept==TRUE){
      # DataFrame with coefficients and significativity with intercept
      Coeff.results <- data.frame(
        name = c("Intercept" ,nameVarX),
        Estimate = result_param,
        Std.Error = t(Rob_Std),
        zvalue = t(result_param/Rob_Std),
        Pr = (1-pnorm(q=abs(t(result_param/Rob_Std))))*2,
        stringsAsFactors = FALSE
      )
    }else{
      # DataFrame with coefficients and significativity
      Coeff.results <- data.frame(
        name = c(nameVarX),
        Estimate = result_param,
        Std.Error = t(Rob_Std),
        zvalue = t(result_param/Rob_Std),
        Pr = (1-pnorm(q=abs(t(result_param/Rob_Std))))*2,
        stringsAsFactors = FALSE
      )
    }

    results <- list(Estimation=Coeff.results,
                    AIC=AIC, BIC=BIC, R2=R2, LogLik=loglikelihood, VarCov=VCM)

    # ----------------------------------- Modele 2 -----------------------------------

  } else if (type_model==2) {

    # Initialization
    coeff_initialize <- as.vector(rep(0,ncol(Var_X_Reg)))
    coeff_initialize <- solve(t(Var_X_Reg)%*%(Var_X_Reg))%*%t(Var_X_Reg)%*%(Var_Y_Reg)

    if (Intercept==TRUE) {
      coeff_initialize[1] <- (coeff_initialize[1]-0.4)/0.25
      for (compteur in 2:ncol(Var_X_Reg)) {
        coeff_initialize[compteur] <- coeff_initialize[compteur]/0.25
      }
    }else{
      for (compteur in 1:ncol(Var_X_Reg)) {
        coeff_initialize[compteur] <- coeff_initialize[compteur]/0.25
      }
    }


    F_mod2 <- function(par) {

      # Sample Size
      T_size_function <- length(Var_Y_Reg)

      # Number of Coeff
      last <- length(par)

      # Vector of parameters except the lagged index
      beta <- as.vector(rep(0,last))
      beta <- par[1:(last)]

      # =================================
      # === Construction of the index ===
      # =================================

      # Initialization of the vector of probability
      expf <- as.vector(rep(0,T_size_function))

      for (compteur in 1:T_size_function) {
        expf[compteur] <- exp(Var_X_Reg[compteur,] %*% beta)/(1+exp(Var_X_Reg[compteur,] %*% beta))
      }

      prob <- as.vector(rep(0,T_size_function))
      for (compteur in 1:T_size_function) {
        prob[compteur] <- min(max(expf[compteur],0.0000001),0.999999)
      }

      # Vector of individual log-likelihood
      lik <- Var_Y_Reg*log(prob)+(1-Var_Y_Reg)*log(1-prob)

      # log-likelihood maximization function
      return(-sum(lik))
    }

    results <- optim(coeff_initialize, F_mod2, gr = NULL, method = "Nelder-Mead", control = list(maxit = 500000, factr = TRUE, abstol=0.00001), hessian = FALSE)

    # Estimated parameters
    Beta <- as.vector(rep(0,ncol(Var_X_Reg)))
    Beta <- results$par[1:ncol(Var_X_Reg)]

    # Approximation Hessian Matrix
    hessc <- hessian(func=F_mod2, x=Beta , "Richardson")

    # Initialization of the vector of index
    T_size <- length(Var_Y_Reg)

    # Initialization of the vector of probability
    expf <- as.vector(rep(0,T_size))

    # Initialization of the vector of density
    pdf <- as.vector(rep(0,T_size))

    for (compteur in 1:T_size) {
      expf[compteur] <- exp(Var_X_Reg[compteur,] %*% Beta)/(1+exp(Var_X_Reg[compteur,] %*% Beta))
      pdf[compteur] <- exp(Var_X_Reg[compteur,] %*% Beta)/((1+exp(Var_X_Reg[compteur,] %*% Beta))^2)
    }

    # The vector of estimated probabilities
    prob <- as.vector(T_size)

    for (compteur in 1:T_size) {
      prob[compteur] <- min(max(expf[compteur],0.0000001),0.999999)
    }

    # Matrix of explicative variables
    X <- matrix(0,nrow=T_size-1,ncol=(ncol(Var_X_Reg)))
    X <- Var_X_Reg[2:T_size,]

    # Matrix of gradient
    vect_gradient <- as.vector(rep(0,T_size-1))
    vect_gradient <- (Var_Y_Reg[2:T_size]-prob[2:T_size])/(prob[2:T_size]*(1-prob[2:T_size]))*pdf[2:T_size]

    gradient <- matrix(0,nrow=(T_size-1),ncol=(ncol(Var_X_Reg)))
    for (compteur in 1:(ncol(Var_X_Reg))) {
      gradient[,compteur] <- vect_gradient
    }

    gradient <- gradient * X

    # Matrix of covariance of gradient
    I=t(gradient)%*%gradient/(T_size-1)

    # Bandwith parameter
    bdw <- floor(4*(T_size/100)^(2/9))

    for (compteur in 1:bdw) {
      u=abs(compteur/bdw)

      if ((0.5<u)&&(u<=1)){
        w=2*(1-abs(u))^3
      }else if((0<u)&&(u<=0.5)){
        w=1-6*abs(u)^2+6*abs(u)^3
      }

      # Matrix of estimated autovariance of gradient
      Covariance <- t(gradient[(1+compteur):(T_size-1),]) %*% gradient[1:(T_size-1-compteur),]/(T_size-compteur)

      # Matrix of asymptotic correction
      I <- I + w * (Covariance + t(Covariance))
    }


    # ===============
    # === Results ===
    # ===============

    # Estimated parameters
    result_param <- as.vector(rep(0,ncol(Var_X_Reg)))
    result_param <- Beta

    ind <- as.vector(rep(0,T_size))
    ind <- Var_X_Reg %*% Beta

    # Asymptotic Matrix of Var-Cov
    V <- solve(hessc)

    # Asymptotic standard errors (non robust)
    Std <- t(diag(sqrt(abs(V))))

    # Robust var-cov
    VCM <- V %*% I %*% V * T_size

    # Robust standard errors
    Rob_Std <- t(diag(sqrt(abs(VCM))))

    # Log-Likelihood
    loglikelihood <- -results$value

    # AIC information criteria
    AIC <- -2*loglikelihood + ncol(Var_X_Reg)

    # BIC information criteria
    BIC <- -2*loglikelihood + (ncol(Var_X_Reg))*log(T_size)

    # R2
    Lc <- sum(Var_Y_Reg*log(mean(Var_Y_Reg))+(1-Var_Y_Reg)*log(1-mean(Var_Y_Reg)))
    R2 <- 1 - ((loglikelihood -ncol(Var_X_Reg))/Lc)


    # initialize the coefficients matrix results
    if (Intercept==TRUE){
      nb_Var_X <- ncol(Var_X_Reg)
      nameVarX <- as.vector(rep(0,nb_Var_X-2))
      nameVarX <- c(1:(nb_Var_X-2))
    }else{
      nb_Var_X <- ncol(Var_X_Reg)
      nameVarX <- as.vector(rep(0,nb_Var_X-1))
      nameVarX <- c(1:(nb_Var_X-1))
    }

    name <- as.vector(rep(0,nb_Var_X))
    Estimate <- as.vector(rep(0,nb_Var_X))
    Std.Error <- as.vector(rep(0,nb_Var_X))
    zvalue <- as.vector(rep(0,nb_Var_X))
    Pr <- as.vector(rep(0,nb_Var_X))


    if (Intercept==TRUE){
      # DataFrame with coefficients and significativity with intercept
      if (Nb_Id > 1){
        Coeff.results <- data.frame(
          name = c("Intercept" , nameVarX[1:(nb_Var_X-2-Nb_Id+1)] ,"Binary_Lag", nameVarX[(nb_Var_X-2-Nb_Id+2):(nb_Var_X-2)]),
          Estimate = result_param,
          Std.Error = t(Rob_Std),
          zvalue = t(result_param/Rob_Std),
          Pr = (1-pnorm(q=abs(t(result_param/Rob_Std))))*2,
          stringsAsFactors = FALSE
        )
      }else{
        Coeff.results <- data.frame(
          name = c("Intercept" , nameVarX ,"Binary_Lag"),
          Estimate = result_param,
          Std.Error = t(Rob_Std),
          zvalue = t(result_param/Rob_Std),
          Pr = (1-pnorm(q=abs(t(result_param/Rob_Std))))*2,
          stringsAsFactors = FALSE
        )
      }
    }else{

      # DataFrame with coefficients and significativity
      if (Nb_Id>1){
        Coeff.results <- data.frame(
          name = c(nameVarX[1:(nb_Var_X-1-Nb_Id+1)] ,"Binary_Lag", nameVarX[(nb_Var_X-1-Nb_Id+2):(nb_Var_X-1)]),
          Estimate = result_param,
          Std.Error = t(Rob_Std),
          zvalue = t(result_param/Rob_Std),
          Pr = (1-pnorm(q=abs(t(result_param/Rob_Std))))*2,
          stringsAsFactors = FALSE
        )
      }else{
        Coeff.results <- data.frame(
          name = c(nameVarX ,"Binary_Lag"),
          Estimate = result_param,
          Std.Error = t(Rob_Std),
          zvalue = t(result_param/Rob_Std),
          Pr = (1-pnorm(q=abs(t(result_param/Rob_Std))))*2,
          stringsAsFactors = FALSE
        )
      }
    }

    results <- list(Estimation=Coeff.results,
                    AIC=AIC, BIC=BIC, R2=R2, LogLik=loglikelihood, VarCov=VCM)

    # ----------------------------------- Modele 3 -----------------------------------

  } else if (type_model==3){

    if ( length(Var_X_Reg) == length(Var_Y_Reg) ) {

      # Initialization of coefficients
      coeff_initialize <- as.vector( rep(0,2) )
      coeff_initialize[1] <- solve( t(Var_X_Reg) %*% (Var_X_Reg) ) %*% t(Var_X_Reg) %*% (Var_Y_Reg)

      # Initialization of the Index
      Pi <- 0
      coeff_initialize[2] <- Pi

    }else{

      # Initialization of coefficients
      coeff_initialize <- as.vector( rep ( 0, ncol(Var_X_Reg) + 1 ) )
      coeff_initialize[1:ncol(Var_X_Reg)] <- solve(t(Var_X_Reg)%*%(Var_X_Reg))%*%t(Var_X_Reg)%*%(Var_Y_Reg)

      # Initialization of the Index
      Pi <- 0
      coeff_initialize[ncol(Var_X_Reg)+1] <- Pi

    }

    if (Intercept==TRUE) {
      coeff_initialize[1] <- (coeff_initialize[1]-0.4)/0.25

      for (compteur in 2:ncol(Var_X_Reg)) {
        coeff_initialize[compteur] <- coeff_initialize[compteur]/0.25
      }
    }else{
      if ( length(Var_X_Reg) == length(Var_X_Reg) ) {
        coeff_initialize <- coeff_initialize/0.25
      }else{
        for (compteur in 1:ncol(Var_X_Reg)) {
          coeff_initialize[compteur] <- coeff_initialize[compteur]/0.25
        }
      }
    }



    F_mod3 <- function(par) {

      # Sample Size
      T_size_function <- length(Var_Y_Reg)

      # Number of Coeff
      last <- length(par)

      # Vector of parameters except the lagged index
      beta <- as.vector(rep(0,last-1))
      beta <- par[1:(last-1)]

      # Parameter of the lagged index (here a logistic transformation)
      alpha <- par[last]/(1+abs(par[last]))

      # =================================
      # === Construction of the index ===
      # =================================

      # Initialization of the vector of the index
      ind <- as.vector(rep(0,T_size_function))

      # Initialization of the vector of probability
      expf <- as.vector(rep(0,T_size_function))

      # Initial value for the index
      mean_Var_X <- as.vector(rep(0,(last-1)))

      for (compteur in 1:(last-1)) {
        if ( length(Var_X_Reg) == length(Var_Y_Reg) ) {
          mean_Var_X[compteur] <- mean(Var_X_Reg)
          p0 <- mean_Var_X * beta/(1-alpha)
          ind[1] <- p0
          expf[1] <- exp(ind[1])/(1+exp(ind[1]))
        } else {
          mean_Var_X[compteur] <- mean(Var_X_Reg[,compteur])
          p0 <- mean_Var_X %*% beta/(1-alpha)
          ind[1] <- p0
          expf[1] <- exp(ind[1])/(1+exp(ind[1]))
        }
      }

      for (compteur in 2:T_size_function) {
        if ( length(Var_X_Reg) == length(Var_Y_Reg) ) {
          ind[compteur] <- alpha * ind[compteur-1] + Var_X_Reg[compteur] * beta
        } else {
          ind[compteur] <- alpha * ind[compteur-1] + Var_X_Reg[compteur,] %*% beta
        }

        expf[compteur] <- exp(ind[compteur])/(1+exp(ind[compteur]))
      }

      prob <- as.vector(rep(0,T_size_function))
      for (compteur in 1:T_size_function) {
        prob[compteur] <- min(max(expf[compteur],0.0000001),0.999999)
      }

      # Vector of individual log-likelihood
      lik <- Var_Y_Reg*log(prob)+(1-Var_Y_Reg)*log(1-prob)

      # We do not consider the initial condition on the index
      estim_lik <- as.vector(rep(0,(T_size_function-1)))
      estim_lik <- lik[2:T_size_function]

      # log-likelihood maximization function
      return(-sum(lik))
    }

    results <- optim(coeff_initialize, F_mod3, gr = NULL, method = "Nelder-Mead", control = list(maxit = 500000, factr = TRUE, abstol=0.00001), hessian = FALSE)

    # Logistic transformation of alpha (inverse)
    if ( length(Var_X_Reg) == length(Var_Y_Reg) ) {
      psi <- results$par[2]

      # Estimated parameter alpha
      alpha <- psi/(1+abs(psi))

      # Intermediate element required by tge analytical gradient
      h <- 0.000001

      # Analytical Gradient for alpha / psi
      d2 <- (((psi+h)/(1+abs(psi+h)))-((psi-h)/(1+abs(psi-h))))/(2*h)

      # Matrix of Taylor development of the logistic function
      C <- diag(2)
      C[2,2] <- d2
    }else{
      psi <- results$par[ncol(Var_X_Reg)+1]

      # Estimated parameter alpha
      alpha <- psi/(1+abs(psi))

      # Intermediate element required by tge analytical gradient
      h <- 0.000001

      # Analytical Gradient for alpha / psi
      d2 <- (((psi+h)/(1+abs(psi+h)))-((psi-h)/(1+abs(psi-h))))/(2*h)

      # Matrix of Taylor development of the logistic function
      C <- diag(ncol(Var_X_Reg)+1)
      C[(ncol(Var_X_Reg)+1),(ncol(Var_X_Reg)+1)] <- d2

    }

    # =======================
    # === Correction term ===
    # =======================

    # Estimated parameters
    if ( length(Var_X_Reg) == length(Var_Y_Reg) ) {
      Beta <- results$par[1]
    }else{
      Beta <- as.vector(rep(0,ncol(Var_X_Reg)))
      Beta <- results$par[1:ncol(Var_X_Reg)]
    }

    # Initialization of the vector of index
    T_size <- length(Var_Y_Reg)
    ind <- as.vector(rep(0,T_size))

    # Initialization of the vector of probability
    expf <- as.vector(rep(0,T_size))

    # Initialization of the vector of density
    pdf <- as.vector(rep(0,T_size))

    # Initial value for the index
    if ( length(Var_X_Reg) == length(Var_Y_Reg) ) {
      last <- 2
      mean_Var_X <- mean(Var_X_Reg)
      p0 <- mean_Var_X * Beta/(1-alpha)
      ind[1] <- p0
      expf[1] <- exp(ind[1])/(1+exp(ind[1]))
      pdf[1] <- exp(ind[1])/((1+exp(ind[1]))^2)

      for (compteur in 2:T_size) {
        ind[compteur] <- alpha * ind[compteur-1] + Var_X_Reg[compteur] * Beta
        expf[compteur] <- exp(ind[compteur])/(1+exp(ind[compteur]))
        pdf[compteur] <- exp(ind[compteur])/((1+exp(ind[compteur]))^2)
      }
    }else{
      last <- ncol(Var_X_Reg)+1
      mean_Var_X <- as.vector(rep(0,(last-1)))

      for (compteur in 1:(last-1)) {
        mean_Var_X[compteur] <- mean(Var_X_Reg[,compteur])
      }

      p0 <- mean_Var_X %*% Beta/(1-alpha)
      ind[1] <- p0
      expf[1] <- exp(ind[1])/(1+exp(ind[1]))
      pdf[1] <- exp(ind[1])/((1+exp(ind[1]))^2)

      for (compteur in 2:T_size) {
        ind[compteur] <- alpha * ind[compteur-1] + Var_X_Reg[compteur,] %*% Beta
        expf[compteur] <- exp(ind[compteur])/(1+exp(ind[compteur]))
        pdf[compteur] <- exp(ind[compteur])/((1+exp(ind[compteur]))^2)
      }
    }

    # ==============================
    # === Robust standard Errors ===
    # ==============================

    # The vector of estimated probabilities
    prob <- as.vector(T_size)
    for (compteur in 1:T_size) {
      prob[compteur] <- min(max(expf[compteur],0.0000001),0.999999)
    }

    # Initial value for the index
    if ( length(Var_X_Reg) == length(Var_Y_Reg) ) {

      # Matrix of explanatory variables (the first of observation is equal to p0)
      X <- matrix(0,nrow=(T_size-1),ncol=2)
      X[,1] <- Var_X_Reg[2:T_size]
      X[,2] <- ind[1:(T_size-1)]

      # Matrix of gradient
      vect_gradient <- as.vector(rep(0,(T_size-1)))
      vect_gradient <- (Var_Y_Reg[2:T_size]-prob[2:T_size])/(prob[2:T_size]*(1-prob[2:T_size]))*pdf[2:T_size]

      gradient <- matrix(0,nrow=(T_size-1),ncol=2)
      for (compteur in 1:2){
        gradient[,compteur] <- vect_gradient
      }
      gradient <- gradient * X

    }else{

      # Matrix of explanatory variables (the first of observation is equal to p0)
      X <- matrix(0,nrow=(T_size-1),ncol=(ncol(Var_X_Reg)+1))
      X[,1:(ncol(Var_X_Reg))] <- Var_X_Reg[2:T_size,]
      X[,ncol(Var_X_Reg)+1] <- ind[1:(T_size-1)]

      # Matrix of gradient
      vect_gradient <- as.vector(rep(0,T_size-1))
      vect_gradient <- (Var_Y_Reg[2:T_size]-prob[2:T_size])/(prob[2:T_size]*(1-prob[2:T_size]))*pdf[2:T_size]

      gradient <- matrix(0,nrow=(T_size-1),ncol=(ncol(Var_X_Reg)+1))

      for (compteur in 1:(ncol(Var_X_Reg)+1)) {
        gradient[,compteur] <- vect_gradient
      }

      gradient <- gradient * X
    }

    # Matrix of covariance of gradient
    I=t(gradient)%*%gradient/(T_size-1)

    # Bandwith parameter
    bdw <- floor(4*(T_size/100)^(2/9))

    for (compteur in 1:bdw) {
      u=abs(compteur/bdw)

      if ((0.5<u)&&(u<=1)) {
        w=2*(1-abs(u))^3
      }else if((0<u)&&(u<=0.5)){
        w=1-6*abs(u)^2+6*abs(u)^3
      }

      # Matrix of estimated autovariance of gradient
      Covariance <- t(gradient[(1+compteur):(T_size-1),]) %*% gradient[1:(T_size-1-compteur),]/(T_size-compteur)

      # Matrix of asymptotic correction
      I <- I + w * (Covariance + t(Covariance))
    }

    # ===============
    # === Results ===
    # ===============

    if ( length(Var_X_Reg) == length(Var_Y_Reg) ) {

      # Parameter alpha
      alpha <- psi / (1+abs(psi))

      # Estimated parameters
      result_param <- as.vector(rep(0,2))
      result_param <- c(Beta, alpha)

      # Hessian matrix
      hessc <- matrix(0,2,2)
      hessc <- hessian(func=F_mod3, x=result_param , "Richardson")

      # Non-robust variance-covariance matrix in the new space
      V <- C %*% solve(hessc) %*% C

      # Asymptotic Matrix of Var-Cov
      V_nonrobust <- V

      # Asymptotic standard errors (non robust)
      Std <- t(diag(sqrt(abs(V_nonrobust))))

      # Robust var-cov
      VCM <- V %*% I %*% V * T_size

      # Robust standard errors
      Rob_Std <- t(diag(sqrt(abs(VCM))))

      # Log-Likelihood
      loglikelihood <- - results$value

      # AIC information criteria
      AIC <- -2*loglikelihood + 2

      # BIC information criteria
      BIC <- -2*loglikelihood + 2*log(T_size)

      # R2
      Lc <- sum(Var_Y_Reg*log(mean(Var_Y_Reg))+(1-Var_Y_Reg)*log(1-mean(Var_Y_Reg)))
      R2 <- 1 - ((loglikelihood -1)/Lc)

    }else{

      # Parameter alpha
      alpha <- psi / (1+abs(psi))

      # Estimated parameters
      result_param <- as.vector(rep(0,ncol(Var_X_Reg)+1))
      result_param <- c(Beta, alpha)

      # Hessian matrix
      hessc <- matrix(0,(ncol(Var_X_Reg)+1),(ncol(Var_X_Reg)+1))
      hessc <- hessian(func=F_mod3, x=result_param , "Richardson")

      # Non-robust variance-covariance matrix in the new space
      V <- C %*% solve(hessc) %*% C

      # Asymptotic Matrix of Var-Cov
      V_nonrobust <- V

      # Asymptotic standard errors (non robust)
      Std <- t(diag(sqrt(abs(V_nonrobust))))

      # Robust var-cov
      VCM <- V %*% I %*% V * T_size

      # Robust standard errors
      Rob_Std <- t(diag(sqrt(abs(VCM))))

      # Log-Likelihood
      loglikelihood <- - results$value

      # AIC information criteria
      AIC <- -2*loglikelihood + ncol(Var_X_Reg)+1

      # BIC information criteria
      BIC <- -2*loglikelihood + (ncol(Var_X_Reg)+1)*log(T_size)

      # R2
      Lc <- sum(Var_Y_Reg*log(mean(Var_Y_Reg))+(1-Var_Y_Reg)*log(1-mean(Var_Y_Reg)))
      R2 <- 1 - ((loglikelihood -ncol(Var_X_Reg))/Lc)
    }

    # initialize the coefficients matrix results
    if (Intercept==TRUE){
      nb_Var_X <- ncol(Var_X_Reg)
      nameVarX <- as.vector(rep(0,nb_Var_X-1))
      nameVarX <- c(1:(nb_Var_X-1))
    }else{
      if (length(Var_X_Reg)==length(Var_Y_Reg)){
        nb_Var_X <- 1
        nameVarX <- 1
      }else{
        nb_Var_X <- ncol(Var_X_Reg)
        nameVarX <- as.vector(rep(0,nb_Var_X))
        nameVarX <- c(1:(nb_Var_X))
      }
    }

    name <- as.vector(rep(0,nb_Var_X+1))
    Estimate <- as.vector(rep(0,nb_Var_X+1))
    Std.Error <- as.vector(rep(0,nb_Var_X+1))
    zvalue <- as.vector(rep(0,nb_Var_X+1))
    Pr <- as.vector(rep(0,nb_Var_X+1))


    if (Intercept==TRUE){
      # DataFrame with coefficients and significativity with intercept
      if (Nb_Id>1){
        Coeff.results <- data.frame(
          name = c("Intercept" , nameVarX[1:(nb_Var_X-1-Nb_Id+1)] ,"Index_Lag", nameVarX[(nb_Var_X-1-Nb_Id+2):(nb_Var_X-1)]),
          Estimate = result_param,
          Std.Error = t(Rob_Std),
          zvalue = t(result_param/Rob_Std),
          Pr = (1-pnorm(q=abs(t(result_param/Rob_Std))))*2,
          stringsAsFactors = FALSE
        )
      }else{
        Coeff.results <- data.frame(
          name = c("Intercept" , nameVarX ,"Index_Lag"),
          Estimate = result_param,
          Std.Error = t(Rob_Std),
          zvalue = t(result_param/Rob_Std),
          Pr = (1-pnorm(q=abs(t(result_param/Rob_Std))))*2,
          stringsAsFactors = FALSE
        )
      }
    }else{
      # DataFrame with coefficients and significativity
      if (Nb_Id>1){
        Coeff.results <- data.frame(
          name = c(nameVarX[1:(nb_Var_X-Nb_Id+1)] ,"Index_Lag", nameVarX[(nb_Var_X-Nb_Id+2):(nb_Var_X)]),
          Estimate = result_param,
          Std.Error = t(Rob_Std),
          zvalue = t(result_param/Rob_Std),
          Pr = (1-pnorm(q=abs(t(result_param/Rob_Std))))*2,
          stringsAsFactors = FALSE
        )
      }else{
        Coeff.results <- data.frame(
          name = c(nameVarX ,"Index_Lag"),
          Estimate = result_param,
          Std.Error = t(Rob_Std),
          zvalue = t(result_param/Rob_Std),
          Pr = (1-pnorm(q=abs(t(result_param/Rob_Std))))*2,
          stringsAsFactors = FALSE
        )
      }
    }

    results <- list(Estimation=Coeff.results,
                    AIC=AIC, BIC=BIC, R2=R2, LogLik=loglikelihood, VarCov=VCM)

    # ----------------------------------- Modele 4 -----------------------------------
  } else if (type_model==4) {

    # Initialization of coefficients
    coeff_initialize <- as.vector(rep(0,ncol(Var_X_Reg)+1))
    coeff_initialize <- solve(t(Var_X_Reg)%*%(Var_X_Reg))%*%t(Var_X_Reg)%*%(Var_Y_Reg)

    # Initialisation of the index
    Pi <- 0

    coeff_initialize[ncol(Var_X_Reg)+1] <- Pi

    if (Intercept==1) {
      coeff_initialize[1] <- (coeff_initialize[1]-0.4)/0.25
      for (compteur in 2:ncol(Var_X_Reg)) {
        coeff_initialize[compteur] <- coeff_initialize[compteur]/0.25
      }
    }else{
      for (compteur in 1:ncol(Var_X_Reg)) {
        coeff_initialize[compteur] <- coeff_initialize[compteur]/0.25
      }
    }

    F_mod4 <- function(par) {

      # Sample Size
      T_size_function <- length(Var_Y_Reg)

      # Number of Coeff
      last <- length(par)

      # Vector of parameters except the lagged index
      beta <- as.vector(rep(0,last-1))
      beta <- par[1:(last-1)]

      # Parameter of the lagged index (here a logistic transformation)
      alpha <- par[last]/(1+abs(par[last]))

      # =================================
      # === Construction of the index ===
      # =================================

      # Initialization of the vector of the index
      ind <- as.vector(rep(0,T_size_function))

      # Initialization of the vector of probability
      expf <- as.vector(rep(0,T_size_function))

      # Initial value for the index
      mean_Var_X <- as.vector(rep(0,(last-1)))

      for (compteur in 1:(last-1)) {
        mean_Var_X[compteur] <- mean(Var_X_Reg[,compteur])
      }


      p0 <- mean_Var_X %*% beta/(1-alpha)

      ind[1] <- p0
      expf[1] <- exp(ind[1])/(1+exp(ind[1]))

      for (compteur in 2:T_size_function) {
        ind[compteur] <- alpha * ind[compteur-1] + Var_X_Reg[compteur,] %*% beta
        expf[compteur] <- exp(ind[compteur])/(1+exp(ind[compteur]))
      }

      prob <- as.vector(rep(0,T_size_function))
      for (compteur in 1:T_size_function) {
        prob[compteur] <- min(max(expf[compteur],0.0000001),0.999999)
      }

      # Vector of individual log-likelihood
      lik <- Var_Y_Reg*log(prob)+(1-Var_Y_Reg)*log(1-prob)

      # We do not consider the initial condition on the index
      estim_lik <- as.vector(rep(0,(T_size_function-1)))
      estim_lik <- lik[2:T_size_function]

      # log-likelihood maximization function
      return(-sum(lik))
    }

    results <- optim(coeff_initialize, F_mod4, gr = NULL, method = "Nelder-Mead", control = list(maxit = 500000, factr = TRUE, abstol=0.00001), hessian = FALSE)

    # Logistic transformation of alpha (inverse)
    psi <- results$par[ncol(Var_X_Reg)+1]

    # Estimated parameter alpha
    alpha <- psi/(1+abs(psi))

    # Intermediate element required by tge analytical gradient
    h <- 0.000001

    # Analytical Gradient for alpha / psi
    d2 <- (((psi+h)/(1+abs(psi+h)))-((psi-h)/(1+abs(psi-h))))/(2*h)

    # Matrix of Taylor development of the logistic function
    C <- diag(ncol(Var_X_Reg)+1)
    C[(ncol(Var_X_Reg)+1),(ncol(Var_X_Reg)+1)] <- d2


    # =======================
    # === Correction term ===
    # =======================

    # Estimated parameters
    Beta <- as.vector(rep(0,ncol(Var_X_Reg)))
    Beta <- results$par[1:ncol(Var_X_Reg)]

    # Initialization of the vector of index
    T_size <- length(Var_Y_Reg)
    ind <- as.vector(rep(0,T_size))

    # Initialization of the vector of probability
    expf <- as.vector(rep(0,T_size))

    # Initialisation of the vector of density
    pdf <- as.vector(rep(0,T_size))

    # Initial value for the index
    last <- ncol(Var_X_Reg)+1
    mean_Var_X <- as.vector(rep(0,(last-1)))
    for (compteur in 1:(last-1)) {
      mean_Var_X[compteur] <- mean(Var_X_Reg[,compteur])
    }

    p0 <- mean_Var_X %*% Beta/(1-alpha)

    ind[1] <- p0

    expf[1] <- exp(ind[1])/(1+exp(ind[1]))

    pdf[1] <- exp(ind[1])/((1+exp(ind[1]))^2)

    for (compteur in 2:T_size){
      ind[compteur] <- alpha * ind[compteur-1] + Var_X_Reg[compteur,] %*% Beta
      expf[compteur] <- exp(ind[compteur])/(1+exp(ind[compteur]))
      pdf[compteur] <- exp(ind[compteur])/((1+exp(ind[compteur]))^2)
    }

    # ==============================
    # === Robust standard Errors ===
    # ==============================

    # The vector of estimated probabilities
    prob <- as.vector(rep(0,T_size))
    for (compteur in 1:T_size){
      prob[compteur] <- min(max(expf[compteur],0.0000001),0.999999)
    }


    # Matrix of explanatory variables (the first of observation is equal to p0)
    X <- matrix(0,nrow=(T_size-1),ncol=(ncol(Var_X_Reg)+1))
    X[,1:(ncol(Var_X_Reg))] <- Var_X_Reg[2:T_size,]
    X[,ncol(Var_X_Reg)+1] <- ind[1:(T_size-1)]

    # Matrix of gradient
    vect_gradient <- as.vector(rep(0,T_size-1))
    vect_gradient <- (Var_Y_Reg[2:T_size]-prob[2:T_size])/(prob[2:T_size]*(1-prob[2:T_size]))*pdf[2:T_size]

    gradient <- matrix(0,nrow=(T_size-1),ncol=(ncol(Var_X_Reg)+1))
    for (compteur in 1:(ncol(Var_X_Reg)+1)) {
      gradient[,compteur] <- vect_gradient
    }

    gradient <- gradient * X

    # Matrix of covariance of gradient
    I=t(gradient)%*%gradient/(T_size-1)

    # Bandwith parameter
    bdw <- floor(4*(T_size/100)^(2/9))


    for (compteur in 1:bdw) {
      u=abs(compteur/bdw)

      if ((0.5<u)&&(u<=1)) {
        w=2*(1-abs(u))^3
      }else if((0<u)&&(u<=0.5)){
        w=1-6*abs(u)^2+6*abs(u)^3
      }

      # Matrix of estimated autovariance of gradient
      Covariance <- t(gradient[(1+compteur):(T_size-1),]) %*% gradient[1:(T_size-1-compteur),]/(T_size-compteur)

      # Matrix of asymptotic correction
      I <- I + w * (Covariance + t(Covariance))
    }

    # ===============
    # === Results ===
    # ===============


    # Parameter alpha
    alpha <- psi / (1+abs(psi))

    # Estimated parameters
    result_param <- as.vector(rep(0,ncol(Var_X_Reg)+1))
    result_param <- c(Beta, alpha)

    # Hessian matrix
    hessc <- matrix(0,(ncol(Var_X_Reg)+1),(ncol(Var_X_Reg)+1))
    hessc <- hessian(func=F_mod4, x=result_param , "Richardson")

    # Non-robust variance-covariance matrix in the new space
    V <- C %*% solve(hessc) %*% C

    # Asymptotic Matrix of Var-Cov
    V_nonrobust <- V

    # Asymptotic standard errors (non robust)
    Std <- t(diag(sqrt(abs(V_nonrobust))))

    # Robust var-cov
    VCM <- V %*% I %*% V * T_size

    # Robust standard errors
    Rob_Std <- t(diag(sqrt(abs(VCM))))

    # Log-Likelihood
    loglikelihood <- -results$value

    # AIC information criteria
    AIC <- -2*loglikelihood + ncol(Var_X_Reg)+1

    # BIC information criteria
    BIC <- -2*loglikelihood + (ncol(Var_X_Reg)+1)*log(T_size)

    # R2
    Lc <- sum(Var_Y_Reg*log(mean(Var_Y_Reg))+(1-Var_Y_Reg)*log(1-mean(Var_Y_Reg)))
    R2 <- 1 - ((loglikelihood -ncol(Var_X_Reg))/Lc)

    # initialize the coefficients matrix results
    if (Intercept==TRUE){
      nb_Var_X <- ncol(Var_X_Reg)
      nameVarX <- as.vector(rep(0,nb_Var_X-2))
      nameVarX <- c(1:(nb_Var_X-2))
    }else{
      nb_Var_X <- ncol(Var_X_Reg)
      nameVarX <- as.vector(rep(0,nb_Var_X-1))
      nameVarX <- c(1:(nb_Var_X-1))
    }

    name <- as.vector(rep(0,nb_Var_X+1))
    Estimate <- as.vector(rep(0,nb_Var_X+1))
    Std.Error <- as.vector(rep(0,nb_Var_X+1))
    zvalue <- as.vector(rep(0,nb_Var_X+1))
    Pr <- as.vector(rep(0,nb_Var_X+1))

    if (Intercept==TRUE){
      # DataFrame with coefficients and significativity with intercept
      if (Nb_Id>1){
        Coeff.results <- data.frame(
          name = c("Intercept" , nameVarX[1:(nb_Var_X-2-Nb_Id+1)],"Binary_Lag" ,"Index_Lag", nameVarX[(nb_Var_X-2-Nb_Id+2):(nb_Var_X-2)]),
          Estimate = result_param,
          Std.Error = t(Rob_Std),
          zvalue = t(result_param/Rob_Std),
          Pr = (1-pnorm(q=abs(t(result_param/Rob_Std))))*2,
          stringsAsFactors = FALSE
        )
      }else{
        Coeff.results <- data.frame(
          name = c("Intercept" , nameVarX , "Binary_Lag", "Index_Lag"),
          Estimate = result_param,
          Std.Error = t(Rob_Std),
          zvalue = t(result_param/Rob_Std),
          Pr = (1-pnorm(q=abs(t(result_param/Rob_Std))))*2,
          stringsAsFactors = FALSE
        )
      }
    }else{
      # DataFrame with coefficients and significativity
      if (Nb_Id>1){
        Coeff.results <- data.frame(
          name = c(nameVarX[1:(nb_Var_X-1-Nb_Id+1)],"Binary_Lag" ,"Index_Lag", nameVarX[(nb_Var_X-1-Nb_Id+2):(nb_Var_X-1)]),
          Estimate = result_param,
          Std.Error = t(Rob_Std),
          zvalue = t(result_param/Rob_Std),
          Pr = (1-pnorm(q=abs(t(result_param/Rob_Std))))*2,
          stringsAsFactors = FALSE
        )
      }else{
        Coeff.results <- data.frame(
          name = c(nameVarX , "Binary_Lag", "Index_Lag"),
          Estimate = result_param,
          Std.Error = t(Rob_Std),
          zvalue = t(result_param/Rob_Std),
          Pr = (1-pnorm(q=abs(t(result_param/Rob_Std))))*2,
          stringsAsFactors = FALSE
        )
      }
    }

    results <- list(Estimation=Coeff.results,
                    AIC=AIC, BIC=BIC, R2=R2, LogLik=loglikelihood, VarCov=VCM)
  }

  return(results)
}



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%% CODE FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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



# ================================================================
# ================================================================
# ================================================================


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



# ================================================================
# ================================================================
# ================================================================


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%% Fuction EWS: seuil NSR %%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Input:
# 	- Var_Proba: vector of calculated probabilities
# 	- Dicho_Y: Binary variable
#	  - cutoff_interval
# Output:
# threshold with NSR method

EWS_NSR_Criterion <- function(Var_Proba, Dicho_Y, cutoff_interval){

  nb_period <- length(Var_Proba)

  # initialization Matrix results
  Matrix_results <- matrix(data=0, nrow= (floor(1/cutoff_interval)+1) , ncol= 2)
  counter_results <- 1

  # test of each threshold with optim_interval parameter
  for (optim_target in seq(0 , 1 , by = cutoff_interval ))
  {
    # initialization counter for matching and missmatching binary variable
    counter_matching <- 0
    counter_missmatching <- 0

    # loop to compute matching and missmatching binary variable
    for (counter_interval in 1:nb_period){
      if (Var_Proba[counter_interval] >= optim_target){
        dummy_dicho <- 1
      }else{
        dummy_dicho <- 0
      }
      if ((Dicho_Y[counter_interval] == 1) && (dummy_dicho == 1)){
        counter_matching <- counter_matching + 1
      }else if ((Dicho_Y[counter_interval] == 0) && (dummy_dicho == 1)){
        counter_missmatching <- counter_missmatching + 1
      }else{
        counter_matching <- counter_matching
      }
    }

    # recovery of results
    if (counter_matching != 0){
      ratio_tampon <- counter_missmatching /counter_matching
      Matrix_results[counter_results, 1] <- optim_target
      Matrix_results[counter_results, 2] <- ratio_tampon
    }else{
      Matrix_results[counter_results, 1] <- NA
      Matrix_results[counter_results, 2] <- NA
    }
    counter_results <- counter_results + 1
  }

  optim_cutoff_NSR <- Matrix_results[order(Matrix_results[,2]),][1,1]

  return(optim_cutoff_NSR)
}



# ================================================================
# ================================================================
# ================================================================


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%% Function EWS: AM Criteria %%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Input:
# 	- Var_Proba: vector of calculated probabilities
# 	- Dicho_Y: Binary variable
#	  - cutoff_interval
# Output:
# threshold with AM method

EWS_AM_Criterion <- function(Var_Proba, Dicho_Y, cutoff_interval){

  nb_period <- length(Dicho_Y)

  # number of 1 in Dicho_Y
  counter_1 <- sum(Dicho_Y)

  # number of 0 in Dicho_Y
  counter_0 <- nb_period - counter_1

  # initialization Matrix results
  Matrix_results <- matrix(data=0, nrow= (floor(1/cutoff_interval)+1) , ncol= 2)
  counter_results <- 1

  # Loop in order to test each threshold with optim_interval parameter
  for (optim_target in seq(0,1,cutoff_interval) ){

    # Initilization of counter to calculate sensibility and specificity
    counter_sensibility <- 0
    counter_specificity <- 0

    # Calculate sensibility and specificity counters
    for (counter_interval in 1:nb_period){
      if ( (Var_Proba[counter_interval] >= optim_target) && (Dicho_Y[counter_interval] == 1) )
      {
        counter_sensibility <- counter_sensibility + 1
      }else if ( (Var_Proba[counter_interval] < optim_target) && (Dicho_Y[counter_interval] == 0) ){
        counter_specificity <- counter_specificity + 1
      }else{
      }
    }

    # Calculate sensibility and specificity
    Sensibility <- counter_sensibility / counter_1
    Specificity <- counter_specificity / counter_0


    # recovery of results
    ratio_tampon <- abs(Sensibility + Specificity - 1)
    Matrix_results[counter_results, 1] <- optim_target
    Matrix_results[counter_results, 2] <- ratio_tampon
    counter_results <- counter_results + 1
  }

  optim_cutoff_AM <- Matrix_results[order(Matrix_results[,2], decreasing = TRUE),][1,1]

  return(optim_cutoff_AM)
}



# ================================================================
# ================================================================
# ================================================================


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%% Function EWS: seuil CSA %%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Input:
# 	- Var_Proba: vector of calculated probabilities
# 	- Dicho_Y: Binary variable
#	  - cutoff_interval
# Output:
# threshold with CSA method

EWS_CSA_Criterion <- function(Var_Proba, Dicho_Y , cutoff_interval){

  nb_period <- length(Dicho_Y)

  # number of 1 in Dicho_Y
  counter_1 <- sum(Dicho_Y)

  # number of 0 in Dicho_Y
  counter_0 <- nb_period - counter_1

  # initialization Matrix results
  Matrix_results <- matrix(data=0, nrow= (floor(1/cutoff_interval)+1) , ncol= 2)
  counter_results <- 1

  # Initilization of counter to calculate sensibility and specificity
  for (optim_target in seq(0,1,cutoff_interval) )
  {
    # Initilization of sensibility and specificity counters
    counter_sensibility <- 0
    counter_specificity <- 0

    # Calculate sensibility and specificity counters
    for (counter_interval in 1:nb_period){
      if ( (Var_Proba[counter_interval] >= optim_target) && (Dicho_Y[counter_interval]==1) ){
        counter_sensibility <- counter_sensibility + 1
      }else if ((Var_Proba[counter_interval] < optim_target) && (Dicho_Y[counter_interval]==0) ){
        counter_specificity <- counter_specificity + 1
      }else{
      }
    }

    # Calculate sensibility and specifity
    Sensibility <- counter_sensibility / counter_1
    Specificity <- counter_specificity / counter_0

    # recovery of results
    ratio_tampon <- abs(Sensibility - Specificity)
    Matrix_results[counter_results, 1] <- optim_target
    Matrix_results[counter_results, 2] <- ratio_tampon
    counter_results <- counter_results + 1
  }

  optim_cutoff_CSA <- Matrix_results[order(Matrix_results[,2]),][1,1]

  return(optim_cutoff_CSA)
}



# ================================================================
# ================================================================
# ================================================================


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%% Vector error %%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Input:
# 	- Dicho_Y: Binary variable
# 	- Exp_X: Matrix of explanatory variables
#	  - Intercept: Boolean variable that is equal to TRUE (with intercept) or FALSE (without)
# 	- Nb_Id: number of individuals
# 	- Lag: number of lag for the logistic estimation
# 	- type_model:
#			1: static
#			2: dynamic with the lag binary variable
#			3: dynamiC with the lag index variable
#			4: dynamiC with both lag binary and lag index variable
# Output:
# Vector of estimation errors


Vector_Error <- function(Dicho_Y, Exp_X, Intercept, Nb_Id, Lag, type_model)
{
  # Estimation
  logistic_results <- Logistic_Estimation(Dicho_Y, Exp_X, Intercept, Nb_Id, Lag, type_model)

  # Error Calculation
  Results <- as.vector(rep(0,length(logistic_results$prob)))			# Initialization
  Results <- Dicho_Y[1:length(logistic_results$prob)] - logistic_results$prob 		# Estimation

  return(Results)
}



# ================================================================
# ================================================================
# ================================================================


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%% Logistic Estimation Function - 4 models %%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Input:
# 	- Dicho_Y: Binary variable
# 	- Exp_X: Matrix of explanatory variables
#	  - Intercept: Boolean variable that is equal to TRUE (with intercept) or FALSE (without)
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
                    AIC=AIC, BIC=BIC, R2=R2, index=ind[1:T_size], prob = prob, LogLik=loglikelihood, VarCov=VCM)

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
                    AIC=AIC, BIC=BIC, R2=R2, index=ind[1:T_size], prob = prob, LogLik=loglikelihood, VarCov=VCM)

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
                    AIC=AIC, BIC=BIC, R2=R2, index=ind[1:T_size], prob = prob, LogLik=loglikelihood, VarCov=VCM)


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
                    AIC=AIC, BIC=BIC, R2=R2, index=ind[1:T_size], prob = prob, LogLik=loglikelihood, VarCov=VCM)

  }

  return(results)
}



# ================================================================
# ================================================================
# ================================================================


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%% Function BlockBootstrapp %%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Input:
# 	- Dicho_Y: Binary variable
# 	- Exp_X: Matrix of explanatory variables
#	  - Intercept: Boolean variable that is equal to TRUE (with intercept) or FALSE (without)
#	  - n_simul: number of simulations for the block bootstrapp
# Output:
# Matrix with bootstrapp series

BlockBootstrapp <- function(Dicho_Y, Exp_X, Intercept, n_simul)
{

  # optimal block size
  size_block <- floor(length(Dicho_Y)^(1/5))

  if (length(Dicho_Y)==length(Exp_X) && Intercept==TRUE){
    Exp_X <- matrix(data=c(rep(1,length(Dicho_Y)), Exp_X), ncol=2, nrow=length(Dicho_Y))
  }else if (length(Dicho_Y)!=length(Exp_X) && Intercept==TRUE){
    Exp_X <- matrix(data=c(rep(1,length(Dicho_Y)), Exp_X), ncol=(1+ncol(Exp_X)), nrow=length(Dicho_Y))
  }

  # number of colomns for simulation matrix
  if (length(Dicho_Y)==length(Exp_X))
  {
    nvalue <- 2 			# 1 explanatory variable + 1 binary variable
  }else
  {
    nvalue <- ncol(Exp_X)+1        # n explanatory variables + 1 binary variable
  }

  # Initialization of matrix results
  matrix_results <- matrix(data=0, ncol= (nvalue*n_simul), nrow=length(Dicho_Y))

  for (compteur_simul in 1:n_simul)
  {
    # block position
    block_position <- sample(1:(length(Dicho_Y)-size_block+1),1)

    # block recovery
    block_value <- matrix(data=0, ncol= nvalue, nrow= size_block )# initialisation de la taille du block

    if (length(Dicho_Y)==length(Exp_X))
    {
      block_value[,1] <- Dicho_Y[block_position:(block_position+size_block-1)]
      block_value[,2] <- Exp_X[block_position:(block_position+size_block-1)]
    } else {
      block_value[,1] <- Dicho_Y[block_position:(block_position+size_block-1)]
      block_value[,2:nvalue] <- Exp_X[block_position:(block_position+size_block-1),]
    }

    # Recovery of results
    if (length(Dicho_Y)==length(Exp_X))
    {
      matrix_results[1:(length(Dicho_Y)-size_block),(1+(compteur_simul-1)*nvalue)] <- Dicho_Y[(size_block+1):(length(Dicho_Y))]
      matrix_results[1:(length(Dicho_Y)-size_block),(compteur_simul*nvalue)] <- Exp_X[(size_block+1):(length(Dicho_Y))]
      matrix_results[(length(Dicho_Y)-size_block+1):length(Dicho_Y),(1+(compteur_simul-1)*nvalue):(compteur_simul*nvalue)] <- block_value
    } else {
      matrix_results[1:(length(Dicho_Y)-size_block),(1+(compteur_simul-1)*nvalue)] <- Dicho_Y[(size_block+1):(length(Dicho_Y))]
      matrix_results[1:(length(Dicho_Y)-size_block),(2+(compteur_simul-1)*nvalue):(compteur_simul*nvalue)] <- Exp_X[(size_block+1):(length(Dicho_Y)),]
      matrix_results[(length(Dicho_Y)-size_block+1):length(Dicho_Y),(1+(compteur_simul-1)*nvalue):(compteur_simul*nvalue)] <- block_value
    }
  }
  return(matrix_results)
}



# ================================================================
# ================================================================
# ================================================================


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%% GIRF function %%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Input:
# 	- Dicho_Y: Binary variable
# 	- Exp_X: Matrix of explanatory variables
#	  - Lag: number of lag to calculate the GIRF
#	  - Int: Boolean variable that is equal to TRUE (with intercept) or FALSE (without)
# 	- t_mod:
#			1: static
#			2: dynamic with the lag binary variable
#			3: dynamiC with the lag index variable
#			4: dynamiC with both lag binary and lag index variable
#	  - horizon: horizon target for the GIRF analysis
#	  - shock_size: size of the shock
#	  - OC: threshold to determine the value of the dichotomous variable as a function of the index level
# Output:
# Matrix with:
# 			- column 1: horizon
# 			- column 2: index
# 			- column 3: index with shock
# 			- column 4: probability associated to the index
# 			- column 5: probability associated to the index with shock
# 			- column 6: binary variable associated to the index
# 			- column 7: binary variable associated to the index with shock
GIRF_Dicho <- function(Dicho_Y, Exp_X, Lag, Int, t_mod, horizon, shock_size, OC)
  {

  # Initialization
  matrix_results <- matrix(data=0,nrow=(horizon+1), ncol=7)

  for (compteur_horizon in 0:horizon)
  {
    if (compteur_horizon==0)
    {
      # index at time 0
      results_forecast <- Logistic_Estimation(Dicho_Y, Exp_X, Int, 1, Lag, t_mod)

      # estimated coefficients
      Coeff_estimated <- as.vector(rep(0,length(results_forecast$Estimation[,2]))) # Initialization
      Coeff_estimated <- results_forecast$Estimation[,2]

      # last explanatory variables
      Last_Exp_X <- as.vector(rep(0,length(Coeff_estimated)))

      # Data Processing
      if (Int==TRUE && length(Dicho_Y)==length(Exp_X)){
        Exp_X <- matrix(data=c(rep(1,length(Dicho_Y)),Exp_X), nrow=length(Dicho_Y) , ncol=2 )
      } else if (Int==TRUE && length(Dicho_Y)!=length(Exp_X)){
        Exp_X <- matrix(data=c(rep(1,length(Dicho_Y)),Exp_X), nrow=length(Dicho_Y) , ncol=(1+ncol(Exp_X)) )
      }

      if (t_mod==1){
        # model 1: Exp_X
        if (length(Dicho_Y)==length(Exp_X))
        {
          Last_Exp_X <- c(Exp_X[(length(Dicho_Y)-1)])
        } else {
          Last_Exp_X <- c(Exp_X[(length(Dicho_Y)-1),])
        }
      } else if (t_mod==2) {
        # model 2: Exp_X + Dicho_Y
        if (length(Dicho_Y)==length(Exp_X))
        {
          Last_Exp_X <- c(Exp_X[(length(Dicho_Y)-1)], Dicho_Y[length(Dicho_Y)])
        } else {
          Last_Exp_X <- c(Exp_X[(length(Dicho_Y)-1),], Dicho_Y[length(Dicho_Y)])
        }
      } else if (t_mod==3) {
        # model 3: Exp_X + Index
        Last_Index <- results_forecast$index[length(results_forecast$index)]	# Recover last index
        if (length(Dicho_Y)==length(Exp_X))
        {
          Last_Exp_X <- c(Exp_X[(length(Dicho_Y)-1)], Last_Index)
        } else {
          Last_Exp_X <- c(Exp_X[(length(Dicho_Y)-1),], Last_Index)
        }
      } else if (t_mod==4) {
        # model 4: Exp_X + Index + Dicho_Y
        Last_Index <- results_forecast$index[length(results_forecast$index)]	# Recover last index
        if (length(Dicho_Y)==length(Exp_X))
        {
          Last_Exp_X <- c(Exp_X[(length(Dicho_Y)-1)], Dicho_Y[length(Dicho_Y)], Last_Index)
        } else {
          Last_Exp_X <- c(Exp_X[(length(Dicho_Y)-1),], Dicho_Y[length(Dicho_Y)], Last_Index)
        }
      }

      # Recovery of the index and index with shock
      index <- results_forecast$index[length(results_forecast$index)]
      proba <- exp(index)/(1+exp(index))
      index_shock <- index + shock_size
      proba_shock <- exp(index_shock)/(1+exp(index_shock))

      # Calculation of the dichotomous variable as a function of the threshold and the index
      if (proba > OC){
        Var_Dicho <- 1
      }else{
        Var_Dicho <- 0
      }

      # Calculation of the dichotomous variable as a function of the threshold and the index with shock
      if (proba_shock > OC){
        Var_Dicho_shock <- 1
      }else{
        Var_Dicho_shock <- 0
      }
    } else {

      # horizon > 0

      if (t_mod==1)
      {

        # Index recovery and associated probability
        index <- results_forecast$index[length(results_forecast$index)]
        proba <- exp(index)/(1+exp(index))

        # Index with shock recovery and associated probability
        index_shock <- results_forecast$index[length(results_forecast$index)]
        proba_shock <- exp(index_shock)/(1+exp(index_shock))

        # Calculation of the dichotomous variable as a function of the threshold and the index
        if (proba > OC){
          Var_Dicho <- 1
        }else{
          Var_Dicho <- 0
        }

        # Calculation of the dichotomous variable as a function of the threshold and the index with shock
        if (proba_shock > OC){
          Var_Dicho_shock <- 1
        }else{
          Var_Dicho_shock <- 0
        }

      } else if (t_mod==2)
      {

        # Calculation of the new index
        Last_Exp_X <- c(Last_Exp_X[1:(length(Last_Exp_X)-1)],Var_Dicho)	# last values of X and last binary
        index <- Coeff_estimated %*% Last_Exp_X
        proba <- exp(index)/(1+exp(index))

        # Calculation of the new index with shock
        Last_Exp_X <- c(Last_Exp_X[1:(length(Last_Exp_X)-1)],Var_Dicho_shock)	# last values of X and last binary with shock
        index_shock <- Coeff_estimated %*% Last_Exp_X
        proba_shock <- exp(index_shock)/(1+exp(index_shock))

        # Calculation of the dichotomous variable as a function of the threshold and the index
        if (proba > OC){
          Var_Dicho <- 1
        }else{
          Var_Dicho <- 0
        }

        # Calculation of the dichotomous variable as a function of the threshold and the index with shock
        if (proba_shock > OC){
          Var_Dicho_shock <- 1
        }else{
          Var_Dicho_shock <- 0
        }


      } else if (t_mod==3)
      {

        # Calculation of the new index
        Last_Exp_X <- c(Last_Exp_X[1:(length(Last_Exp_X)-1)],index)	# last values of X and last index
        index <- Coeff_estimated %*% Last_Exp_X
        proba <- exp(index)/(1+exp(index))

        # Calculation of the new index with shock
        Last_Exp_X <- c(Last_Exp_X[1:(length(Last_Exp_X)-1)],index_shock)	# last values of X and last index with shock
        index_shock <- Coeff_estimated %*% Last_Exp_X
        proba_shock <- exp(index_shock)/(1+exp(index_shock))

        # Calculation of the dichotomous variable as a function of the threshold and the index
        if (proba > OC){
          Var_Dicho <- 1
        }else{
          Var_Dicho <- 0
        }

        # Calculation of the dichotomous variable as a function of the threshold and the index with shock
        if (proba_shock > OC){
          Var_Dicho_shock <- 1
        }else{
          Var_Dicho_shock <- 0
        }


      } else if (t_mod==4)
      {

        # Calculation of the new index
        Last_Exp_X <- c(Last_Exp_X[1:(length(Last_Exp_X)-2)],Var_Dicho, index)	# last values of X, last binary and last index
        index <- Coeff_estimated %*% Last_Exp_X
        proba <- exp(index)/(1+exp(index))

        # Calculation of the new index with shock
        Last_Exp_X <- c(Last_Exp_X[1:(length(Last_Exp_X)-2)], Var_Dicho_shock, index_shock)	# last values of X, last binary with shock and last index with shock
        index_shock <- Coeff_estimated %*% Last_Exp_X
        proba_shock <- exp(index_shock)/(1+exp(index_shock))

        # Calculation of the dichotomous variable as a function of the threshold and the index
        if (proba > OC){
          Var_Dicho <- 1
        }else{
          Var_Dicho <- 0
        }

        # Calculation of the dichotomous variable as a function of the threshold and the index with shock
        if (proba_shock > OC){
          Var_Dicho_shock <- 1
        }else{
          Var_Dicho_shock <- 0
        }

      }

    }

    matrix_results[compteur_horizon + 1,1] <- compteur_horizon
    matrix_results[compteur_horizon + 1,2] <- index
    matrix_results[compteur_horizon + 1,3] <- index_shock
    matrix_results[compteur_horizon + 1,4] <- proba
    matrix_results[compteur_horizon + 1,5] <- proba_shock
    matrix_results[compteur_horizon + 1,6] <- Var_Dicho
    matrix_results[compteur_horizon + 1,7] <- Var_Dicho_shock
  }

  return(matrix_results)
}



# ================================================================
# ================================================================
# ================================================================


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%% GIRF pour Intervale de confiance %%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Input:
# 	- Dicho_Y: Binary variable
# 	- Exp_X: Matrix of explanatory variables
#	  - Int: Boolean variable that is equal to TRUE (with intercept) or FALSE (without)
# 	- Lag: number of lag for the logistic estimation
# 	- t_mod:
#			1: static
#			2: dynamic with the lag binary variable
#			3: dynamiC with the lag index variable
#			4: dynamiC with both lag binary and lag index variable
#	  - n_simul: number of simulations for the bootstrap
#	  - centile_shock: percentile of the shock from the estimation errors
#	  - horizon: horizon
#	  - OC: either a value or the name of the optimal cut-off / threshold ("NSR", "CSA", "AM")
#
# Output:
# Matrix where for each simulation there are 7 columns with:
# 			- column 1: horizon
# 			- column 2: index
# 			- column 3: index with shock
# 			- column 6: probability associated to the index
# 			- column 7: probability associated to the index with shock
# 			- column 6: binary variable associated to the index
# 			- column 7: binary variable associated to the index with shock

Simul_GIRF <- function(Dicho_Y, Exp_X, Int, Lag, t_mod, n_simul, centile_shock, horizon, OC)
{

  # Initialization of the results matrix
  matrix_results <- matrix(data=0,ncol=(7*n_simul), nrow=(horizon+1))

  # number of values for the bootstrapp and results
  if (length(Dicho_Y)==length(Exp_X) && Int==FALSE)
  {
    nvalue <- 2 			# 1 explanatory variable + 1 binary variable
  }else if (length(Dicho_Y)==length(Exp_X) && Int==TRUE){
    nvalue <- 3 	    # 1 explanatory variable + 1 binary variable + 1 Intercept
  }else if (length(Dicho_Y)!=length(Exp_X) && Int==FALSE){
    nvalue <- ncol(Exp_X)+1 	# n explanatory variables + 1 binary variable
  }else{
    nvalue <- ncol(Exp_X)+2 	# n explanatory variables + 1 binary variable + 1 Intercept
  }

  # Block Bootstrap estimation
  matrice_bootstrap <- matrix(data=0,nrow=length(Dicho_Y),ncol= n_simul*nvalue)
  matrice_bootstrap <- BlockBootstrapp(Dicho_Y, Exp_X, Int, n_simul)

  # Estimation of coefficients and errors for each simulation
  for (compteur_simul in 1:n_simul)
  {
    # Vector of binary variable
    Dicho_Y_bootstrap <- matrix(data=0,ncol=1,nrow=length(Dicho_Y))
    Dicho_Y_bootstrap <- matrice_bootstrap[,1+(compteur_simul-1)*nvalue]

    # Matrix of explanatory variables
    Exp_X_bootstrap <- matrix(data=0,ncol=(nvalue-1),nrow=length(Dicho_Y))
    Exp_X_bootstrap <- matrice_bootstrap[,(2+(compteur_simul-1)*nvalue):(compteur_simul*nvalue)]

    if (OC=="NSR") {

      # Regression
      Results_serie_bootstrap <- Logistic_Estimation(Dicho_Y_bootstrap, Exp_X_bootstrap, FALSE, 1, Lag, t_mod)

      # Recovery of estimated probabilities
      vecteur_proba_bootstrap <- as.vector(rep(0,length(Dicho_Y_bootstrap)-Lag))
      vecteur_proba_bootstrap <- Results_serie_bootstrap$prob

      # Recovery of the binary variable
      vecteur_binary_bootstrap <- as.vector(rep(0,length(Dicho_Y_bootstrap)-Lag))
      vecteur_binary_bootstrap <- Dicho_Y_bootstrap[(1+Lag):length(Dicho_Y_bootstrap)]

      pas_optim <- 0.0001 # intervals for OC (optimal cut-off) estimation

      # NSR calculation
      threshold_estimated <- EWS_NSR_Criterion(vecteur_proba_bootstrap, vecteur_binary_bootstrap, pas_optim)

    } else if (OC=="CSA") {

      # Regression
      Results_serie_bootstrap <- Logistic_Estimation(Dicho_Y_bootstrap, Exp_X_bootstrap, FALSE, 1, Lag, t_mod)

      # Recovery of estimated probabilities
      vecteur_proba_bootstrap <- as.vector(rep(0,length(Dicho_Y_bootstrap)-Lag))
      vecteur_proba_bootstrap <- Results_serie_bootstrap$prob

      # Recovery of the binary variable
      vecteur_binary_bootstrap <- as.vector(rep(0,length(Dicho_Y_bootstrap)-Lag))
      vecteur_binary_bootstrap <- Dicho_Y_bootstrap[(1+Lag):length(Dicho_Y_bootstrap)]

      pas_optim <- 0.0001 # intervals for OC (optimal cut-off) estimation

      # CSA calculation
      threshold_estimated <- EWS_CSA_Criterion(vecteur_proba_bootstrap, vecteur_binary_bootstrap, pas_optim)

    } else if (OC=="AM") {

      # Regression
      Results_serie_bootstrap <- Logistic_Estimation(Dicho_Y_bootstrap, Exp_X_bootstrap, FALSE, 1, Lag, t_mod)

      # Recovery of estimated probabilities
      vecteur_proba_bootstrap <- as.vector(rep(0,length(Dicho_Y_bootstrap)-Lag))
      vecteur_proba_bootstrap <- Results_serie_bootstrap$prob

      # Recovery of the binary variable
      vecteur_binary_bootstrap <- as.vector(rep(0,length(Dicho_Y_bootstrap)-Lag))
      vecteur_binary_bootstrap <- Dicho_Y_bootstrap[(1+Lag):length(Dicho_Y_bootstrap)]

      pas_optim <- 0.0001 # intervals for OC (optimal cut-off) estimation

      # AM calculation
      threshold_estimated <- EWS_AM_Criterion(vecteur_proba_bootstrap, vecteur_binary_bootstrap, pas_optim)

    } else {
    }

    # Error Calculation
    Residuals_bootstrap <- Vector_Error(Dicho_Y_bootstrap, Exp_X_bootstrap, FALSE, 1, Lag, t_mod)

    # Recovery of the shock in the estimation error vector
    size_shock_bootstrap <- quantile(Residuals_bootstrap, centile_shock)

    # Calculation of the response function
    matrix_results[,(1+(7*(compteur_simul-1))):(7*compteur_simul)] <- GIRF_Dicho(Dicho_Y_bootstrap, Exp_X_bootstrap, Lag, FALSE, t_mod, horizon, size_shock_bootstrap, threshold_estimated)

  }

  return(matrix_results)
}



# ================================================================
# ================================================================
# ================================================================


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%% GIRF Index IC %%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Input:
#   - results_simulation_GIRF: matrix output of the Simulation_GIRF function
# 	- CI_bounds: size of the confidence intervals
#   - n_simul: number of simulations
# Output:
# List with:
#       - Simulation_CI: the index values that belong in the CI for each horizon
#       - values_CI: Index values with lower bound, average index, and upper bound for each horizon

GIRF_Index_CI <- function(results_simul_GIRF, CI_bounds, n_simul, horizon_forecast){

  # Index recovery
  storage_index <- matrix(data=0, nrow=(horizon_forecast+1), ncol=n_simul)
  for (counter_simul in 1:n_simul)
  {
    storage_index[,counter_simul] <- results_simul_GIRF[,(3+(7*(counter_simul-1)))] - results_simul_GIRF[,(2+(7*(counter_simul-1)))]
  }

  # store index for each horizon
  for (counter_forecast in 1:(horizon_forecast+1))
  {
    storage_index[counter_forecast,] <-sort(storage_index[counter_forecast,])
  }

  # Remove values outside of CI
  simul_inf <- ceiling( ((1-CI_bounds)/2) * n_simul ) 	# simulation number of the lower bound
  simul_sup <- n_simul - simul_inf + 1 	# simulation number of the upper bound
  result_CI <- matrix(data=0,nrow=(horizon_forecast+1), ncol=(simul_sup - simul_inf + 1)) 	# Initialization
  result_CI <- storage_index[,simul_inf:simul_sup]

  # Index average for each horizon
  mean_result <- as.vector(horizon_forecast+1)
  for (compteur in 1:(horizon_forecast+1))
  {
    mean_result[compteur] <- mean(storage_index[compteur,simul_inf:simul_sup])
  }

  # Matrix with lower bound, average, and upper bound
  result_Graph <- matrix(data=0,nrow=(horizon_forecast+1), ncol=3)
  result_Graph[,1] <- storage_index[,simul_inf]
  result_Graph[,2] <- mean_result
  result_Graph[,3] <- storage_index[,simul_sup]

  results <- list(Simulation_CI=result_CI, values_CI=result_Graph)

  return(results)
}


# ================================================================
# ================================================================
# ================================================================


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%% GIRF Index IC %%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Input:
#   - results_simulation_GIRF: matrix output of the Simulation_GIRF function
# 	- CI_bounds: size of the confidence intervals
#   - n_simul: number of simulations
# Output:
# List with:
#       - horizon
#       - Simulation_CI_proba_shock: the proba_shock values that belong in the CI for each horizon
#       - Simulation_CI_proba: the proba values that belong in the CI for each horizon
#       - CI_proba_shock= proba_shock values with lower bound, average index, and upper bound for each horizon,
#       - CI_proba= proba values with lower bound, average index, and upper bound for each horizon,

GIRF_Proba_CI <- function( results_simul_GIRF, CI_bounds, n_simul, horizon_forecast){


  # Proba recovery
  storage_proba_shock <- matrix(data=0, nrow=(horizon_forecast+1), ncol=n_simul, horizon_forecast)
  storage_proba <- matrix(data=0, nrow=(horizon_forecast+1), ncol=n_simul)
  for (counter_simul in 1:n_simul)
  {
    storage_proba_shock[,counter_simul] <- results_simul_GIRF[,(5+(7*(counter_simul-1)))]
    storage_proba[,counter_simul] <- results_simul_GIRF[,(4+(7*(counter_simul-1)))]
  }

  # store proba for each horizon
  for (counter_horizon in 1:(horizon_forecast+1))
  {
    storage_proba_shock[counter_horizon,] <-sort(storage_proba_shock[counter_horizon,])
    storage_proba[counter_horizon,] <-sort(storage_proba[counter_horizon,])
  }

  # Remove values outside of CI
  simul_inf <- ceiling( ((1-CI_bounds)/2) * n_simul ) 	# simulation number of the lower bound
  simul_sup <- n_simul - simul_inf + 1 	# simulation number of the upper bound
  result_CI_proba_shock <- matrix(data=0,nrow=(horizon_forecast+1), ncol=(simul_sup - simul_inf + 1)) 	# Initialization
  result_CI_proba_shock <- storage_proba_shock[,simul_inf:simul_sup]
  result_CI_proba <- matrix(data=0,nrow=(horizon_forecast+1), ncol=(simul_sup - simul_inf + 1)) 	# Initialization
  result_CI_proba <-  storage_proba[,simul_inf:simul_sup]

  # Averages proba for each horizon
  mean_result_proba_shock <- as.vector(horizon_forecast+1)
  mean_result_proba <- as.vector(horizon_forecast+1)
  for (counter_horizon in 1:(horizon_forecast+1))
  {
    mean_result_proba_shock[counter_horizon] <- mean(result_CI_proba_shock[counter_horizon,])
    mean_result_proba[counter_horizon] <- mean(result_CI_proba[counter_horizon,])
  }

  # Matrix with lower bound, average, and upper bound for proba shock and proba without shock
  result_proba_shock <- matrix(data=0,nrow=(horizon_forecast+1), ncol=3)
  result_proba_shock[,1] <- storage_proba_shock[,simul_inf]
  result_proba_shock[,2] <-mean_result_proba_shock
  result_proba_shock[,3] <- storage_proba_shock[,simul_sup]
  result_proba <- matrix(data=0,nrow=(horizon_forecast+1), ncol=3)
  result_proba[,1] <- storage_proba[,simul_inf]
  result_proba[,2] <- mean_result_proba
  result_proba[,3] <- storage_proba[,simul_sup]

  horizon_vect <- as.vector(rep(0,horizon_forecast+1))
  horizon_vect <- c(0:horizon_forecast)

  results <- list(horizon = horizon_vect, Simulation_CI_proba_shock=result_CI_proba_shock, Simulation_CI_proba=result_CI_proba,
                  CI_proba_shock=result_proba_shock, CI_proba=result_proba)
  return(results)
}


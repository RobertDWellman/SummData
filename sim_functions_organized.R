# Functions to get summary statistics ----
get.summstat.survival <- function(E,Y,X,B,A,prescription.mode=seq(30,trunc,by=30),
                         my.presc.K=1,tie.method="efron",interact=FALSE){
  if (is.null(E)){
    stop("E is required for survival data")
  }

  glm.ctrl <- glm.control(epsilon = 1e-8, maxit = 25, trace = FALSE)
  cox.ctrl <- survival::coxph.control(eps=1e-09, toler.chol=.Machine$double.eps^0.75,
                                      iter.max=20, toler.inf=sqrt(1e-09), outer.max=10)

  ### (Correlated) binary covariates B
  B <- as.matrix(B)
  n.B <- ncol(B)

  #### Marginal mean, i.e. prevalence of each binary covariate
  P <- apply(as.matrix(B),2,function(XX){as.numeric(table(XX)[[2]]/sum(table(XX)))})
  #### Common probability, i.e. number of one's shared by each pair of binary variables divided by sample size
  Common.P <- t(as.matrix(B))%*%as.matrix(B)/nrow(B)

  ### Categorical variables A
  #### Intercept and coefficients from multinomial regression for categorical variables, fitted on all binary covariates B.
  #### Number of categorical
  A <- as.matrix(A)
  n.A <- ncol(A)

  ## Create list containing distributions of levels for
  ## each categorical variable
  A.dist <- apply(A,2,function(XX){as.numeric(table(XX)/sum(table(XX)))})

  ### P.orddata and Corr.orddata
  P.ord = append(lapply(P,FUN=function(x){c(1-x,x)}),
                 A.dist)

  if (is.numeric(A)) {
    Corr.ord = cor(cbind(B,A))
  } else {
    Corr.ord = cor(cbind(B,as.factor(A)))
  }

  ### Categorical variables A
  #### Intercept and coefficients from multinomial regression for categorical variables, fitted on all binary covariates B.
  #### Number of categorical
  #### A -> model.matrix type "A.indicator"
  # 9.24: catlabs may be redundant
  if(n.A>0)
  {
    A.indicator <- NULL
    catlabs <- NULL
    for(i in 1:n.A)
    {
      dummy <- NULL
      levels <- sort(unique(A[,i]))
      for(level in levels[-1]) {
        catlabs <- c(catlabs, paste0('cat',i,'_','level',level))
        dummy <- cbind(dummy, A[,i]==level)
      }
      A.indicator <- cbind(A.indicator, dummy)
    }
    colnames(A.indicator) <- catlabs
  }


  ####################################
  #### Chain of coef for B
  coef.chain <- vector("list", n.B)
  names(coef.chain) <- colnames(B)
  coef.chain[[1]] <- mean(B[,1])
  for (i in 2:n.B) {
    coef.chain[[i]] <- tryCatch({
      coef(glm(B[,i] ~ B[, 1:(i-1)], family = "binomial"))
    }, error = function(e) {
      message("Chain regression failed for B[", i, "]: ", e$message)
      NA
    })
  }

  #### Coef for relationship between binary and categorical confounders

  coef.AonB <- vector("list", n.A)
  names(coef.AonB) <- colnames(A)

  # For the first categorical variable, regress on all B variables
  coef.AonB[[1]] <- tryCatch({
    coef(nnet::multinom(A[, 1] ~ B, trace = FALSE))
  }, error = function(e) {
    message("Model fitting failed for A[1]: ", e$message)
    NA
  })

  # For subsequent categorical variables, regress on the indicators of previous categorical variables and B
  if (n.A > 1) {
    for (i in 2:n.A) {
      # Get the corresponding indicator columns for all previous categorical variables
      start_index <- 1
      end_index <- sum(sapply(1:(i - 1), function(k) length(unique(A[, k])) - 1))

      # Extract the relevant indicator columns
      previous_A_indicators <- A.indicator[, start_index:end_index, drop = FALSE]

      # Create a formula using previous indicators and B as predictors
      formula <- as.formula(paste("A[, i] ~ . + B"))  # Use . to represent all predictors

      # Fit the multinomial logistic regression model, catch errors due to mismatched rows or other issues
      coef.AonB[[i]] <- tryCatch({
        model_data <- data.frame(A = A[, i], previous_A_indicators, B)
        model <- nnet::multinom(A ~ ., data = model_data, trace = FALSE)
        coef(model)
      }, warning = function(w) {
        message("Warning during model fitting for A[", i, "]: ", w$message)
        NA
      }, error = function(e) {
        message("Model fitting failed for A[", i, "]: ", e$message)
        NA
      })
    }
  }


  ####################################

  ### Exposure variable X|B,A (Propensity Score Model Coefficients)
  #### Intercept and coefficients from logistic regression
  ps.fit <- glm.fit(cbind(1,B,A.indicator), X, control=glm.ctrl, family=binomial())
  class(ps.fit) <- 'glm'
  coef.XonZ <- coef(ps.fit)
  ps.by.x <- cbind(X, fitted(ps.fit))
  ps.vcov <- vcov(ps.fit)
  probs <- predict(ps.fit, type = "response")
  ps.c <- mean(sample(probs[X == 1L], 1000000L, TRUE) > sample(probs[X == 0L], 1000000L, TRUE))


  ####  Censoring distribution -- Weibull shape and scale parameters estimated through parametric survival regression:
  #### Note: use only censoring data to estimate proportion among the censored.
  #### Unique censoring distribution shape
  #### Note: If one decide to add extra distribution properties, (1) and (2) should fit model using data eliminating the unique observations,
  #### e.g. observations with E <- the top two most frequent presc pattern bumps

  ## Proportion with a prescription bump amongst those censored
  P.presc <- NULL
  # P.presc <- sapply(prescription.mode, function(pm) mean(E[Y==0] == pm)) # 9.24: changed below

  for( i in 1:length(prescription.mode))
  {
    P.presc<-cbind(P.presc,mean(E[Y==0]==prescription.mode[i]))
  }

  P.presc <- data.frame(P.presc)
  names(P.presc) <- as.character(prescription.mode)

  #### CENSORING DISTRIBUTION ####

  # Covariate matrix
  if (interact) {
    interact.X.covs <- cbind(B,A.indicator)*X
    colnames(interact.X.covs) <- paste("X",colnames(interact.X.covs),sep="_")
    Z <- cbind(X, cbind(B,A.indicator,interact.X.covs))
  } else {
    Z <- cbind(X, B, A.indicator)
  }


  # Simple Censoring No Bumps
  simple.cens <- survival::survreg(survival::Surv(E, 1-Y) ~ 1, dist="weib")
  simple.coef.cens <- coef(simple.cens)
  simple.scale.cens <- simple.cens$scale
  simple.vcov.cens <- vcov(simple.cens)
  names(simple.coef.cens) <- c("Intercept")

  # Covariate Censoring no Bumps
  rownames(Z) <- paste(1:dim(Z)[[1]])
  cov.cens <- survival::survreg(survival::Surv(E, 1-Y) ~ Z, dist="weib")
  cov.coef.cens <- coef(cov.cens)
  cov.scale.cens <- cov.cens$scale
  cov.vcov.cens <- vcov(cov.cens)
  names(cov.coef.cens) <- c("Intercept",colnames(Z))

  #Simple Censoring With Bumps
  # Select only the top my.presc.K bumps
  prescription.mode.topK <- prescription.mode[tail(order(unlist(P.presc)),my.presc.K)]
  P.presc.topK <- P.presc[tail(order(unlist(P.presc)),my.presc.K)]

  ## Bump Censoring Types
  ind<- (!E%in%prescription.mode.topK)

  simplebump.cens <- survival::survreg(survival::Surv(E[ind], 1-Y[ind]) ~ 1, dist="weib")
  simplebump.coef.cens <- coef(simplebump.cens)
  simplebump.scale.cens <- simplebump.cens$scale
  names(simplebump.coef.cens) <- c("Intercept")
  simplebump.vcov.cens <- vcov(simplebump.cens)

  # Covariate Censoring with Bumps
  covbump.cens <- survival::survreg(survival::Surv(E[ind], 1-Y[ind]) ~ Z[ind,], dist="weib")
  covbump.coef.cens <- coef(covbump.cens)
  covbump.scale.cens <- covbump.cens$scale
  names(covbump.coef.cens) <- c("Intercept",colnames(Z))
  covbump.vcov.cens <- vcov(covbump.cens)

  ### EVENT DISTRIBUTION
  #### Weibull shape and scale parameters estimated through parametric survival regression,
  adj.event <- survival::survreg(survival::Surv(E, Y) ~ Z, dist="weib")

  if(sum(is.nan(coef(adj.event)))>0) {
    adj.event <- survival::survreg(survival::Surv(E, Y) ~ Z, dist="exp")
    adj.coef.event <- coef(adj.event)
    adj.scale.event <- 1
  } else {
    adj.coef.event <- coef(adj.event)
    adj.scale.event <- adj.event$scale
  }
  names(adj.coef.event) <- c("Intercept",colnames(Z))
  adj.vcov.event <- vcov(adj.event)

  #### Hazard ratio estimated through Cox PH model, w/adjustment for treatment X and covariates Z
  require(survival)
  # 9.24: add the cox.ctrl
  cox.adjusted <- coxph(Surv(E, Y) ~ . , data=as.data.frame(model.matrix(~Z)[,-1]),
                        method=tie.method,x=TRUE,control=cox.ctrl)
  # class(cox.adjusted) <- "coxph"
  cox.coef.adjusted <- coef(cox.adjusted)
  names(cox.coef.adjusted) <- colnames(Z)
  cox.vcov <- vcov(cox.adjusted)

  ## Cox censoring model for plasmode simulation
  cox.adjusted.cens <-  coxph(Surv(E, 1-Y) ~ . , data=as.data.frame(model.matrix(~Z)[,-1]),
                              method=tie.method,x=TRUE,control=cox.ctrl)
  # class(cox.adjusted.cens) <- "coxph"

  #### Additional information
  control.rate <- sum(Y[X==0])/sum(E[X==0])
  compare.rate <- sum(Y[X==1])/sum(E[X==1])
  control.events <- sum(Y[X==0])
  compare.events <- sum(Y[X==1])
  N.X <- table(X)
  P.time <- c(sum(E[X==0]), sum(E[X==1]))/table(X)

  # return survival outcome version
  return(list(n=length(Y),N.X=N.X, P.time=P.time,
              control.events=control.events,compare.events=compare.events,
              control.rate=control.rate,compare.rate=compare.rate,
              P=P,Common.P=Common.P,P.ord=P.ord,
              Corr.ord=Corr.ord,coef.XonZ=coef.XonZ,Coef.bin=coef.chain, Coef.cat=coef.AonB,
              propensity=ps.by.x,
              propensity.vcov=ps.vcov,P.presc=P.presc, cStat=ps.c,
              prescription.mode=prescription.mode,P.presc.topK=P.presc.topK,
              prescription.mode.topK=prescription.mode.topK,
              simple.coef.cens=simple.coef.cens,
              simple.scale.cens=simple.scale.cens, simple.vcov.cens=simple.vcov.cens,
              simplebump.coef.cens=simplebump.coef.cens,
              simplebump.scale.cens=simplebump.scale.cens,
              simplebump.vcov.cens=simplebump.vcov.cens,
              cov.coef.cens=cov.coef.cens,cov.scale.cens=cov.scale.cens,
              cov.vcov.cens=cov.vcov.cens,
              covbump.coef.cens=covbump.coef.cens,
              covbump.scale.cens=covbump.scale.cens,
              covbump.vcov.cens=covbump.vcov.cens,
              adj.coef.event=adj.coef.event,adj.scale.event=adj.scale.event,
              adj.vcov.event=adj.vcov.event,cox.coef.adjusted=cox.coef.adjusted,cox.vcov=cox.vcov,
              cox.fit.event=cox.adjusted,cox.fit.cens=cox.adjusted.cens
  ))
}

get.summstat.binary <- function(Y,X,B,A){
  glm.ctrl <- glm.control(epsilon = 1e-8, maxit = 25, trace = FALSE)

  ### (Correlated) binary covariates B
  B <- as.matrix(B)
  n.B <- ncol(B)

  #### Marginal mean, i.e. prevalence of each binary covariate
  P <- apply(as.matrix(B),2,function(XX){as.numeric(table(XX)[[2]]/sum(table(XX)))})
  #### Common probability, i.e. number of one's shared by each pair of binary variables divided by sample size
  Common.P <- t(as.matrix(B))%*%as.matrix(B)/nrow(B)

  ## Categorical variables A
  #### Intercept and coefficients from multinomial regression for categorical variables, fitted on all binary covariates B.
  #### Number of categorical
  A <- as.matrix(A)
  n.A <- ncol(A)

  ## Create list containing distributions of levels for
  ## each categorical variable
  A.dist <- apply(A,2,function(XX){as.numeric(table(XX)/sum(table(XX)))})

  ### P.orddata and Corr.orddata
  P.ord = append(lapply(P,FUN=function(x){c(1-x,x)}),
                 A.dist)

  if (is.numeric(A)) {
    Corr.ord = cor(cbind(B,A))
  } else {
    Corr.ord = cor(cbind(B,as.factor(A)))
  }

  ### Categorical variables A
  #### Intercept and coefficients from multinomial regression for categorical variables, fitted on all binary covariates B.
  #### Number of categorical
  #### A -> model.matrix type "A.indicator"
  if(n.A>0)
  {
    A.indicator <- NULL
    catlabs <- NULL
    for(i in 1:n.A)
    {
      dummy <- NULL
      levels <- sort(unique(A[,i]))
      for(level in levels[-1]) {
        catlabs <- c(catlabs, paste0('cat',i,'_','level',level))
        dummy <- cbind(dummy, A[,i]==level)
      }
      A.indicator <- cbind(A.indicator, dummy)
    }
    colnames(A.indicator) <- catlabs
  }

  ####################################
  #### Chain of coef for B
  coef.chain <- vector("list", n.B)
  names(coef.chain) <- colnames(B)
  coef.chain[[1]] <- mean(B[,1])
  for (i in 2:n.B) {
    coef.chain[[i]] <- tryCatch({
      coef(glm(B[,i] ~ B[, 1:(i-1)], family = "binomial"))
    }, error = function(e) {
      message("Chain regression failed for B[", i, "]: ", e$message)
      NA
    })
  }

  #### Coef for relationship between binary and categorical confounders

  coef.AonB <- vector("list", n.A)
  names(coef.AonB) <- colnames(A)

  # For the first categorical variable, regress on all B variables
  coef.AonB[[1]] <- tryCatch({
    coef(nnet::multinom(A[, 1] ~ B, trace = FALSE))
  }, error = function(e) {
    message("Model fitting failed for A[1]: ", e$message)
    NA
  })

  # For subsequent categorical variables, regress on the indicators of previous categorical variables and B
  if (n.A > 1) {
    for (i in 2:n.A) {
      # Get the corresponding indicator columns for all previous categorical variables
      start_index <- 1
      end_index <- sum(sapply(1:(i - 1), function(k) length(unique(A[, k])) - 1))

      # Extract the relevant indicator columns
      previous_A_indicators <- A.indicator[, start_index:end_index, drop = FALSE]

      # Create a formula using previous indicators and B as predictors
      formula <- as.formula(paste("A[, i] ~ . + B"))  # Use . to represent all predictors

      # Fit the multinomial logistic regression model, catch errors due to mismatched rows or other issues
      coef.AonB[[i]] <- tryCatch({
        model_data <- data.frame(A = A[, i], previous_A_indicators, B)
        model <- nnet::multinom(A ~ ., data = model_data, trace = FALSE)
        coef(model)
      }, warning = function(w) {
        message("Warning during model fitting for A[", i, "]: ", w$message)
        NA
      }, error = function(e) {
        message("Model fitting failed for A[", i, "]: ", e$message)
        NA
      })
    }
  }


  ####################################

  ### Exposure variable X|B,A (Propensity Score Model Coefficients)
  #### Intercept and coefficients from logistic regression
  ps.fit <- glm.fit(cbind(1,B,A.indicator), X, control=glm.ctrl, family=binomial())
  class(ps.fit) <- 'glm'
  coef.XonZ <- coef(ps.fit)
  ps.by.x <- cbind(X, fitted(ps.fit))
  ps.vcov <- vcov(ps.fit)
  probs <- predict(ps.fit, type = "response")
  ps.c <- mean(sample(probs[X == 1L], 1000000L, TRUE) > sample(probs[X == 0L], 1000000L, TRUE))

  #### Matthew added: Binary outcome
  #### Logistic regression for the event distribution

  logit.fit <- glm.fit(cbind(1,X,B,A.indicator), Y, control=glm.ctrl, family=binomial())
  class(logit.fit) <- 'glm'
  coef.Yon1 <- coef(logit.fit)[1]
  coef.YonX <- coef(logit.fit)[2]
  coef.YonZ <- coef(logit.fit)[-c(1,2)]

  control.events <- sum(Y[X==0])
  compare.events <- sum(Y[X==1])
  N.X <- table(X)

  # return binary outcome version

  return(list(n=length(Y),N.X=N.X,
              control.events=control.events,compare.events=compare.events,
              P=P,Common.P=Common.P,P.ord=P.ord,
              Corr.ord=Corr.ord,coef.XonZ=coef.XonZ,Coef.bin=coef.chain, Coef.cat=coef.AonB,
              propensity=ps.by.x,
              propensity.vcov=ps.vcov,cStat=ps.c,
              coef.Yon1=coef.Yon1,coef.YonX=coef.YonX,coef.YonZ=coef.YonZ  ## Matthew added
  ))
}
#####


# Hidden functions used to generate data ----
.ordgendata <- function(n, sigma, quants.norm){
  retval = mvtnorm::rmvnorm(n = n, sigma = sigma)
  for (i in 1:ncol(sigma)) {
    retval[, i] = cut(x = retval[, i], breaks = c(-1/0, quants.norm[[i]]),
                      right = FALSE)
  }
  retval - 1
}


.gencov.ord <- function(n, P.ord, Quant.norm, Corr.norm, coef.XonZ){

  Z <- .ordgendata(n, sigma=Corr.norm, quants.norm=Quant.norm)
  n.B <- length(which(unlist(lapply(P.ord,FUN=function(x){length(x)})==2)))
  n.A <- ncol(Z) - n.B

  if (n.B > 0 ) {
    B <- Z[, 1:n.B, drop=FALSE]
    colnames(B) <- paste0('bin', 1:n.B)
  }

  if(n.A>0)
  {
    A <- Z[,(n.B+1):ncol(Z)]
    A.indicator <- NULL
    catlabs <- NULL
    for(i in 1:n.A)
    {
      dummy <- NULL
      levels <- sort(unique(A[,i]))
      for(level in levels[-1]) {
        catlabs <- c(catlabs, paste0('cat',i,'_','level',level))
        dummy <- cbind(dummy, A[,i]==level)
      }
      A.indicator <- cbind(A.indicator, dummy)
    }
    colnames(A.indicator) <- catlabs
  }

  ### Confounders Z
  Z.model.data <- as.matrix(cbind(B, A.indicator))

  ### Exposure variable X
  X <- rbinom(n,size=1,
              p=1 / (1 + exp(-(coef.XonZ[1] +
                                 Z.model.data %*% coef.XonZ[2:length(coef.XonZ)]))))

  return(list(X=X, B=B, A.indicator=A.indicator))
}


.gencov <- function(n, P, Common.P, coef.AonB, coef.XonZ){
  B <- bindata::rmvbin(n, margprob=P, commonprob=Common.P)

  ### Categorical variables A
  ## Obtain correct predicted probabilities from multinomial coefficients by person
  P.A <- exp(cbind(rep(1,n),B)%*%t(coef.AonB[[1]]) )
  P.A <- cbind(rep(1,n),P.A)
  P.A <- P.A/apply(P.A,1,sum) # gives same result as fitted(multinom(age ~ z + x))

  # Generate Confounders A
  A.indicator <- t(apply(P.A,1,FUN=function(y){rmultinom(n=1,size=1,prob=y)}))
  #A <- factor( as.character(c(1,2,3,4)%*%t(A.indicator)) )
  A.indicator <- A.indicator[,-1]

  # For multiple categorical variables
  n.A <- length(coef.AonB)
  if(n.A>1)
  {
    for(i in 2:length(coef.AonB))
    {
      P.A2 <- exp( cbind(rep(1,n),A.indicator,B)%*%t(coef.AonB[[i]]) )
      # P.A2 <- exp( cbind(rep(1,n),B)%*%t(coef.AonB[[i]]) )
      P.A2 <- cbind(rep(1,n),P.A2)
      P.A2 <- P.A2/apply(P.A2,1,sum) # gives same result as fitted(multinom(age ~ z + x))

      A.indicator.tmp <- t(apply(P.A2,1,FUN=function(y){rmultinom(n=1,size=1,prob=y)}))
      #A <- cbind(A, factor( as.character(c(1,2,3,4)%*%t(A.indicator.tmp)) ))
      A.indicator <- cbind(A.indicator, A.indicator.tmp[,-1])
    }
  }

  ### Confounders Z
  Z.model.data <- as.matrix(cbind(B, A.indicator))

  ### Exposure variable X
  X <- rbinom(n,size=1,
              p=1 / (1 + exp(-(coef.XonZ[1] + Z.model.data %*% coef.XonZ[2:length(coef.XonZ)]))))

  return(list(X=X, B=B, A.indicator=A.indicator))
}


.gencov.chain <- function(n, coef.chain, coef.AonB, coef.XonZ,names.cat=NULL){
  n.B = length(coef.chain)
  B = rbinom(n=n,size=1,prob=coef.chain[[1]])
  for(i in 2:n.B){
    P.B = (1+exp(- (cbind(rep(1,n),B)%*%(coef.chain[[i]])) ))^(-1)
    B = cbind(B, rbinom(n=n,size=1,prob=P.B))
  }
  colnames(B) <- names(coef.chain)
  ### Categorical variables A
  ## Obtain correct predicted probabilities from multinomial coefficients by person
  P.A <- exp(cbind(rep(1,n),B)%*%t(coef.AonB[[1]]) )
  P.A <- cbind(rep(1,n),P.A)
  P.A <- P.A/apply(P.A,1,sum) # gives same result as fitted(multinom(age ~ z + x))

  # Generate Confounders A
  A.indicator <- t(apply(P.A,1,FUN=function(y){rmultinom(n=1,size=1,prob=y)}))
  #A <- factor( as.character(c(1,2,3,4)%*%t(A.indicator)) )
  A.indicator <- A.indicator[,-1]

  # For multiple categorical variables
  n.A <- length(coef.AonB)
  if(n.A>1)
  {
    for(i in 2:length(coef.AonB))
    {
      P.A2 <- exp( cbind(rep(1,n),A.indicator,B)%*%t(coef.AonB[[i]]) ) # Matthew: Q: exp( cbind(rep(1,n),B)%*%t(coef.AonB[[i]]) ) Question: want to regress A[2] on A[1]?
      P.A2 <- cbind(rep(1,n),P.A2)
      P.A2 <- P.A2/apply(P.A2,1,sum) # gives same result as fitted(multinom(age ~ z + x))

      A.indicator.tmp <- t(apply(P.A2,1,FUN=function(y){rmultinom(n=1,size=1,prob=y)}))
      #A <- cbind(A, factor( as.character(c(1,2,3,4)%*%t(A.indicator.tmp)) ))
      A.indicator <- cbind(A.indicator, A.indicator.tmp[,-1])
    }
  }

  colnames(A.indicator) <- names.cat

  ### Confounders Z
  Z.model.data <- as.matrix(cbind(B, A.indicator))

  ### Exposure variable X
  X <- rbinom(n,size=1,
              p=1 / (1 + exp(-(coef.XonZ[1] + Z.model.data %*% coef.XonZ[2:length(coef.XonZ)]))))


  return(list(X=X, B=B, A.indicator=A.indicator))
}


.gendata.survival <- function(logHR.X.site=NULL, n, P, Common.P, coef.AonB=NULL, coef.XonZ,
                              coef.cens, scale.cens,
                              coef.event, scale.event,
                              censtype="simple", trunc=366,
                              P.presc.topK=NULL, prescription.mode.topK=NULL,
                              method=1, Corr.norm=NULL, Quant.norm=NULL, P.ord=NULL,
                              coef.chain=NULL, user.data=NULL,noX=FALSE,
                              strat.var.cens=NULL,strat.var.event=NULL)
{
  if (is.null(coef.cens) || is.null(scale.cens) || is.null(coef.event) || is.null(scale.event)){
    stop("For survival outcome, coef.cens, scale.cens, coef.event, and scale.event must be specified")
  }
  if (method == 4 && is.null(user.data)) stop("User data must be provided for method 4")

  ### set up specified HR of X
  if (!noX && !is.null(logHR.X.site)) {
    coef.event[grep("X",names(coef.event))] <- -logHR.X.site*scale.event
  }

  ### Generate covariates via specified method
  exppluscovs <- switch(method,
                        `1` = .gencov.ord(n, P.ord, Quant.norm, Corr.norm, coef.XonZ),
                        `2` = .gencov(n, P, Common.P, coef.AonB, coef.XonZ),
                        `3` = .gencov.chain(n, coef.chain, coef.AonB, coef.XonZ),
                        `4` = {
                          exppluscovs <- list(B=NULL, A.indicator=NULL)
                          exppluscovs$A.indicator <- user.data[,-1]
                          exppluscovs$X <- if (noX) NULL else user.data[, 1]
                          exppluscovs
                        },
                        stop("Invalid method specified"))

  ### Confounders Z
  Z.model.data <- cbind(exppluscovs$B, exppluscovs$A.indicator)

  ### Exposure variable X
  X <- exppluscovs$X

  ### Censoring
  censorT<-NULL
  u <- runif(n)

  if(censtype%in%c("simple","simplebump"))
  {
    censorT <- ceiling((-log(u) * exp(coef.cens / scale.cens)) ^ (scale.cens))
    censorT.1 <- ceiling((-log(runif(n)) * exp(coef.cens / scale.cens)) ^ (scale.cens))
    censorT.0 <- ceiling((-log(runif(n)) * exp(coef.cens / scale.cens)) ^ (scale.cens))
  }
  if(censtype%in%c("cov","covbump"))
  {
    if (length(scale.cens)==1) {

      censorT <- ceiling((-log(u)*exp( cbind(1,X,Z.model.data) %*% coef.cens/scale.cens ))^(scale.cens))
      if (!is.null(X)) {
        censorT.1 <- ceiling((-log(runif(n))*exp( cbind(rep(1,n),1,Z.model.data) %*% coef.cens/scale.cens ))^(scale.cens))
        censorT.0 <- ceiling((-log(runif(n))*exp( cbind(rep(1,n),0,Z.model.data) %*% coef.cens/scale.cens ))^(scale.cens))
      } else {
        censorT.1 <- ceiling((-log(runif(n))*exp( cbind(rep(1,n),Z.model.data) %*% coef.cens/scale.cens ))^(scale.cens))
        censorT.0 <- ceiling((-log(runif(n))*exp( cbind(rep(1,n),Z.model.data) %*% coef.cens/scale.cens ))^(scale.cens))
      }

    } else {

      censorT <- rep(0,times=n)
      censorT[strat.var.cens==0] <- ceiling((-log(u[strat.var.cens==0])*exp(
        cbind(1,X,Z.model.data)[strat.var.cens==0,] %*%
          coef.cens/scale.cens[[1]] ))^(scale.cens[[1]]))
      censorT[strat.var.cens==1] <- ceiling((-log(u[strat.var.cens==1])*exp(
        cbind(1,X,Z.model.data)[strat.var.cens==1,] %*%
          coef.cens/scale.cens[[2]] ))^(scale.cens[[2]]))
    }

  }
  if(is.null(censorT)) stop("censoring model is misspecified")

  ### Event time
  ue <- runif(n)

  # TODO: 9.24 adj.coef.event has a na column, need to find why. KPNC, cat_level2012
  if (length(scale.event)==1) {
    eventT <- ceiling((-log(ue) * exp( cbind(1,X,Z.model.data) %*% coef.event / scale.event )) ^ (scale.event))
    if (!is.null(X)) {
      eventT.1 <- ceiling((-log(runif(n)) * exp( cbind(1,1,Z.model.data) %*% coef.event / scale.event )) ^ (scale.event))
      eventT.0 <- ceiling((-log(runif(n)) * exp( cbind(1,0,Z.model.data) %*% coef.event / scale.event )) ^ (scale.event))
    } else {
      eventT.1 <- ceiling((-log(runif(n)) * exp( cbind(1,Z.model.data) %*% coef.event / scale.event )) ^ (scale.event))
      eventT.0 <- ceiling((-log(runif(n)) * exp( cbind(1,Z.model.data) %*% coef.event / scale.event )) ^ (scale.event))
    }

  } else {
    eventT <- rep(0,times=n)
    eventT[strat.var.event==0] <- ceiling((-log(ue[strat.var.event==0])*exp(
      cbind(1,X,Z.model.data)[strat.var.event==0,] %*%
        coef.event/scale.event[[1]] ))^(scale.event[[1]]))
    eventT[strat.var.event==1] <- ceiling((-log(ue[strat.var.event==1])*exp(
      cbind(1,X,Z.model.data)[strat.var.event==1,] %*%
        coef.event/scale.event[[2]] ))^(scale.event[[2]]))
  }


  ### Observed time and event indicator
  if(censtype%in%c("simplebump","covbump"))
  {

    trunc.presc <- apply(rmultinom(n,1,c(P.presc.topK,'NA' = 1-sum(P.presc.topK))),
                         2,function(XX) {c(prescription.mode.topK,NA)[which(XX==1)]})
    CensorTbump <- censorT
    CensorTbump[!is.na(trunc.presc)] <- trunc.presc[!is.na(trunc.presc)]

    E <- apply(cbind(eventT,trunc,CensorTbump),1,min,na.rm=T)
    Y   <- as.numeric(ifelse(E == eventT,1,0))

    CensorTbump.1 <- censorT.1
    CensorTbump.1[!is.na(trunc.presc)] <- trunc.presc[!is.na(trunc.presc)]
    CensorTbump.0 <- censorT.0
    CensorTbump.0[!is.na(trunc.presc)] <- trunc.presc[!is.na(trunc.presc)]

    E.1 <- apply(cbind(eventT.1,trunc,CensorTbump.1),1,min,na.rm=T)
    E.0 <- apply(cbind(eventT.0,trunc,CensorTbump.0),1,min,na.rm=T)

    Y.1 <- as.numeric(ifelse(E.1 == eventT.1,1,0))
    Y.0 <- as.numeric(ifelse(E.0 == eventT.0,1,0))

  }
  else
  {
    E <- apply(cbind(eventT,trunc,censorT),1,min,na.rm=T)
    Y <- as.numeric(ifelse(E == eventT,1,0))

    E.1 <- apply(cbind(eventT.1,trunc,censorT.1),1,min,na.rm=T)
    E.0 <- apply(cbind(eventT.0,trunc,censorT.0),1,min,na.rm=T)

    Y.1 <- as.numeric(ifelse(E.1 == eventT.1,1,0))
    Y.0 <- as.numeric(ifelse(E.0 == eventT.0,1,0))
  }


  ## Marginal sample needs to be same size as original
  marg.x1 <- sample(1:n,size=sum(X))

  return(list(marginal.dat=rbind(cbind(x=1, y=Y.1, obst=E.1)[marg.x1,],cbind(x=0, y=Y.0, obst=E.0)[-marg.x1,]),
              B=exppluscovs$B,A.indicator=exppluscovs$A.indicator,
              X=exppluscovs$X,Y=Y,E=E)
    )

}


.gendata.binary <- function(n, P, Common.P, coef.AonB=NULL, coef.XonZ,
                            method=1, Corr.norm=NULL, Quant.norm=NULL, P.ord=NULL,
                            coef.chain=NULL, user.data=NULL,noX=FALSE,
                            coef.Yon1, coef.YonX, coef.YonZ, set.coef.YonX=NULL)
{
  # test 10.8
  # n=SS$n
  # P=SS$P
  # Common.P=SS$Common.P
  # coef.XonZ=SS$coef.XonZ
  # coef.chain=SS$Coef.bin
  # coef.AonB=SS$Coef.cat
  # method=method
  # Corr.norm=SS$Corr.norm
  # Quant.norm=SS$Quants.norm
  # P.ord=SS$P.ord
  # coef.Yon1=SS$coef.Yon1
  # coef.YonX=SS$coef.YonX
  # coef.YonZ=SS$coef.YonZ


  if (is.null(coef.Yon1) || is.null(coef.YonX) || is.null(coef.YonZ)) {
    stop("For binary outcome, coef.Yon1, coef.YonX, and coef.YonZ must be specified")
  }
  if (method == 4 && is.null(user.data)) stop("User data must be provided for method 4")

  ### Set custom coef.YonX if provided
  if (!is.null(set.coef.YonX)) {
    coef.YonX <- set.coef.YonX
  }

  ### Generate covariates via specified method

  exppluscovs <- switch(method, # 9.24: try this new code
                        `1` = .gencov.ord(n, P.ord, Quant.norm, Corr.norm, coef.XonZ),
                        `2` = .gencov(n, P, Common.P, coef.AonB, coef.XonZ),
                        `3` = .gencov.chain(n, coef.chain, coef.AonB, coef.XonZ),
                        `4` = {
                          exppluscovs <- list(B=NULL,A.indicator=NULL)
                          exppluscovs$A.indicator <- user.data[,-c(1)]
                          exppluscovs$X <- if(noX) NULL else user.data[,1]
                          exppluscovs
                        },
                        stop("Invalid method specified")
  )

  ### Confounders Z
  Z.model.data <- cbind(exppluscovs$B, exppluscovs$A.indicator)

  ### Exposure variable X
  X <- exppluscovs$X

  #### Matthew added: Binary outcome
  Y <- rbinom(n,size=1,
              p=1 / (1 + exp(-(coef.Yon1 + X * coef.YonX + Z.model.data %*% coef.YonZ))))
  Y.1 <- rbinom(n,size=1,
                p=1 / (1 + exp(-(coef.Yon1 + coef.YonX + Z.model.data %*% coef.YonZ))))
  Y.0 <- rbinom(n,size=1,
                p=1 / (1 + exp(-(coef.Yon1 + Z.model.data %*% coef.YonZ))))

  ## Marginal sample needs to be same size as original
  marg.x1 <- sample(1:n,size=sum(X))

  return(list(marginal.dat=rbind(cbind(x=1, y=Y.1)[marg.x1,],cbind(x=0, y=Y.0)[-marg.x1,]),
              B=exppluscovs$B,A.indicator=exppluscovs$A.indicator,
              X=exppluscovs$X,Y=Y)
  )
}
#####


# Functions for generating data ----
generate.data.survival <- function(Summ.Stat,censtype="simple", trunc=365,method=1, set.logHR.X=NULL){
  n.sites<-length(Summ.Stat)

  ## Data.Site: simulated data across sites
  Data.Simulated <- NULL
  Data.Marg <- NULL

  ## Data.Pool.forPS: pooled data for est pooled PS
  Data.Pool.forPS <- NULL

  ## Loop through sites to generate site specific data
  for(i in 1:n.sites)
  {
    ## read in summary statistics SS
    SS <- Summ.Stat[[i]]
    ## get common parameters
    if(is.null(set.logHR.X)) {logHR.X.site <- SS$logHR.X} else {logHR.X.site <- set.logHR.X}
    n <- SS$n
    P <- SS$P
    Common.P <- SS$Common.P
    coef.XonZ <- SS$coef.XonZ
    coef.event <- c(SS$intercept, SS$adj.coef.event[-1])
    scale.event <- SS$adj.scale.event
    coef.chain <- SS$Coef.bin
    coef.AonB <- SS$Coef.cat
    Corr.norm <- SS$Corr.norm
    Quant.norm <- SS$Quants.norm
    P.ord <- SS$P.ord

    # set coef.cens and scale.cens for different methods
    if (censtype == "simple") {
      coef.cens <- SS$simple.coef.cens
      scale.cens <- SS$simple.scale.cens
    } else if (censtype == "simplebump") {
      coef.cens <- SS$simplebump.coef.cens
      scale.cens <- SS$simplebump.scale.cens
    } else if (censtype == "cov") {
      coef.cens <- SS$cov.coef.cens
      scale.cens <- SS$cov.scale.cens
    } else if (censtype == "covbump") {
      coef.cens <- SS$covbump.coef.cens
      scale.cens <- SS$covbump.scale.cens
    } else {
      stop("Censoring type is misspecified")
    }

    P.presc.topK <- if (censtype %in% c("simplebump", "covbump")) SS$P.presc.topK else NULL
    prescription.mode.topK <- if (censtype %in% c("simplebump", "covbump")) SS$prescription.mode.topK else NULL

    ## generate site specific data
    DS <- .gendata.survival(
      logHR.X.site = logHR.X.site,
      n = n, P = P, Common.P = Common.P, coef.XonZ = coef.XonZ,
      coef.cens = coef.cens, scale.cens = scale.cens,
      coef.event = coef.event, scale.event = scale.event,
      censtype = censtype, trunc = trunc, coef.chain = coef.chain, coef.AonB = coef.AonB,
      P.presc.topK = P.presc.topK, prescription.mode.topK = prescription.mode.topK,
      method = method, Corr.norm = Corr.norm, Quant.norm = Quant.norm, P.ord = P.ord
    )

    # TODO: method 4 need to be tested

    if(is.null(DS)) stop("Censoring type is misspecified")


    ## save site specific data
    if (n.sites > 1) {
      Data.Simulated <- rbind(Data.Simulated, data.frame(B=DS$B,A=DS$A.indicator,X=DS$X,Y=DS$Y,E=DS$E,site=i))
      Data.Marg <- rbind(Data.Marg,data.frame(site=i,DS$marginal.dat))
    } else {
      Data.Simulated <- data.frame(B=DS$B, A=DS$A.indicator, X=DS$X,Y=DS$Y,E=DS$E)
      Data.Marg <- data.frame(DS$marginal.dat)
    }

  }

  return(list(Data.Simulated=Data.Simulated,Data.Marginal=Data.Marg))
}


generate.data.binary <- function(Summ.Stat,method=1, set.coef.YonX=NULL){
  n.sites<-length(Summ.Stat)

  ## Data.Site: simulated data across sites
  Data.Simulated <- NULL
  Data.Marg <- NULL

  ## Data.Pool.forPS: pooled data for est pooled PS
  Data.Pool.forPS <- NULL

  ## Loop through sites to generate site specific data
  for(i in 1:n.sites)
  {
    ## read in summary statistics SS
    SS <- Summ.Stat[[i]]
    ## generate site specific data

    # TODO: method 4 working?

    DS <- .gendata.binary(n=SS$n, P=SS$P, Common.P=SS$Common.P, coef.XonZ=SS$coef.XonZ,
                         coef.chain=SS$Coef.bin, coef.AonB=SS$Coef.cat,
                         method=method, Corr.norm=SS$Corr.norm, Quant.norm=SS$Quants.norm, P.ord=SS$P.ord,
                         coef.Yon1=SS$coef.Yon1, coef.YonX=SS$coef.YonX, coef.YonZ=SS$coef.YonZ,
                         set.coef.YonX=set.coef.YonX)

    ## save site specific data
    if (n.sites > 1) {
      Data.Simulated <- rbind(Data.Simulated, data.frame(B=DS$B,A=DS$A.indicator,X=DS$X,Y=DS$Y,site=i))
      Data.Marg <- rbind(Data.Marg,data.frame(site=i,DS$marginal.dat))
    } else {
      Data.Simulated <- data.frame(B=DS$B, A=DS$A.indicator, X=DS$X,Y=DS$Y)
      Data.Marg <- data.frame(DS$marginal.dat)
    }

  }

  return(list(Data.Simulated=Data.Simulated,Data.Marginal=Data.Marg))
}


# TODO: (1) debug, make these functions able to replicate previous plot (survival+binary): done
# (2) check all the 4 methods in the survival functions can work: done, except for method 4 which require user data
# (3) Fixed the issue that categorical variables did not regress on previous categorical variables: done
# (4) think about how to edit the Hazard ratio, probably in the .gendata.survival function: done
# (5) add notations critical for the package

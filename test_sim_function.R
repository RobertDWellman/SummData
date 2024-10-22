rm(list=ls())
library(foreach)
library(doParallel)

source('~/Yuexiang Peng/UW/Research/Jennifer Nelson/Cook/summer project/simu_0817/Final Code/functions.R')
source('D:/OneDrive - UW/Documents/GitHub/sim.realistic.data/R/sim_functions_organized.R')
datdir <- "~/Yuexiang Peng/UW/Research/Jennifer Nelson/Cook/summer project/simu_0817/Data/"
resdir <- "~/Yuexiang Peng/UW/Research/Jennifer Nelson/Cook/summer project/simu_0817/Result/"

# source('G:/CTRHS/Sentinel/Y6_Task_Orders/Methods/Cook-SurvivalII/Programming/Matthew/Final Code/functions_original.R')
# datdir <- "G:/CTRHS/Sentinel/Y7_Task_Orders_2015/Big sim/Data/datatoGHC/datatoGHC/Angioedema/"
# resdir <- '//groups/data/CTRHS/Sentinel/Y6_Task_Orders/Methods/Cook-SurvivalII/Programming/Matthew/Data/'


require(survival)
glm_control <- glm.control(epsilon = 1e-8, maxit = 25, trace = FALSE)
cox.ctrl <- survival::coxph.control(eps=1e-09, toler.chol=.Machine$double.eps^0.75,
                                    iter.max=10, toler.inf=sqrt(1e-09), outer.max=10)
sites <- c("AEOS","KPNC","HUOS","HCOS","HPHC")
names(sites) <- rep("site", length(sites))

## Type for censoring time, to fit Weibull model:
#### censtype = "simple", "simplebump", "cov", "covbump"
censtype <- "cov"

# ## Relax cox fitting settings to run faster
# cox.ctrl <- survival::coxph.control(eps=1e-09, toler.chol=.Machine$double.eps^0.75,
#                                     iter.max=20, toler.inf=sqrt(1e-09), outer.max=10)
tie.method <- "efron" ##"exact"

## Maximum follow-up Time (i.e. 1 year = 366 days)
trunc <- 365 # will be used when generating tab1

## Compute proportionate sample sizes using summary stats from full data
Summ.Stat <- readRDS(paste0(resdir,"angio_summary_20160507_cov_chain.rds"))

Sample.size.site <- do.call("cbind", lapply(Summ.Stat, function(XX) {
  n.curr <- NULL
  n.curr <- matrix(XX$n)
  colnames(n.curr) <- XX$site
  n.curr
}))

N.sim <- 150000
N.sim.site <- ceiling(N.sim*(Sample.size.site/sum(Sample.size.site)))
N.sim.site[[length(sites)]] <- N.sim - sum(N.sim.site[1:(length(sites)-1)])

## Run summary statistics at each site
begin <- proc.time()
cl <- makeCluster(length(sites))
registerDoParallel(cl)
seed1 <- 596841 # change to see if we avoid creating na. original:596840
clusterSetRNGStream(cl = cl, iseed = seed1) #Multiple streams of seeds
Summ.Stat <- NULL
system.time(
  Summ.Stat <- foreach(site=sites,size=N.sim.site) %dopar% {
    # #test
    # site <- "KPNC"
    # size <- 19275

    Y <- E <- X <- comor <- ipvisits <- age_cat <- white <- sex <- dat <- NULL
    SS.list <- NULL
    dat <- readRDS(paste0(datdir,paste0(site,".rds")))

    names(dat) <- tolower(names(dat))
    dat <- dat[dat$year < 2013,]
    dat <- dat[sample(1:dim(dat)[[1]],size),] # Matthew: why do we need resample? I think is to make the sample size smaller. If not resample, create NA at adj.coef.event for unknown reason

    Y <- dat$event_astreated
    E <- dat$followuptime_astreated
    Y[E>365] <- 0
    E[E>365] <- 365

    ## Pattern of Prescriptions fills for bump detection
    prescription.mode <- as.numeric(names(tail(sort(table(E[Y==0])/length(E[Y==0])),3)))

    ## number of top frequent bumps to pick up
    my.presc.K <- length(prescription.mode)

    X <- dat$exposure

    ## Binary covariates must be coded as 0/1
    comor <- cut(dat$comorbidscore, breaks = c(-2.5,0,100), labels=FALSE) - 1
    ipvisits <- cut(dat$numip, breaks = c(-1,0,100), labels=FALSE) - 1
    age_cat <- cut(dat$age, breaks = c(0,44,54,64,150), labels=FALSE)
    ed1plus <- dat$numed > 0
    sex <- dat$sex == 'F'

    B <- cbind(ipvisits,ed1plus,comor,sex)
    C <- cbind(age_cat, year=dat$year)

    # test for get.summstat.survival function
    SS.list <- get.summstat.survival(E=E,Y=Y,X=X,B=B,A=C,prescription.mode,
                            my.presc.K,tie.method)

    # test for get.summstat.binary function
    # SS.list <- get.summstat.binary(Y=Y,X=X,B=B,A=C)

    SS.list <- append(SS.list, site, after=0)

    # norm.spec <- .ordtonorm(probs=SS.list$P.ord, Cor=SS.list$Corr.ord)


    append(SS.list,
           list(
             # Corr.norm=norm.spec$corr.norm,
             # Quants.norm=norm.spec$quants.norm,
             # logHR.X=SS.list$cox.coef.adjusted[[1]],
             # intercept=SS.list$adj.coef.event[[1]],
             dat.boot=data.frame(X,Y,E,B,C))
           )

  }
)
stopCluster(cl)
proc.time() - begin
saveRDS(Summ.Stat,paste0(resdir,"angio_summary_cov_chain_survival_241020_2.rds")) ## Matthew
# saveRDS(Summ.Stat,paste0(resdir,"angio_summary_cov_chain_binary_240924.rds"))


## Sample the same dataset that was used for summary statistics.
## This will be the basis for boostrap sampling.

Summ.Stat <- readRDS(paste0(resdir,"angio_summary_cov_chain_survival_241020_2.rds"))
# Summ.Stat <- readRDS(paste0(resdir,"angio_summary_cov_chain_binary_240924.rds"))

names(Summ.Stat) <- sites
dat.all <- do.call(rbind, lapply(1:5, function(XX) cbind(Summ.Stat[[XX]]$dat.boot,site=XX))) # Matthew: still need numeric

B <- with(dat.all, cbind(ipvisits,ed1plus,comor,sex))
C <- with(dat.all, cbind(age_cat,year))
S <- as.numeric(dat.all$site) # Matthew: generate NA here

C.indicator <- NULL
for (i in 1:ncol(C)) {
  for (j in sort(unique(C[,i]))[-c(1)]) {
    C.indicator <- cbind(C.indicator, as.integer(C[,i]==j))
  }
}

S.indicator <- NULL
for (j in sort(unique(S))[-c(1)]) {
  S.indicator <- cbind(S.indicator, as.integer(S==j))
}


dat.boot <- cbind(B,C.indicator,X=dat.all$X,Y=dat.all$Y,
                  E=dat.all$E,site=as.numeric(dat.all$site))

# dat.boot <- cbind(B,C.indicator,X=dat.all$X,Y=dat.all$Y,
#                   E=dat.all$E,site=dat.all$site)

prop.labels <- c("Intercept","1+ IP Visits","1+ ED Visits", "Comorbid 1+","Female",
                 "Age 45-54","Age 55-64","Age 65+","2009","2010","2011","2012",
                 "Site 2","Site 3","Site 4","Site 5")

out.labels <- c("ACEI","1+ IP Visits", "1+ ED Visits", "Comorbid 1+","Female",
                "Age 45-54","Age 55-64","Age 65+","2009","2010","2011","2012",
                "Site 2","Site 3","Site 4","Site 5")


## Summary of dataset used in example
ps.boot <- glm(X~.,data=as.data.frame(cbind(dat.boot[,1:12],S.indicator)))
ps.boot.coef <- ps.boot$coefficients
ps.boot.se <- sqrt(diag(vcov(ps.boot)))

ps.boot.tab <- round(cbind(ps.boot.coef,ps.boot.se),4)
rownames(ps.boot.tab) <- prop.labels

ps.boot.site <- NULL
ps.boot.site.or <- NULL
for (i in 1:length(sites)) {
  temp <- temp.coef <- temp.se <- temp.tab <- index <- NULL
  index <- dat.boot[,'site'] == i
  temp <- glm(X~.,data=as.data.frame(cbind(dat.boot[index,1:12])))
  temp.coef <- temp$coefficients
  temp.or <- exp(temp.coef)
  temp.se <- sqrt(diag(vcov(temp)))
  temp.or.se <- sqrt(temp.or*temp.or*temp.se*temp.se)
  temp.tab <- round(cbind(temp.coef,temp.se),4)
  temp.tab.or <- round(cbind(temp.or,temp.or.se),4)
  ps.boot.site <-cbind(ps.boot.site,temp.tab)
  ps.boot.site.or <-cbind(ps.boot.site.or,temp.tab.or)
}
rownames(ps.boot.site) <- prop.labels[1:12]
rownames(ps.boot.site.or) <- prop.labels[1:12]

# het <- unlist(foreach(i=Summ.Stat) %do% exp(i$cox.coef.adjusted[[1]]))# Matthew: No cox.coef.adjusted for binary. delete the 2 rows
# logHR <- list(1.0,1.0,1.5,1.5,2.0,2.0,het,het)

# # parallel computing version ----
# samp.size <- dim(dat.boot)[[1]]
# n.sim <- 5
# #seed1<-83745 bump
# #seed1<-379 simple
# seed1 <-566391
# n_cores <- detectCores()
# cl<-makeCluster(n_cores) #Don't need to load parallel
#
# registerDoParallel(cl, cores = n_cores )
# clusterSetRNGStream(cl = cl, iseed = seed1) #Multiple streams of seeds
# res<-NULL
#
# start <- proc.time()[3]
# res <- foreach(c=1:n_cores,
#                .verbose=TRUE,
#                .errorhandling=c('pass'),.packages=c("foreach","survival")) %dopar% {
#                  # # test
#                  # c=1
#                  # .verbose=TRUE
#                  # .errorhandling=c('pass')
#                  # .packages=c("foreach","survival")
#                  sim.estimates <- sapply(1:n.sim,function(XX) NULL)
#                  PS.beta.site.norm <- sapply(1:length(sites),function(XX) NULL)
#                  PS.beta.site.chain <- sapply(1:length(sites),function(XX) NULL)
#                  PS.beta.site.boot <- sapply(1:length(sites),function(XX) NULL)
#                  PS.beta.norm <- NULL
#                  PS.beta.chain <- NULL
#                  PS.beta.boot <- NULL
#                  cox.site.norm  <- sapply(1:length(sites),function(XX) NULL)
#                  cox.site.chain <- sapply(1:length(sites),function(XX) NULL)
#                  cox.site.boot  <- sapply(1:length(sites),function(XX) NULL)
#
#                  cox.pooled.norm  <- NULL
#                  cox.pooled.chain <- NULL
#                  cox.pooled.boot  <- NULL
#                  common.prob.norm <- sapply(1:n.sim,function(XX) NULL)
#                  common.prob.chain <- sapply(1:n.sim,function(XX) NULL)
#                  common.prob.boot <- sapply(1:n.sim,function(XX) NULL)
#
#                  for (sim in 1:n.sim) {
#                    sim.dat.norm <- NULL
#                    sim.dat.chain <- NULL
#                    sim.dat.boot <- NULL
#                    PS.list.norm <- NULL
#                    PS.list.chain <- NULL
#                    PS.list.boot <- NULL
#                    cox.norm <- NULL
#                    cox.chain <- NULL
#                    cox.boot <- NULL
#
#                    # test survival
#                    ## simulate data using multivariate normal
#                    sim.dat.norm <- generate.data.survival(Summ.Stat, censtype=censtype, trunc=366)$Data.Simulated # Matthew:
#                    ## simulate data using covariate chain
#                    sim.dat.chain <- generate.data.survival(Summ.Stat, censtype=censtype, trunc=366,method=3)$Data.Simulated
#
#                    # # test binary
#                    # ## simulate data using multivariate normal
#                    # sim.dat.norm <- generate.data.binary(Summ.Stat, censtype=censtype)$Data.Simulated # Matthew:
#                    # ## simulate data using covariate chain
#                    # sim.dat.chain <- generate.data.binary(Summ.Stat, censtype=censtype,method=3)$Data.Simulated
#
#
#                    sim.dat.boot <- NULL
#                    sim.dat.boot0 <- NULL
#                    for (site in 1:length(sites)) {
#                      site.indx <- dat.boot[,"site"]==site
#                      sim.dat.boot0 <- as.data.frame(dat.boot[site.indx,][sample(1:sum(site.indx),replace=TRUE),])
#                      covs <- sim.dat.boot0[,c(1:11)]
#                      x <- sim.dat.boot0[,c(12)]
#                      y <- sim.dat.boot0[,c(13)]
#                      e <- sim.dat.boot0[,c(14)]
#                      s <- sim.dat.boot0[,c(15)]
#                      names(covs) <- colnames(Summ.Stat[[site]]$cox.fit.event$x)[-1]
#                      sim.dat.boot <- rbind(sim.dat.boot,cbind(covs,X=x,Y=y,E=e,site=s))
#
#                      site.indx <- NULL
#                      sim.dat.boot0 <- NULL
#                    }
#
#                    cov.cols <- 1:11
#
#                    common.prob.norm[[sim]] <- (t(as.matrix(sim.dat.norm[,cov.cols]))%*%as.matrix(sim.dat.norm[,cov.cols]))/samp.size
#                    common.prob.chain[[sim]] <- (t(as.matrix(sim.dat.chain[,cov.cols]))%*%as.matrix(sim.dat.chain[,cov.cols]))/samp.size
#                    common.prob.boot[[sim]] <- (t(as.matrix(sim.dat.boot[,cov.cols]))%*%as.matrix(sim.dat.boot[,cov.cols]))/samp.size
#
#
#                    PS.list.norm <- comp.pscore(X=sim.dat.norm$X,Z=sim.dat.norm[,cov.cols],S=sim.dat.norm$site)
#                    PS.list.chain <- comp.pscore(X=sim.dat.chain$X,Z=sim.dat.chain[,cov.cols],S=sim.dat.chain$site)
#                    PS.list.boot <- comp.pscore(X=sim.dat.boot[,'X'],Z=sim.dat.boot[,cov.cols],S=sim.dat.boot[,'site'])
#
#                    PS.beta.site.norm <- foreach(i=PS.beta.site.norm,j=PS.list.norm[[4]]) %do% cbind(i,j)
#                    PS.beta.site.chain <- foreach(i=PS.beta.site.chain,j=PS.list.chain[[4]]) %do% cbind(i,j)
#                    PS.beta.site.boot <- foreach(i=PS.beta.site.boot,j=PS.list.boot[[4]]) %do% cbind(i,j)
#                    PS.beta.norm <- cbind(PS.beta.norm,PS.list.norm[[2]])
#                    PS.beta.chain <- cbind(PS.beta.chain,PS.list.chain[[2]])
#                    PS.beta.boot <- cbind(PS.beta.boot,PS.list.boot[[2]])
#
#                    # Matthew: need a binary outcome version. still name as cox for now, change it later
#                    # test survival version
#
#                    # # test why Y has na
#                    # any(is.na(sim.dat.norm$Y))
#                    # sum(is.na(sim.dat.norm$Y))/length(sim.dat.norm$Y)
#
#                    # TODO: Warning message: Ran out of iterations and did not converge
#                    cox.norm <- cox.est(X=sim.dat.norm$X,Y=sim.dat.norm$Y,E=sim.dat.norm$E,
#                                        Z=sim.dat.norm[,cov.cols],S=sim.dat.norm$site)
#                    cox.chain <- cox.est(X=sim.dat.chain$X,Y=sim.dat.chain$Y,E=sim.dat.chain$E,
#                                         Z=sim.dat.chain[,cov.cols],S=sim.dat.chain$site)
#                    cox.boot <- cox.est(X=sim.dat.boot[,'X'],Y=sim.dat.boot[,'Y'],E=sim.dat.boot[,'E'],
#                                        Z=sim.dat.boot[,cov.cols],S=sim.dat.boot[,'site'])
#
#                    # # test binary version
#                    # cox.norm <- logistic.est(X=sim.dat.norm$X,Y=sim.dat.norm$Y,
#                    #                          Z=sim.dat.norm[,cov.cols],S=sim.dat.norm$site)
#                    # cox.chain <- logistic.est(X=sim.dat.chain$X,Y=sim.dat.chain$Y,
#                    #                           Z=sim.dat.chain[,cov.cols],S=sim.dat.chain$site)
#                    # cox.boot <- logistic.est(X=sim.dat.boot[,'X'],Y=sim.dat.boot[,'Y'],
#                    #                          Z=sim.dat.boot[,cov.cols],S=sim.dat.boot[,'site'])
#
#                    for (k in 1:5) {
#                      if (!is.null(cox.norm$Site.Specific[[k]])) {
#                        cox.site.norm[[k]]  <-  cbind(cox.site.norm[[k]], cox.norm$Site.Specific[[k]])
#                      } else {
#                        cox.site.norm[[k]]  <-  cbind(cox.site.norm[[k]], rep(NA,times=length(cov.cols)+1))
#                      }
#                      if (!is.null(cox.chain$Site.Specific[[k]])) {
#                        cox.site.chain[[k]]  <-  cbind(cox.site.chain[[k]], cox.chain$Site.Specific[[k]])
#                      } else {
#                        cox.site.chain[[k]]  <-  cbind(cox.site.chain[[k]], rep(NA,times=length(cov.cols)+1))
#                      }
#                      if (!is.null(cox.boot$Site.Specific[[k]])) {
#                        cox.site.boot[[k]]  <-  cbind(cox.site.boot[[k]], cox.boot$Site.Specific[[k]])
#                      } else {
#                        cox.site.boot[[k]]  <-  cbind(cox.site.boot[[k]], rep(NA,times=length(cov.cols)+1))
#                      }
#                    }
#
#
#                    cox.pooled.norm  <- cbind(cox.pooled.norm,cox.norm$Pooled)
#                    cox.pooled.chain <- cbind(cox.pooled.chain,cox.chain$Pooled)
#                    cox.pooled.boot  <- cbind(cox.pooled.boot,cox.boot$Pooled)
#
#                  }
#                  list(PS.beta.site.norm=PS.beta.site.norm,
#                       PS.beta.site.chain=PS.beta.site.chain,
#                       PS.beta.site.boot=PS.beta.site.boot,
#                       cox.site.norm=cox.site.norm,cox.site.chain=cox.site.chain,
#                       cox.site.boot=cox.site.boot, PS.beta.norm=PS.beta.norm,
#                       PS.beta.chain=PS.beta.chain,PS.beta.boot=PS.beta.boot,
#                       cox.pooled.norm=cox.pooled.norm,cox.pooled.chain=cox.pooled.chain,
#                       cox.pooled.boot=cox.pooled.boot,common.prob.norm=common.prob.norm,
#                       common.prob.chain=common.prob.chain,common.prob.boot=common.prob.boot)
#                }
#
# stopCluster(cl)
# end <- proc.time()[3]
# end-start
# #####

# no parallel version ----
samp.size <- dim(dat.boot)[[1]]
sim_num = 3
n.sim <- sim_num
# seed1<-83745 bump
# seed1<-379 simple
seed1 <- 566391
n_cores <- detectCores()

set.seed(seed1)
res <- NULL

start <- proc.time()[3]

res <- list()

for (c in 1:n_cores) {
  sim.estimates <- sapply(1:n.sim, function(XX) NULL)
  PS.beta.site.norm <- sapply(1:length(sites), function(XX) NULL)
  PS.beta.site.chain <- sapply(1:length(sites), function(XX) NULL)
  PS.beta.site.boot <- sapply(1:length(sites), function(XX) NULL)
  PS.beta.norm <- NULL
  PS.beta.chain <- NULL
  PS.beta.boot <- NULL
  cox.site.norm <- sapply(1:length(sites), function(XX) NULL)
  cox.site.chain <- sapply(1:length(sites), function(XX) NULL)
  cox.site.boot <- sapply(1:length(sites), function(XX) NULL)

  cox.pooled.norm <- NULL
  cox.pooled.chain <- NULL
  cox.pooled.boot <- NULL
  common.prob.norm <- sapply(1:n.sim, function(XX) NULL)
  common.prob.chain <- sapply(1:n.sim, function(XX) NULL)
  common.prob.boot <- sapply(1:n.sim, function(XX) NULL)

  for (sim in 1:n.sim) {
    sim.dat.norm <- generate.data.survival(Summ.Stat, censtype = censtype, trunc = 366)$Data.Simulated
    sim.dat.chain <- generate.data.survival(Summ.Stat, censtype = censtype, trunc = 366, method = 3)$Data.Simulated # 3

    sim.dat.boot <- NULL
    sim.dat.boot0 <- NULL
    for (site in 1:length(sites)) {
      site.indx <- dat.boot[, "site"] == site
      sim.dat.boot0 <- as.data.frame(dat.boot[site.indx, ][sample(1:sum(site.indx), replace = TRUE), ])
      covs <- sim.dat.boot0[, c(1:11)]
      x <- sim.dat.boot0[, c(12)]
      y <- sim.dat.boot0[, c(13)]
      e <- sim.dat.boot0[, c(14)]
      s <- sim.dat.boot0[, c(15)]
      names(covs) <- colnames(Summ.Stat[[site]]$cox.fit.event$x)[-1]
      sim.dat.boot <- rbind(sim.dat.boot, cbind(covs, X = x, Y = y, E = e, site = s))

      site.indx <- NULL
      sim.dat.boot0 <- NULL
    }

    cov.cols <- 1:11

    common.prob.norm[[sim]] <- (t(as.matrix(sim.dat.norm[, cov.cols])) %*% as.matrix(sim.dat.norm[, cov.cols])) / samp.size
    common.prob.chain[[sim]] <- (t(as.matrix(sim.dat.chain[, cov.cols])) %*% as.matrix(sim.dat.chain[, cov.cols])) / samp.size
    common.prob.boot[[sim]] <- (t(as.matrix(sim.dat.boot[, cov.cols])) %*% as.matrix(sim.dat.boot[, cov.cols])) / samp.size

    PS.list.norm <- comp.pscore(X = sim.dat.norm$X, Z = sim.dat.norm[, cov.cols], S = sim.dat.norm$site)
    PS.list.chain <- comp.pscore(X = sim.dat.chain$X, Z = sim.dat.chain[, cov.cols], S = sim.dat.chain$site)
    PS.list.boot <- comp.pscore(X = sim.dat.boot[, 'X'], Z = sim.dat.boot[, cov.cols], S = sim.dat.boot[, 'site'])

    PS.beta.site.norm <- mapply(function(i, j) cbind(i, j), PS.beta.site.norm, PS.list.norm[[4]], SIMPLIFY = FALSE)
    PS.beta.site.chain <- mapply(function(i, j) cbind(i, j), PS.beta.site.chain, PS.list.chain[[4]], SIMPLIFY = FALSE)
    PS.beta.site.boot <- mapply(function(i, j) cbind(i, j), PS.beta.site.boot, PS.list.boot[[4]], SIMPLIFY = FALSE)
    PS.beta.norm <- cbind(PS.beta.norm, PS.list.norm[[2]])
    PS.beta.chain <- cbind(PS.beta.chain, PS.list.chain[[2]])
    PS.beta.boot <- cbind(PS.beta.boot, PS.list.boot[[2]])

    cox.norm <- cox.est(X = sim.dat.norm$X, Y = sim.dat.norm$Y, E = sim.dat.norm$E,
                        Z = sim.dat.norm[, cov.cols], S = sim.dat.norm$site)
    cox.chain <- cox.est(X = sim.dat.chain$X, Y = sim.dat.chain$Y, E = sim.dat.chain$E,
                         Z = sim.dat.chain[, cov.cols], S = sim.dat.chain$site)
    cox.boot <- cox.est(X = sim.dat.boot[, 'X'], Y = sim.dat.boot[, 'Y'], E = sim.dat.boot[, 'E'],
                        Z = sim.dat.boot[, cov.cols], S = sim.dat.boot[, 'site'])

    for (k in 1:5) {
      if (!is.null(cox.norm$Site.Specific[[k]])) {
        cox.site.norm[[k]] <- cbind(cox.site.norm[[k]], cox.norm$Site.Specific[[k]])
      } else {
        cox.site.norm[[k]] <- cbind(cox.site.norm[[k]], rep(NA, times = length(cov.cols) + 1))
      }
      if (!is.null(cox.chain$Site.Specific[[k]])) {
        cox.site.chain[[k]] <- cbind(cox.site.chain[[k]], cox.chain$Site.Specific[[k]])
      } else {
        cox.site.chain[[k]] <- cbind(cox.site.chain[[k]], rep(NA, times = length(cov.cols) + 1))
      }
      if (!is.null(cox.boot$Site.Specific[[k]])) {
        cox.site.boot[[k]] <- cbind(cox.site.boot[[k]], cox.boot$Site.Specific[[k]])
      } else {
        cox.site.boot[[k]] <- cbind(cox.site.boot[[k]], rep(NA, times = length(cov.cols) + 1))
      }
    }

    cox.pooled.norm <- cbind(cox.pooled.norm, cox.norm$Pooled)
    cox.pooled.chain <- cbind(cox.pooled.chain, cox.chain$Pooled)
    cox.pooled.boot <- cbind(cox.pooled.boot, cox.boot$Pooled)
  }

  res[[c]] <- list(
    PS.beta.site.norm = PS.beta.site.norm,
    PS.beta.site.chain = PS.beta.site.chain,
    PS.beta.site.boot = PS.beta.site.boot,
    cox.site.norm = cox.site.norm,
    cox.site.chain = cox.site.chain,
    cox.site.boot = cox.site.boot,
    PS.beta.norm = PS.beta.norm,
    PS.beta.chain = PS.beta.chain,
    PS.beta.boot = PS.beta.boot,
    cox.pooled.norm = cox.pooled.norm,
    cox.pooled.chain = cox.pooled.chain,
    cox.pooled.boot = cox.pooled.boot,
    common.prob.norm = common.prob.norm,
    common.prob.chain = common.prob.chain,
    common.prob.boot = common.prob.boot
  )
}

end <- proc.time()[3]
end - start
#####

saveRDS(res,paste0(resdir,"angio_datasim_cov_bootstarp_survival_241020_2.rds"))
# saveRDS(res,paste0(resdir,"angio_datasim_cov_2016_0525.rds")) #check1

# res <- readRDS(paste0(resdir,"angio_datasim_simple_2016_0525.rds")) # Matthew: where does this data come from? what is the difference between this and the previous data?
#
# res <- readRDS(paste0(resdir,"angio_datasim_simplebump_2016_0525.rds")) # Matthew: use the simple one first


# Plots -----
# Matthew: run below functions first -----
comb.coef <- function(margin=1, m1,m2,m3) {
  if (margin==1) {
    dim.curr <- dim(m1)[[1]]
  } else {
    dim.curr <- dim(m1)[[2]]
  }

  which.coef <- rep(1:dim.curr,time=3)
  which.meth <- rep(1:3,each=dim.curr)

  if (margin==1) {
    return(t(rbind(m1,m2,m3)[order(which.coef,which.meth),]))
  } else {
    return(cbind(m1,m2,m3)[,order(which.coef,which.meth)])
  }
}


interleave <- function(m1,m2){
  ord1 <- 2*(1:dim(m1)[[2]])-1
  ord2 <- 2*(1:dim(m2)[[2]])
  cbind(m1,m2)[,order(c(ord1,ord2))]
}


res <- readRDS(paste0(resdir,"angio_datasim_cov_bootstarp_survival_241020_2.rds"))


carr.norm <- array(NA,dim=c(11,11,sim_num))
carr.chain <- array(NA,dim=c(11,11,sim_num))
carr.boot <- array(NA,dim=c(11,11,sim_num))

cprob.norm <- 0
cprob.chain <- 0
cprob.boot <- 0
for (i in 1:16) { # TODO: may change to the core number
  for (j in 1:sim_num) {
    carr.norm[,,j] <- res[[i]]$common.prob.norm[[j]]
    carr.chain[,,j] <- res[[i]]$common.prob.chain[[j]]
    carr.boot[,,j] <- res[[i]]$common.prob.boot[[j]]
    cprob.norm <- cprob.norm + res[[i]]$common.prob.norm[[j]]
    cprob.chain <- cprob.chain + res[[i]]$common.prob.chain[[j]]
    cprob.boot <- cprob.boot + res[[i]]$common.prob.boot[[j]]
  }
}

round(cprob.norm/5000,6)
round(cprob.chain/5000,6)
round(cprob.boot/5000, 6)

apply(carr.norm,c(1,2),sd)
apply(carr.chain,c(1,2),sd)
apply(carr.boot,c(1,2),sd)


PS.coef.summ <- do.call(cbind,lapply(Summ.Stat,function(XX) XX$coef.XonZ))
rownames(PS.coef.summ) <- c("Intercept","1+ IP Visits","1+ ED Visits",
                            "1+ Comorbid Score","Sex","Age 45-54","Age 55-64",
                            "Age 65+","2009","2010","2011","2012")

Out.coef.summ <- do.call(cbind,lapply(Summ.Stat, # Matthew: TODO
                                      function(XX) -XX$adj.coef.event/XX$adj.scale.event))
Out.coef.summ.cox <- do.call(cbind,lapply(Summ.Stat, # Matthew: TODO
                                          function(XX) XX$cox.coef.adjusted))
rownames(Out.coef.summ) <- c("Intercept","X", "1+ IP Visits","1+ ED Visits", # Matthew: TODO
                             "1+ Comorbid Score","Sex","Age 45-54","Age 55-64",
                             "Age 65+","2009","2010","2011","2012") # Matthew: Add "X" here.

# Table.descr <- do.call("rbind", # Matthew: TODO
#                        lapply(Summ.Stat,
#                               function(XX) {
#                                 temp <- cbind(N=XX$N.X,
#                                               P_Time=round(XX$P.time,1),
#                                               Events=c(XX$control.events,
#                                                        XX$compare.events),
#                                               "Rates"=round(c(1000*365.25*XX$control.rate,
#                                                               1000*365.25*XX$compare.rate),3))
#                                 data.frame(Drug=c('BB','ACEI'), temp)
#                               }))
#
# site.col <- rep("", times=2*length(c("Site 1","Site 2","Site 3","Site 4","Site 5")))
# site.col[c(1,3,5,7,9)] <-  toupper(c("Site 1","Site 2","Site 3","Site 4","Site 5"))
# Table.descr <- data.frame(Site=site.col,Table.descr)
#
# Table.descr_new <- as.data.frame(lapply(Table.descr, as.character),
#                                  stringsAsFactors = FALSE)

Cat.dist <-  lapply(Summ.Stat, function(XX) {
  temp<- as.matrix(unlist(XX$P.ord))
  colnames(temp) <- toupper(XX[[1]])
  temp.frame <- data.frame('id'=trimws(as.character(rownames(temp))),
                           round(temp,2)*100)
  temp.frame[order(trimws(temp.frame$id)),]
})


merged.data.frame <- Reduce(function(...) merge(..., by='id', all=TRUE), Cat.dist)
table.proportions <- merged.data.frame[order(trimws(merged.data.frame$id)),]

# Converted.weib <- lapply(Summ.Stat, function(XX) {  # Matthew: TODO
#   if ("Log(scale)"%in%rownames(XX$adj.vcov.event)) {
#     temp.grad <- c(-1/(XX$adj.scale.event), XX$adj.coef.event[["X"]]/XX$adj.scale.event)
#     temp.vcov <- XX$adj.vcov.event[c("ZX","Log(scale)"),c("ZX","Log(scale)")]
#     vcov.new <- temp.grad %*% temp.vcov %*% temp.grad
#     temp.est <- -XX$adj.coef.event[["X"]]/(XX$adj.scale.event)
#   } else {
#     temp.vcov <- XX$adj.vcov.event[c("ZX"),c("ZX")]
#     vcov.new <- XX$adj.vcov.event[c("ZX"),c("ZX")]
#     temp.est <- -XX$adj.coef.event[["X"]]/(XX$adj.scale.event)
#   }
#   exp(c(temp.est, temp.est +
#           c(qnorm(0.025),qnorm(0.975))*sqrt(vcov.new)))
# })
# # names(Converted.weib) <- site_mask # Matthew: not nessesary?
# Converted.weib <- do.call(rbind, Converted.weib)
#
# Cox.estimates <- lapply(Summ.Stat, function(XX) { # Matthew: TODO
#   vname <- substr(XX[[1]],1,4)
#   temp.vcov <- XX$cox.vcov[1,1]
#   temp.est <- XX$cox.coef.adjusted[[1]]
#
#   exp(c(temp.est, temp.est +
#           c(qnorm(0.025),qnorm(0.975))*sqrt(temp.vcov)))
# })
# merged.cox.coefs <- do.call(rbind, Cox.estimates)

interleave <- function(m1,m2)
{
  ord1 <- 2*(1:dim(m1)[[1]])-1
  ord2 <- 2*(1:dim(m2)[[1]])
  rbind(m1,m2)[order(c(ord1,ord2)),]
}
main.effect <- interleave(Converted.weib,merged.cox.coefs)

main.effect <- main.effect[-c(10:9),]

site_mask <- sites
op <- par(mar=c(5.1,6.1,4.1,1.1),adj=.5)
y.lab.points <- seq(7,1,by=-2) + 0.5
lab.points <- c(7.75,
                7.25,5.75,5.25,3.75,3.25,1.75,1.25)
plot(main.effect[,1],lab.points, axes=FALSE, pch=15,
     col=c("tomato","steelblue"), ylab='',
     xlab = "Estimated HR and 95% CI",
     xlim=c(min(main.effect[,2],na.rm=TRUE)-0.5,
            max(main.effect[,3],na.rm=TRUE)+.5), ylim = c(0,9.5))
axis(side=2, lwd=0,at=y.lab.points, labels=trimws(toupper(site_mask[-5])), las=1,
     cex.axis=0.80, hadj=0.5)
axis(side=1,cex.axis=0.85,at=c(1:12))
cols <- rep(c("tomato","steelblue"),times=9)
abline(v=1,lty=2)
for (i in 8:1) {
  lines(c(main.effect[9-i,2],main.effect[9-i,3]),
        c(lab.points[[9-i]],lab.points[[9-i]]),
        col= cols[[9-i]], lwd=2)
}
points(8.5,7,pch=15,col="tomato")
points(8.5,6.5,pch=15,col="steelblue")
text(8.6,7,"Weibull",pos=4,cex=0.9)
text(8.6,6.5,"Cox",pos=4,cex=0.9)
par(op)


op <- par(mar=c(4,3,3,1))
hist(res[[1]]$cox.pooled.boot[1,],breaks=20)
rbind(apply(res[[4]]$cox.pooled.norm,1,function(X) c(mean(X),sd(X))),
      apply(res[[4]]$cox.pooled.chain,1,function(X) c(mean(X),sd(X))),
      apply(res[[4]]$cox.pooled.boot,1,function(X) c(mean(X),sd(X))))
par(op)

p1 <- t(do.call(cbind, lapply(res, function(l) l$PS.beta.norm)))
p2 <- t(do.call(cbind, lapply(res, function(l) l$PS.beta.chain)))
p3 <- t(do.call(cbind, lapply(res, function(l) l$PS.beta.boot)))

o1 <- t(do.call(cbind, lapply(res, function(l) l$cox.pooled.norm)))
o2 <- t(do.call(cbind, lapply(res, function(l) l$cox.pooled.chain)))
o3 <- t(do.call(cbind, lapply(res, function(l) l$cox.pooled.boot)))

o1 <- t(do.call(cbind, lapply(res, function(l) l$PS.beta.site.norm[[1]])))
o2 <- t(do.call(cbind, lapply(res, function(l) l$PS.beta.site.chain[[1]])))
o3 <- t(do.call(cbind, lapply(res, function(l) l$PS.beta.site.boot[[1]])))

op <- par(mar=c(4,3,3,1),mfrow=c(3,1))
boxplot(p1)
boxplot(p2)
boxplot(p3)
par(op)

op <- par(mar=c(4,3,3,1),mfrow=c(3,1))
boxplot(o1)
boxplot(o2)
boxplot(o3)
par(op)

test1 <- comb.coef(margin=2,p1,p2,p3)
test <- comb.coef(margin=2,o1,o2,o3)

apply(test,2,function(X) c(mean(X),sd(X)))
# Matthew: this is for figure 2

op <- par(mar=c(3,3,3,1))
tiff(paste0(resdir,"pooled_propensity_package_version_1020_2.tif"),width=10,height=7,units='in',res=300)
boxplot(test1, ylim = c(-0.8,0.6), col=c(rgb(1,0,0,0.4),rgb(0,1,0,0.4),rgb(0,0,1,0.4)),axes=FALSE)
# dev.off()
axis(side=2,at=seq(-0.8,0.6,by = 0.2),cex.axis=0.85)
axis(side=1,at=seq(2,47,by=3),labels=FALSE,line= 0.5)
text(seq(2,47,by=3), -0.9, labels = prop.labels, srt = 45, adj = c(1.1,1.1)
     ,xpd = TRUE, cex=0.75)
points(20,-.4,pch=15,col=rgb(1,0,0,0.4))
points(20,-.45,pch=15,col=rgb(0,1,0,0.4))
points(20,-.5,pch=15,col=rgb(0,0,1,0.4))
text(20.5,-.4,"Multivariate Normals",pos=4,cex=0.75)
text(20.5,-.45,"Chains of Regressions",pos=4,cex=0.75)
text(20.5,-.5,"Bootstrap",pos=4,cex=0.75)
dev.off()

########################################
par(op)

o1 <- cbind(res[[1]]$cox.pooled.norm,res[[2]]$cox.pooled.norm)
o2 <- cbind(res[[1]]$cox.pooled.chain,res[[2]]$cox.pooled.chain)
o3 <- cbind(res[[1]]$cox.pooled.boot,res[[2]]$cox.pooled.boot)

test <- comb.coef(margin=1,o1,o2,o3)

apply(test,2,function(X) c(median(X),sd(X)))
op <- par(mar=c(5,4,4,1))

tiff(paste0(resdir,"Cox_coeff_package_version_1020_2.tif"),width=10,height=7,units='in',res=300)
# tiff(paste0(resdir,"Logit_coeff_package_version.tif"),width=10,height=7,units='in',res=300)

boxplot(test,ylim = c(-2,2), col=c(rgb(1,0,0,0.4),rgb(0,1,0,0.4),rgb(0,0,1,0.4)),axes=FALSE)
axis(side=2,cex.axis=0.85)
axis(side=1,at=seq(2,47,by=3),labels=FALSE)
text(seq(2,47, by=3), -2.3, labels = out.labels, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.85)
points(10,-1.3,pch=15,col=rgb(1,0,0,0.4))
points(10,-1.4,pch=15,col=rgb(0,1,0,0.4))
points(10,-1.5,pch=15,col=rgb(0,0,1,0.4))
text(10.5,-1.3,"Multivariate Normals",pos=4,cex=0.7)
text(10.5,-1.4,"Chain of Regressions",pos=4,cex=0.7)
text(10.5,-1.5,"Plasmode Simulation",pos=4,cex=0.7)
dev.off()
par(op)



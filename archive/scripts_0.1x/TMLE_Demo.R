# Project: Targeted Learning for Real World Data
# Author: Rachael Phillips (rachaelvphillips@berkeley.edu), Susan Gruber (sgruber@putnamds.com)
# June 29, 2020
# This file shows how to analyze a dataset with and without missing outcome data
# to obtain estimates of the ATE using an unadjusted estimate and 
# all variants of matching, IPW, and TMLE investigated for Task 1.
# It also shows how to obtain a bootstrap variance estimate.

 ################################################################################
#####################  Task 1 - Demo  #######################
################################################################################

library(sandwich)
library(MatchIt)
library(tmle)
# for super learning, also requires superlearner, glmnet, xgboost, gam, dbarts, but those will be 
# loaded automatically

#####  truncate values of x to lie between min and max specified by bds)
bound <- function(x, bds){
	x[x > max(bds)] <- max(bds)
	x[x < min(bds)] <- min(bds)
	return(x)
}


###############################################################################
# Create dataset with outcome Y, treatment A, and covariates W1, W2, W3
# Treatment is not randomized.
# Also create a second outcome, Y2, where some outcomes are missing
################################################################################

##### Calculate P(Y=1|A,W)
calc_pYA <- function(A, W,intercept = -1.2,  coef_A = -1.2, coefs =  c(.4, .4, -.1)){
  plogis(intercept + coef_A* A +  W %*% coefs)
}

##### One time only, evaluate the true ATE be calculating
### the mean counterfactual outcome under both levels of treatment

# # n <- 10^7
# W1 <- rnorm(n)
# W2 <- rbinom(n, 1, .4)
# W3 <- rbinom(n, 1,  .3)
# Q1 <- calc_pYA(A=rep(1,n), W = cbind(W1, W2, W3))
# Q0 <- calc_pYA(A=rep(0,n), W = cbind(W1, W2, W3))
# psi0 <- mean(Q1-Q0)
# -0.1636357
#########################
psi0 <- -0.1636357

# Simulate data for analysis
set.seed(100)
n <- 1000
W1 <- rnorm(n)
W2 <- rbinom(n, 1, .4)
W3 <- rbinom(n, 1,  .3)
pscore <- plogis(-.7 + .3*W1 - .2*W2 - .3*W3)
A <- rbinom(n, 1, pscore)
pYA <- calc_pYA(A=A, W = cbind(W1, W2, W3))
Y <- rbinom(n, 1, pYA)

# Set 15% of outcomes to missing
# conditional probability outcome is observed (mean = 0.85)
pDelta1 <- plogis(1.67 - .2* A + .5*W3 + .4*W1*W2)
Delta <- rbinom(n, 1, pDelta1)
Y.miss <- Y
Y.miss[Delta == 0] <- NA


############################### Unadjusted #####################################
# In an unadjusted analyses with a binary outcome
# the coefficient in a linear regression of Y on A is exactly equal 
# to the difference in the counterfactual outcomes obtained using a logistic regression model.
# This is true when weights are incorporated into the regression, but NOT when covariates other
# than the treatment are included in the model

m <- lm(Y ~ A)
est.unadj <- coef(m)[2]
se.unadj <- sqrt(vcov(m)[2,2])
# Identical to the following
# m2 <- glm(Y~A, family = "binomial")
# EY1 <- predict(m2, newdata = data.frame(A=1), type = "response")
# EY0 <- predict(m2, newdata = data.frame(A=0), type = "response")
# unadj <- EY1 - EY0

# unadjusted complete case analysis when outcomes are missing
m2 <- lm(Y.miss ~ A)
est.unadj.miss <- coef(m2)[2]
se.unadj.miss <- sqrt(vcov(m2)[2,2])

################################## Matching ####################################
# Full matching on logitPS scale without regard to missingness in the outcome

 pscore_formula <- as.formula("A ~ W1 + W2 + W3")
 matched <- suppressWarnings(MatchIt::matchit(pscore_formula, 
                                                 data = data.frame(A, W1, W2, W3),
                                                 method = "full", 
                                                 replace = TRUE))
 mdat <- MatchIt::match.data(matched)
 m_unadj <- lm(Y ~ A, weights = mdat$weights)
 m_adj <- glm(Y ~ A + W1 + W2 + W3, weights = mdat$weights, family = "quasibinomial")
 EY1 <- predict(m_adj, newdata = data.frame(A=1, W1, W2, W3), type = "response")
 EY0 <- predict(m_adj, newdata = data.frame(A=0, W1, W2, W3), type = "response")
 est.match.unadj <- coef(m_unadj)[2]
 se.match.unadj <- sqrt(vcovHC(m_unadj, type = "HC0")[2,2])
 est.match.adj <- mean(EY1-EY0)
 # analytic variance estimate not available for the adjusted matching estimate of the ATE parameter 

# complete cases analysis after matching, without imputing the missing outcome
 m_unadj.miss <- lm(Y.miss ~ A, weights = mdat$weights, subset = Delta == 1)
 m_adj.miss <- glm(Y.miss ~ A + W1 + W2 + W3, weights = mdat$weights, family = "quasibinomial",
 								subset = Delta == 1)
 EY1.miss <- predict(m_adj, newdata = data.frame(A=1, W1, W2, W3), type = "response")
 EY0.miss <- predict(m_adj, newdata = data.frame(A=0, W1, W2, W3), type = "response")
 est.match.unadj.miss <- coef(m_unadj.miss)[2]
 se.match.unadj.miss <- sqrt(vcovHC(m_unadj.miss, type = "HC0")[2,2])
 est.match.adj.miss <- mean(EY1.miss-EY0.miss)


#################################### IPW #######################################
# estimate PS and P(Delta = 1 | A, W) using main terms logistic regression
# (this is correct for the PS, but misspecified for estimating the missingness probabilities
# investigate impact of bounding the IP weight at 4 different levels.
IPWbounds <- c(1e6, 100, 40, 5/sqrt(n)/log(n))
est.iptw <- est.iptcw <- rep(NA,  length(IPWbounds))
se.iptw <- se.iptcw <- rep(NA,  length(IPWbounds))
m <- glm(A ~ W1 + W2 +W3, family = "binomial")
ps <- predict(m, type = "response")
wt.stab <- mean(A) * A/ ps + mean(1-A) * (1-A)/(1-ps)
m.D <- glm(Delta ~ A + W1 + W2 + W3, family = "binomial")
pDelta <- predict(m.D, type = "response")
wt.stab.miss <- mean(Delta) * Delta / pDelta * wt.stab
for (i in 1:length(IPWbounds)){
	m.iptw <- lm(Y ~ A, weight = bound(wt.stab, c(0,IPWbounds[i])))
	est.iptw[i] <- coef(m.iptw)[2]
	se.iptw[i] <- sqrt(vcovHC(m.iptw, type = "HC0")[2,2])
	# with missing outcomes
	m.iptcw <- lm(Y.miss ~ A, weight = bound(wt.stab.miss, c(0,IPWbounds[i])), subset = Delta == 1)
	est.iptcw[i] <- coef(m.iptcw)[2]
	se.iptcw[i] <- sqrt(vcovHC(m.iptcw, type = "HC0")[2,2])
}

########### TMLE - no prescreening of W before SL estimation of PS #############
# investigate impact of bounding the PS at 4 different levels.
# Four tmle variants
# 1. default SL libraries, no pre-screening of covariates
# 2. default SL libraries, retain covariates associated with Y
# 3. default SL libraries, retain covariates associated with residuals
# 4. demonstrating custom SL libraries for Q, g, Delta
PSbounds <- 1/IPWbounds
est.tmle1 <- est.tmle2 <- est.tmle3 <- est.tmle4 <- rep(NA, length(PSbounds))
se.tmle1 <- se.tmle2 <- se.tmle3 <- se.tmle4 <- rep(NA, length(PSbounds))
CI.tmle1 <- CI.tmle2 <- CI.tmle3 <- CI.tmle4 <- matrix(NA, nrow = 2, ncol = length(PSbounds))

est.tmle1.miss <- est.tmle2.miss <- est.tmle3.miss <- est.tmle4.miss <- rep(NA, length(PSbounds))
se.tmle1.miss <- se.tmle2.miss <- se.tmle3.miss <- se.tmle4.miss <- rep(NA, length(PSbounds))
CI.tmle1.miss <- CI.tmle2.miss <- CI.tmle3.miss <- CI.tmle4.miss <- matrix(NA, nrow = 2, ncol = length(PSbounds))

for(i in 1:length(PSbounds)){
	 tmle1 <- tmle(Y, A, W = cbind(W1, W2, W3),
                       gbound = PSbounds[i], family = "binomial",
                       prescreenW.g = FALSE)
     tmle1.miss <- tmle(Y.miss, A, W = cbind(W1, W2, W3), Delta =Delta,
                       gbound = PSbounds[i], family = "binomial",
                       prescreenW.g = FALSE)
     est.tmle1[i] <- tmle1$est$ATE$psi
     se.tmle1[i] <- sqrt(tmle1$est$ATE$var.psi)
     CI.tmle1[,i] <- tmle1$est$ATE$CI
     est.tmle1.miss[i] <- tmle1.miss$est$ATE$psi
     se.tmle1.miss[i] <- sqrt(tmle1.miss$est$ATE$var.psi)
     CI.tmle1.miss[,i] <- tmle1.miss$est$ATE$CI
     
     # For efficiency, pass in SL-based estimates of Q from previous call to tmle
      tmle2 <- tmle(Y, A, W = cbind(W1, W2, W3), Q = tmle1$Qinit$Q,
                       gbound = PSbounds[i], family = "binomial")
      tmle2.miss <- tmle(Y.miss, A, W = cbind(W1, W2, W3), Delta = Delta, Q = tmle1$Qinit$Q,
                       gbound = PSbounds[i], family = "binomial")
	 est.tmle2[i] <- tmle2$est$ATE$psi
     se.tmle2[i] <- sqrt(tmle2$est$ATE$var.psi)
     CI.tmle2[,i] <- tmle2$est$ATE$CI     
     est.tmle2.miss[i] <- tmle2.miss$est$ATE$psi
     se.tmle2.miss[i] <- sqrt(tmle2.miss$est$ATE$var.psi)
     CI.tmle2.miss[,i] <- tmle2.miss$est$ATE$CI                 
                       
      tmle3 <- tmle(Y, A, W = cbind(W1, W2, W3), Q = tmle1$Qinit$Q,
                       gbound = PSbounds[i], family = "binomial", RESID=TRUE) 
   	  tmle3.miss <- tmle(Y.miss, A, W = cbind(W1, W2, W3), Delta = Delta, Q = tmle1$Qinit$Q,
                       gbound = PSbounds[i], family = "binomial", RESID=TRUE) 
      est.tmle3[i] <- tmle3$est$ATE$psi
      se.tmle3[i] <- sqrt(tmle3$est$ATE$var.psi)
      CI.tmle3[,i] <- tmle3$est$ATE$CI   
      est.tmle3.miss[i] <- tmle3.miss$est$ATE$psi
     se.tmle3.miss[i] <- sqrt(tmle3.miss$est$ATE$var.psi)
     CI.tmle3.miss[,i] <- tmle3.miss$est$ATE$CI
                       
      # example using custom SL libraries                 
      tmle4 <- tmle(Y, A, W = cbind(W1, W2, W3),
                       gbound = PSbounds[i], family = "binomial",
                       Q.SL.library = c("SL.glm", "SL.xgboost"),
                       g.SL.library = c("SL.glm", "SL.glmnet"))
       tmle4.miss <- tmle(Y, A, W = cbind(W1, W2, W3), Delta = Delta,
                       gbound = PSbounds[i], family = "binomial",
                       Q.SL.library = c("SL.glm", "SL.xgboost"),
                       g.SL.library = c("SL.glm", "SL.glmnet"),
                      g.Delta.SL.library =  c("SL.glm", "tmle.SL.dbarts.k.5", "SL.glmnet"))                 
       est.tmle4[i] <- tmle4$est$ATE$psi
       se.tmle4[i] <- sqrt(tmle4$est$ATE$var.psi)
       CI.tmle4[,i] <- tmle4$est$ATE$CI     
       est.tmle4.miss[i] <- tmle4.miss$est$ATE$psi
       se.tmle4.miss[i] <- sqrt(tmle4.miss$est$ATE$var.psi)
       CI.tmle4.miss[,i] <- tmle4.miss$est$ATE$CI              
 }
 
 
################################ IPW - SL weights ##############################
# We'll re-use the SL-based estimates for the PS and PDelta1 from the previous call to TMLE,
# instead of re-estimating them

est.iptw.SL <- est.iptcw.SL <- rep(NA,  length(IPWbounds))
se.iptw.SL <- se.iptcw.SL <- rep(NA,  length(IPWbounds))

ps <- tmle1$g$g1W
pDelta <-  tmle1.miss$g.Delta$g1W[cbind(1:n, A+1)]
wt.stab <- mean(A) * A/ ps + mean(1-A) * (1-A)/(1-ps)
wt.stab.miss <- mean(Delta) * Delta / pDelta * wt.stab
for (i in 1:length(IPWbounds)){
	m.iptw <- lm(Y ~ A, weight = bound(wt.stab, c(0,IPWbounds[i])))
	est.iptw.SL[i] <- coef(m.iptw)[2]
	se.iptw.SL[i] <- sqrt(vcovHC(m.iptw, type = "HC0")[2,2])
	# with missing outcomes
	m.iptcw <- lm(Y.miss ~ A, weight = bound(wt.stab.miss, c(0,IPWbounds[i])), subset = Delta == 1)
	est.iptcw.SL[i] <- coef(m.iptcw)[2]
	se.iptcw.SL[i] <- sqrt(vcovHC(m.iptcw, type = "HC0")[2,2])
}

# Results

res <- cbind(est = c(est.unadj, est.match.unadj, est.match.adj, 
					est.iptw, est.iptw.SL,
					est.tmle1, est.tmle2, est.tmle3, est.tmle4),  
					bias = c(est.unadj, est.match.unadj, est.match.adj, 
					est.iptw, est.iptw.SL, est.tmle1, est.tmle2, est.tmle3, est.tmle4) - psi0,  
					se = c(se.unadj, se.match.unadj, NA, se.iptw, se.iptw.SL,
					se.tmle1, se.tmle2, se.tmle3, se.tmle4))
res.missing <- cbind(est = c(est.unadj.miss, est.match.unadj.miss, est.match.adj.miss, 
					est.iptcw, est.iptcw.SL, est.tmle1.miss, est.tmle2.miss, est.tmle3.miss, est.tmle4.miss),  
					bias = c(est.unadj.miss, est.match.unadj.miss, est.match.adj.miss, 
					est.iptcw, est.iptcw.SL, est.tmle1.miss, est.tmle2.miss, est.tmle3.miss, est.tmle4.miss)-psi0,  
					se = c(se.unadj.miss, se.match.unadj.miss, NA, se.iptcw, se.iptcw.SL,
					se.tmle1.miss, se.tmle2.miss, se.tmle3.miss, se.tmle4.miss))
rownames(res) <- rownames(res.missing) <- c("unadj", "match.unadj", "match.adj",
		"ipw.1mill", "ipw.100", "ipw.40", "ipw.rootn",
		"ipwSL.1mill", "ipwSL.100", "ipwSL.40", "ipwSL.rootn",
		"tmle1.mill", "tmle1.01", "tmle1.025", "tmle1.rootn",
		"tmle2.mill", "tmle2.01", "tmle2.025", "tmle2.rootn",
		"tmle3.mill", "tmle3.01", "tmle3.025", "tmle3.rootn",
		"tmle4.mill", "tmle4.01", "tmle4.025", "tmle4.rootn")
round(res, 3)	
# #                est  bias    se
# unadj       -0.144 0.019 0.028
# match.unadj -0.154 0.010 0.030
# match.adj   -0.147 0.017    NA
# ipw.1mill   -0.159 0.005 0.025
# ipw.100     -0.159 0.005 0.025
# ipw.40      -0.159 0.005 0.025
# ipw.rootn   -0.144 0.019 0.025
# ipwSL.1mill -0.159 0.005 0.025
# ipwSL.100   -0.159 0.005 0.025
# ipwSL.40    -0.159 0.005 0.025
# ipwSL.rootn -0.144 0.019 0.025
# tmle1.mill  -0.158 0.006 0.025
# tmle1.01    -0.159 0.004 0.025
# tmle1.025   -0.157 0.006 0.025
# tmle1.rootn -0.157 0.007 0.013
# tmle2.mill  -0.158 0.006 0.025
# tmle2.01    -0.159 0.004 0.025
# tmle2.025   -0.158 0.006 0.025
# tmle2.rootn -0.157 0.007 0.013
# tmle3.mill  -0.154 0.010 0.026
# tmle3.01    -0.159 0.004 0.024
# tmle3.025   -0.153 0.011 0.025
# tmle3.rootn -0.157 0.007 0.013
# tmle4.mill  -0.159 0.005 0.025
# tmle4.01    -0.158 0.005 0.025
# tmle4.025   -0.158 0.006 0.025
# tmle4.rootn -0.158 0.006 0.013	


round(res.missing, 3)	
               # est  bias    se
# unadj       -0.139 0.025 0.031
# match.unadj -0.144 0.019 0.033
# match.adj   -0.147 0.017    NA
# ipw.1mill   -0.152 0.012 0.028
# ipw.100     -0.152 0.012 0.028
# ipw.40      -0.152 0.012 0.028
# ipw.rootn   -0.139 0.025 0.028
# ipwSL.1mill -0.152 0.012 0.028
# ipwSL.100   -0.152 0.012 0.028
# ipwSL.40    -0.152 0.012 0.028
# ipwSL.rootn -0.139 0.025 0.028
# tmle1.mill  -0.152 0.012 0.027
# tmle1.01    -0.153 0.010 0.027
# tmle1.025   -0.152 0.012 0.027
# tmle1.rootn -0.151 0.012 0.012
# tmle2.mill  -0.153 0.011 0.027
# tmle2.01    -0.155 0.009 0.027
# tmle2.025   -0.152 0.012 0.027
# tmle2.rootn -0.151 0.013 0.012
# tmle3.mill  -0.148 0.016 0.028
# tmle3.01    -0.150 0.013 0.028
# tmle3.025   -0.147 0.017 0.028
# tmle3.rootn -0.151 0.013 0.012
# tmle4.mill  -0.152 0.011 0.027
# tmle4.01    -0.152 0.012 0.027
# tmle4.025   -0.152 0.012 0.027
# tmle4.rootn -0.152 0.012 0.012
	
#!/usr/bin/env Rscript
# Minimal smoke test: pure-R matching + manual TMLE + 2 bootstrap draws
set.seed(42)
n     <- 400
TRUNC <- 0.05
COVS  <- c("age","sex","biomarker","comorbidity","ckd")

age <- rnorm(n,55,10); sex <- rbinom(n,1,.55); bio <- rnorm(n)
com <- sample(0:2,n,TRUE,c(.5,.3,.2)); ckd <- rbinom(n,1,.12)
lp  <- -0.5 + 1.5*(0.03*(age-55)+0.8*sex+0.6*bio+0.5*ckd+0.3*com)
A   <- rbinom(n,1,plogis(lp))
Y   <- rbinom(n,1,plogis(-2.5+0.3*sex+(-0.333)*A))
dat <- data.frame(age=age, sex=sex, biomarker=bio, comorbidity=com,
                  ckd=ckd, treatment=A, event_24=Y)
cat("data generated: n=",n,"treated=",sum(A),"\n")

# ---- IPTW ----
g_mod <- glm(treatment~age+sex+biomarker+comorbidity+ckd, data=dat, family=binomial())
ps    <- pmin(pmax(predict(g_mod, type="response"), TRUNC), 1-TRUNC)
pA    <- mean(A)
w     <- ifelse(A==1, pA/ps, (1-pA)/(1-ps))
mu1   <- weighted.mean(Y[A==1], w[A==1])
mu0   <- weighted.mean(Y[A==0], w[A==0])
iptw_est <- mu1 - mu0
cat("IPTW ATE =", round(iptw_est, 4), "\n")

# ---- Pure-R 1:1 matching ----
logps <- qlogis(ps)
cal   <- 0.2 * sd(logps)
ti    <- which(A==1); ci <- which(A==0)
uc    <- logical(length(ci)); mc <- integer(length(ti))
for (i in seq_along(ti)) {
  d    <- abs(logps[ti[i]] - logps[ci]); d[uc] <- Inf
  best <- which.min(d)
  if (d[best] <= cal) { mc[i] <- ci[best]; uc[best] <- TRUE
  } else mc[i] <- NA_integer_
}
valid <- !is.na(mc)
idx   <- c(ti[valid], mc[valid])
md    <- dat[idx, ]
cat("matched n =", nrow(md), "\n")

# ---- Manual TMLE on matched data ----
Am <- md[["treatment"]]; Ym <- md[["event_24"]]
Wm <- md[, COVS, drop=FALSE]
gm <- pmin(pmax(predict(glm(Am~., data=cbind(Am=Am,Wm), family=binomial()),
                         type="response"), TRUNC), 1-TRUNC)
qm <- glm(Ym~., data=cbind(Ym=Ym,Am=Am,Wm), family=binomial())
Q1 <- predict(qm, newdata=cbind(Am=1,Wm), type="response")
Q0 <- predict(qm, newdata=cbind(Am=0,Wm), type="response")
QAW <- pmin(pmax(Am*Q1+(1-Am)*Q0, 1e-6), 1-1e-6)
H   <- Am/gm - (1-Am)/(1-gm)
eps <- coef(glm(Ym~H+offset(qlogis(QAW)), family=binomial()))["H"]
Q1s <- plogis(qlogis(pmin(pmax(Q1,1e-6),1-1e-6)) + eps/gm)
Q0s <- plogis(qlogis(pmin(pmax(Q0,1e-6),1-1e-6)) - eps/(1-gm))
ate <- mean(Q1s) - mean(Q0s)
eif <- (Am/gm-(1-Am)/(1-gm))*(Ym-QAW) + (Q1s-Q0s) - ate
se  <- sqrt(var(eif)/nrow(md))
cat("Match-TMLE ATE =", round(ate,4), " se =", round(se,4), "\n")

# ---- 2 bootstrap draws ----
for (b in 1:2) {
  idx2 <- sample.int(n, n, replace=TRUE)
  d2   <- dat[idx2, ]
  A2   <- d2[["treatment"]]; W2 <- d2[, COVS, drop=FALSE]
  g2   <- pmin(pmax(predict(glm(A2~., data=cbind(A2=A2,W2), family=binomial()),
                              type="response"), TRUNC), 1-TRUNC)
  lp2  <- qlogis(g2); cal2 <- 0.2*sd(lp2)
  ti2  <- which(A2==1); ci2 <- which(A2==0)
  uc2  <- logical(length(ci2)); mc2 <- integer(length(ti2))
  for (i in seq_along(ti2)) {
    dd   <- abs(lp2[ti2[i]] - lp2[ci2]); dd[uc2] <- Inf
    bb   <- which.min(dd)
    if (dd[bb] <= cal2) { mc2[i] <- ci2[bb]; uc2[bb] <- TRUE
    } else mc2[i] <- NA_integer_
  }
  v2  <- !is.na(mc2)
  md2 <- d2[c(ti2[v2], mc2[v2]), ]
  Am2 <- md2[["treatment"]]; Ym2 <- md2[["event_24"]]; Wm2 <- md2[, COVS, drop=FALSE]
  qm2 <- glm(Ym2~., data=cbind(Ym2=Ym2,Am2=Am2,Wm2), family=binomial())
  b_ate <- mean(predict(qm2, newdata=cbind(Am2=1,Wm2), type="response")) -
           mean(predict(qm2, newdata=cbind(Am2=0,Wm2), type="response"))
  cat("  boot", b, "ATE =", round(b_ate,4), "\n")
}

cat("SMOKE_PASS\n")

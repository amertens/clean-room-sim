# ============================================================
# Reviewer item #2 (SMOKE): does the bootstrap repair the marginal-overlap
# undercoverage of TMLE and cross-fitted TMLE?
# ------------------------------------------------------------
# The main simulation showed plain TMLE / TMLE_CF are anti-conservative under
# marginal overlap (IF-based SE/SD ~ 0.79, coverage ~ 0.87). The existing
# bootstrap-variance study only covered IPTW and Match_TMLE. Here we add TMLE
# and TMLE_CF and compare influence-function vs nonparametric-bootstrap CIs,
# focused on the marginal-overlap scenario.
#
# Usage: Rscript sandbox/variance/var_smoke.R [n_mc] [b_boot] [overlap]
#   defaults: n_mc = 40, b_boot = 150, overlap = 1.5 (marginal)
# ============================================================
args    <- commandArgs(trailingOnly = TRUE)
N_MC    <- if (length(args) >= 1) as.integer(args[[1]]) else 40L
B_BOOT  <- if (length(args) >= 2) as.integer(args[[2]]) else 150L
OVERLAP <- if (length(args) >= 3) as.numeric(args[[3]]) else 1.5
N_OBS   <- 2000L; TRUNC <- 0.05; SEED <- 2026L
COVARS  <- c("age","sex","biomarker","comorbidity","ckd")
set.seed(SEED)

generate_data <- function(n, overlap_strength) {
  age <- rnorm(n,55,10); sex <- rbinom(n,1,.55); biomarker <- rnorm(n)
  comorbidity <- sample(0:2,n,TRUE,c(.5,.3,.2)); ckd <- rbinom(n,1,.12)
  lp_trt <- -0.5 + overlap_strength*(0.03*(age-55)+0.8*sex+0.6*biomarker+0.5*ckd+0.3*comorbidity)
  treatment <- rbinom(n,1,plogis(lp_trt))
  lp_out <- -2.5 + 0.015*(age-55)+0.3*sex+0.2*biomarker+0.6*ckd+0.25*comorbidity+(-0.05/0.15)*treatment
  data.frame(age=age,sex=sex,biomarker=biomarker,comorbidity=comorbidity,ckd=ckd,
             treatment=treatment, event_24=rbinom(n,1,plogis(lp_out)))
}
compute_truth <- function(n) {
  set.seed(SEED); age<-rnorm(n,55,10); sex<-rbinom(n,1,.55); bio<-rnorm(n)
  com<-sample(0:2,n,T,c(.5,.3,.2)); ckd<-rbinom(n,1,.12)
  lp<- -2.5+0.015*(age-55)+0.3*sex+0.2*bio+0.6*ckd+0.25*com
  mean(plogis(lp+(-0.05/0.15)))-mean(plogis(lp))
}

.target <- function(Y,A,g,Q1,Q0) {
  Q_AW <- pmin(pmax(A*Q1+(1-A)*Q0,1e-6),1-1e-6)
  H <- A/g-(1-A)/(1-g)
  eps <- tryCatch(coef(glm(Y~H+offset(qlogis(Q_AW)),family=binomial()))["H"],
                  error=function(e) 0)
  Q1s <- plogis(qlogis(pmin(pmax(Q1,1e-6),1-1e-6))+eps/g)
  Q0s <- plogis(qlogis(pmin(pmax(Q0,1e-6),1-1e-6))-eps/(1-g))
  ate <- mean(Q1s)-mean(Q0s)
  eif <- H*(Y-Q_AW)+(Q1s-Q0s)-ate
  list(estimate=ate, se=sqrt(var(eif)/length(Y)))
}
# Plain TMLE (GLM nuisances, full sample)
.tmle_fit <- function(dat) {
  Y<-dat$event_24; A<-dat$treatment; W<-dat[,COVARS,drop=FALSE]
  gm<-glm(A~., data=cbind(A=A,W), family=binomial())
  g <- pmin(pmax(predict(gm,type="response"),TRUNC),1-TRUNC)
  qm<-glm(Y~., data=cbind(Y=Y,A=A,W), family=binomial())
  Q1<-predict(qm,newdata=cbind(A=1,W),type="response"); Q0<-predict(qm,newdata=cbind(A=0,W),type="response")
  r<-.target(Y,A,g,Q1,Q0); list(estimate=r$estimate,se=r$se,
    ci_lower=r$estimate-1.96*r$se, ci_upper=r$estimate+1.96*r$se)
}
# Cross-fitted TMLE (2-fold, GLM nuisances, pooled targeting)
.tmle_cf_fit <- function(dat, V=2L) {
  Y<-dat$event_24; A<-dat$treatment; W<-dat[,COVARS,drop=FALSE]; n<-length(Y)
  folds<-sample(rep(seq_len(V),length.out=n))
  g<-Q1<-Q0<-numeric(n)
  for (v in seq_len(V)) {
    tr<-folds!=v; te<-folds==v
    gm<-glm(A~.,data=cbind(A=A,W)[tr,,drop=FALSE],family=binomial())
    g[te]<-predict(gm,newdata=W[te,,drop=FALSE],type="response")
    qm<-glm(Y~.,data=cbind(Y=Y,A=A,W)[tr,,drop=FALSE],family=binomial())
    Q1[te]<-predict(qm,newdata=cbind(A=1,W[te,,drop=FALSE]),type="response")
    Q0[te]<-predict(qm,newdata=cbind(A=0,W[te,,drop=FALSE]),type="response")
  }
  g<-pmin(pmax(g,TRUNC),1-TRUNC)
  r<-.target(Y,A,g,Q1,Q0); list(estimate=r$estimate,se=r$se,
    ci_lower=r$estimate-1.96*r$se, ci_upper=r$estimate+1.96*r$se)
}
.boot_se <- function(dat, fn, B=B_BOOT, seed=1L) {
  set.seed(seed); n<-nrow(dat)
  bs<-vapply(seq_len(B),function(b){
    idx<-sample.int(n,n,replace=TRUE)
    tryCatch(fn(dat[idx,,drop=FALSE])$estimate,error=function(e) NA_real_)},numeric(1))
  bs<-bs[is.finite(bs)]; if(length(bs)<10) return(list(se=NA,ci=c(NA,NA)))
  list(se=sd(bs), ci=unname(quantile(bs,c(.025,.975))))
}

truth <- compute_truth(500000L)
cat(sprintf("True RD = %.5f  | overlap=%.2f  N_MC=%d  B=%d\n",truth,OVERLAP,N_MC,B_BOOT))

ests <- c("TMLE","TMLE_CF"); fns <- list(TMLE=.tmle_fit, TMLE_CF=.tmle_cf_fit)
rows <- list()
for (i in seq_len(N_MC)) {
  dat <- generate_data(N_OBS, OVERLAP)
  for (m in ests) {
    f  <- tryCatch(fns[[m]](dat), error=function(e) NULL); if (is.null(f)) next
    bt <- tryCatch(.boot_se(dat, fns[[m]], seed=i), error=function(e) NULL)
    rows[[length(rows)+1]] <- data.frame(rep=i, method=m,
      estimate=f$estimate, se_if=f$se,
      se_boot=if(!is.null(bt)) bt$se else NA,
      cov_if = as.integer(f$ci_lower<=truth & truth<=f$ci_upper),
      cov_boot = if(!is.null(bt)&&all(is.finite(bt$ci))) as.integer(bt$ci[1]<=truth & truth<=bt$ci[2]) else NA)
  }
  if (i %% 10 == 0) { cat(sprintf("  rep %d/%d\n",i,N_MC)); flush(stdout()) }
}
R <- do.call(rbind, rows)
summ <- do.call(rbind, lapply(ests, function(m){
  d<-R[R$method==m,]; emp<-sd(d$estimate)
  data.frame(method=m, n=nrow(d), bias=round(mean(d$estimate)-truth,5),
    emp_sd=round(emp,4),
    se_if=round(mean(d$se_if),4), ratio_if=round(mean(d$se_if)/emp,3),
    cov_if=round(mean(d$cov_if),3),
    se_boot=round(mean(d$se_boot,na.rm=TRUE),4),
    ratio_boot=round(mean(d$se_boot,na.rm=TRUE)/emp,3),
    cov_boot=round(mean(d$cov_boot,na.rm=TRUE),3))
}))
cat("\n=== Variance smoke (marginal overlap): IF vs bootstrap ===\n")
print(summ, row.names=FALSE)
cat("\nVAR_SMOKE_DONE\n")

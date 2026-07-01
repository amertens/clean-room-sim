# Quick integration check for Scenario D: parse run_simulation.R, and verify the
# misspecified DGP's truth + that crude is biased (confounding present).
ok <- tryCatch({ invisible(parse("run_simulation.R")); TRUE },
               error = function(e) { cat("PARSE ERROR:", conditionMessage(e), "\n"); FALSE })
cat("run_simulation.R parse:", if (ok) "OK" else "FAILED", "\n")

# Replicate the integrated misspec DGP (must match generate_data/compute_truth).
gen <- function(n, eff = -0.05, seed = 1) {
  set.seed(seed)
  age<-rnorm(n,55,10); sex<-rbinom(n,1,.55); bio<-rnorm(n); com<-sample(0:2,n,T,c(.5,.3,.2)); ckd<-rbinom(n,1,.12)
  a<-(age-55)/10
  ps<-plogis(-0.1+0.5*a+1.0*sex*bio+0.5*ckd+0.2*com); A<-rbinom(n,1,ps)
  lpy<- -0.6+0.4*a+0.7*a^2+1.0*sex*bio+0.6*ckd+0.3*com + eff/0.15*A*(1+0.4*sex)
  list(Y=rbinom(n,1,plogis(lpy)), A=A, ps=ps)
}
big<-3e5; set.seed(7)
age<-rnorm(big,55,10); sex<-rbinom(big,1,.55); bio<-rnorm(big); com<-sample(0:2,big,T,c(.5,.3,.2)); ckd<-rbinom(big,1,.12)
a<-(age-55)/10; base<- -0.6+0.4*a+0.7*a^2+1.0*sex*bio+0.6*ckd+0.3*com
truth<-mean(plogis(base+(-0.05/0.15)*(1+0.4*sex)))-mean(plogis(base))
d<-gen(2000, seed=42)
crude<-mean(d$Y[d$A==1])-mean(d$Y[d$A==0])
cat(sprintf("Scenario D truth RD = %.5f\n", truth))
cat(sprintf("PS range [%.3f, %.3f]  frac>.95=%.1f%%\n", min(d$ps), max(d$ps), 100*mean(d$ps>.95)))
cat(sprintf("crude RD (n=2000) = %.5f  -> crude bias = %+.5f (confounding present if large)\n",
            crude, crude-truth))
cat("SCENARIO_D_CHECK_DONE\n")

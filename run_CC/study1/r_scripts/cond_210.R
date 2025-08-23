
rm(list=ls())

library(rstan)  # for model fitting
source('fx.R')  # for simulation

options(mc.cores=2) # 2 chains per model


# Specify constants (to be replaced) --------------------------------------
COND  = 210
CONST = c(70,10,1,0.16,0.4,0.16,0.2,0,0,0,0.1,0,0,0,0,0,1,1,1,100,2,3000,1000,1,1,0.48)
names(CONST) = c("N", "J", "G00", "TAU0", "G10", "TAU1", "G20", "TAU2", 
                 "G30", "TAU3", "RHO01", "RHO02", "RHO12", "X1_mu", 
                 "X2_mu", "X3_mu", "X1_sd", "X2_sd", "X3_sd", "S", 
                 "nchain","niter","burn","thin","marg","SIG2")
CONST = as.data.frame(t(CONST))

  
# Compile models ----------------------------------------------------
# Fit initial model to dummy data, then just call update below to avoid recompiling time.
# This is OS-specific, so has to be run once on the cluster before each simulation cell.
# Fit five models:
# A. True model
# B. The (incorrect) absence of a random effect (TAU1 = 0)
# C. The (incorrect) presence of a random effect (TAU2 > 0) 
# D. The (incorrect) absence of a fixed effect (G20 = 0)
# E. The (incorrect) presence of a fixed effect (G30 != 0)

formA = y ~ 1 + X1 + X2 + (X1|g)   
formB = y ~ 1 + X1 + X2 + (1|g)
formC = y ~ 1 + X1 + X2 + (X1+X2|g)
formD = y ~ 1 + X1 + (X1|g)   
formE = y ~ 1 + X1 + X2 + X3 + (X1|g)

DM    = expand.grid(1:10,1:10) # design matrix (dummy N & J)
DM$X1 = rnorm(nrow(DM)) # mean and sd don't matter here, just setting up model structure
DM$X2 = rnorm(nrow(DM))
DM$X3 = rnorm(nrow(DM))
DM$y  = rnorm(nrow(DM)) # same for the outcome
colnames(DM) = c('obs','g','X1','X2','X3', 'y')

options(mc.cores=CONST$nchain)

modA_0 = stan('mlmmarg.stan',data=mod.mat.mlm(formA,DM),chains=CONST$nchain, iter=CONST$niter, warmup=CONST$burn, thin=CONST$thin)
modB_0 = stan('mlmmarg.stan',data=mod.mat.mlm(formB,DM),chains=CONST$nchain, iter=CONST$niter, warmup=CONST$burn, thin=CONST$thin)
modC_0 = stan('mlmmarg.stan',data=mod.mat.mlm(formC,DM),chains=CONST$nchain, iter=CONST$niter, warmup=CONST$burn, thin=CONST$thin)
modD_0 = stan('mlmmarg.stan',data=mod.mat.mlm(formD,DM),chains=CONST$nchain, iter=CONST$niter, warmup=CONST$burn, thin=CONST$thin)
modE_0 = stan('mlmmarg.stan',data=mod.mat.mlm(formE,DM),chains=CONST$nchain, iter=CONST$niter, warmup=CONST$burn, thin=CONST$thin)

# Generate dataset --------------------------------------------------------

# Simulate from true model
g  = rep(1:CONST$N, each=CONST$J)
X1 = rnorm(n=CONST$N*CONST$J, mean=CONST$X1_mu, sd=CONST$X1_sd)
X2 = rnorm(n=CONST$N*CONST$J, mean=CONST$X2_mu, sd=CONST$X2_sd)
X3 = rnorm(n=CONST$N*CONST$J, mean=CONST$X3_mu, sd=CONST$X3_sd)
dm = data.frame(X1=X1,X2=X2,X3=X3,g=g)

# Store pars in list
pars = list(GAMMA=c(CONST$G00, CONST$G10, CONST$G20), 
            TAU=matrix(c(CONST$TAU0, 
                         CONST$RHO01*sqrt(CONST$TAU0)*sqrt(CONST$TAU1), 
                         CONST$RHO01*sqrt(CONST$TAU0)*sqrt(CONST$TAU1), 
                         CONST$TAU1), 2, 2), 
            SIGMA=CONST$SIG)
names(pars$GAMMA) = c('(Intercept)','X1','X2')

# simulate from true model S times
y = simulate_lmer(y ~ X1 + X2 + (X1|g), DM=dm, pars=pars, nsim = CONST$S)

# Fit each model to each dataset and compute fit indices --------------------------------------------------

out   = list() # output to save

for(s in 1:CONST$S) {
  
  cat('\n=========================fitting sample ', s,'=========================\n')
  
  dat = dm
  dat$y = y[,s]
  
  # fit
  cat('fitting models...\n')
  suppressMessages(
    {
      modA = stan('mlmmarg.stan',data=mod.mat.mlm(formA,dat),chains=CONST$nchain,iter=CONST$niter, warmup=CONST$burn, thin=CONST$thin,fit=modA_0)
      modB = stan('mlmmarg.stan',data=mod.mat.mlm(formB,dat),chains=CONST$nchain,iter=CONST$niter, warmup=CONST$burn, thin=CONST$thin,fit=modB_0)
      modC = stan('mlmmarg.stan',data=mod.mat.mlm(formC,dat),chains=CONST$nchain,iter=CONST$niter, warmup=CONST$burn, thin=CONST$thin,fit=modC_0)
      modD = stan('mlmmarg.stan',data=mod.mat.mlm(formD,dat),chains=CONST$nchain,iter=CONST$niter, warmup=CONST$burn, thin=CONST$thin,fit=modD_0)
      modE = stan('mlmmarg.stan',data=mod.mat.mlm(formE,dat),chains=CONST$nchain,iter=CONST$niter, warmup=CONST$burn, thin=CONST$thin,fit=modE_0)
    }
  )
  
  # Compute fit indices for each model 
  cat('computing fit indices...\n')
  
  suppressWarnings(
    {
      dicA  = DIC_marg(modA, mod.mat.mlm(formA,dat))
      waicA = waic(extract_log_lik(modA))$estimates
      looA  = loo(extract_log_lik(modA))$estimates
      
      dicB  = DIC_marg(modB, mod.mat.mlm(formB,dat))
      waicB = waic(extract_log_lik(modB))$estimates
      looB  = loo(extract_log_lik(modB))$estimates
      
      dicC  = DIC_marg(modC, mod.mat.mlm(formC,dat))
      waicC = waic(extract_log_lik(modC))$estimates
      looC  = loo(extract_log_lik(modC))$estimates
      
      dicD  = DIC_marg(modD, mod.mat.mlm(formD,dat))
      waicD = waic(extract_log_lik(modD))$estimates
      looD  = loo(extract_log_lik(modD))$estimates
      
      dicE  = DIC_marg(modE, mod.mat.mlm(formE,dat))
      waicE = waic(extract_log_lik(modE))$estimates
      looE  = loo(extract_log_lik(modE))$estimates
    }
  )
  
  # Save output
  # (saving the whole models would take too much memory)
  # (same with the pll for now)
  cat('storing output...\n')
  
  out[[s]] = list(cond=COND, iter=s, marg=CONST$marg,
             data  = dm,
             # pll = list(A=extract_log_lik(modA),B=extract_log_lik(modB),C=extract_log_lik(modC),D=extract_log_lik(modD),E=extract_log_lik(modE)),
             fixef = list(A=fixef_marg(modA),B=fixef_marg(modB),C=fixef_marg(modC),D=fixef_marg(modD),E=fixef_marg(modE)),
             ranef = list(A=ranef_marg(modA),B=ranef_marg(modB),C=ranef_marg(modC),D=ranef_marg(modD),E=ranef_marg(modE)),
             diag  = list(A=diag_marg(modA),B=diag_marg(modB),C=diag_marg(modC),D=diag_marg(modD),E=diag_marg(modE)),
             dic   = list(A=dicA,B=dicB,C=dicC,D=dicD,E=dicE),
             waic  = list(A=waicA,B=waicB,C=waicC,D=waicD,E=waicE), 
             loo   = list(A=looA,B=looB,C=looC,D=looD,E=looE)
             )
  
} # end of sample loop 

# save output to disk again at end
fname = paste0('/home/sdevine/scratch/cond_',COND,'.rds')
saveRDS(out, fname)

# Sys.sleep(500) # give ample time to save to avoid data loss



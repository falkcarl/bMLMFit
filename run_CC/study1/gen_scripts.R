
rm(list=ls()) # clear workspace

source('fx.R') # for error variance computation

# Specify constants -----------------------------------------------------

set.seed(2023)

# Set constants
J      = c(10,50,100) # cluster sizes to try (level 2)
N      = c(20,50,70)  # observation levels to try (level 1)
G00    = 1            # fixed intercept
TAU0   = c(.04,.16)   # random intercept variance
G10    = c(.2, .4)    # fixed slope of X1
TAU1   = c(.04,.16)   # random slope of X1
G20    = c(.2,.4)     # fixed slope of X2
TAU2   = 0            # random slope of X2
G30    = 0            # fixed slope of X3
TAU3   = 0            # random slope of X3
RHO01  = .1           # corr between TAU0 and TAU1
RHO02  = 0            # corr between TAU0 and TAU2
RHO12  = 0            # corr between TAU1 and TAU2
X1_mu  = 0            # mean for X1
X2_mu  = 0            # mean for X2
X3_mu  = 0            # mean for X3
X1_sd  = 1            # std. dev. for X1
X2_sd  = 1            # std. dev. for X2
X3_sd  = 1            # std. dev. for X3
S      = 100          # number of iterations per sample per model
nchain = 2            # number of chains per MCMC sample per model
niter  = 3000         # number of MCMC samples per model
burn   = 1000         # number of burn-in MCMC samples per model
thin   = 1            # amount of thinning per MCMC sample per model
marg   = c(0,1)       # marginal or conditional estimation method (1 = marginal, 0 = conditional)

# specify and save experimental conditions
conds = expand.grid(N=N, J=J,G00=G00,TAU0=TAU0,G10=G10,
                    TAU1=TAU1,G20=G20,TAU2=TAU2,G30=G30,
                    TAU3=TAU3,RHO01=RHO01,RHO02=RHO02,RHO12=RHO12,
                    X1_mu=X1_mu,X2_mu=X2_mu,X3_mu=X3_mu,X1_sd=X1_sd,
                    X2_sd=X2_sd,X3_sd=X3_sd,S=S,
                    nchain=nchain,niter=niter,burn=burn,thin=thin,
                    marg=marg)

# add error variance, computed per row such that total variance = 1
conds$SIG2 = 
  sapply(split(conds, seq(nrow(conds))), function(x) {
  get_vareY(x$G00, x$G10, x$G20, x$X1_sd, x$X2_sd, 
            x$RHO12, x$TAU0, x$TAU1, x$RHO01)
  })
write.csv(conds, 'conds.csv')

# compute proportion of variance explained... 
prop_var =
  sapply(split(conds, seq(nrow(conds))), function(x) {
    get_vareY(x$G00, x$G10, x$G20, x$X1_sd, x$X2_sd,
              x$RHO12, x$TAU0, x$TAU1, x$RHO01, return_prop_var = T)
  })

tmp = cbind(conds, as.data.frame(t(prop_var)))
tmp = apply(tmp,2,as.numeric)
write.csv(tmp, 'conds_with_var.csv')


# Generate R scripts ------------------------------------------------------

#' Creates an R script for each condition
#' This is done from a template script, `template_r_script.R`

# Load blank template bash script
Bscript = paste(readLines("template_bash_script.sh"), collapse="\n")

# Specify which fits have already finished, if any, to skip these
fit_conds = list.files('out/')
fit_conds = as.vector(sapply(fit_conds, function(x) 
  as.numeric(gsub('.rds', '', strsplit(x,'_')[[1]][2]))))

# Iterate over conditions
for(i in seq(nrow(conds))) {

  # if(i %in% fit_conds) next 
  
  # load template R script
  Rscript = ifelse(conds$marg[i], 
                   paste(readLines("template_marginal.R"), collapse="\n"),
                   paste(readLines("template_conditional.R"), collapse="\n")
                   )
  
  # replace constants with proper constants for this iteration
  new_r_script = gsub('__CONST__', paste(conds[i,],collapse = ','), Rscript)
  new_r_script = gsub('__COND__', i, new_r_script)

  # write new R script to file
  r_fname = paste0('r_scripts/cond_',i,'.R')
  cat(new_r_script,file=r_fname,sep='')
  
  # replace __JOBNAME__ and __RFILE__ in bash script
  new_sh_script = gsub('__JOBNAME__', paste0('bMLMfit_cond_',i), Bscript)
  new_sh_script = gsub('__RFILE__', r_fname, new_sh_script)
  
  if((i-1)%%10==0) {
    # email every 10 conditions
    new_sh_script = gsub('__MAIL__',  'ALL', new_sh_script)
  } else {
    new_sh_script = gsub('__MAIL__',  'NONE', new_sh_script)
  }
  
  # write bash to file
  cat(new_sh_script,file=paste0('bash_scripts/cond_',i,'.sh'),sep='')
  
}


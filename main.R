rm(list=ls())

library(reshape2) # for wide -> long data transform
alpha   = scales::alpha
tab_mod = sjPlot::tab_model
ci      = function(x,a=qnorm(.975)) sd(x)/sqrt(length(x)) * a 

# Set constants -----------------------------------------------------------

TRUE_MOD   = 1 # i.e., "model A"
FIT_ICS    = c('DIC','WAIC','LOO')
IVS        = c('N','J','TAU0','TAU1','G10','G20')
MODELS     = 1:5
MODEL_COLS = alpha(RColorBrewer::brewer.pal(length(MODELS), 'Dark2'),.5) # for plotting

# Load summary data file --------------------------------------------------

dat = readRDS('data/fits_summary.rds')

# Summarize diagnostics ---------------------------------------------------
## Rhats ------------------------------------------------------------------

rhats    = sapply(dat$max_rhats, function(x) x)
rhat_p   = mean(rhats > 1.05)

pdf('figs/rhats_hists.pdf', 6, 6)

hist(rhats[1,], col = MODEL_COLS[1], xlim=c(.99,1.10), 
     xlab=expression(hat(R)), main='', ylim=c(0,12500))
for(i in 2:nrow(rhats)) hist(rhats[i,], col=MODEL_COLS[i], add=T)
abline(v=1.05, lty=2)

legend('topright', bty='n', fill=MODEL_COLS, legend=LETTERS[1:nrow(rhats)], title='Model')
legend('right',bty='n', legend=paste0('P(>1.05) = ',round(rhat_p,4)))

dev.off()



# ESS ----------------------------------------------------------------------

ess   = sapply(dat$min_ess, function(x) x)
ess_p = mean(ess < 100)

pdf('figs/ess_hists.pdf', 6, 6)

hist(ess[1,], col = MODEL_COLS[1], breaks = 5,
     xlab = 'ESS', main='', ylim=c(0,10000))
for(i in 2:nrow(ess)) hist(ess[i,], col=MODEL_COLS[i], add=T, breaks = 5)
abline(v=100, lty=2)

legend('topright', bty='n', fill=MODEL_COLS, legend=LETTERS[1:nrow(ess)], title='Model')
legend('right',bty='n', legend=paste0('P(<100) = ',round(ess_p,4)))

dev.off()


# IC Performance -------------------------------------------------------------

## Model accurate model selection -----------

# compute accuracy
for(i in FIT_ICS){
  dat[,paste0('acc_',i)] = as.numeric(dat[,paste0('win_',i)]==TRUE_MOD)
}

# melt df
melt_dat = melt(dat, id.vars = IVS,
                 measure.vars = paste0('acc_',FIT_ICS), 
                 variable.name = 'ic',
                 value.name = 'accuracy')
melt_dat$ic = sapply(melt_dat$ic, function(x) 
  strsplit(as.character(x), split = '_')[[1]][2])

# make fit index a factor
melt_dat$ic = factor(melt_dat$ic, levels = c('DIC','WAIC','LOO'))

# center predictors
melt_dat[,IVS] = apply(melt_dat[,IVS],2,function(x) x - median(x))

# fit
acc_mod = glm(accuracy ~ ic*N*J*TAU0*TAU1*G10*G20,
              data=melt_dat,
              family='binomial')

capture.output(summary(acc_mod), file = 'out/acc_model_summary.txt')



## Visualize performance ---------------------------------------------------

#' create indexable array with proportions of model selection per: 
#'    1. Fit index
#'    2. Model selected
#'    3. N
#'    4. J
#'    5. Tau0
#'    6. Tau1
#'    7. gamma10
#'    8. gamma20


mod_select = array(NA,dim = c(length(FIT_ICS), 
                              length(MODELS), 
                              length(unique(dat$N)), 
                              length(unique(dat$J)), 
                              length(unique(dat$TAU0)),
                              length(unique(dat$TAU1)),
                              length(unique(dat$G10)),
                              length(unique(dat$G20))),
                   dimnames = list(FIT_ICS,
                                   MODELS,
                                   unique(dat$N),
                                   unique(dat$J),
                                   unique(dat$TAU0),
                                   unique(dat$TAU1),
                                   unique(dat$G10),
                                   unique(dat$G20))
)

for(i in FIT_ICS) {
  
  win_vec = dat[,paste0('win_',i)]
  
  for(m in MODELS) {
    mod_select[i,m,,,,,,] = tapply(win_vec, 
                                   list(dat$N,
                                        dat$J, 
                                        dat$TAU0, 
                                        dat$TAU1,
                                        dat$G10, 
                                        dat$G20), 
                                   function(x) mean(x==m))
    
  }
}

## Accurate model selection ------------------------------------------------

### ICs --------
pdf('figs/pcorrect_ic_barplot.pdf', 6, 6)

pcor = apply(mod_select[,TRUE_MOD,,,,,,], 1, mean)
ecor = apply(mod_select[,TRUE_MOD,,,,,,], 1, ci)

bp = barplot(pcor, ylim=c(0,1), ylab='P(Correct)')
arrows(bp, pcor-ecor, bp, pcor+ecor,length=0)

dev.off()

### ICs X N --------
pdf('figs/pcorrect_ic_n_barplot.pdf', 6, 6)

pcor = apply(mod_select[,TRUE_MOD,,,,,,], c(2,1), mean)
ecor = apply(mod_select[,TRUE_MOD,,,,,,], c(2,1), ci)

bp = barplot(pcor, beside=T,ylim=c(0,1), ylab='P(Correct)', 
             legend.text = T,
             args.legend = list(x='topleft',bty='n',
                                title='Cluster Size'))
arrows(bp, pcor-ecor, bp, pcor+ecor,length=0)

dev.off()

### ICs X J --------
pdf('figs/pcorrect_ic_j_barplot.pdf', 6, 6)

pcor = apply(mod_select[,TRUE_MOD,,,,,,], c(3,1), mean)
ecor = apply(mod_select[,TRUE_MOD,,,,,,], c(3,1), ci)

bp = barplot(pcor, beside=T,ylim=c(0,1), ylab='P(Correct)', 
             legend.text = T,
             args.legend = list(x='topleft',bty='n',
                                title='Observation Size'))
arrows(bp, pcor-ecor, bp, pcor+ecor,length=0)

dev.off()

### ICs X N X J --------
pdf('figs/pcorrect_ic_n_j_barplot.pdf', 8, 4)

pcor = apply(mod_select[,TRUE_MOD,,,,,,], c(3,2,1), mean)
ecor = apply(mod_select[,TRUE_MOD,,,,,,], c(3,2,1), ci)

layout(matrix(1:3,1,3))

for(i in 1:3) {
  bp = barplot(pcor[,,i], beside=T,ylim=c(0,1), ylab='P(Correct)', 
               legend.text = T, main = dimnames(pcor)[[3]][i],
               args.legend = list(x='topleft',bty='n',
                                  title='Observation Size'))
  arrows(bp, pcor[,,i]-ecor[,,i], bp, pcor[,,i]+ecor[,,i],length=0)
  
}

dev.off()


### ICs X N X J X Tau0 --------

pcor = apply(mod_select[,TRUE_MOD,,,,,,], c(3,2,4,1), mean)
ecor = apply(mod_select[,TRUE_MOD,,,,,,], c(3,2,4,1), ci)


for(i in 1:3) {
  pdf(paste0('figs/pcorrect_ic_n_j_tau0_',FIT_ICS[i],'_barplot.pdf'), 8, 4)
  layout(matrix(1:2,1,2))
  
  for(j in 1:2) {
    bp = barplot(pcor[,,j,i], beside=T,ylim=c(0,1), 
                 xlab = 'Cluster Size', ylab='P(Correct)', 
                 main = bquote(tau[0] ~ ' = ' ~ .(dimnames(pcor)[[3]][j])),
                 legend.text = T, 
                 args.legend = list(x='topleft',bty='n',
                                    title='Observation Size',
                                    ncol=3))
    arrows(bp, pcor[,,j,i]-ecor[,,j,i], bp, pcor[,,j,i]+ecor[,,j,i],length=0)
    
  } # end of tau0 loop
  dev.off()
  
} # end of fit index

### ICs X N X J X Tau1 --------

pcor = apply(mod_select[,TRUE_MOD,,,,,,], c(3,2,5,1), mean)
ecor = apply(mod_select[,TRUE_MOD,,,,,,], c(3,2,5,1), ci)

for(i in 1:3) {
  pdf(paste0('figs/pcorrect_ic_n_j_tau1_',FIT_ICS[i],'_barplot.pdf'), 8, 4)
  layout(matrix(1:2,1,2))
  
  for(j in 1:2) {
    bp = barplot(pcor[,,j,i], beside=T,ylim=c(0,1), 
                 xlab = 'Cluster Size', ylab='P(Correct)', 
                 main = bquote(tau[1] ~ ' = ' ~ .(dimnames(pcor)[[3]][j])),
                 legend.text = T, 
                 args.legend = list(x='topleft',bty='n',
                                    title='Observation Size',
                                    ncol=3))
    arrows(bp, pcor[,,j,i]-ecor[,,j,i], bp, pcor[,,j,i]+ecor[,,j,i],length=0)
    
  } # end of tau0 loop
  dev.off()
  
} # end of fit index


### ICs X N X J X gam1 --------

pcor = apply(mod_select[,TRUE_MOD,,,,,,], c(3,2,6,1), mean)
ecor = apply(mod_select[,TRUE_MOD,,,,,,], c(3,2,6,1), ci)

for(i in 1:3) {
  pdf(paste0('figs/pcorrect_ic_n_j_gam1_',FIT_ICS[i],'_barplot.pdf'), 8, 4)
  layout(matrix(1:2,1,2))
  
  for(j in 1:2) {
    bp = barplot(pcor[,,j,i], beside=T,ylim=c(0,1), 
                 xlab = 'Cluster Size', ylab='P(Correct)', 
                 main = bquote(gamma[1] ~ ' = ' ~ .(dimnames(pcor)[[3]][j])),
                 legend.text = T, 
                 args.legend = list(x='topleft',bty='n',
                                    title='Observation Size',
                                    ncol=3))
    arrows(bp, pcor[,,j,i]-ecor[,,j,i], bp, pcor[,,j,i]+ecor[,,j,i],length=0)
    
  } # end of tau0 loop
  dev.off()
  
} # end of fit index


### ICs X N X J X gam2 --------

pcor = apply(mod_select[,TRUE_MOD,,,,,,], c(3,2,7,1), mean)
ecor = apply(mod_select[,TRUE_MOD,,,,,,], c(3,2,7,1), ci)

for(i in 1:3) {
  pdf(paste0('figs/pcorrect_ic_n_j_gam2_',FIT_ICS[i],'_barplot.pdf'), 8, 4)
  layout(matrix(1:2,1,2))
  
  for(j in 1:2) {
    bp = barplot(pcor[,,j,i], beside=T,ylim=c(0,1), 
                 xlab = 'Cluster Size', ylab='P(Correct)', 
                 main = bquote(gamma[2] ~ ' = ' ~ .(dimnames(pcor)[[3]][j])),
                 legend.text = T, 
                 args.legend = list(x='topleft',bty='n',
                                    title='Observation Size',
                                    ncol=3))
    arrows(bp, pcor[,,j,i]-ecor[,,j,i], bp, pcor[,,j,i]+ecor[,,j,i],length=0)
    
  } # end of tau0 loop
  dev.off()
  
} # end of fit index


### Overall model selection -------------------------------------------------

pdf('figs/overall_mod_select_barplot.pdf', 6, 6)

pchoose = apply(mod_select[,,,,,,,], c(1,2), mean)
echoose = apply(mod_select[,,,,,,,], c(1,2), ci)

bp = barplot(pchoose, beside=T, ylim=c(0,1), 
             xlab='Chosen Model', ylab='P(Chosen)', 
             names.arg = LETTERS[1:5], 
             legend.text = T,
             args.legend = 
               list(x='topleft',
                    bty='n'))
arrows(bp,pchoose-echoose,bp,pchoose+echoose,length=0)
abline(h=1/length(MODELS),lty=2)

dev.off()


# Create table with proportions -------------------------------------------

prop_table = melt(mod_select, 
                  varnames = c('IC','Model',IVS),
                  value.name = 'p_choose')
write.csv(prop_table,'out/prop_table.csv')





rm(list=ls())

library(reshape2) # for wide -> long data transform
alpha   = scales::alpha
tab_mod = sjPlot::tab_model
ci      = function(x,a=qnorm(.975),na.rm=T) sd(x,na.rm=na.rm)/sqrt(length(x)) * a 
ci_binom = function(x, bound=1) prop.test(sum(x), length(x))$conf.int[bound]

options(max.print=999999) ## for storing regression output


# Set constants -----------------------------------------------------------

TRUE_MOD   = 1 # i.e., "model A"
PARS       = c('g0','g1','g2','t0','t1')
FIT_ICS    = c('DIC','WAIC','LOO')
IVS        = c('marg','N','J','TAU0','TAU1','G10','G20')
MODELS     = 1:5
MODEL_COLS = alpha(RColorBrewer::brewer.pal(length(MODELS), 'Dark2'),.5) # for plotting

# Load summary data file --------------------------------------------------

dat = readRDS('data/fits_summary.rds')

# table(dat[!duplicated(dat$cond), 'marg'])

# Summarize diagnostics ---------------------------------------------------
## Rhats ------------------------------------------------------------------

rhats_cond = sapply(dat[dat$marg==0,'max_rhats'], function(x) x)
rhats_marg = sapply(dat[dat$marg==1,'max_rhats'], function(x) x)
all_rhats  = cbind(rhats_cond, rhats_marg)
rhat_p   = mean(all_rhats > 1.05)

### overall ---------

pdf('figs/diag/rhats_hists.pdf', 6, 6)

hist(all_rhats[1,], col = MODEL_COLS[1], xlim=c(.99,1.10), 
     xlab=expression('Max'~hat(R)), main='', ylim=c(0,30000))
for(i in 2:nrow(all_rhats)) hist(all_rhats[i,], col=MODEL_COLS[i], add=T)
abline(v=1.05, lty=2)

legend('topright', bty='n', fill=MODEL_COLS, legend=LETTERS[1:nrow(all_rhats)], title='Model')
legend('right',bty='n', legend=paste0('P(>1.05) = ',round(rhat_p,4)))

dev.off()

### by estimation strategy -----------

pdf('figs/diag/rhats_hists_by_strat.pdf', 12, 6)
layout(matrix(1:2,1,2))

### cond
hist(rhats_cond[1,], col = MODEL_COLS[1], xlim=c(.99,1.10), 
     xlab=expression('Max'~hat(R)), ylim=c(0,12500), 
     main="Conditional Estimation")
for(i in 2:nrow(rhats_cond)) hist(rhats_cond[i,], col=MODEL_COLS[i], add=T)
abline(v=1.05, lty=2)
legend('topright', bty='n', fill=MODEL_COLS, legend=LETTERS[1:nrow(rhats_cond)], title='Model')

### marg
hist(rhats_marg[1,], col = MODEL_COLS[1], xlim=c(.99,1.05), 
     xlab=expression('Max'~hat(R)), ylim=c(0,12500), 
     main="Marginal Estimation")
for(i in 2:nrow(rhats_marg)) hist(rhats_marg[i,], col=MODEL_COLS[i], add=T)

dev.off()


## ESS ----------------------------------------------------------------------

ess_cond = sapply(dat[dat$marg==0,'min_ess'], function(x) x)
ess_marg = sapply(dat[dat$marg==1,'min_ess'], function(x) x)
ess_all  = sapply(dat$min_ess, function(x) x)
ess_p    = mean(ess_all < 100)

### overall -------------------
pdf('figs/diag/ess_hists.pdf', 6, 6)

hist(ess_all[1,], col = MODEL_COLS[1], breaks = 10,
     xlab = 'Min ESS', main='', ylim=c(0,30000), xlim=c(0,6000))
for(i in 2:nrow(ess_all)) hist(ess_all[i,], col=MODEL_COLS[i], add=T, breaks = 10)
abline(v=100, lty=2)
legend('topright', bty='n', fill=MODEL_COLS, legend=LETTERS[1:nrow(ess_all)], title='Model')
legend('right',bty='n', legend=paste0('P(<100) = ',round(ess_p,4)))

dev.off()

### by estimation strategy -------------------

pdf('figs/diag/ess_hists_by_strat.pdf', 12, 6)
layout(matrix(1:2,1,2))

### cond
hist(ess_cond[1,], col = MODEL_COLS[1], breaks = 10,
     xlab = 'Min ESS', main='Conditional Estimation', ylim=c(0,8000), xlim=c(0,6000))
for(i in 2:nrow(ess_cond)) hist(ess_cond[i,], col=MODEL_COLS[i], add=T, breaks = 15)
abline(v=100, lty=2)
legend('topright', bty='n', fill=MODEL_COLS, legend=LETTERS[1:nrow(ess_cond)], title='Model')
mean(ess_cond < 100)
legend('right',bty='n', legend=paste0('P(<100) = ',round(mean(ess_cond<100),4)))


### marg
hist(ess_marg[1,], col = MODEL_COLS[1], breaks = 10,
     xlab = 'Min ESS', main='Marginal Estimation', ylim=c(0,10000), xlim=c(0,6000))
for(i in 2:nrow(ess_marg)) hist(ess_marg[i,], col=MODEL_COLS[i], add=T, breaks = 15)
abline(v=100, lty=2)
legend('topright', bty='n', fill=MODEL_COLS, legend=LETTERS[1:nrow(ess_marg)], title='Model')
legend('right',bty='n', legend=paste0('P(<100) = ',round(mean(ess_marg<100),4)))

dev.off()

# Quantify parameter bias ------------------------------------------

#' create indexable array with proportions of model selection per: 
#'    1. Parameter (5: g00, g10, g20, tau0, tau1)
#'    2. Model
#'    3. N
#'    4. J
#'    5. Tau0
#'    6. Tau1
#'    7. gamma10
#'    8. gamma20
#'    9. Estimation strategy

bias_mat = array(NA,dim = c(length(PARS), 
                            length(MODELS), 
                            length(unique(dat$N)), 
                            length(unique(dat$J)), 
                            length(unique(dat$TAU0)),
                            length(unique(dat$TAU1)),
                            length(unique(dat$G10)),
                            length(unique(dat$G20)), 
                            length(unique(dat$marg))),
                 dimnames = list(PARS,
                                 MODELS,
                                 unique(dat$N),
                                 unique(dat$J),
                                 unique(dat$TAU0),
                                 unique(dat$TAU1),
                                 unique(dat$G10),
                                 unique(dat$G20), 
                                 unique(dat$marg))
)

for(i in PARS) {
  
  d_list = dat[,paste0('d_',i,'s')]
  
  for(m in MODELS) {
    d_mod = sapply(d_list, function(x) x[[m]])
    bias_mat[i,m,,,,,,,] = tapply(d_mod, 
                                  list(dat$N,
                                       dat$J, 
                                       dat$TAU0, 
                                       dat$TAU1,
                                       dat$G10, 
                                       dat$G20, 
                                       dat$marg), 
                                  mean)
    
  }
}

#' Convenience function to visualize bias rates across a condition
#'
#' @param fname file name to save figues
#' @param cond_idx index corresponding to colours of bars
#' @param leg_title legend title
#' @param main_title plot title
#' @param leg_pos legend position
#' @param est_idx which estimation stategy to use (defaults to average across strats)

bias_barplot_2way = function(fname,
                             parnames,
                             leg_title='', 
                             legend_text='',
                             main_title='',
                             leg_pos='topright',
                             est_idx=1:2
) 
{
  
  ## Compute average error
  md = apply(bias_mat[parnames,,,,,,,,est_idx], c(1,2), mean, na.rm=T)
  ed = apply(bias_mat[parnames,,,,,,,,est_idx], c(1,2), ci)
  
  
  ## Visualize
  pdf(fname, 6, 6)
  
  ylim = c(-.04, .04)#range(pretty(c(md+ed+.01, md-ed-.01)))
  b = barplot(md, beside=T, ylim=ylim, xpd=F, 
              xlab='Model', ylab=expression('Avg. '~Delta),
              legend.text = T, main=main_title,
              args.legend = list(x=leg_pos,bty='n', cex=1.5,
                                 legend=legend_text))
  arrows(b, md-ed, b, md+ed, length=0)
  abline(h=0,lty=1)
  dev.off()
  
}

#' Convenience function to visualize bias rates across two conditions
#'
#' @param fname file name to save figues
#' @param cond_idx index corresponding to colours of bars
#' @param leg_title legend title
#' @param group_title plot title
#' @param leg_pos legend position
#' @param est_idx which estimation stategy to use (defaults to average across strats)
bias_barplot_3way = function(fname, 
                             cond_idx,
                             leg_title='',
                             group_title='',
                             leg_pos='bottomright',
                             est_idx=1:2) 
{
  
  ## Compute error
  md = apply(bias_mat[,,,,,,,,est_idx], c(cond_idx,2,1), mean)
  ed = apply(bias_mat[,,,,,,,,est_idx], c(cond_idx,2,1), ci)
  
  ## Visualize
  pdf(fname, 12, 6)
  layout(matrix(1:6,2,3,byrow=T))
  
  titles = list('g0'=expression(gamma['00']),
                'g1'=expression(gamma[10]),
                'g2'=expression(gamma[20]),
                't0'=expression(tau[0]^2),
                't1'=expression(tau[1]^2))
  
  ylim = range(pretty(c(md+ed, md-ed)))
  
  ## loop through pars for each plot
  for(p in (PARS)) {
    b = barplot(md[,,p], beside=T, ylim=ylim, xpd=F, 
                xlab='Model', ylab=expression('Avg. '~Delta),
                main=titles[[p]], 
                legend.text = T, 
                args.legend = list(x=leg_pos,bty='n', cex=1, 
                                   title=leg_title))
    arrows(b, md[,,p]-ed[,,p], b, md[,,p]+ed[,,p], length=0)
  }
  mtext(group_title, side=3, line=-1.25, outer=TRUE, cex=1, font=3)
  
  dev.off()
  
}

## model X gamma ----------------------------------------------------------

# Avg. error in fixefs array. Dims: 
# 1: models (1-5)
# 2: coef (g0,g1,g2)

### overall
bias_barplot_2way('figs/bias/delta_fixef_barplot.pdf', 
                  c('g0','g1','g2'), 
                  legend_text = c(expression(gamma['00']),
                                  expression(gamma[10]),
                                  expression(gamma[20])))

### conditional
bias_barplot_2way('figs/bias/delta_fixef_barplot_conditional.pdf', 
                  c('g0','g1','g2'), 
                  legend_text = c(expression(gamma['00']),
                                  expression(gamma[10]),
                                  expression(gamma[20])), 
                  est_idx = '0', main_title="Conditional Estimation")

### marginal
bias_barplot_2way('figs/bias/delta_fixef_barplot_marginal.pdf', 
                  c('g0','g1','g2'), 
                  legend_text = c(expression(gamma['00']),
                                  expression(gamma[10]),
                                  expression(gamma[20])), 
                  est_idx = '1', main_title="Marginal Estimation")



## model X tau ----------------------------------------------------------

# Avg. error in fixefs array. Dims: 
# 1: models (1-5)
# 2: coef (g0,g1,g2)

### overall
bias_barplot_2way('figs/bias/delta_ranef_barplot.pdf', 
                  c('t0','t1'), leg_pos = 'bottomleft',
                  legend_text = c(expression(tau['0']),
                                  expression(tau['1'])))

### conditional
bias_barplot_2way('figs/bias/delta_ranef_barplot_conditional.pdf', 
                  c('t0','t1'), leg_pos = 'bottomleft',
                  legend_text = c(expression(tau['0']),
                                  expression(tau['1'])), 
                  est_idx = '0', main_title="Conditional Estimation")

### marginal
bias_barplot_2way('figs/bias/delta_ranef_barplot_marginal.pdf', 
                  c('t0','t1'), leg_pos = 'bottomleft',
                  legend_text = c(expression(tau['0']),
                                  expression(tau['1'])), 
                  est_idx = '1', main_title="Marginal Estimation")


## model X par X n ----------------------------------------------------------

### overall
bias_barplot_3way('figs/bias/delta_n_barplot.pdf', 3, leg_title = 'Cluster Size')

### conditional
bias_barplot_3way('figs/bias/delta_n_barplot_conditional.pdf', 3, leg_title = 'Cluster Size',
                  est_idx = '0',
                  group_title = 'Conditional Estimation')

### marginal
bias_barplot_3way('figs/bias/delta_n_barplot_marginal.pdf', 3, leg_title = 'Cluster Size',
                  est_idx = '1',
                  group_title = 'Marginal Estimation')


## model X par X j ----------------------------------------------------------

### overall
bias_barplot_3way('figs/bias/delta_j_barplot.pdf', 4, leg_title = 'Observation Size')

### conditional
bias_barplot_3way('figs/bias/delta_j_barplot_conditional.pdf', 4, leg_title = 'Observation Size',
                  est_idx = '0',
                  group_title = 'Conditional Estimation')

### marginal
bias_barplot_3way('figs/bias/delta_j_barplot_marginal.pdf', 4, leg_title = 'Observation Size',
                  est_idx = '1',
                  group_title = 'Marginal Estimation')

## model X par X t0 ----------------------------------------------------------

### overall
bias_barplot_3way('figs/bias/delta_t0_barplot.pdf', 5, leg_title = expression(tau[0]^2))

### conditional
bias_barplot_3way('figs/bias/delta_t0_barplot_conditional.pdf', 5, leg_title = expression(tau[0]^2),
                  est_idx = '0',
                  group_title = 'Conditional Estimation')

### marginal
bias_barplot_3way('figs/bias/delta_t0_barplot_marginal.pdf', 5, leg_title = expression(tau[0]^2),
                  est_idx = '1',
                  group_title = 'Marginal Estimation')


## model X par X t1 ----------------------------------------------------------

### overall
bias_barplot_3way('figs/bias/delta_t1_barplot.pdf', 6, leg_title = expression(tau[1]^2))

### conditional
bias_barplot_3way('figs/bias/delta_t1_barplot_conditional.pdf', 6, leg_title = expression(tau[1]^2),
                  est_idx = '0',
                  group_title = 'Conditional Estimation')

### marginal
bias_barplot_3way('figs/bias/delta_t1_barplot_marginal.pdf', 6, leg_title = expression(tau[1]^2),
                  est_idx = '1',
                  group_title = 'Marginal Estimation')


## model X par X g0 ----------------------------------------------------------

### overall
bias_barplot_3way('figs/bias/delta_gam0_barplot.pdf', 7, leg_title = expression(gamma['00']))

### conditional
bias_barplot_3way('figs/bias/delta_gam0_barplot_conditional.pdf', 7, leg_title = expression(gamma['00']),
                  est_idx = '0',
                  group_title = 'Conditional Estimation')

### marginal
bias_barplot_3way('figs/bias/delta_gam0_barplot_marginal.pdf', 7, leg_title = expression(gamma['00']),
                  est_idx = '1',
                  group_title = 'Marginal Estimation')


## model X par X g1 ----------------------------------------------------------

### overall
bias_barplot_3way('figs/bias/delta_gam1_barplot.pdf', 8, leg_title = expression(gamma['10']))

### conditional
bias_barplot_3way('figs/bias/delta_gam1_barplot_conditional.pdf', 8, leg_title = expression(gamma['10']),
                  est_idx = '0',
                  group_title = 'Conditional Estimation')

### marginal
bias_barplot_3way('figs/bias/delta_gam1_barplot_marginal.pdf', 8, leg_title = expression(gamma['10']),
                  est_idx = '1',
                  group_title = 'Marginal Estimation')



# IC Performance without Uncertainty -------------------------------------------------------------

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
melt_dat[,IVS] = apply(melt_dat[,IVS],2,function(x) x - median(unique(x)))

# fit
# acc_mod = glm(accuracy ~ ic*marg*N*J*TAU0*TAU1*G10*G20,
#               data=melt_dat,
#               family='binomial')
# 
# capture.output(summary(acc_mod), file = 'out/acc_model_summary.txt')



## Visualize performance ---------------------------------------------------


## Melt
melt_dat = melt(dat, 
                id.vars = IVS,
                measure.vars = paste0('acc_',FIT_ICS), 
                variable.name = 'ic',
                value.name = 'accuracy')
melt_dat$ic = sapply(melt_dat$ic, function(x) 
  strsplit(as.character(x), split = '_')[[1]][2])
melt_dat$ic = factor(melt_dat$ic, levels = c('DIC','WAIC','LOO'))



melt_dat$marg = ifelse(melt_dat$marg==0, 'Conditional', 'Marginal')


### ICs --------

#### overall -------------- 

pdf('figs/ic_performance/wo_uncertainty/pcorrect_ic_barplot.pdf', 6, 6)

## Compute accuaracy and error
pcor = tapply(melt_dat$accuracy, melt_dat$ic, mean)
elb  = tapply(melt_dat$accuracy, melt_dat$ic, ci_binom, bound=1)
eub  = tapply(melt_dat$accuracy, melt_dat$ic, ci_binom, bound=2)

bp = barplot(pcor, ylim=c(0,1.1), ylab='P(Correct)')
arrows(bp, elb, bp, eub,length=0)

dev.off()

#### by estimation strategy -------------- 

pdf('figs/ic_performance/wo_uncertainty/pcorrect_ic_barplot_by_strat.pdf', 6, 6)

pcor = tapply(melt_dat$accuracy, list(melt_dat$marg, melt_dat$ic), mean)
elb  = tapply(melt_dat$accuracy, list(melt_dat$marg, melt_dat$ic), ci_binom, bound=1)
eub  = tapply(melt_dat$accuracy, list(melt_dat$marg, melt_dat$ic), ci_binom, bound=2)

bp = barplot(pcor, beside=T, ylim=c(0,1.1), ylab='P(Correct)', 
             legend.text = T, 
             args.legend = list(x='topleft', bty='n', cex=1.5, title='Estimation Approach'))
arrows(bp, elb, bp, eub,length=0)

dev.off()


### ICs X N --------

pcor = tapply(melt_dat$accuracy, list(melt_dat$N, melt_dat$ic, melt_dat$marg), mean)
elb  = tapply(melt_dat$accuracy, list(melt_dat$N, melt_dat$ic, melt_dat$marg), ci_binom, bound=1)
eub  = tapply(melt_dat$accuracy, list(melt_dat$N, melt_dat$ic, melt_dat$marg), ci_binom, bound=2)


### overall
pdf('figs/ic_performance/wo_uncertainty/pcorrect_ic_n_barplot.pdf', 6, 6)

bp = barplot(apply(pcor, 1:2, mean), beside=T, ylim=c(0,1.1), ylab='P(Correct)', 
             legend.text = T, 
             args.legend = list(x='topleft', bty='n', cex=1.5, title='Cluster Size',ncol=3))
arrows(bp, apply(elb, 1:2, mean), bp, apply(eub, 1:2, mean),length=0)

dev.off()

### conditional
pdf('figs/ic_performance/wo_uncertainty/pcorrect_ic_n_barplot_conditional.pdf', 6, 6)

bp = barplot(pcor[,,1], beside=T, ylim=c(0,1.1), ylab='P(Correct)', 
             legend.text = T, 
             args.legend = list(x='topleft', bty='n', cex=1.5, title='Cluster Size',ncol=3))
arrows(bp, elb[,,1], bp, eub[,,1],length=0)

dev.off()

### marginal
pdf('figs/ic_performance/wo_uncertainty/pcorrect_ic_n_barplot_marginal.pdf', 6, 6)

bp = barplot(pcor[,,2], beside=T, ylim=c(0,1.1), ylab='P(Correct)', 
             legend.text = T, 
             args.legend = list(x='topleft', bty='n', cex=1.5, title='Cluster Size',ncol=3))
arrows(bp, elb[,,2], bp, eub[,,2],length=0)

dev.off()



### ICs X J --------

pcor = tapply(melt_dat$accuracy, list(melt_dat$J, melt_dat$ic, melt_dat$marg), mean)
elb  = tapply(melt_dat$accuracy, list(melt_dat$J, melt_dat$ic, melt_dat$marg), ci_binom, bound=1)
eub  = tapply(melt_dat$accuracy, list(melt_dat$J, melt_dat$ic, melt_dat$marg), ci_binom, bound=2)


### overall
pdf('figs/ic_performance/wo_uncertainty/pcorrect_ic_j_barplot.pdf', 6, 6)

bp = barplot(apply(pcor, 1:2, mean), beside=T, ylim=c(0,1.1), ylab='P(Correct)', 
             legend.text = T, 
             args.legend = list(x='topleft', bty='n', cex=1.5, title='Cluster Size',ncol=3))
arrows(bp, apply(elb, 1:2, mean), bp, apply(eub, 1:2, mean),length=0)

dev.off()

### conditional
pdf('figs/ic_performance/wo_uncertainty/pcorrect_ic_j_barplot_conditional.pdf', 6, 6)

bp = barplot(pcor[,,1], beside=T, ylim=c(0,1.1), ylab='P(Correct)', 
             legend.text = T, 
             args.legend = list(x='topleft', bty='n', cex=1.5, title='Cluster Size',ncol=3))
arrows(bp, elb[,,1], bp, eub[,,1],length=0)

dev.off()

### marginal
pdf('figs/ic_performance/wo_uncertainty/pcorrect_ic_j_barplot_marginal.pdf', 6, 6)

bp = barplot(pcor[,,2], beside=T, ylim=c(0,1.1), ylab='P(Correct)', 
             legend.text = T, 
             args.legend = list(x='topleft', bty='n', cex=1.5, title='Cluster Size',ncol=3))
arrows(bp, elb[,,2], bp, eub[,,2],length=0)

dev.off()

### ICs X Tau0 --------

pcor = tapply(melt_dat$accuracy, list(melt_dat$TAU0, melt_dat$ic, melt_dat$marg), mean)
elb  = tapply(melt_dat$accuracy, list(melt_dat$TAU0, melt_dat$ic, melt_dat$marg), ci_binom, bound=1)
eub  = tapply(melt_dat$accuracy, list(melt_dat$TAU0, melt_dat$ic, melt_dat$marg), ci_binom, bound=2)


### overall
pdf('figs/ic_performance/wo_uncertainty/pcorrect_ic_t0_barplot.pdf', 6, 6)

bp = barplot(apply(pcor, 1:2, mean), beside=T, ylim=c(0,1.1), ylab='P(Correct)', 
             legend.text = T, 
             args.legend = list(x='topleft', bty='n', cex=1.5, title=expression(tau[0]),ncol=1))
arrows(bp, apply(elb, 1:2, mean), bp, apply(eub, 1:2, mean),length=0)

dev.off()

### conditional
pdf('figs/ic_performance/wo_uncertainty/pcorrect_ic_t0_barplot_conditional.pdf', 6, 6)

bp = barplot(pcor[,,1], beside=T, ylim=c(0,1.1), ylab='P(Correct)', 
             legend.text = T, main='Conditional Estimation',
             args.legend = list(x='topleft', bty='n', cex=1.5, title=expression(tau[0]),ncol=1))
arrows(bp, elb[,,1], bp, eub[,,1],length=0)

dev.off()

### marginal
pdf('figs/ic_performance/wo_uncertainty/pcorrect_ic_t0_barplot_marginal.pdf', 6, 6)

bp = barplot(pcor[,,2], beside=T, ylim=c(0,1.1), ylab='P(Correct)', 
             legend.text = T,  main='Marginal Estimation',
             args.legend = list(x='topleft', bty='n', cex=1.5, title=expression(tau[0]),ncol=1))
arrows(bp, elb[,,2], bp, eub[,,2],length=0)

dev.off()

### ICs X Tau1 --------

pcor = tapply(melt_dat$accuracy, list(melt_dat$TAU1, melt_dat$ic, melt_dat$marg), mean)
elb  = tapply(melt_dat$accuracy, list(melt_dat$TAU1, melt_dat$ic, melt_dat$marg), ci_binom, bound=1)
eub  = tapply(melt_dat$accuracy, list(melt_dat$TAU1, melt_dat$ic, melt_dat$marg), ci_binom, bound=2)


### overall
pdf('figs/ic_performance/wo_uncertainty/pcorrect_ic_t1_barplot.pdf', 6, 6)

bp = barplot(apply(pcor, 1:2, mean), beside=T, ylim=c(0,1.1), ylab='P(Correct)', 
             legend.text = T, 
             args.legend = list(x='topleft', bty='n', cex=1.5, title=expression(tau[1]),ncol=1))
arrows(bp, apply(elb, 1:2, mean), bp, apply(eub, 1:2, mean),length=0)

dev.off()

### conditional
pdf('figs/ic_performance/wo_uncertainty/pcorrect_ic_t1_barplot_conditional.pdf', 6, 6)

bp = barplot(pcor[,,1], beside=T, ylim=c(0,1.1), ylab='P(Correct)', 
             legend.text = T,  main='Conditional Estimation',
             args.legend = list(x='topleft', bty='n', cex=1.5, title=expression(tau[1]),ncol=1))
arrows(bp, elb[,,1], bp, eub[,,1],length=0)

dev.off()

### marginal
pdf('figs/ic_performance/wo_uncertainty/pcorrect_ic_t1_barplot_marginal.pdf', 6, 6)

bp = barplot(pcor[,,2], beside=T, ylim=c(0,1.1), ylab='P(Correct)', 
             legend.text = T,  main='Marginal Estimation',
             args.legend = list(x='topleft', bty='n', cex=1.5, title=expression(tau[1]),ncol=1))
arrows(bp, elb[,,2], bp, eub[,,2],length=0)

dev.off()

### ICs X gam1 --------

pcor = tapply(melt_dat$accuracy, list(melt_dat$G10, melt_dat$ic, melt_dat$marg), mean)
elb  = tapply(melt_dat$accuracy, list(melt_dat$G10, melt_dat$ic, melt_dat$marg), ci_binom, bound=1)
eub  = tapply(melt_dat$accuracy, list(melt_dat$G10, melt_dat$ic, melt_dat$marg), ci_binom, bound=2)


### overall
pdf('figs/ic_performance/wo_uncertainty/pcorrect_ic_g1_barplot.pdf', 6, 6)

bp = barplot(apply(pcor, 1:2, mean), beside=T, ylim=c(0,1.1), ylab='P(Correct)', 
             legend.text = T, 
             args.legend = list(x='topleft', bty='n', cex=1.5, title=expression(gamma[1]),ncol=1))
arrows(bp, apply(elb, 1:2, mean), bp, apply(eub, 1:2, mean),length=0)

dev.off()

### conditional
pdf('figs/ic_performance/wo_uncertainty/pcorrect_ic_g1_barplot_conditional.pdf', 6, 6)

bp = barplot(pcor[,,1], beside=T, ylim=c(0,1.1), ylab='P(Correct)', 
             legend.text = T,  main='Conditional Estimation',
             args.legend = list(x='topleft', bty='n', cex=1.5, title=expression(gamma[1]),ncol=1))
arrows(bp, elb[,,1], bp, eub[,,1],length=0)

dev.off()

### marginal
pdf('figs/ic_performance/wo_uncertainty/pcorrect_ic_g1_barplot_marginal.pdf', 6, 6)

bp = barplot(pcor[,,2], beside=T, ylim=c(0,1.1), ylab='P(Correct)', 
             legend.text = T,  main='Marginal Estimation',
             args.legend = list(x='topleft', bty='n', cex=1.5, title=expression(gamma[1]),ncol=1))
arrows(bp, elb[,,2], bp, eub[,,2],length=0)

dev.off()

### ICs X gam2 --------

pcor = tapply(melt_dat$accuracy, list(melt_dat$G20, melt_dat$ic, melt_dat$marg), mean)
elb  = tapply(melt_dat$accuracy, list(melt_dat$G20, melt_dat$ic, melt_dat$marg), ci_binom, bound=1)
eub  = tapply(melt_dat$accuracy, list(melt_dat$G20, melt_dat$ic, melt_dat$marg), ci_binom, bound=2)

### overall
pdf('figs/ic_performance/wo_uncertainty/pcorrect_ic_g2_barplot.pdf', 6, 6)

bp = barplot(apply(pcor, 1:2, mean), beside=T, ylim=c(0,1.1), ylab='P(Correct)', 
             legend.text = T, 
             args.legend = list(x='topleft', bty='n', cex=1.5, title=expression(gamma[2]),ncol=1))
arrows(bp, apply(elb, 1:2, mean), bp, apply(eub, 1:2, mean),length=0)

dev.off()

### conditional
pdf('figs/ic_performance/wo_uncertainty/pcorrect_ic_g2_barplot_conditional.pdf', 6, 6)

bp = barplot(pcor[,,1], beside=T, ylim=c(0,1.1), ylab='P(Correct)', 
             legend.text = T,  main='Conditional Estimation',
             args.legend = list(x='topleft', bty='n', cex=1.5, title=expression(gamma[2]),ncol=1))
arrows(bp, elb[,,1], bp, eub[,,1],length=0)

dev.off()

### marginal
pdf('figs/ic_performance/wo_uncertainty/pcorrect_ic_g2_barplot_marginal.pdf', 6, 6)

bp = barplot(pcor[,,2], beside=T, ylim=c(0,1.1), ylab='P(Correct)', 
             legend.text = T,  main='Marginal Estimation',
             args.legend = list(x='topleft', bty='n', cex=1.5, title=expression(gamma[2]),ncol=1))
arrows(bp, elb[,,2], bp, eub[,,2],length=0)

dev.off()



### ICs X N X J --------

pcor = tapply(melt_dat$accuracy, list(melt_dat$J, melt_dat$N, melt_dat$ic, melt_dat$marg), mean)
elb  = tapply(melt_dat$accuracy, list(melt_dat$J, melt_dat$N, melt_dat$ic, melt_dat$marg), ci_binom, bound=1)
eub  = tapply(melt_dat$accuracy, list(melt_dat$J, melt_dat$N, melt_dat$ic, melt_dat$marg), ci_binom, bound=2)

### overall
pdf('figs/ic_performance/wo_uncertainty/pcorrect_ic_n_j_barplot.pdf', 8, 4)
par(mar=c(5, 4, 6, 2) + 0.1)
layout(matrix(1:3,1,3))

for(n in as.character(unique(melt_dat$N))) {
  bp = barplot(apply(pcor, 1:3, mean)[,n,], beside=T, ylim=c(0,1.1),
               ylab='P(Correct)', xlab="",
               legend.text = T, 
               main=paste0("Cluster Size = ",n),
               args.legend = list(x='topleft', bty='n', title='Observation Size',ncol=3))
  arrows(bp, apply(elb, 1:3, mean)[,n,], bp, apply(eub, 1:3, mean)[,n,],length=0)
  
}

dev.off()

### conditional
pdf('figs/ic_performance/wo_uncertainty/pcorrect_ic_n_j_barplot_conditional.pdf', 8, 4)
par(mar=c(5, 4, 6, 2) + 0.1)
layout(matrix(1:3,1,3))

for(n in as.character(unique(melt_dat$N))) {
  bp = barplot(apply(pcor, 1:4, mean)[,n,,1], beside=T, ylim=c(0,1.1),
               ylab='P(Correct)', xlab="",
               legend.text = T, 
               main=paste0("Cluster Size = ",n),
               args.legend = list(x='topleft', bty='n', title='Observation Size',ncol=3))
  arrows(bp, apply(elb, 1:4, mean)[,n,,1], bp, apply(eub, 1:4, mean)[,n,,1],length=0)
  
}

mtext('Conditional Estimation', side=3, line=-1.25, outer=TRUE, cex=1, font=3)

dev.off()

### marginal
pdf('figs/ic_performance/wo_uncertainty/pcorrect_ic_n_j_barplot_marginal.pdf', 8, 4)
par(mar=c(5, 4, 6, 2) + 0.1)
layout(matrix(1:3,1,3))

for(n in as.character(unique(melt_dat$N))) {
  bp = barplot(apply(pcor, 1:4, mean)[,n,,2], beside=T, ylim=c(0,1.1),
               ylab='P(Correct)', xlab="",
               legend.text = T, 
               main=paste0("Cluster Size = ",n),
               args.legend = list(x='topleft', bty='n', title='Observation Size',ncol=3))
  arrows(bp, apply(elb, 1:4, mean)[,n,,2], bp, apply(eub, 1:4, mean)[,n,,2],length=0)
  
}

mtext('Marginal Estimation', side=3, line=-1.25, outer=TRUE, cex=1, font=3)

dev.off()


## Overall model selection -------------------------------------------------


#' create indexable array with proportions of model selection per: 
#'    1. Fit index
#'    2. Model selected
#'    3. N
#'    4. J
#'    5. Tau0
#'    6. Tau1
#'    7. gamma10
#'    8. gamma20
#'    9. Estimation strategy

mod_select = array(NA,dim = c(length(FIT_ICS), 
                              length(MODELS), 
                              length(unique(dat$N)), 
                              length(unique(dat$J)), 
                              length(unique(dat$TAU0)),
                              length(unique(dat$TAU1)),
                              length(unique(dat$G10)),
                              length(unique(dat$G20)),
                              length(unique(dat$marg))),
                   dimnames = list(FIT_ICS,
                                   MODELS,
                                   unique(dat$N),
                                   unique(dat$J),
                                   unique(dat$TAU0),
                                   unique(dat$TAU1),
                                   unique(dat$G10),
                                   unique(dat$G20), 
                                   unique(dat$marg))
)

for(i in FIT_ICS) {
  
  win_vec = dat[,paste0('win_',i)]
  
  for(m in MODELS) {
    mod_select[i,m,,,,,,,] = tapply(win_vec, 
                                    list(dat$N,
                                         dat$J, 
                                         dat$TAU0, 
                                         dat$TAU1,
                                         dat$G10, 
                                         dat$G20, 
                                         dat$marg), 
                                    function(x) mean(x==m, na.rm=T))
    
  }
}

#' Convenience function to visualize model selection rates across ICs
#'
#' @param fname file name to save figues
#' @param leg_title legend title
#' @param main_title plot title
#' @param leg_pos legend position
#' @param ic_idx which ICs to use
#' @param est_idx which estimation stategy to use (defaults to average across strats)
IC_barplot_modselect = function(fname,
                                cond_idx,
                                main_title='',
                                leg_pos='topleft',
                                ic_idx=1,est_idx=1:2) 
{
  
  ## Compute choice rates
  pchoose = apply(mod_select[,,,,,,,,est_idx], c(ic_idx,2), mean)
  echoose = apply(mod_select[,,,,,,,,est_idx], c(ic_idx,2), ci)
  
  ## Visualize
  pdf(fname, 6, 6)
  
  bp = barplot(pchoose, beside=T, ylim=c(0,1), 
               xlab='Chosen Model', ylab='P(Chosen)', 
               main=main_title,
               names.arg = LETTERS[1:5], 
               legend.text = T,
               args.legend = 
                 list(x=leg_pos,
                      bty='n'))
  arrows(bp,pchoose-echoose,bp,pchoose+echoose,length=0)
  abline(h=1/length(MODELS),lty=2)
  
  dev.off()
  
}

### overall
IC_barplot_modselect('figs/ic_performance/wo_uncertainty/overall_mod_select_barplot.pdf')

### conditional
IC_barplot_modselect('figs/ic_performance/wo_uncertainty/overall_mod_select_barplot_conditional.pdf', 
                     est_idx = '0', main_title = 'Conditional Estimation')

### marginal
IC_barplot_modselect('figs/ic_performance/wo_uncertainty/overall_mod_select_barplot_marginal.pdf', 
                     est_idx = '1', main_title = 'Marginal Estimation')


### Create table with proportions -------------------------------------------

prop_table = melt(mod_select, 
                  varnames = c('IC','Model',IVS),
                  value.name = 'p_choose')
write.csv(prop_table,'out/prop_table.csv')


# IC Performance With Uncertainty -------------------------------------------------------------

## Model accurate model selection -----------

# compute accuracy
for(i in FIT_ICS){
  dat[,paste0('acc_1se_',i)] = as.numeric(dat[,paste0('win_1se_',i)]==TRUE_MOD)
}

# melt df
melt_dat = melt(dat, id.vars = IVS,
                measure.vars = paste0('acc_1se_',FIT_ICS), 
                variable.name = 'ic',
                value.name = 'accuracy')
melt_dat$ic = sapply(melt_dat$ic, function(x) 
  strsplit(as.character(x), split = '_')[[1]][3])

# make fit index a factor
melt_dat$ic = factor(melt_dat$ic, levels = c('DIC','WAIC','LOO'))

# center predictors
melt_dat[,IVS] = apply(melt_dat[,IVS],2,function(x) x - median(x))

# fit
# acc_mod_1se = glm(accuracy ~ ic*marg*N*J*TAU0*TAU1*G10*G20,
#               data=melt_dat,
#               family='binomial')
# 
# capture.output(summary(acc_mod_1se), file = 'out/acc_model_1se_summary.txt')


## Visualize performance ---------------------------------------------------


## Melt
melt_dat = melt(dat, 
                id.vars = IVS,
                measure.vars = paste0('acc_1se_',FIT_ICS), 
                variable.name = 'ic',
                value.name = 'accuracy')
melt_dat$ic = sapply(melt_dat$ic, function(x) 
  strsplit(as.character(x), split = '_')[[1]][3])
melt_dat$ic = factor(melt_dat$ic, levels = c('DIC','WAIC','LOO'))

melt_dat$marg = ifelse(melt_dat$marg==0, 'Conditional', 'Marginal')


### ICs --------

#### overall -------------- 

pdf('figs/ic_performance/w_uncertainty/pcorrect_ic_barplot.pdf', 6, 6)

## Compute accuaracy and error
pcor = tapply(melt_dat$accuracy, melt_dat$ic, mean)
elb  = tapply(melt_dat$accuracy, melt_dat$ic, ci_binom, bound=1)
eub  = tapply(melt_dat$accuracy, melt_dat$ic, ci_binom, bound=2)

bp = barplot(pcor, ylim=c(0,1.1), ylab='P(Correct)')
arrows(bp, elb, bp, eub,length=0)

dev.off()

#### by estimation strategy -------------- 

pdf('figs/ic_performance/w_uncertainty/pcorrect_ic_barplot_by_strat.pdf', 6, 6)

pcor = tapply(melt_dat$accuracy, list(melt_dat$marg, melt_dat$ic), mean)
elb  = tapply(melt_dat$accuracy, list(melt_dat$marg, melt_dat$ic), ci_binom, bound=1)
eub  = tapply(melt_dat$accuracy, list(melt_dat$marg, melt_dat$ic), ci_binom, bound=2)

bp = barplot(pcor, beside=T, ylim=c(0,1.1), ylab='P(Correct)', 
             legend.text = T, 
             args.legend = list(x='topleft', bty='n', 1, title='Estimation Approach'))
arrows(bp, elb, bp, eub,length=0)

dev.off()


### ICs X N --------

pcor = tapply(melt_dat$accuracy, list(melt_dat$N, melt_dat$ic, melt_dat$marg), mean)
elb  = tapply(melt_dat$accuracy, list(melt_dat$N, melt_dat$ic, melt_dat$marg), ci_binom, bound=1)
eub  = tapply(melt_dat$accuracy, list(melt_dat$N, melt_dat$ic, melt_dat$marg), ci_binom, bound=2)


### overall
pdf('figs/ic_performance/w_uncertainty/pcorrect_ic_n_barplot.pdf', 6, 6)

bp = barplot(apply(pcor, 1:2, mean), beside=T, ylim=c(0,1.1), ylab='P(Correct)', 
             legend.text = T, 
             args.legend = list(x='topleft', bty='n', 1, title='Cluster Size',ncol=3))
arrows(bp, apply(elb, 1:2, mean), bp, apply(eub, 1:2, mean),length=0)

dev.off()

### conditional
pdf('figs/ic_performance/w_uncertainty/pcorrect_ic_n_barplot_conditional.pdf', 6, 6)

bp = barplot(pcor[,,1], beside=T, ylim=c(0,1.1), ylab='P(Correct)', 
             legend.text = T, 
             args.legend = list(x='topleft', bty='n', 1, title='Cluster Size',ncol=3))
arrows(bp, elb[,,1], bp, eub[,,1],length=0)

dev.off()

### marginal
pdf('figs/ic_performance/w_uncertainty/pcorrect_ic_n_barplot_marginal.pdf', 6, 6)

bp = barplot(pcor[,,2], beside=T, ylim=c(0,1.1), ylab='P(Correct)', 
             legend.text = T, 
             args.legend = list(x='topleft', bty='n', 1, title='Cluster Size',ncol=3))
arrows(bp, elb[,,2], bp, eub[,,2],length=0)

dev.off()

### ICs X J --------

pcor = tapply(melt_dat$accuracy, list(melt_dat$J, melt_dat$ic, melt_dat$marg), mean)
elb  = tapply(melt_dat$accuracy, list(melt_dat$J, melt_dat$ic, melt_dat$marg), ci_binom, bound=1)
eub  = tapply(melt_dat$accuracy, list(melt_dat$J, melt_dat$ic, melt_dat$marg), ci_binom, bound=2)


### overall
pdf('figs/ic_performance/w_uncertainty/pcorrect_ic_j_barplot.pdf', 6, 6)

bp = barplot(apply(pcor, 1:2, mean), beside=T, ylim=c(0,1.1), ylab='P(Correct)', 
             legend.text = T, 
             args.legend = list(x='topleft', bty='n', 1, title='Cluster Size',ncol=3))
arrows(bp, apply(elb, 1:2, mean), bp, apply(eub, 1:2, mean),length=0)

dev.off()

### conditional
pdf('figs/ic_performance/w_uncertainty/pcorrect_ic_j_barplot_conditional.pdf', 6, 6)

bp = barplot(pcor[,,1], beside=T, ylim=c(0,1.1), ylab='P(Correct)', 
             legend.text = T, 
             args.legend = list(x='topleft', bty='n', 1, title='Cluster Size',ncol=3))
arrows(bp, elb[,,1], bp, eub[,,1],length=0)

dev.off()

### marginal
pdf('figs/ic_performance/w_uncertainty/pcorrect_ic_j_barplot_marginal.pdf', 6, 6)

bp = barplot(pcor[,,2], beside=T, ylim=c(0,1.1), ylab='P(Correct)', 
             legend.text = T, 
             args.legend = list(x='topleft', bty='n', 1, title='Cluster Size',ncol=3))
arrows(bp, elb[,,2], bp, eub[,,2],length=0)

dev.off()

### ICs X Tau0 --------

pcor = tapply(melt_dat$accuracy, list(melt_dat$TAU0, melt_dat$ic, melt_dat$marg), mean)
elb  = tapply(melt_dat$accuracy, list(melt_dat$TAU0, melt_dat$ic, melt_dat$marg), ci_binom, bound=1)
eub  = tapply(melt_dat$accuracy, list(melt_dat$TAU0, melt_dat$ic, melt_dat$marg), ci_binom, bound=2)


### overall
pdf('figs/ic_performance/w_uncertainty/pcorrect_ic_t0_barplot.pdf', 6, 6)

bp = barplot(apply(pcor, 1:2, mean), beside=T, ylim=c(0,1.1), ylab='P(Correct)', 
             legend.text = T, 
             args.legend = list(x='topleft', bty='n', 1, title=expression(tau[0]),ncol=1))
arrows(bp, apply(elb, 1:2, mean), bp, apply(eub, 1:2, mean),length=0)

dev.off()

### conditional
pdf('figs/ic_performance/w_uncertainty/pcorrect_ic_t0_barplot_conditional.pdf', 6, 6)

bp = barplot(pcor[,,1], beside=T, ylim=c(0,1.1), ylab='P(Correct)', 
             legend.text = T,  main='Conditional Estimation',
             args.legend = list(x='topleft', bty='n', 1, title=expression(tau[0]),ncol=1))
arrows(bp, elb[,,1], bp, eub[,,1],length=0)

dev.off()

### marginal
pdf('figs/ic_performance/w_uncertainty/pcorrect_ic_t0_barplot_marginal.pdf', 6, 6)

bp = barplot(pcor[,,2], beside=T, ylim=c(0,1.1), ylab='P(Correct)', 
             legend.text = T,  main='Marginal Estimation',
             args.legend = list(x='topleft', bty='n', 1, title=expression(tau[0]),ncol=1))
arrows(bp, elb[,,2], bp, eub[,,2],length=0)

dev.off()

### ICs X Tau1 --------

pcor = tapply(melt_dat$accuracy, list(melt_dat$TAU1, melt_dat$ic, melt_dat$marg), mean)
elb  = tapply(melt_dat$accuracy, list(melt_dat$TAU1, melt_dat$ic, melt_dat$marg), ci_binom, bound=1)
eub  = tapply(melt_dat$accuracy, list(melt_dat$TAU1, melt_dat$ic, melt_dat$marg), ci_binom, bound=2)


### overall
pdf('figs/ic_performance/w_uncertainty/pcorrect_ic_t1_barplot.pdf', 6, 6)

bp = barplot(apply(pcor, 1:2, mean), beside=T, ylim=c(0,1.1), ylab='P(Correct)', 
             legend.text = T, 
             args.legend = list(x='topleft', bty='n', 1, title=expression(tau[1]),ncol=1))
arrows(bp, apply(elb, 1:2, mean), bp, apply(eub, 1:2, mean),length=0)

dev.off()

### conditional
pdf('figs/ic_performance/w_uncertainty/pcorrect_ic_t1_barplot_conditional.pdf', 6, 6)

bp = barplot(pcor[,,1], beside=T, ylim=c(0,1.1), ylab='P(Correct)', 
             legend.text = T,  main='Conditional Estimation',
             args.legend = list(x='topleft', bty='n', 1, title=expression(tau[1]),ncol=1))
arrows(bp, elb[,,1], bp, eub[,,1],length=0)

dev.off()

### marginal
pdf('figs/ic_performance/w_uncertainty/pcorrect_ic_t1_barplot_marginal.pdf', 6, 6)

bp = barplot(pcor[,,2], beside=T, ylim=c(0,1.1), ylab='P(Correct)', 
             legend.text = T,  main='Marginal Estimation',
             args.legend = list(x='topleft', bty='n', 1, title=expression(tau[1]),ncol=1))
arrows(bp, elb[,,2], bp, eub[,,2],length=0)

dev.off()

### ICs X gam1 --------

pcor = tapply(melt_dat$accuracy, list(melt_dat$G10, melt_dat$ic, melt_dat$marg), mean)
elb  = tapply(melt_dat$accuracy, list(melt_dat$G10, melt_dat$ic, melt_dat$marg), ci_binom, bound=1)
eub  = tapply(melt_dat$accuracy, list(melt_dat$G10, melt_dat$ic, melt_dat$marg), ci_binom, bound=2)


### overall
pdf('figs/ic_performance/w_uncertainty/pcorrect_ic_g1_barplot.pdf', 6, 6)

bp = barplot(apply(pcor, 1:2, mean), beside=T, ylim=c(0,1.1), ylab='P(Correct)', 
             legend.text = T, 
             args.legend = list(x='topleft', bty='n', 1, title=expression(gamma[1]),ncol=1))
arrows(bp, apply(elb, 1:2, mean), bp, apply(eub, 1:2, mean),length=0)

dev.off()

### conditional
pdf('figs/ic_performance/w_uncertainty/pcorrect_ic_g1_barplot_conditional.pdf', 6, 6)

bp = barplot(pcor[,,1], beside=T, ylim=c(0,1.1), ylab='P(Correct)', 
             legend.text = T,  main='Conditional Estimation',
             args.legend = list(x='topleft', bty='n', 1, title=expression(gamma[1]),ncol=1))
arrows(bp, elb[,,1], bp, eub[,,1],length=0)

dev.off()

### marginal
pdf('figs/ic_performance/w_uncertainty/pcorrect_ic_g1_barplot_marginal.pdf', 6, 6)

bp = barplot(pcor[,,2], beside=T, ylim=c(0,1.1), ylab='P(Correct)', 
             legend.text = T,  main='Marginal Estimation',
             args.legend = list(x='topleft', bty='n', 1, title=expression(gamma[1]),ncol=1))
arrows(bp, elb[,,2], bp, eub[,,2],length=0)

dev.off()

### ICs X gam2 --------

pcor = tapply(melt_dat$accuracy, list(melt_dat$G20, melt_dat$ic, melt_dat$marg), mean)
elb  = tapply(melt_dat$accuracy, list(melt_dat$G20, melt_dat$ic, melt_dat$marg), ci_binom, bound=1)
eub  = tapply(melt_dat$accuracy, list(melt_dat$G20, melt_dat$ic, melt_dat$marg), ci_binom, bound=2)

### overall
pdf('figs/ic_performance/w_uncertainty/pcorrect_ic_g2_barplot.pdf', 6, 6)

bp = barplot(apply(pcor, 1:2, mean), beside=T, ylim=c(0,1.1), ylab='P(Correct)', 
             legend.text = T, 
             args.legend = list(x='topleft', bty='n', 1, title=expression(gamma[2]),ncol=1))
arrows(bp, apply(elb, 1:2, mean), bp, apply(eub, 1:2, mean),length=0)

dev.off()

### conditional
pdf('figs/ic_performance/w_uncertainty/pcorrect_ic_g2_barplot_conditional.pdf', 6, 6)

bp = barplot(pcor[,,1], beside=T, ylim=c(0,1.1), ylab='P(Correct)', 
             legend.text = T,  main='Conditional Estimation',
             args.legend = list(x='topleft', bty='n', 1, title=expression(gamma[2]),ncol=1))
arrows(bp, elb[,,1], bp, eub[,,1],length=0)

dev.off()

### marginal
pdf('figs/ic_performance/w_uncertainty/pcorrect_ic_g2_barplot_marginal.pdf', 6, 6)

bp = barplot(pcor[,,2], beside=T, ylim=c(0,1.1), ylab='P(Correct)', 
             legend.text = T, main='Marginal Estimation',
             args.legend = list(x='topleft', bty='n', 1, title=expression(gamma[2]),ncol=1))
arrows(bp, elb[,,2], bp, eub[,,2],length=0)

dev.off()



### ICs X N X J --------

pcor = tapply(melt_dat$accuracy, list(melt_dat$J, melt_dat$N, melt_dat$ic, melt_dat$marg), mean)
elb  = tapply(melt_dat$accuracy, list(melt_dat$J, melt_dat$N, melt_dat$ic, melt_dat$marg), ci_binom, bound=1)
eub  = tapply(melt_dat$accuracy, list(melt_dat$J, melt_dat$N, melt_dat$ic, melt_dat$marg), ci_binom, bound=2)

### overall
pdf('figs/ic_performance/w_uncertainty/pcorrect_ic_n_j_barplot.pdf', 8, 4)
par(mar=c(5, 4, 6, 2) + 0.1)
layout(matrix(1:3,1,3))

for(n in as.character(unique(melt_dat$N))) {
  bp = barplot(apply(pcor, 1:3, mean)[,n,], beside=T, ylim=c(0,1.1),
               ylab='P(Correct)', xlab="",
               legend.text = T, 
               main=paste0("Cluster Size = ",n),
               args.legend = list(x='topleft', bty='n', title='Observation Size',ncol=3))
  arrows(bp, apply(elb, 1:3, mean)[,n,], bp, apply(eub, 1:3, mean)[,n,],length=0)
  
}

dev.off()

### conditional
pdf('figs/ic_performance/w_uncertainty/pcorrect_ic_n_j_barplot_conditional.pdf', 8, 4)
par(mar=c(5, 4, 6, 2) + 0.1)
layout(matrix(1:3,1,3))

for(n in as.character(unique(melt_dat$N))) {
  bp = barplot(apply(pcor, 1:4, mean)[,n,,1], beside=T, ylim=c(0,1.1),
               ylab='P(Correct)', xlab="",
               legend.text = T, 
               main=paste0("Cluster Size = ",n),
               args.legend = list(x='topleft', bty='n', title='Observation Size',ncol=3))
  arrows(bp, apply(elb, 1:4, mean)[,n,,1], bp, apply(eub, 1:4, mean)[,n,,1],length=0)
  
}

mtext('Conditional Estimation', side=3, line=-1.25, outer=TRUE, cex=1, font=3)

dev.off()

### marginal
pdf('figs/ic_performance/w_uncertainty/pcorrect_ic_n_j_barplot_marginal.pdf', 8, 4)
par(mar=c(5, 4, 6, 2) + 0.1)
layout(matrix(1:3,1,3))

for(n in as.character(unique(melt_dat$N))) {
  bp = barplot(apply(pcor, 1:4, mean)[,n,,2], beside=T, ylim=c(0,1.1),
               ylab='P(Correct)', xlab="",
               legend.text = T, 
               main=paste0("Cluster Size = ",n),
               args.legend = list(x='topleft', bty='n', title='Observation Size',ncol=3))
  arrows(bp, apply(elb, 1:4, mean)[,n,,2], bp, apply(eub, 1:4, mean)[,n,,2],length=0)
  
}

mtext('Marginal Estimation', side=3, line=-1.25, outer=TRUE, cex=1, font=3)

dev.off()


## Overall model selection -------------------------------------------------

#' create indexable array with proportions of model selection per: 
#'    1. Fit index
#'    2. Model selected
#'    3. N
#'    4. J
#'    5. Tau0
#'    6. Tau1
#'    7. gamma10
#'    8. gamma20
#' but this time, incorporate the following uncertainty metrics: 
#'    DIC : at least 4 points difference
#'    WAIC: 1 SE difference
#'    LOO : 1 SE difference
#' these are compared between top-2 models and 
#' simpler model is chosen in the case of no robust difference
#' (see data/parse_output_nested.R for implementation)

# specify model complexity per model
# (ignoring random effect estimates per cluster)

mod_select = array(NA,dim = c(length(FIT_ICS),
                              length(MODELS),
                              length(unique(dat$N)),
                              length(unique(dat$J)),
                              length(unique(dat$TAU0)),
                              length(unique(dat$TAU1)),
                              length(unique(dat$G10)),
                              length(unique(dat$G20)), 
                              length(unique(dat$marg))),
                   dimnames = list(FIT_ICS,
                                   MODELS,
                                   unique(dat$N),
                                   unique(dat$J),
                                   unique(dat$TAU0),
                                   unique(dat$TAU1),
                                   unique(dat$G10),
                                   unique(dat$G20), 
                                   unique(dat$marg))
)

for(i in FIT_ICS) {
  
  win_vec = dat[,paste0('win_1se_',i)]
  
  for(m in MODELS) {
    mod_select[i,m,,,,,,,] = tapply(win_vec,
                                    list(dat$N,
                                         dat$J,
                                         dat$TAU0,
                                         dat$TAU1,
                                         dat$G10,
                                         dat$G20, 
                                         dat$marg),
                                    function(x) mean(x==m))
    
  }
}



### overall
IC_barplot_modselect('figs/ic_performance/w_uncertainty/overall_mod_select_barplot.pdf', 
                     leg_pos = 'topright')

### conditional
IC_barplot_modselect('figs/ic_performance/w_uncertainty/overall_mod_select_barplot_conditional.pdf', 
                     est_idx = '0', main_title = 'Conditional Estimation', 
                     leg_pos = 'topright')

### marginal
IC_barplot_modselect('figs/ic_performance/w_uncertainty/overall_mod_select_barplot_marginal.pdf', 
                     est_idx = '1', main_title = 'Marginal Estimation', 
                     leg_pos = 'topright')


### Create table with proportions -------------------------------------------

prop_table = melt(mod_select, 
                  varnames = c('IC','Model',IVS),
                  value.name = 'p_choose')
write.csv(prop_table,'out/prop_table_1se.csv')


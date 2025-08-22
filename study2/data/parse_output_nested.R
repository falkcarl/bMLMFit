
rm(list=ls())

select_model = function(vals, ses=NULL, u=0, ic='DIC', 
                        n_pars= c(7, 5, 10, 6, 8) ) {
  #' Select a model based on fit index and uncertainty
  #' 
  #' A, 6 pars:   gamma00, gamma10, gamma20, tau0, tau1, tau01, sigma
  #' B, 5 pars:   gamma00, gamma10, gamma20, tau0, sigma
  #' C, 10 pars:  gamma00, gamma10, gamma20, tau0, tau1, tau2, tau01, tau02, tau12, sigma
  #' D, 6 pars:   gamma00, gamma10, tau0, tau1, tau01, sigma
  #' E, 8 pars:   gamma00, gamma10, gamma20, gamma30, tau0, tau1, tau01, sigma
  
  if(u==0) {
    return(as.numeric(which.min(vals)))
  }
  
  # pick best two
  m1 = which.min(vals)
  m2 = which.min(vals[-m1])
  m2 = which(names(vals)==names(m2)) # fix index
  
  # check if these differ based on uncertainty
  # (d==TRUE, if difference is robust)
  if(ic=='DIC') {
    d = (vals[m2] - vals[m1]) >= u
  } else {
    ci_1 = c(vals[m1]-u*ses[m1], vals[m1]+u*ses[m1])
    ci_2 = c(vals[m2]-u*ses[m2], vals[m2]+u*ses[m2])
    d = !(ci_1[1] <= ci_2[2]) && (ci_1[2] >= ci_2[1])
  }
  
  # pick best and/or most concise model
  if(d) {
    win = m1
  } else {
    # choose simplest model if CIs overlap
    # if both models are equally simple, 
      # choose model with best fitting point estimate
    win = ifelse(n_pars[m2] < n_pars[m1], m2, m1)
  }
  win
}

# Load files and prepare df -----------------------------------------------

rds_files = list.files('raw_data/',pattern='.rds',full.names = T)
n_samples = 100
conds = read.csv('conds_with_var.csv')[,c('N','J','marg','TAU0','G10',
                                          'TAU1','G20','SIG2')]
colnames = c('cond','iter',
             'marg',
             names(conds), 
             'g0s',
             'g1s',
             'g2s',
             't0s',
             't1s',
             'd_g0s',
             'd_g1s',
             'd_g2s',
             'd_t0s',
             'd_t1s',
             'max_rhats',
             'min_ess',
             'DIC',
             'WAIC',
             'LOO',
             'WAIC_se',
             'LOO_se',
             'DIC_elpd',
             'WAIC_elpd',
             'LOO_elpd',
             'DIC_p', 
             'WAIC_p',
             'LOO_p',
             paste0('win_',c('DIC','WAIC','LOO')),
             paste0('win_1se_',c('DIC','WAIC','LOO'))
             )

# preallocate
df  = matrix(NA, length(rds_files)*n_samples, length(colnames), 
             dimnames = list(NULL, colnames))
df = as.data.frame(df)

# iterate -----------------------------------------------------------------

row = 1
pb = txtProgressBar(1,nrow(df))

for(f in rds_files) {
  
  setTxtProgressBar(pb,row)
  
  x = readRDS(f)

  
  conds_x   = conds[x[[1]]$cond,]
  cond_num  = x[[1]]$cond
  cond_marg = x[[1]]$marg
  
  ## different estimation strategies have slightly different naming conventions
  ests_name = ifelse(cond_marg, 'Estimates','Estimate')

  for(i in 1:length(x)) {
    
    # mean gammas per model
    g0s = sapply(x[[i]]$fixef, function(j) j['Intercept',ests_name])
    g1s = sapply(x[[i]]$fixef, function(j) j['X1',ests_name])
    g2s = sapply(x[[i]]$fixef, function(j) {
      if('X2' %in% rownames(j)) return(j['X2',ests_name])
      return(NA)
    })
    
    # difference from true value
    d_g0s = 1-g0s
    d_g1s = conds_x[['G10']] - g1s
    d_g2s = conds_x[['G20']] - g2s
    
    # mean taus per model
    ## different estimation strategy have different naming convention
    if(cond_marg) {
      ## marginal
      t0s = sapply(x[[i]]$ranef, function(j) j$g$tau[1,'Estimates'])
      t1s = sapply(x[[i]]$ranef, function(j) {
        if('Tau[2,2]' %in% rownames(j$g$tau)) return(j$g$tau['Tau[2,2]','Estimates'])
        return(NA)
      })
      
    } else {
      ## conditional
      t0s = sapply(x[[i]]$ranef, function(j) j$g$sd['Intercept',ests_name])
      t1s = sapply(x[[i]]$ranef, function(j) {
        if('X1' %in% rownames(j$g$sd)) return(j$g$sd['X1','Estimate'])
        return(NA)
      })
      
      ## square to put on variance scale
      t0s = t0s^2
      t1s = t1s^2
    }
    
    # difference from true value
    d_t0s = conds_x[['TAU0']] - t0s
    d_t1s = conds_x[['TAU1']] - t1s
    
    # point estimates of ic
    dics  = sapply(x[[i]]$dic,function(j) j[['DIC']])
    waics = sapply(x[[i]]$waic,function(j) j['waic','Estimate'])
    loos  = sapply(x[[i]]$loo,function(j) j['looic','Estimate'])
    
    # standard errors
    waics_e = sapply(x[[i]]$waic,function(j) j['waic','SE'])
    loos_e  = sapply(x[[i]]$loo,function(j)  j['looic','SE'])
    
    # elpd
    dics_elpd  = sapply(x[[i]]$dic,function(j) j[['elpd']])
    waics_elpd = sapply(x[[i]]$waic,function(j) j['elpd_waic','Estimate'])
    loos_elpd  = sapply(x[[i]]$loo,function(j) j['elpd_loo','Estimate'])
    
    # num pars
    dics_p  = sapply(x[[i]]$dic,function(j) j[['pDIC']])
    waics_p = sapply(x[[i]]$waic,function(j) j['p_waic','Estimate'])
    loos_p  = sapply(x[[i]]$loo,function(j) j['p_loo','Estimate'])
    
    # lowest estimate
    winning_dic  = select_model(dics)
    winning_waic = select_model(waics)
    winning_loos = select_model(loos)
    
    # lowest estimate with uncertainty incorporated
    winning_1se_dic  = select_model(dics, u = 4, ic = 'DIC')
    winning_1se_waic = select_model(waics, ses = waics_e, u = 1, ic = 'WAIC')
    winning_1se_loos = select_model(loos,  ses = loos_e,  u = 1, ic = 'LOO')
    
    # summary stats on diagnostics
    max_rhat = as.vector(sapply(x[[i]]$diag, function(j) max(j[,'Rhat'], na.rm=T)))
    min_ess  = as.vector(sapply(x[[i]]$diag, function(j) min(j[,'ESS'], na.rm=T)))
      
    # save
    df[row,] = c(cond_num, 
                 i,
                 cond_marg,
                 unlist(conds_x), 
                 list(list(g0s)),
                 list(list(g1s)),
                 list(list(g2s)),
                 list(list(t0s)),
                 list(list(t1s)),
                 list(list(d_g0s)),
                 list(list(d_g1s)),
                 list(list(d_g2s)),
                 list(list(d_t0s)),
                 list(list(d_t1s)),
                 list(list(max_rhat)),
                 list(list(min_ess)),
                 list(list(dics)),
                 list(list(waics)),
                 list(list(loos)),
                 list(list(waics_e)),
                 list(list(loos_e)),
                 list(list(dics_elpd)),
                 list(list(waics_elpd)),
                 list(list(loos_elpd)),
                 list(list(dics_p)), 
                 list(list(waics_p)), 
                 list(list(loos_p)),
                 winning_dic, 
                 winning_waic, 
                 winning_loos,
                 winning_1se_dic, 
                 winning_1se_waic, 
                 winning_1se_loos
                 )
    
    row = row+1
    
  } # end of iter loop
} # end of file loop



# save --------------------------------------------------------------------


saveRDS(df, 'fits_summary.rds')




rm(list=ls())

rds_files = list.files('raw_data/',pattern='.rds',full.names = T)
n_samples = 100

conds = read.csv('conds_with_var.csv')[,c('N','J','TAU0','G10',
                                          'TAU1','G20','SIG2')]

colnames = c('cond','iter',
             names(conds), 
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
             paste0('win_',c('DIC','WAIC','LOO'))
             )
             

# preallocate
df  = matrix(NA, length(rds_files)*n_samples, length(colnames), 
             dimnames = list(NULL, colnames))
df = as.data.frame(df)

row = 1
pb = txtProgressBar(1,nrow(df))

for(f in rds_files) {
  
  setTxtProgressBar(pb,row)
  
  x = readRDS(f)
  
  conds_x  = conds[x[[1]]$cond,]
  cond_num = x[[1]]$cond
  
  for(i in 1:length(x)) {
    
    # point estimates
    dics  = sapply(x[[i]]$dic,function(j) j[['DIC']])
    waics = sapply(x[[i]]$waic,function(j) j['waic','Estimate'])
    loos  = sapply(x[[i]]$loo,function(j) j['looic','Estimate'])
    
    # standard errors
    waics_e = sapply(x[[i]]$waic,function(j) j['waic','SE'])
    loos_e = sapply(x[[i]]$loo,function(j)  j['looic','SE'])
    
    # elpd
    dics_elpd  = sapply(x[[i]]$dic,function(j) j[['elpd']])
    waics_elpd = sapply(x[[i]]$waic,function(j) j['elpd_waic','Estimate'])
    loos_elpd  = sapply(x[[i]]$loo,function(j) j['elpd_loo','Estimate'])
    
    # lowest estimate
    winning_dic  = as.numeric(which.min(dics))
    winning_waic = as.numeric(which.min(waics))
    winning_loos = as.numeric(which.min(loos))
    
    # summary stats on diagnostics
    max_rhat = as.vector(sapply(x[[i]]$diag, function(j) max(j$Rhat)))
    min_ess  = as.vector(sapply(x[[i]]$diag, function(j) min(j$ESS)))
      
    # save
    df[row,] = c(cond_num, 
                 i,
                 unlist(conds_x), 
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
                 winning_dic, 
                 winning_waic, 
                 winning_loos
    )
    
    row = row+1
    
  } # end of iter loop
} # end of file loop

# save output
saveRDS(df, 'fits_summary.rds')


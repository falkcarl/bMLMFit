
rm(list=ls())

rds_files = list.files('out',pattern='.rds',full.names = T)
n_samples = 100

conds = read.csv('conds_with_var.csv')[,c('N','J','TAU0','G10',
                                          'TAU1','G20','SIG2')]

colnames = c('cond','iter','marg',
             names(conds), 
             paste0('DIC_',LETTERS[1:5]), 
             paste0('WAIC_',LETTERS[1:5]),
             paste0('LOO_',LETTERS[1:5]), 
             paste0('WAIC_se_',LETTERS[1:5]),
             paste0('LOO_se_',LETTERS[1:5]),
             paste0('DIC_elpd',LETTERS[1:5]), 
             paste0('WAIC_elpd',LETTERS[1:5]),
             paste0('LOO_elpd',LETTERS[1:5]), 
             paste0('win_',c('DIC','WAIC','LOO')))

# preallocate
df  = matrix(NA, length(rds_files)*n_samples, length(colnames), 
             dimnames = list(NULL, colnames))
row = 1

err_sum = 0

pb = txtProgressBar(1,nrow(df))

for(f in rds_files) {
  
  setTxtProgressBar(pb,row)
  
  x = tryCatch({
    readRDS(f)
  }, error = function(e) {
    cat('error on file', f, '\n')
    return(999)
  })
  
  if(is.numeric(x)) {
    err_sum = err_sum+1
    next
  } 
  
  conds_x  = conds[x[[1]]$cond,]
  cond_num = x[[1]]$cond
  cond_marg = x[[1]]$marg
  
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
    
    # save
    df[row,] = c(cond_num, 
                 i,
                 cond_marg,
                 as.vector(unlist(conds_x)), 
                 as.vector(dics),
                 as.vector(waics),
                 as.vector(loos),
                 as.vector(waics_e),
                 as.vector(loos_e),
                 as.vector(dics_elpd),
                 as.vector(waics_elpd),
                 as.vector(loos_elpd),
                 winning_dic, 
                 winning_waic, 
                 winning_loos
                 )
    
    row = row+1
    
  } # end of iter loop
} # end of file loop

# save output
df = as.data.frame(df)
# write.csv(df, '../../analysis/data/fits_summary.csv')

# table(df[!duplicated(df$cond),'marg'])

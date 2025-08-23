rm(list=ls())

files = list.files('scratch/', full.names = T)
cond_nums = unique(unname(sapply(files, function(x) strsplit(x,'_')[[1]][2])))

shorted = c()

for(cn in cond_nums) {
  
  cond_f = files[grepl(cn, files)]
  out = list()
  for(f in cond_f) {
    x = readRDS(f)
    out = append(out, x)
  }
  
  cat('condition',cn, ':', length(out), 'samples\n')
  if(length(out) < 100) shorted = c(shorted, cn)
  
  ## update iteration number 
  for(i in 1:length(out)) {
    out[[i]]$iter = i
  }
  ## sapply(out, function(x) x$iter)
  
  ## save
  saveRDS(out, file=paste0('out/cond_',cn,'.rds'))
  
}

## test
# x = readRDS('out/cond_212.rds')
# length(x)
# sapply(x, function(x) x$iter)


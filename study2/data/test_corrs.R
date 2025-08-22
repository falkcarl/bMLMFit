

rm(list=ls())

files = list.files('raw_data/', pattern = '.rds', full.names = T)

out = matrix(NA, 
             length(files)*100, 7, 
             dimnames = list(NULL, c('cond','s','N','J','rX1X2', 'rX1X3', 'rX2X3')))
row = 1

pb = txtProgressBar(1, nrow(out))
for(f in files) {
  
  setTxtProgressBar(pb, row)
  x = readRDS(f)
  
  for(i in 1:length(x)) {
    
    dat     = x[[i]]$data
    n       = length(unique(dat$g))
    j       = unique(table(dat$g))
    rx1x2   = cor(dat$X1, dat$X2)
    rx1x3   = cor(dat$X1, dat$X3)
    rx2x3   = cor(dat$X2, dat$X3)
    
    out[row, ] = c(x[[i]]$cond, i, n, j,rx1x2,rx1x3,rx2x3)
    row = row +1
    
  } ## end of sample loop
} ## end of file loop

out = as.data.frame(out)

# -------------------------------------------------------------------------



for(r in c('rX1X2', 'rX1X3', 'rX2X3')){
  
  pdf(paste0(r, '.pdf'), 10, 10)
  layout(matrix(1:9, 3, 3, byrow = T))
  
  for(n in unique(out$N)) {
    for(j in unique(out$J)) {
      sub = out[out$N==n & out$J==j,r]
      if(length(sub)==0) next
      hist(sub, 
           xlim=c(-1,1), ylim=c(0,100), 
           main=paste0('N = ',n, '\nJ = ', j), 
           xlab=paste0('Cor(',gsub('r','',r), ')'))
      abline(v=0.4, lty=2)
    }
  }
  
  dev.off()
}


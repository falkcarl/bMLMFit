
# Run all batch files in bash_files/
## only run marginal estimation for now to fix bug
files = list.files('bash_scripts',full.names = T)

to_rerun = read.csv('to_rerun.csv', header = F)[,1]
to_rerun = as.numeric(unlist(regmatches(to_rerun, gregexpr("[0-9]+", to_rerun))))

for(i in seq(files)) {
  if(as.numeric(unlist(regmatches(files[i], gregexpr("[0-9]+", files[i])))) %in% to_rerun) {
    ## only rerun those that need longer run times
    system(paste("sbatch ", files[i],sep=""))
    Sys.sleep(2) # give some time between run calls
  }
}

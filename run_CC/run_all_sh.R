
# Run all batch files in bash_files/
files = list.files('bash_scripts',full.names = T)

for(i in seq(files)) {
  system(paste("sbatch ", files[i],sep=""))
  Sys.sleep(2) # give some time between run calls
}

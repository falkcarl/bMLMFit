#!/bin/bash -login

## walltime HH:MM:SS
#SBATCH -t 168:00:00

## nodes
#SBATCH --cpus-per-task=2

## mem
#SBATCH --mem-per-cpu=10G

## name of job
#SBATCH --job-name=bMLMfit_cond_208

## account to use
#SBATCH --account=def-cffalk

## email notification
#SBATCH --mail-user=seandamiandevine@gmail.com
#SBATCH --mail-type=NONE

## correct version of R
module load r/4.4.0 gcc 

## cd to the directory from which sbatch was run from
cd $SLURM_SUBMIT_DIR

## execute job
R CMD BATCH --no_save r_scripts/cond_208.R console_out/bMLMfit_cond_208_out.txt

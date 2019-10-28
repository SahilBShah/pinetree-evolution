#!/bin/bash
#
# Simple SLURM script for submitting multiple serial
# commands (e.g. parametric studies) using a script wrapper
# to launch the commands.
#
# To use, change this job script to accommodate 
# running your serial application(s) in your WORKDIR 
# directory (usually the directory of submission).
# Edit the commands file to specify the executions
# for the launcher (paramrun) to perform.
#-------------------------------------------------------
#-------------------------------------------------------
# 
#         <------ Setup Parameters ------>
#
#SBATCH -J launcher_test
#SBATCH -N 1
#SBATCH -n 65             #use site recommended # of cores
#SBATCH -p normal
#SBATCH -o launcher_test.o%j
#SBATCH -e launcher_test.e%j
#SBATCH -t 20:00:00
#SBATCH -A A-bio7   #uncomment and insert acct name if necessary
#------------------------------------------------------

#                         # USING SLURM; plugins defines SLURM env. vars.
export LAUNCHER_RMI=SLURM
export LAUNCHER_PLUGIN_DIR=$LAUNCHER_DIR/plugins


#                         # JOB_FILE is a list of executions to run

export LAUNCHER_JOB_FILE=`pwd`/tacc_grant_rerun
export LAUNCHER_SCHED=interleaved
export LAUNCHER_WORKDIR=`pwd`
 
$LAUNCHER_DIR/paramrun    # will run the executions in the LAUNCHER_JOB_FILE file
                          # "JOB" is a misnomer--these are not slurm jobs
                          # Each line in the commands file is an execution.

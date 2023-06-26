#! /bin/tcsh -f

# This script wakes up in a directory containing pre-formatted .svm files (with or without added decoys) having activity in the first column
# The last element in the chromosome is the name of the descriptor set to employ

source $GACONF/common.pars
set mode=denovoDockpareto
set pid=$$
source $GACONF/argproc.cmd

touch $pidlist/$pid >& /dev/null
if (!($?runit)) then
  if ($?SLURM_JOB_ID) then
    set tmpdir=/scratch/job.$SLURM_JOB_ID
    if (!(-d $tmpdir)) mkdir $tmpdir
    set runit=$SLURM_JOB_ID.$pid
  else
    set runit=$pid
  endif
else
  if ($?SLURM_JOB_ID) then
    set runit=$SLURM_JOB_ID.$runit
  endif
endif


# Get a current model number
source $GACONF/mktempdir.cmd

# Install EVO environment into the new dir
ln -s $EVOSRC/* $dir/
$GACONF/getChromo.csh dir=$dir mode=$mode attempt=$attempt 
# The chromosome contains the SMILES
set old=$PWD
cd $dir
$GACONF/smidockPLANTS.csh ligfile=chromo recdir=$recdir dock=dock
if (-z dock1/score) then
  echo `cat chromo` = $dir 0 -9999 -9999 -9999  > $old/done_so_far.temp$pid
else
  echo `cat chromo` = $dir 0 `awk '{printf "%.2f %.2f %.2f\n",-1*$1,-1*$3,-1*$2}' dock1/score`  > $old/done_so_far.temp$pid
endif
cd $old

rm -r $dir $pidlist/$pid >& /dev/null

#

#! /bin/tcsh

set par=8  # By default, try to use 8 cores for descriptor preprocessing
set minvar=0.02
source $GACONF/argproc.cmd
set parfile=$GACONF/$mode.rng

rm no_go >& /dev/null
if (!($?workdir)) then
  echo FATAL - need to specify a working directory name workdir=path_to_working_directory. It will be created unless you specified cont=y| tee -a no_go
  exit
endif

if ((-e $workdir) && !(-d $workdir))  then
  echo FATAL - cannot create a directory $workdir, because a file of same name exists| tee -a no_go
  echo remove it, or use workdir=other_location on command line| tee -a no_go
  exit
endif

if (!(-e $parfile) || (-z $parfile)) then
  echo FATAL - $parfile is not a valid definition of the GA problem space
  exit
endif

if ((-d $workdir) && ($cont:s/y// == $cont)) then
  echo FATAL - there already is a workdir $workdir, but you have not used the cont=yes or cont=y flag| tee -a no_go
  exit
else if (!(-d $workdir) && !($cont:s/y// == $cont)) then
  echo FATAL - you tried to contine a run in $workdir, but there is no such directory | tee -a no_go
  exit
endif

if (!($cont:s/y// == $cont)) then
  touch $workdir/best_pop $workdir/done_so_far
  exit
endif

if (!(-d $data_dir)) then
  echo WARNING - you have not defined a valid data_dir - it is assumed that your run does not concern model building based on descriptor files to import
  echo This will just create an empty $workdir and if applicable copy $parfile into it
  mkdir $workdir
  cp $parfile $workdir >& /dev/null
  touch $workdir/best_pop $workdir/done_so_far
  exit
endif


# Set flag "run in svm mode"
if (($mode == $mode:s/GTM//) && ($mode == $mode:s/BTM//) && ($mode == $mode:s/NTM//) && ($mode == $mode:s/SOM//)) then
  set svmode=1
else
  set svmode=no #By default, do not enter property into svm files at preproc stage
endif

if (!(-d $workdir)) then

  set propfile=(`ls -1 $data_dir/*.$mode`)
  if (!($#propfile == 1)) then
    echo FATAL - $data_dir must contain one and only one .$mode file with activities in the only column of that text file - no title line| tee -a no_go
    exit
  endif

  set descfiles=(`ls -1 $data_dir/*.svm`)
  if ($#descfiles == 0) then
    echo FATAL: $data_dir should also contain some .svm descriptor files - none found| tee -a no_go
    exit
  endif

  # Fix min-max of property range
  if (-e $data_dir/proprange.dat) then
   set proprange=(`cat $data_dir/proprange.dat`)
  else
    set proprange=(`awk 'BEGIN {min=+99999;max=-min} {if ($1>max) max=$1;if ($1<min) min=$1} END {r=max-min;print min-0.1*r,max+0.1*r}' $propfile`)
  endif



  # Prepare preliminary stuff in $workdir
  mkdir $workdir >& /dev/null
  touch $workdir/best_pop $workdir/done_so_far

  awk -v i=1 '{nlin ++;s+=$i;s2+=$i*$i} END {print sqrt(s2/nlin-s*s/nlin/nlin)}' $propfile > $workdir/data.stdev
  if (!($mode == $mode:s/class//)) then
    awk -v i=1 '{if (!($i in seen)) {n++;seen[$i]=1}} END {print 1.0/n}' $propfile > $workdir/minimal_useful_fitlevel
  else
    echo 0.0 > $workdir/minimal_useful_fitlevel
  endif
  cp $data_dir/proprange.dat $workdir >& /dev/null
  set descopt=descriptors
  foreach svm ($data_dir/*.svm)
    if (!($decoy_dir == nodecoys) && !(-e $decoy_dir/$svm:t)) then
      echo FATAL - you have specified a decoy set, but the equivalent decoy file $decoy_dir/$svm:t of $svm does not exist| tee -a no_go
      rm -r $workdir
      exit
    endif

    if ($par<2) then
      $GACONF/prepProcDesc.csh svm=$svm workdir=$workdir data_dir=$data_dir propfile=$propfile svmode=$svmode proprange=$proprange[1] decoy_dir=$decoy_dir prune=$prune minvar=$minvar
    else
      echo $GACONF/prepProcDesc.csh svm=$svm workdir=$workdir data_dir=$data_dir propfile=$propfile svmode=$svmode proprange=$proprange[1] decoy_dir=$decoy_dir prune=$prune minvar=$minvar >> $workdir/par.cmd
    endif

  end

  if ($par>1) $GACONF/$MACHTYPE/parallel -j $par -a $workdir/par.cmd

  foreach var (.orig.EDprops _pruned.orig.EDprops)
    foreach f (`ls -1 $workdir/*$var`)
      set descopt=($descopt $f:t:r:r)
    end
  end

  echo $descopt > $workdir/$parfile:t
  cat $parfile >>  $workdir/$parfile:t

  # In GTM mode, check whether mapsets are enabled
  if ($svmode == no) then
    set mapsets=(mapsets `awk '{n=split($1,w,"[+]");for (i=1;i<=n;i++) {if (!(w[i] in seen)) {print w[i];seen[w[i]]=1}}}' $workdir/*.orig.svm`)
    if ($#mapsets > 1) echo $mapsets >>  $workdir/$parfile:t
  endif

endif


# Exit

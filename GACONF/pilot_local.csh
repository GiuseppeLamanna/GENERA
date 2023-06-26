#! /bin/tcsh

# This script pilots force field fitting based on a genetic algorithm in the parameters space defined by $parfile
#set nnodes=50
set nnodes=1
set deploy=local
source $GACONF/common.pars
#set nc=1
set nc = `nproc`      # use all CPUs by default #Careful, modified!

# Specify data sources and optionally Override default pars
# by passing par=value as command line option when calling this script
# Use cont=yes to restart an unexpectedly terminated run (due to an effing power cut, for example)
source $GACONF/argproc.cmd
if (!($?job)) then
  set job=$GACONF/local_$mode.csh #The piloted script running the actual SVM challenge at picked configuration
endif
if (($chromdiv == chromdiv.awk) && (-e $GACONF/chromdiv_$mode.awk)) set chromdiv=chromdiv_$mode.awk

# Prepare workdir
if (!(-e $workdir/$mode.rng) || (-z $workdir/$mode.rng)) cp $GACONF/$mode.rng $workdir
if  ( ($?continue) || ($?use_chromo) )  set cont=yes
$GACONF/prepWorkDir.csh workdir=$workdir data_dir=$data_dir decoy_dir=$decoy_dir lowscale=$lowscale upscale=$upscale cont=$cont prune=$prune mode=$mode par=$nc deploy=$deploy
if (-e no_go) exit

# Make sure that variables pointing to directories or files are reported as absolute pathes. Otherwise (not beginning with a /) append $PWD
# This is important because scripts may change directories, and would get lost if relative path names were used
foreach var (use_chromo target recdir)
  sed s/VARNAME/$var/g $GACONF/absolutize_pathdata.template > my.src
  source my.src
end
rm my.src
cd $workdir
rm -r $pidlist >& /dev/null
mkdir $pidlist
if (!(-e $mode.prob) && (-x $GACONF/prepare_$mode.rng.csh)) then
  echo Running $GACONF/prepare_$mode.rng.csh ...
  $GACONF/prepare_$mode.rng.csh $argv
endif
if (!(-e $mode.rng) || (-z $mode.rng)) then
  echo FATAL - failed to find or create $mode.rng
  exit
endif
rm stop_now new_top>& /dev/null
if ($?wait) exit

# Detect the nr. of parameters to optimize: note that in the done_so_far report, these npar parameters will be listed next to the
# subdir name in which the model building attempt happened, followed by the fitting and cross-validation scores. The last column of that
# file is the cross-validation score, at position npars + 4 because an equal sign separates parameters and output on the line
@ keycol = 0
if (!($?minaccept))  set minaccept=-99999

# To rebuild models based on previously checked setup schemes
if ($?use_chromo) then
  @ nreb = `cat $use_chromo|wc -l`
  rm rebuild.cmd >& /dev/null
  foreach k (`seq 1 $nreb` )
    echo "$job $argv runit=REBUILD.$k use_chromo=$use_chromo" >> rebuild.cmd
  end
  echo Submitted $nreb model rebuilding jobs, which will be found in directories suffixed as .REBUILD.k when finished, with k=1...$nreb
  $GACONF/mparallel.csh cmdfile=rebuild.cmd clean=y
  rm -r $pidlist >& /dev/null
  exit
endif

@ safetrials = 2 * $maxconfigs
rm -r pilot.cmd >& /dev/null
foreach k (`seq 1 $safetrials`)
  echo "$job $argv  >& /dev/null" >> pilot.cmd
end
$GACONF/mparallel.csh nslots=$nnodes cmdfile=pilot.cmd deploy=$deploy nowait=y clean=y


set topFit=-99999
@ latest = 0
touch best_pop done_so_far
rm new_top >& /dev/null
sleep 60

while ((`cat done_so_far|wc -l` < $maxconfigs) && !(-e stop_now) && (`echo $topFit $goodFit| awk '{if ($1>=$2) {print 1} else {print 0}}'` == 0) && ($latest < $noprog))

  set news=(`ls -1 done_so_far.temp*`)
  @ kwait = 0
  while (($kwait < 3000) && ($#news < $minNovel))
    sleep 5
    @ kwait ++
    set news=(`ls -1 done_so_far.temp*`)
  end
  if ($keycol == 0) then
    # First time you see results, set keycol as the last field in done_so_far
    @ keycol = `awk '{print NF}' done_so_far $news| sort -nr| head -1`
  endif
  awk -v attempt=$attempt -v ncol=$keycol 'NF==ncol {l=$0;sub(attempt"[.][^ ]*","",l);if (!(l in seen)) {seen[l]=1;print $0}}' done_so_far $news > dsf.temp
  mv dsf.temp done_so_far
  #@ nok = `awk -v minaccept=$minaccept -v minpop=$minpop '{print $0 > "done_so_far";if (1*$NF>1*minaccept) n++} END {if (n<2*minpop) n=2*minpop;print n}' dsf.temp`
  rm $news >& /dev/null
  if ($mode:s/pareto// == $mode) then
    # Standard one-objective in last col
    awk -v ncol=$keycol 'NF==ncol' done_so_far| sort -gr -k$keycol | uniq > sorted_so_far
    awk -v fract=$faccept -v minpop=$minpop -v topfit=$topFit -v minaccept=$minaccept -f $GACONF/$chromdiv  sorted_so_far > bp.temp
    @ latest = `cat sorted_so_far| wc -l` - `awk '{print $2}' new_top`
    set topFit=`awk '{print $1}' new_top`
    mv sorted_so_far done_so_far
    mv bp.temp best_pop
    @ focus = $keycol
  else  if (`cat done_so_far| wc -l` > $minpop) then
    cp best_pop best_pop.old >& /dev/null
    @ nb = `cat best_pop.old|wc -l`
    if (($dom>0) && ($nb > $maxbreed)) then
      # Too many in best_pop: reduce $dom
      @ dom --
    else if ($nb < $minbreed) then
      # Too few in best_pop: increase dom and therefore undo the filter.awk
      rm best_pop.old filter.awk >& /dev/null
      @ dom ++
    endif
    echo $dom > dominance
    if ((-e best_pop.old) && !(-z best_pop.old)) then
      awk '{for (i=NF;$i==1*$i;i--) {if ($i<min[i] || min[i]=="") min[i]=$i; if ($i>max[i] || max[i]=="") max[i]=$i};i++} END {if (i=="" || i<2) {print " {print $0}"} else {printf "($%d>=%s",i,min[i];for (k=i+1;k<=NF;k++) printf " && $%d>=%s",k,min[k]; printf " ) || $%d>%s",i,max[i]; for (k=i+1;k<=NF;k++) printf " || $%d > %s",k,max[k];print " {print $0}"}}' best_pop.old > filter.awk
      awk -f ./filter.awk done_so_far| $GACONF/$MACHTYPE/pareto -i /dev/stdin -o best_pop.new -d $dom
    else
      $GACONF/$MACHTYPE/pareto -i done_so_far -o best_pop.new -d $dom
    endif
    mv best_pop.new best_pop
    set topFit = 0
    set goodFit = 1
    if (`diff -q best_pop best_pop.old|& wc -l` > 0) then
      @ latest = 0
      echo new_front@: `cat done_so_far|wc -l` > new_top
    else
      @ latest ++
    endif
    @ focus = 0
  endif #end if pareto mode
  # Migrate if appropriate
  if (($?migrantDest) && (-e best_pop) && !(-z best_pop))  $GACONF/migrate.csh $argv focus=$focus
end  # loop over total explored configs

echo ALL done
if (-e terminator.cmd) source terminator.cmd
echo Remaining processes `ls -1 $pidlist/*` killed
foreach p (`ls -1 $pidlist/*`)
  if (`ps --pid $p:t --no-header| wc -l` == 1) kill -TERM $p:t
end
echo Remaining processes `ls -1 $pidlist/*` killed
rm -r $pidlist
mkdir loosers
foreach d ( $attempt.* )
  if (!(`fgrep $d best_pop | wc -l` == 1)) then
    mv $d/stat loosers/$d.stat
    mv $d/extval loosers/$d.extval >& /dev/null
    rm -f -r $d
  endif
end
gzip loosers/*

#

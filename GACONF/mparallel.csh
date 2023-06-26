#! /bin/tcsh -f

#Takes a command file and dispatches it on a cluster using qsub commands
#set nslots=10  One may provide a default for nslots, but it's better not to and take the nr of lines in cmdfile as the default
#However, if one specifies nnodes=NNN on command line, it will submit only nnodes qsub jobs, each containing several commands

set template=$GACONF/template.sh
set cmdfile=none
set deploy=local
set wtime=60
set pid=$$
source $GACONF/argproc.cmd

if (!(-e $cmdfile) || (-z $cmdfile)) then
  echo FATAL - this needs a file of commands to be passed as cmdfile=CCCCC
  exit
endif

if (!($?nslots)) then
  set nslots=`cat $cmdfile|wc -l`
endif

if ($deploy == local) then
  if (`nproc` < $nslots) set nslots=`nproc`
  $GACONF/$MACHTYPE/parallel -j $nslots -a $cmdfile >& /dev/null &
  set childpid=$!
  if ($?nowait) then
    if ($?clean) then 
      echo kill -TERM $childpid > terminator.cmd
      echo kill -TERM $childpid >> terminator.cmd
      sed s/PARENTPID/$childpid/ $GACONF/killtree.template >> terminator.cmd
    endif
    exit
  else
    wait
  endif
else if ($deploy == slurm) then
  source $GACONF/slurm.oppars
  mkdir mpar.$pid
  foreach k (`seq 1 $nslots`)
    awk -v k=$k -v n=$nslots 'NR-n*int(NR/n)==k-1' $cmdfile>mpar.$pid/cmd$k
    sbatch  -A $account -p $queue -N 1 -n $nc $GACONF/$MACHTYPE/parallel -j $nc -a mpar.$pid/cmd$k|awk '{print $NF}' >> mpar.$pid/joblist
  end
  if (!($?nowait)) then
    @ running = $nslots
    while ($running > 0)
      sleep $wtime
      @ running = `squeue -u $USER | awk -v jfile=mpar.$pid/joblist 'BEGIN {while(getline<jfile>0) job[$1]=1} $1 in job && ($5=="R" || $5=="PD") {n++} END {print n}'`
    end
    if ($?clean) then
      foreach j (`cat mpar.$pid/joblist`)
        rm slurm-$j.out >& /dev/null
      end
      rm -r mpar.$pid
    endif
  else
    if ($?clean) then
      echo scancel `cat mpar.$pid/joblist` > terminator.cmd
      foreach j (`cat mpar.$pid/joblist`)
        echo rm slurm-$j.out >> terminator.cmd
      end
      echo rm -r mpar.$pid >> terminator.cmd
    endif
  endif


else if ($deploy == qsub) then
  mkdir mpar.$pid
  foreach k (`seq 1 $nslots`)
    awk -v k=$k -v n=$nslots 'NR-n*int(NR/n)==k-1' $cmdfile > mpar.$pid/cmd$k
  end
  qsub -t 1-$nslots":"1 -V $template mpar.$pid/cmd|awk '{print int($3)}' >> mpar.$pid/joblist
  if ($?nowait) then
    echo Preparing to exit - creating terminator.cmd
    rm terminator.cmd >& /dev/null
    foreach j ( `cat mpar.$pid/joblist`)
      echo qdel $j >> terminator.cmd
    end
    if ($?clean) then
      foreach j (`cat mpar.$pid/joblist`)
        echo "rm *.o$j *.e$j" >> terminator.cmd
      end
      echo "rm -r mpar."$pid >> terminator.cmd
    endif
    exit
  endif
  @ running = $nslots
  while ($running > 0)
    sleep $wtime
    @ running = `qstat -u $USER | awk -v jfile=mpar.$pid/joblist 'BEGIN {while(getline<jfile>0) job[$1]=1} $1 in job && (toupper($5)=="R" || toupper($5)=="W") {n++} END {print 1*n}'`
  end
  if ($?clean) then
     foreach j (`cat mpar.$pid/joblist`)
      rm *.o$j *.e$j >& /dev/null
    end
    rm -r mpar.$pid
  endif
endif

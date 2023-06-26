#! /bin/tcsh

@ focus = 0
source $GACONF/common.pars
source $GACONF/argproc.cmd
if (`awk -v f=$migfrq 'BEGIN {srand();if (rand()<f) {print 1} else {print 0}}'` == 0) exit
rm candidate >& /dev/null
awk -v f=$focus '{for (i=NF;(i>=f) && ($i==1*$i);i--) if (1*$i>max[i] || max[i]=="") {max[i]=$i;top[i]=$0;sub("[ ]+=.*","",top[i])}} END {for (i in top) pick[top[i]]=1;while(getline<"expats">0) if ($0 in pick) delete pick[$0];for ( p in pick) print p}' best_pop| head -1 > candidate
if (!(-z candidate)) then
  cat candidate >> expats
  scp candidate $migrantDest/migrant >& /dev/null &
  #set childpid=$!
  #sleep 2  # Give the scp 2 seconds to complete
  #kill $childpid >& /dev/null
  #rm candidate
endif

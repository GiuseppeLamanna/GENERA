#! /bin/tcsh -f

@ ntrials = 20
source $GACONF/argproc.cmd
cd $dir

if (!(-e best_pop)) then
  echo FATAL-best_pop-missing > chromo
  exit
endif

awk -v n=`cat best_pop| wc -l| $GACONF/$MACHTYPE/rand` 'NR==n {print $1}' best_pop > ref.smi
python3 $SCRIPTS/removeStereochemistry.py ref.smi ref.smi
@ len = `awk '{print length($1)}' ref.smi`
@ t = 0
cat ref.smi >> config_history
ls -lt
while ( (!(-e chromo) || (-z chromo)) && ($t < $ntrials))
  @ t ++

  # Generate a smiles
  rm chromo >& /dev/null
  #S=5 al massimo
  set mutlist=`awk -v l=$len 'BEGIN {srand(PROCINFO["pid"]);nmax=int(l/2);if (nmax>5) nmax=5;mutlist="";while (mutlist=="") {delete mut;nmut=1+int(rand()*nmax);for (i=1;i<=nmut;i++) mut[int(1+rand()*l)]=1;for (m in mut) mutlist=mutlist""m","};sub(",$","",mutlist);print mutlist}'`
  
  bash ./run `cat ref.smi` smiles 1 "[$mutlist]"
  echo $mutlist
  cat smiles.csv > ref.multi
  if (!(-z ref.multi)) then
    # First standardize, then check uniqueness

    python3 $SCRIPTS/standardize.py ref.multi ref.multi
    awk ' {print $NF}' ref.multi > chromo
    # Check if unique & recognizable as SMILES
    awk 'BEGIN {while (getline<"ref.multi">0) if ($NF!="O") this[$NF]=1} {if ($1 in this) {delete this; exit}} END {for (s in this) print s}' config_history > chromo

    # Check if free of bad groups
    if (!(-z chromo)) then
	python3 $SCRIPTS/removeBadGroups.py chromo $SCRIPTS/bad_groups good.smi
	awk 'NR>1' good.smi > chromo
    endif
  endif

end

# If there is no chromo at this point, too bad

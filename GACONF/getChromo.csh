#! /bin/tcsh

@ npop = 1
set sex=X
set sorted=n
source $GACONF/argproc.cmd

@ ndf = `awk '{if ($1=="end") ref=NR} END {print NR-ref}' $mode.rng`

if ((-e migrant) && !(-z migrant)) then
  mv migrant $dir/chromo
  if (!(`cat $dir/chromo|wc -l` == 1) || !(`awk 'NR==1 {print NF}' $dir/chromo` == $ndf)) then
    rm $dir/chromo
  else
    cat $dir/chromo >> expats   #Mark this as expat, to make sure it won't be reexported
  endif
else if (($?use_chromo) && ($?line)) then
    awk -v line=$line 'line==1*line && NR==line' $use_chromo > $dir/chromo
else if ($?use_chromo) then
    # Implicitly assume first line
    head -1  $use_chromo >  $dir/chromo
endif


if (!(-e $dir/chromo) || (-z $dir/chromo)) then
  # In Darwin we trust
  if ($mode:s/denovo// == $mode) then
    sed 's/=.*//' done_so_far | awk -v ndf=$ndf 'NF==ndf' > $dir/config_history
  else
    awk '{print $1}' done_so_far > $dir/config_history
  #find $dir/../ -name "chromo" -exec awk -v ndf=$ndf 'NF==ndf' {} + >> $dir/config_history  #This is bad because it finds ANY hidden chromos
  if ($mode:s/SEL// == $mode) then
    awk -v ndf=$ndf 'NF==ndf' $dir/../*/chromo >> $dir/config_history
  else
    awk -v ndf=$ndf 'NF==ndf' $dir/../*/selchromo >> $dir/config_history
  endif
  if ($mode == $mode:s/denovo//) then
    sed 's/=.*//' best_pop | awk -f $GACONF/make_children$sex.awk -v parfile=$mode.rng -v npop=$npop -v histfile=$dir/config_history -v sorted=$sorted > $dir/chromo
  else
    cp best_pop $dir
    $GACONF/denovo.csh dir=$dir
  endif
endif
rm $dir/config_history >& /dev/null

# Create acolumn version of the line chromosome and check
awk -f $GACONF/mat_transp.awk $dir/chromo > $dir/chromoT
if ( !( `cat $dir/chromoT| wc -l` == $ndf ) ) then
  echo FATAL - faulty chromosome...
  #rm -f -r $dir
  exit
endif

#end

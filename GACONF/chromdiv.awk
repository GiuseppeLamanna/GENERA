BEGIN {
    if (minaccept=="") minaccept=-99999; #do not exclude any items because of insufficinet fitness
    if (keep=="") keep=9999; #keep all diverse pushed thru the filter
    if (recsep=="") recsep="="; #process chromosome-fitness line until a recsep is found
    if (maxocc=="") maxocc=1;   #A valueappearing at least maxocc times in the selection does no longer count as "novel"
    if (minnew=="") minnew=3;   #Require at least 3 novel values not present in the chromosome before
    kept=0;
    if (maxfit=="") maxfit=-999999;
    if (fract=="") fract=0.7;
    if (minpop=="") minpop=50;
    nlin=0;
    # Stop reading through the pile after 10K solutions, unless minpop is not very large
    nmax=2*minpop;
    if (nmax<10000) nmax=10000;
}

($0~"=" && NR<nmax) {
    if (nlin==0) {
        if (keycol=="") keycol=NF;
        if (keycol>0 && keycol<=NF) maxfit=$keycol;
        ncols=NF;
    }
    if (NF==ncols) {
        nlin ++;
        lastpos[$keycol]=nlin;  #Mark the latest line on which this fitness score has been first encountered
        found=0;
        line="";
        for (i=1;(i<=NF && $i!=recsep) ;i++) {
            line=line":"$i;
            if (1*nocc[i":"$i]<maxocc) found ++;
        }
        if (found>=minnew && !(line in saved) && ($NF*1>minaccept)) {
            for (i=1;(i<=NF && $i!=recsep);i++) {
                nocc[i":"$i]++;
            }
            kept ++;
            saved[line]=1;
            chromo[kept]=$0;
            fit[kept]=$keycol;
        }
    }
}

END {
    if (topfit=="" || maxfit>topfit) print maxfit,NR > "new_top";
    if (nlin>minpop && kept>1) {
        # Get the fitness cutoff corresponding to better than fract
        for (f in lastpos) {
            ferr=1-1*lastpos[f]/nlin-fract;
            ferr=ferr*ferr;
            if (best=="" || ferr<1*best) {
                best=ferr;
                cutoff=f;
            }
        }
        for (k=1;k<=kept;k++) if (1*fit[k]>=1*cutoff) print chromo[k];
    }
}

BEGIN {
    if (npop=="") npop=1;
    if (parfile == "") parfile=ENVIRON["FFDIR"]"/param.rng";
    if (histfile == "") histfile="config_history";
    if (coupled_runs == "") coupled_runs="coupled_runs.lst";
    if (config2fit_sep=="") config2fit_sep="=";
    if (min_for_sex=="") min_for_sex=100;
    stderr="/dev/stderr";

    
    nvars=0;
    while (getline<parfile>0) {
	nvars ++;
	mode[nvars]=$NF
	if ($NF=="R") {
            #Range mode - min max prec
	    min[nvars]=$2;
	    max[nvars]=$3;
	    prec[nvars]="%"$4"f";
	    np=$4;
	    sub("[0-9][.]","",np);
	    nopts[nvars]=int((max[nvars]-min[nvars])*exp(np*log(10))+0.5)+1;
	} else if ($NF=="E") {
            #Enumeration mode
	    nopts[nvars]=NF-2;
	    for (i=2;i<=NF-1;i++) {opt[nvars,i-1]=$i};
	    prec[nvars]="%s";
	} else if ($NF=="N") {
            #Open-ended numerical
	    init[nvars]=$2;
            step[nvars]=$3;
	    maxmv[nvars]=$4;
	    min[nvars]=-99999;
	    if ($5==1*$5) min[nvars]=$5;
	    max[nvars]=+99999;
	    if ($6==1*$6) max[nvars]=$6;
	    nn=split(step[nvars],w,"[.]");
	    if (nn==2) {
		prec[nvars]="%"(1+length(w[1])+length(w[2]))"."length(w[2])"f";
	    } else if (nn==1) {
		prec[nvars]="%"length(step)"d";
	    }
            nopt[nvars]=int((max[nvars]-min[nvars])/step[nvars]);
        } else {
	   print "FATAL - option "$NF" not supported in rng file" >"/dev/stderr";
           exit;
        }
    }
    close(parfile);
    
    nhist=0;
    while (getline<histfile>0) {
	nhist ++;
	l=$0;
	sub(config2fit_sep".*$","",l);
	found[l]=nhist;
    }
    close(histfile);
    
    while (getline<coupled_runs>0) {
	histf=$1;
	while (getline<histf>0) {
	    nhist ++;
	    l=$0;
	    sub(config2fit_sep".*$","",l);
	    found[l]=nhist;
	}
	close(histf);
    }
    close(coupled_runs);
    
    srand(PROCINFO["pid"]);
}


{
    inchromo="";
    for (i=1;i<=nvars;i++) {
	x[NR,i]=$i;
	fstr="%s "prec[i];
	inchromo=sprintf(fstr,inchromo,$i);
    }
    found[inchromo]=1;
}


END {
    
    # For N variables, occasionally let the departure point move to one of the chromosomes
    ovr_init=-NR+int(2*NR*rand()+0.999);
    if (ovr_init>0) for (i=1;i<=nvars;i++) init[i]=x[ovr_init,i];
 
    new=0;
    ntrials=0;
    while (new<npop) {
	
	# Draw a random number and decide whether it will be a mutation or a cross_over
	
	ntrials ++;
	delete val;
	rand_threshold=(NR-2)/min_for_sex;
        if (rand_threshold>0.8) rand_threshold=0.8;
	if (rand()>rand_threshold) { 
            # Thre aint enough heredity to play GA yet, do random	    
	    # Get a fully random chromosome
	    print "Fully Random chromosome..." > stderr;
	    for (j=1;j<=nvars;j++) {
		if (mode[j]=="R") {
		  val[j]=min[j]+(max[j]-min[j])*rand();
		 } else if (mode[j]=="E"){
		    jj=1+int(nopts[j]*rand());
		    val[j]=opt[j,jj];
		} else {
                    val[j]="";
                    while (val[j]=="") {
			s=int(maxmv[j]*rand());
			if (rand()<0.5) s=-s;
                        val[j]=init[j]+s*step[j];
		        if (val[j]>max[j] || val[j]<min[j]) val[j]="";
                    }
                }
	    }

	} else if (rand()>0.40) {
	    
	    print "Attempting Cross-Over.." > stderr;
	    # It's a cross-over
	    p1=1;
	    p2=1;
	    while (p1==p2) {
		p1=int(1+NR*rand()*rand());
		p2=int(1+NR*rand());
	    }
	    
	    pos=2+int((nvars-2)*rand());
	    
	    for (j=1;j<=pos;j++) {
		val[j]=x[p1,j];
	    }
	    for (j=pos+1;j<=nvars;j++) {
		val[j]=x[p2,j];
	    }
	    
	} else {
	    
	    print "Attempting Mutation.." > stderr;
	    # Mutation!
	    p1=int(1+NR*rand());
	    pos=1+int(nvars*rand());
            while (nopts[pos] == 1) pos=1+int(nvars*rand());
	    
	    for (j=1;j<=nvars;j++) {
		val[j]=x[p1,j];
		if (j == pos) while(val[j]==x[p1,j]) {
  		  if (mode[j]=="R") {
		    val[j]=min[j]+(max[j]-min[j])*rand();
		   } else if (mode[j]=="E"){
		      jj=1+int(nopts[j]*rand());
		      val[j]=opt[j,jj];
		  } else {
		      s=1+int(maxmv[j]*rand());
		      if (rand()<0.5) s=-s;
		      val[j]=x[p1,j]+s*step[j];
		      if (val[j]>max[j] || val[j]<min[j]) val[j]=x[p1,j];
                  }
               }
	    }
	} 

	
	
	fstr=prec[1];
	chromo=sprintf(fstr,val[1]);;
	for (j=2;j<=nvars;j++) {
	    fstr="%s "prec[j];
	    chromo=sprintf(fstr,chromo,val[j]);
	}
	
	gsub("  *"," ",chromo);
	if (!(chromo in found) || (ntrials>10000*npop)) {
	    new ++;
	    print chromo;
	    print chromo >> histfile;
	    print new" new chromosomes were output" > stderr;
	    nhist ++;
	    found[chromo]=nhist;
	} else {
	    print chromo" discarded because already seen" > stderr;
	}
    }
}




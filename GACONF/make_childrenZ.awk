function num_compare(i1,v1,i2,v2) {
  if(1.0*v1<1.0*v2) {
    return 1;
  } else {
    return -1;
  }
}

BEGIN {
  srand(PROCINFO["pid"]);
  if (npop=="") npop=1;
  if (parfile == "") parfile=ENVIRON["FFDIR"]"/param.rng";
  if (histfile == "") histfile="config_history";
  if (coupled_runs == "") coupled_runs="coupled_runs.lst";
  if (config2fit_sep=="") config2fit_sep="=";
  if (min_for_sex=="") min_for_sex=100;
  if (lfile == "") lfile="active_size"
  if (maxopt == "") maxopt=100; #Allow (roughly) at most 100 DISTINCT options in the bag
  if (maxdup == "") maxdup=10 # Allow at most 10 repeats of the most popular item in the bag
  stderr="/dev/stderr";

  ntypes=0;
  while (getline<parfile>0 && $1!="end") {
    ntypes ++;
    mode[ntypes]=$1;
    if ($1=="B") {
      prec[ntypes]="%s";
      mi=$2;
      ma=$3; #$2 is min, $3 is max
      no=$3-$2+1; #Nr of options in bias file
      bfile=$4;
      nno=0;
      level0="";
      if (bfile in used) {
        # This bias file already employed to create option list for type used[bfile]
        # Current option list may be a subset of that one, within the range ($2,$3)
        parent=used[bfile];
        for (k=1;k<=nopts[parent];k++) if (opt[parent,k]>=mi &&  opt[parent,k]<=ma) {nno++;opt[ntypes,nno]=opt[parent,k]}
        nopts[ntypes]=nno;
      } else {
        while (getline<bfile>0) if ($1>=mi && $1<=ma) {
          if (level0 == "") level0=maxdup/$2;
          level=level0*$2;
          #print $0,level > stderr;
          for (k=1;k<=int(level);k++) {nno++;opt[ntypes,nno]=$1}
          if (rand()<level-int(level)) {nno++;opt[ntypes,nno]=$1}
        }
        close(bfile);
        nopts[ntypes]=nno;
        used[bfile]=ntypes;
      }
      print ntypes,nopts[ntypes] >stderr;
    } else if ($1=="I") {
      # Integer range mode - min max {level enanced:level enhanced:level} (min=$2; max=$3)
      prec[ntypes]="%d";
      no=$3-$2+1;
      step[ntypes]=1;
      level=maxopt/no;
      nno=0;
      for (v=$2;v<=$3;v+=step[ntypes]) if (rand()<level) {
        nno ++;
        opt[ntypes,nno]=sprintf(prec[ntypes],v);
      }
      for (i=4;i<=NF;i++) {
        # Biased elements (parameter_value:multiplicity) follow starting from line 4
        split($i,w,":");
        ww=sprintf(prec[ntypes],w[1]);
        level=w[2]*maxopt/no;
        for (k=1;k<=int(level);k++) {nno++;opt[ntypes,nno]=ww}
        if (rand()<level-int(level)) {nno++;opt[ntypes,nno]=ww}
      }
      nopts[ntypes]=nno;
    } else if ($1=="R") {
      # Real range mode - min max prec {level enanced:level enhanced:level} (min=$2; max=$3)
      prec[ntypes]="%"$4"f";
      np=$4;
      sub("[0-9][.]","",np);
      no=int(($3-$2)*exp(np*log(10))+0.5)+1;
      step[ntypes]=sprintf(prec[ntypes],($3-$2)/no);
      level=maxopt/no;
      nno=0;
      for (v=$2;v<=$3;v+=step[ntypes]) if (rand()<level) {
        nno ++;
        opt[ntypes,nno]=sprintf(prec[ntypes],v);
      }
      for (i=5;i<=NF;i++) {
        # Biased elements (parameter_value:multiplicity) follow starting from line 5
        split($i,w,":");
        ww=sprintf(prec[ntypes],w[1]);
        level=w[2]*maxopt/no;
        for (k=1;k<=int(level);k++) {nno++;opt[ntypes,nno]=ww}
        if (rand()<level-int(level)) {nno++;opt[ntypes,nno]=ww}
      }
      nopts[ntypes]=nno;
    } else if ($1=="E") {
      #Enumeration mode - all on line are options
      no=NF-1;
      prec[ntypes]="%s";
      nno=0;
      for (i=2;i<=NF;i++) {
        level=maxopt/no;
        ww=$i;
        if (split($i,w,":") == 2) {
          level=w[2]*maxopt/no;
          ww=w[1];
        }
        for (k=1;k<=int(level);k++) {nno++;opt[ntypes,nno]=ww}
        if (rand()<level-int(level)) {nno++;opt[ntypes,nno]=ww}
      }
      nopts[ntypes]=nno;
    } else if ($1=="N") {
      #Open-ended numerical
      init[ntypes]=$2;
      step[ntypes]=$3;
      maxmv[ntypes]=$4;
      min[ntypes]=-99999;
      if ($5==1*$5) min[ntypes]=$5;
      max[ntypes]=+99999;
      if ($6==1*$6) max[ntypes]=$6;
      nn=split(step[ntypes],w,"[.]");
      if (nn==2) {
        prec[ntypes]="%"(1+length(w[1])+length(w[2]))"."length(w[2])"f";
      } else if (nn==1) {
        prec[ntypes]="%"length(step[ntypes])"d";
      }
      nopts[ntypes]=int((max[ntypes]-min[ntypes])/step[ntypes]);
    } else {
      print "FATAL - option "$NF" not supported in rng file" >"/dev/stderr";
      exit;
    }
  }

  nvars=0;
  while (getline<parfile>0 ) {
    nvars ++;
    type[nvars]=$2
  }
  close(parfile);

  #limited active subset of chromosome specified (must be beginning)
  getline<lfile;
  if ($1==int($1)) {
    actvars=int($1);
  } else {
    actvars=nvars;
  }
  close(lfile);

  nhist=0;
  while (getline<histfile>0) {
    nhist ++;
    fstr=prec[type[1]];
    l=$1;
    for (j=2;j<=actvars;j++) l+=":"$j;
    found[l]=nhist;
  }
  close(histfile);

  while (getline<coupled_runs>0) {
    histf=$1;
    while (getline<histf>0) {
      nhist ++;
      l=$1;
      for (j=2;j<=actvars;j++) l+=":"$j;
      found[l]=nhist;
    }
    close(histf);
  }
  close(coupled_runs);

}


{
  l=$1;
  for (j=2;j<=actvars;j++) l+=":"$j;
  found[l]=NR;
  for (i=1;i<=nvars;i++) x[NR,i]=$i;
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
      for (jj=1;jj<=nvars;jj++) {
        j=type[jj];
        if (mode[j]!="N") {
          val[jj]=opt[j,1+int(nopts[j]*rand())];
        } else {
          val[jj]="";
          while (val[jj]=="") {
            s=int(maxmv[j]*rand());
            if (rand()<0.5) s=-s;
            val[jj]=init[j]+s*step[j];
            if (val[jj]>max[j] || val[jj]<min[j]) val[jj]="";
          }
        }
      }

    } else if (rand()>0.40 && actvars>3) {

      print "Attempting Cross-Over.." > stderr;
      # It's a cross-over
      p1=1;
      p2=1;
      while (p1==p2) {
        p1=int(1+NR*rand()*rand());
        p2=int(1+NR*rand());
      }

      pos=2+int((actvars-2)*rand());

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
      pos=1+int(actvars*rand());
      while (nopts[type[pos]] == 1) pos=1+int(nvars*rand());

      for (j=1;j<=nvars;j++) {
        val[j]=x[p1,j];
        if (j == pos) while(val[j]==x[p1,j]) {
          if (mode[type[j]]!="N") {
            val[j]=opt[type[j],1+int(nopts[type[j]]*rand())];
          } else {
            s=1+int(maxmv[type[j]]*rand());
            if (rand()<0.5) s=-s;
            val[j]=x[p1,j]+s*step[type[j]];
            if (val[j]>max[type[j]] || val[j]<min[type[j]]) val[j]=x[p1,j];
          }
        }
      }
    }

    if (sorted=="y") {
      print "Sorting.." > stderr;
      asort(val,sval,"num_compare");
      for (j=1;j<=actvars;j++) val[j]=sval[j];
    }

    fstr=prec[type[1]];
    chromo=sprintf(fstr,val[1]);
    for (j=2;j<=actvars;j++) {
      fstr="%s:"prec[type[j]];
      chromo=sprintf(fstr,chromo,val[j]); 
    }

    gsub("  *"," ",chromo);
    if (!(chromo in found) || (ntrials>10000*npop)) {
      nhist ++;
      found[chromo]=nhist;
      new ++;
      fstr=prec[type[1]];
      chromo=sprintf(fstr,val[1]);;
      for (j=2;j<=actvars;j++) {
        fstr="%s "prec[type[j]];
        chromo=sprintf(fstr,chromo,val[j]);
      }
      print chromo;
      print chromo >> histfile;
      print new" new chromosomes were output" > stderr;
    } else {
      print chromo" discarded because already seen" > stderr;
    }
  } #END big while
}

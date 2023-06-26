BEGIN {
  if (npop=="") npop=1;
  if (parfile == "") parfile=ENVIRON["FFDIR"]"/param.rng";
  if (lfile == "") lfile="active_size"
  if (histfile == "") histfile="config_history";
  if (coupled_runs == "") coupled_runs="coupled_runs.lst";
  if (config2fit_sep=="") config2fit_sep="=";
  stderr="/dev/stderr";
  ignore_from=-1;
  ignore_if="shit";

  nvars=0;
  while (getline<parfile>0) {
    nvars ++;
    if ($2==1*$2) {
      absmin[nvars]=$2;
      min[nvars]=$3;
      max[nvars]=$4;
      absmax[nvars]=$5;
      prec[nvars]="%"$6"f";
      np=$6;
      sub("[0-9][.]","",np);
      nopts[nvars]=int((max[nvars]-min[nvars])*exp(np*log(10))+0.5)+1
    } else {
      nopts[nvars]=NF-1;
      for (i=2;i<=NF;i++) {opt[nvars,i-1]=$i};
      prec[nvars]="%s";
    }
    if ($(NF-1)=="ignore_remaining") {
      ignore_from=nvars;
      ignore_if=$NF;
      if (nopts[nvars]!="") nopts[nvars]=nopts[nvars]-2; # Zap two last strings from option list
    }
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

  new=0;
  ntrials=0;
  while (new<npop) {

    # Draw a random number and decide whether it will be a mutation or a cross_over

    ntrials ++;
    delete val;
    rand_threshold=(NR-2)/30;
    if (rand_threshold>0.8) rand_threshold=0.8;
    if (rand()>rand_threshold) {
      # Thre aint enough heredity to play GA yet, do random
      # Get a fully random chromosome
      print "Fully Random chromosome..." > stderr;
      for (j=1;j<=nvars;j++) {
        if (j in absmax) {
          if (rand()<0.8) {
            val[j]=min[j]+(max[j]-min[j])*rand();
          } else  {
            val[j]=absmin[j]+(absmax[j]-absmin[j])*rand();
          }
        } else {
          jj=1+int(nopts[j]*rand());
          val[j]=opt[j,jj];
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

      for (j=1;j<=nvars;j++) {
        if (j != pos) {
          val[j]=x[p1,j];
        } else {
          if (j in absmax) {
            # Got a numeric variable
            # Redraw a fully random value within range
            if (rand()<0.8) {
              val[j]=min[j]+(max[j]-min[j])*rand();
            } else  {
              val[j]=absmin[j]+(absmax[j]-absmin[j])*rand();
            }
          } else {
            jj=1+int(nopts[j]*rand());
            val[j]=opt[j,jj];
          }
        }
      }
    }


    if (ignore_from==1*ignore_from && val[ignore_from]==ignore_if) {
      # Reset all the values beyond position ignore_from, since the setup
      # of variable ignore_from makes coming parameters irrelevant
      for (j=ignore_from+1;j<=nvars;j++) {
        if (j in absmax) {
          val[j]=min[j];
        } else {
          val[j]=opt[j,1];
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

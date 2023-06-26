BEGIN {
  if (minact==0) minact=60.;
  if (maxact==0) maxact=90.;
  FS="@";
}
{
  nnf=split($0,a," ");
  if ((NR!=1) && (nnf!=nfold)) {
    print "FATAL: Line",NR,"has only",nnf,"fields instead of",nfold;
    exit;
  }
  for (i=1;i<=nnf;i++) {
    if ((mode=="inhib") && (i>1) && (NR>1)) {
      x[i,NR]=100.-a[i];
      x[i,NR]=(x[i,NR]-minact)/(maxact-minact);
      if (x[i,NR]>1.) {
	x[i,NR]=1;
      } else if (x[i,NR]<0.) {
	x[i,NR]=0;
      }
    } else {
      x[i,NR]=a[i];
    }
  }
  nfold=nnf;
}

END {
  for (i=1;i<=nfold;i++) {
    for (j=1;j<NR;j++) {
      printf("%s\t",x[i,j]);
    }
    printf("%s\n",x[i,NR]);
  }
}

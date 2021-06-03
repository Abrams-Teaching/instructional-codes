/* 1-D Brownian dynamics with replica-exchange
   
   D R E X E L   U N I V E R S I T Y
   Department of Chemical and Biological Engineering
   CHET580 Modern Molecular Simulations
   Cameron F Abrams cfa22@drexel.edu

*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_histogram.h>
#include <string.h>

/* quartic double-well */
void f ( double x, double a, double b, double c, double * e, double * f ) {
  double x2=x*x;
  double x3=x2*x;
  double x4=x2*x2;
  double bb=b-2*a;
  *e = a*x4+bb*x2 + c*x + a;
  *f = -(4*a*x3+2*bb*x + c);
}

void write_f(char * fn, double a, double b, double c, double xmin, double xmax, int n) {
  FILE * fp=fopen(fn,"w");
  double e, fc;
  double x;
  double dx=(xmax-xmin)/n;
  for (x=xmin+0.5*dx;x<=xmax;x+=dx) {
    f(x,a,b,c,&e,&fc);
    fprintf(fp,"% .5lf % .5lf % .5lf\n",x,e,fc);
  }
  fclose(fp);
  printf("Created %s.\n",fn);
}


void my_swap(int i, int j, double * x, int * c) {
  double tmp;
  int oci=c[i];
  int ocj=c[j];
  c[j]=oci;
  c[i]=ocj;
  tmp=x[i];
  x[i]=x[j];
  x[j]=tmp;
}

int main  ( int argc, char * argv [] ) {
  double * x, dx;
  int nrep=1;
  double Tmin=1, Tmax=5;
  int * nswapattempts, * nswapaccepts;
  double * T, * h2, * e, thise, k;
  int Tswap_every=0;
  int i, j, oj;
  int nsteps = 80000000;
  double gam = 5.0, h1;
  double force;
  double a = 0.02;
  double b = -1.0;
  double c = 0.0;
  char * plot_f=NULL;
  int * cfgi_atT;

  int log_every=0;
  FILE ** log_fp;
  char log_fn[255];
  char * log_fnf="rep%d.log";

  FILE * fp;

  gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937);
  unsigned long int Seed = 23410981;

  double dt = 0.0001;  /* Brownian dynamics time-step */

  int hist_n=100;
  int noutside=0;
  double x_min=-10,x_max=10;
  gsl_histogram ** h;

  char pfn[255];

    /* Here we parse the command line arguments */
  for (i=1;i<argc;i++) {
    if (!strcmp(argv[i],"-nrep")) nrep=atoi(argv[++i]);
    else if (!strcmp(argv[i],"-Tswap-every")) Tswap_every=atoi(argv[++i]);
    else if (!strcmp(argv[i],"-Tmin")) Tmin=atof(argv[++i]);
    else if (!strcmp(argv[i],"-Tmax")) Tmax=atof(argv[++i]);
    else if (!strcmp(argv[i],"-abc")) sscanf(argv[++i],"%lf,%lf,%lf",&a,&b,&c);
    else if (!strcmp(argv[i],"-gamma")) gam=atof(argv[++i]);
    else if (!strcmp(argv[i],"-ns")) nsteps = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-seed")) Seed = (unsigned long int)atoi(argv[++i]);
    else if (!strcmp(argv[i],"-hist-n")) hist_n = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-x-min")) x_min=atof(argv[++i]);
    else if (!strcmp(argv[i],"-x-max")) x_max=atof(argv[++i]);
    else if (!strcmp(argv[i],"-log-every")) log_every=atoi(argv[++i]);
    else if (!strcmp(argv[i],"-dt")) dt=atof(argv[++i]);
    else if (!strcmp(argv[i],"-log-fnf")) log_fnf=argv[++i];
    else if (!strcmp(argv[i],"-plot-f")) plot_f=argv[++i];
    else {
      fprintf(stderr,"Error.  Argument '%s' is not recognized.\n",argv[i]);
      exit(-1);
    }
  }
  
  /* allocate replicas */
  x=(double*)malloc(nrep*sizeof(double));
  T=(double*)malloc(nrep*sizeof(double));
  h2=(double*)malloc(nrep*sizeof(double));
  e=(double*)malloc(nrep*sizeof(double));
  nswapattempts=(int*)malloc(nrep*sizeof(int));
  nswapaccepts=(int*)malloc(nrep*sizeof(int));
  cfgi_atT=(int*)malloc(nrep*sizeof(int));
  log_fp=(FILE**)malloc(nrep*sizeof(FILE*));
  if (nrep>1) {
    k=log(Tmax/Tmin)/(nrep-1);
  } else {
    k=0.0;
  }
  h=(gsl_histogram**)malloc(nrep*sizeof(gsl_histogram*));
  for (i=0;i<nrep;i++) {
    x[i]=0.0;
    T[i]=Tmin*exp(k*i);
    cfgi_atT[i]=i;
    e[i]=0.0;
    h[i]=gsl_histogram_alloc(hist_n);
    nswapattempts[j]=0;
    nswapaccepts[j]=0;
    if (gsl_histogram_set_ranges_uniform(h[i],x_min,x_max)) {
      printf("Histogram error.\n");
      exit(-1);
    }
  }

  /* write out the potential so we can plot it later */
  if (plot_f) write_f(plot_f,a,b,c,x_min,x_max,hist_n);

  /* Seed the random number generator */
  gsl_rng_set(r,Seed);

  /* begin integration */
  h1=dt/gam;
  for (i=0;i<nrep;i++) {
    h2[i]=sqrt(6*T[i]*h1);
  }
  noutside=0;
  if (log_fnf && log_every>0) {
    for (i=0;i<nrep;i++) {
      sprintf(log_fn,log_fnf,i);
      log_fp[i]=fopen(log_fn,"w");
      fprintf(log_fp[i],"# Replica %i T = %.5lf\n",i,T[i]);
      fprintf(log_fp[i],"%d % .5f\n",0,x);
      fflush(log_fp[i]);
    }
  }
  for (i=1;i<=nsteps;i++) {
    for (j=0;j<nrep;j++) {
      f(x[j],a,b,c,&thise,&force);
      e[j]=thise;
      dx=h1*force+h2[j]*2*(0.5-gsl_rng_uniform(r));
      x[j]+=dx;
      if (gsl_histogram_increment(h[j],x[j])==GSL_EDOM) noutside++;
    }
    if (Tswap_every > 0 && i%Tswap_every==0) {
      /* do a swap attempt between two neighboring replicas */
      j=(int)gsl_rng_uniform_int(r,nrep-1);
      nswapattempts[j]++;
      oj=j+1;
      if (gsl_rng_uniform(r) < exp((1.0/T[oj]-1.0/T[j])*(e[j]-e[oj]))) {
        /* accept */
        my_swap(j,oj,x,cfgi_atT);
        nswapaccepts[j]++;
      }
    }
    if (log_every>0 && i%log_every==0) {
      for (j=0;j<nrep;j++) {
        //fprintf(stdout,"%d % .5lf %d %d % .5lf\n",j,T[j],nswapattempts[j],nswapaccepts[j],
        //      nswapattempts[j]>0?((double)nswapaccepts[j])/nswapattempts[j]:0.0);
        //fflush(stdout);
        fprintf(log_fp[j],"%d % .5lf %.3le\n",i,x[cfgi_atT[j]],nswapattempts[j]>0?((double)nswapaccepts[j])/nswapattempts[j]:0.0);
        fflush(log_fp[j]);
      }
    }
  }
  for (i=0;i<nrep;i++) {
    fclose(log_fp[i]);
    sprintf(pfn,"p%d.dat",i);
    fp=fopen(pfn,"w");
    fprintf(fp,"# T = %.5lf\n",T[i]);
    fflush(fp);
    fclose(fp);
    fp=fopen(pfn,"a");
    gsl_histogram_fprintf(fp,h[i],"%.5lf","%.8le");
    fclose(fp);
    fprintf(stdout,"Created %s\n",pfn);
  }
}

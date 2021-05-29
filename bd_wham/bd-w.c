/* 1-D Brownian dynamics on a two-well potential with one harmonic restraint
   
   D R E X E L   U N I V E R S I T Y
   Department of Chemical and Biological Engineering
   Cameron F. Abrams cfa22@drexel.edu

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
  *e = a*x4 + bb*x2 + c*x + a;
  *f = -(4*a*x3+2*bb*x + c);
}

/* harmonic window potential centered on x0 */
void w ( double x, double x0, double k, double * e, double * f ) {
    *e = 0.5*k*(x-x0)*(x-x0);
    *f = -k*(x-x0);
}

void write_f(char * fn, double a, double b, double c, double xmin, double xmax, int n) {
  FILE * fp=fopen(fn,"w");
  double e, fc;
  double x;
  double dx=(xmax-xmin)/n;
  fprintf(fp,"#LABELS x V(x) -dV/dx\n");
  for (x=xmin+0.5*dx;x<=xmax;x+=dx) {
    f(x,a,b,c,&e,&fc);
    fprintf(fp,"% .5lf % .5lf % .5lf\n",x,e,fc);
  }
  fclose(fp);
  printf("Created %s.\n",fn);
}

void write_w(char * fn, double xwin, double kwin, double xmin, double xmax, int n) {
  FILE * fp=fopen(fn,"w");
  double e, fc;
  double x;
  double dx=(xmax-xmin)/n;
  fprintf(fp,"#LABELS x W(x) -dW/dx\n");
  for (x=xmin+0.5*dx;x<=xmax;x+=dx) {
    w(x,xwin,kwin,&e,&fc);
    fprintf(fp,"% .5lf % .5lf % .5lf\n",x,e,fc);
  }
  fclose(fp);
  printf("Created %s.\n",fn);
}

int main  ( int argc, char * argv [] ) {
  double x, dx;
  int i;
  int nsteps = 80000000;
  double gam = 5.0, h1, h2;
  double T = 1.0;
  double force, e;
  double a = 0.02;
  double b = -1.0;
  double c = 0.0;
  char * plot_f=NULL;
  char * plot_w=NULL;
  double x_win = 0.0;
  double k_win = 1.0;
  double e_win, f_win;
  int which_win=0;
  char * log_f=NULL;
  int log_every=0;

  gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937);
  unsigned long int Seed = 23410981;

  double dt = 0.0001;  /* Brownian dynamics time-step */

  int hist_n=100;
  int noutside=0;
  double x_min=-10,x_max=10;
  gsl_histogram * h;

  FILE * fp;
  char pfn[255];

    /* Here we parse the command line arguments */
  for (i=1;i<argc;i++) {
    if (!strcmp(argv[i],"-T")) T=atof(argv[++i]);
    else if (!strcmp(argv[i],"-abc")) sscanf(argv[++i],"%lf,%lf,%lf",&a,&b,&c);
    else if (!strcmp(argv[i],"-gamma")) gam=atof(argv[++i]);
    else if (!strcmp(argv[i],"-ns")) nsteps = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-seed")) Seed = (unsigned long int)atoi(argv[++i]);
    else if (!strcmp(argv[i],"-hist-n")) hist_n = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-x-min")) x_min=atof(argv[++i]);
    else if (!strcmp(argv[i],"-x-max")) x_max=atof(argv[++i]);
    else if (!strcmp(argv[i],"-x-win")) x_win=atof(argv[++i]);
    else if (!strcmp(argv[i],"-k-win")) k_win=atof(argv[++i]);
    else if (!strcmp(argv[i],"-which-win")) which_win=atoi(argv[++i]);
    else if (!strcmp(argv[i],"-log-f")) log_f=argv[++i];
    else if (!strcmp(argv[i],"-log-every")) log_every=atoi(argv[++i]);
    else if (!strcmp(argv[i],"-dt")) dt=atof(argv[++i]);
    else if (!strcmp(argv[i],"-plot-f")) plot_f=argv[++i];
    else if (!strcmp(argv[i],"-plot-w")) plot_w=argv[++i];
    else {
      fprintf(stderr,"Error.  Argument '%s' is not recognized.\n",argv[i]);
      exit(-1);
    }
  }

  if (plot_f) write_f(plot_f,a,b,c,x_min,x_max,hist_n);
  if (plot_w) write_w(plot_w,x_win,k_win,x_min,x_max,hist_n);

  /* Allocate the histogram for Delta-U */
  h = gsl_histogram_alloc(hist_n);
  if (gsl_histogram_set_ranges_uniform(h,x_min,x_max)) {
    printf("Histogram error.\n");
    exit(-1);
  }

  /* Seed the random number generator */
  gsl_rng_set(r,Seed);

  /* begin integration */
  h1=dt/gam;
  h2=sqrt(6*T*h1);
  x=x_win;
  noutside=0;
  if (log_f && log_every>0) {
    fp=fopen(log_f,"w");
    fprintf(fp,"%d % .5f\n",0,x);
  }
  for (i=1;i<=nsteps;i++) {
    f(x,a,b,c,&e,&force);
    w(x,x_win,k_win,&e_win,&f_win);
    force+=f_win;
    e+=e_win;
    /* Brownian dynamics */
    dx=h1*force+h2*2*(0.5-gsl_rng_uniform(r));
    x+=dx;
    if (log_every>0 && nsteps%log_every==0) fprintf(fp,"%d % .5f\n",i,x);
    if (gsl_histogram_increment(h,x)==GSL_EDOM) noutside++;
  }
  sprintf(pfn,"p%d.dat",which_win);
  fp=fopen(pfn,"w");
  fprintf(fp,"# gsl-histogram window %d center %.5lf k %.5lf\n",which_win,x_win,k_win);
  fclose(fp);
  fp=fopen(pfn,"a");
  gsl_histogram_fprintf(fp,h,"%.5lf","%.8le");
  fclose(fp);
  fprintf(stdout,"Created %s\n",pfn);
}

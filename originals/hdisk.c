/* 
   Metropolis Monte Carlo simulation of hard disks confined in a circle

   Cameron F. Abrams

   Written for the course CHE T580 Modern Molecular Simulation
   Spring 20-21

   compile using "gcc -o hdisk hdisk.c -lm -lgsl"
   (assumes the GNU Scientific Library is installed)

   runs as "./hdisk -N <number_of_particles> -rho <density> \
                    -R <radius_of_circle> \
                    -nc <numcycles(1e6)>  \
		    -s <seed(?)>"

   You must have the GNU Scientific Library installed; see
   the coursenotes to learn how to do this.

   Drexel University, Department of Chemical Engineering
   Philadelphia
   (c) 2004-2021
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_rng.h>

void out ( FILE * fp, double * rx, double * ry, int n ) {
  int i;
  for (i=0;i<n;i++) {
    fprintf(fp,"%.5lf %.5lf\n",rx[i],ry[i]);
  }
}

/* Initialize particle positions by assigning them
   randomly while avoiding overlaps. */
void init ( double * rx, double * ry, 
	    int n, double R2, double s2, gsl_rng * r ) {
  int i,j,reject;
  double R = sqrt(R2), r2, sx, sy;
  /* nMax is maximum number of insertion trials */
  int nMax = 10000, nTrials=0;
  
  FILE * tmp_fp;

  for (i=0;i<n;i++) {
    reject = 1;
    nTrials = 0;
    while (reject&&nTrials<nMax) {
      reject = 0;
      /* randomly locate particle i in a circle of radius R */
      rx[i] = R*(1.0-2*gsl_rng_uniform(r));
      ry[i] = sqrt(R2-rx[i]*rx[i])*(1.0-2*gsl_rng_uniform(r));
      /* check for overlaps; reject if an overlap is found */
      for (j=0;j<n;j++) {
	if (j!=i) {
	  sx  = rx[i]-rx[j];
	  sy  = ry[i]-ry[j];
	  r2  = sx*sx + sy*sy;
	  if (r2 < s2) reject=1;
	}
      }
      nTrials++;
    }
    if (nTrials==nMax) {
      fprintf(stderr,"# Error -- could not initialize position "
	      "of particle %i in %i trials\n",i,nTrials);
      fprintf(stderr,"# program ends.\n");
      exit(-1);
    }
  }
  tmp_fp=fopen("hdtmpcfg_init","w");
  out(tmp_fp,rx,ry,n);
  fclose(tmp_fp);
}

int main ( int argc, char * argv[] ) {

  double * rx, * ry;
  int N=-1,c,a;
  double R2=16.0, s2=1.0;
  double rho=0.5;
  double dr=0.1,dx,dy,theta,sx,sy,r2;
  int i,j;
  int nCycles = 10;
  int nAcc, reject;
  int noob=0,novl=0;
  int short_out = 0;
  FILE * tmp_fp;

  gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937);
  unsigned long int Seed = 23410981;

  /* Here we parse the command line arguments */
  for (i=1;i<argc;i++) {
    if (!strcmp(argv[i],"-N")) N=atoi(argv[++i]);
    else if (!strcmp(argv[i],"-R")) R2=atof(argv[++i]);
    else if (!strcmp(argv[i],"-rho")) rho=atof(argv[++i]);
    else if (!strcmp(argv[i],"-dr")) dr=atof(argv[++i]);
    else if (!strcmp(argv[i],"-s")) s2=atof(argv[++i]);
    else if (!strcmp(argv[i],"-nc")) nCycles = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-so")) short_out = 1;
    else if (!strcmp(argv[i],"-s")) Seed = (unsigned long)atoi(argv[++i]);
  }

  /* precompute square of particle diameter */
  s2 = s2*s2;

  /* If N was not set by the user, calculate it based on 
     the given values of density and radius */
  if (N==-1) {
    N = (int)(rho*M_PI*R2);
  }
  /* Otherwise, recompute the radius. */
  else {
    R2 = ((double)N)/rho/M_PI;
  }

  if (!short_out) 
    fprintf(stdout,"# R = %.5lf; rho = %.5lf; N = %i; seed = %ul\n",
	                  sqrt(R2),rho,N,Seed);

  /* Seed the random number generator */
  gsl_rng_set(r,Seed);

  /* allocate the position arrays */
  rx = (double*)malloc(N*sizeof(double));
  ry = (double*)malloc(N*sizeof(double));

  /* generate initial positions that fit inside 
     circle with radius R and guarantee no 
     particles overlap. */
  init(rx,ry,N,R2,s2,r);

  nAcc = 0;
  for (c=0;c<nCycles;c++) {
    /* Make N displacement attempts */
    for (a=0;a<N;a++) {
      /* randomly select a particle */
      i=(int)gsl_rng_uniform_int(r,N);
      /* calculate displacement */
      dx = dr*(0.5-gsl_rng_uniform(r));
      dy = dr*(0.5-gsl_rng_uniform(r));
      /* displace particle */
      rx[i]+=dx;
      ry[i]+=dy;

      /* compute new distance to origin */
      r2=rx[i]*rx[i]+ry[i]*ry[i];
      /* reject move if outside circle */
      reject=(r2>R2);
      if (reject) noob++;

      /* check for overlaps with other particles */
      if (!reject) {
	      for (j=0;j<N;j++) {
	        if (j!=i) {
	          sx  = rx[i]-rx[j];
            sy  = ry[i]-ry[j];
            r2  = sx*sx + sy*sy;
            if (r2 < s2) reject=1;
          }
	      }
	      if (reject) novl++;
      }

      /* if move is rejected, undo displacement */
      if (reject) {
        rx[i]-=dx;
        ry[i]-=dy;
      } else nAcc++;
    }
  }
  tmp_fp=fopen("hdtmpcfg","w");
  out(tmp_fp,rx,ry,N);
  fclose(tmp_fp);
  if (short_out) {
    fprintf(stdout,"%.6lf %.6lf %.6lf %.6lf\n",
	    dr,((double)nAcc)/(N*nCycles),
	    ((double)noob)/(N*nCycles-nAcc),
	    ((double)novl)/(N*nCycles-nAcc)
	    );
  }
  else {
    fprintf(stdout,"Results:\n"
	    "Number of Trial Moves:          %i\n"
	    "Maximum Displacement Length:    %.5lf\n"
	    "Acceptance Ratio:               %.5lf\n"
	    "Reject Fraction Out-of-bounds:  %.5lf\n"
	    "Reject Fraction Overlap:        %.5lf\n",
	    N*nCycles,dr,((double)nAcc)/(N*nCycles),
	    ((double)noob)/(N*nCycles-nAcc),
	    ((double)novl)/(N*nCycles-nAcc)
	    );
  }
}

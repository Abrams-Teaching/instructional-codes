/* 
   Metropolis Monte Carlo simulation of a Lennard-Jones fluid
   with Widom test-particle insertion to compute mu_ex using
   the overlapping distribution method

   Cameron F. Abrams

   Written for the course CHE T580, Modern Molecular Simulation
   Spring 20-21

   compile using "gcc -o mclj_od mclj_od.c -lm -lgsl"

   runs as "./mclj_od -N <number_of_particles> -rho <density> \
                    -nc <numcycles(1e6)>  -dr <delta-r> \
		                -s <seed(?)> -ne <#equil.cycles(100)>"

   You must have the GNU Scientific Library installed.

   Drexel University, Department of Chemical Engineering
   Philadelphia
   (c) 2004-2021
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_histogram.h>

/* energy of particle i */
double e_i ( int i, double * rx, double * ry, double * rz, int N, double L,
		 double rc2, int tailcorr, double ecor, 
		 int shift, double ecut, double * vir, int i0 ) {
  int j;
  double dx, dy, dz, r2, r6i;
  double e = 0.0, hL=L/2.0;
  *vir=0.0;
  for (j=i0;j<N;j++) {
    if (i!=j) {
      dx  = rx[i]-rx[j];
      dy  = ry[i]-ry[j];
      dz  = rz[i]-rz[j];
      /* Periodic boundary conditions: Apply the minimum image
        convention; note that this is *not* used to truncate the
        potential as long as there an explicit cutoff. */
      if (dx>hL)       dx-=L;
      else if (dx<-hL) dx+=L;
      if (dy>hL)       dy-=L;
      else if (dy<-hL) dy+=L;
      if (dz>hL)       dz-=L;
      else if (dz<-hL) dz+=L;
      r2 = dx*dx + dy*dy + dz*dz;
      if (r2<rc2) {
        r6i   = 1.0/(r2*r2*r2);
        e    += 4*(r6i*r6i - r6i) - (shift?ecut:0.0);
        *vir += 48*(r6i*r6i-0.5*r6i);
      }
    }
  }
  return e+(tailcorr?ecor:0.0);
}
/* An N^2 algorithm for computing the total energy.  The virial
   is also computed and returned in *vir. */
double total_e ( double * rx, double * ry, double * rz, int N, double L,
		 double rc2, int tailcorr, double ecor, 
		 int shift, double ecut, double * vir ) {
  int i;
  double tvir;
  double e = 0.0;
  *vir=0.0;
  for (i=0;i<N-1;i++) {
    e    += e_i(i,rx,ry,rz,N,L,rc2,tailcorr,ecor,shift,ecut,&tvir,i+1);
    *vir += tvir;
  }
  return e;
}

/* Writes configuration in XYZ format */
void write_xyz(FILE * fp, double * rx, double * ry, double * rz, int n, double L) {
  int i;
  fprintf(fp,"%i\n",n);
  fprintf(fp,"BOX %.5lf %.5lf %.5lf\n",L,L,L);
  for (i=0;i<n;i++) {
    fprintf(fp,"%s %.5lf %.5lf %.5lf\n","Ar",rx[i],ry[i],rz[i]);
  }
}

/* Initialize particle positions by assigning them
   on a cubic grid, then scaling positions 
   to achieve a given box size and thereby, volume,
   and density */
void init ( double * rx, double * ry, double * rz,
	    int n, double L, gsl_rng * r ) {
  int i,ix,iy,iz;
  
  int n3=2;
  /* Find the lowest perfect cube, n3, greater than or equal to the
     number of particles */
  while ((n3*n3*n3)<n) n3++;

  ix=iy=iz=0;
  /* Assign particle positions */
  for (i=0;i<n;i++) {
    rx[i] = ((double)ix+0.5)*L/n3;
    ry[i] = ((double)iy+0.5)*L/n3;
    rz[i] = ((double)iz+0.5)*L/n3;
    ix++;
    if (ix==n3) {
      ix=0;
      iy++;
      if (iy==n3) {
        iy=0;
        iz++;
      }
    }
  }
}

void make_real ( double * rx, double * ry, double * rz, int N, double L, 
	     double rc2, int shift, double ecut,
	     gsl_rng * r, double * e ) {

  int j;
  double dx,dy,dz,hL=L/2.0,r2,r6i;

  (*e)=0.0;
  rx[N]=(gsl_rng_uniform(r)-0.5)*L;
  ry[N]=(gsl_rng_uniform(r)-0.5)*L;
  rz[N]=(gsl_rng_uniform(r)-0.5)*L;

  for (j=0;j<N;j++) {
    dx  = rx[N]-rx[j];
    dy  = ry[N]-ry[j];
    dz  = rz[N]-rz[j];
    if (dx>hL)       dx-=L;
    else if (dx<-hL) dx+=L;
    if (dy>hL)       dy-=L;
    else if (dy<-hL) dy+=L;
    if (dz>hL)       dz-=L;
    else if (dz<-hL) dz+=L;
    r2 = dx*dx + dy*dy + dz*dz;
    if (r2<rc2) {
      r6i   = 1.0/(r2*r2*r2);
      *e    += 4*(r6i*r6i - r6i) - (shift?ecut:0.0);
    }
  }
}

int make_pdf(gsl_histogram * h, int noutside) {
  double norm1=gsl_histogram_sum(h), bs, scale;
  int i;
  for (i=0;i<h->n;i++) {
    bs=h->range[i+1]-h->range[i];
    scale=1.0/(bs*(norm1+noutside));
    h->bin[i]*=scale;
  }
  return 0;
}

enum {XYZ, NONE};
int main ( int argc, char * argv[] ) {

  double * rx, * ry, * rz;
  int N=216,c,a,p;
  double L=0.0;
  double we;
  double beta;
  double rho=0.5, T=1.0, rc2 = 3.5, vir, vir_old, p_sum, pcor, V, E_old;
  double ei_new, ei_old, delta_e, ivir_old, ivir_new, delta_vir, esum, rr3, ecor, ecut;
  double dr=0.2,dx,dy,dz;
  double rxold,ryold,rzold;
  double w_sum, anti_w_sum;
  int i,j;
  int nCycles = 10, nSamp, nEq=1000;
  int nAcc;
  int short_out=0;
  int shift=0;
  int tailcorr=1;
  int prog=0;
  int hist_n=100;
  double du_min=-30,du_max=30;
  gsl_histogram * h;
  int which_sim=0;
  char pfn[255];
  double bf;
  int noutside=0;

  gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937);
  unsigned long int Seed = 23410981;

  FILE * fp;
  char * traj_fn=NULL;
  int traj_out=XYZ;
  int traj_samp=100;

  /* Here we parse the command line arguments */
  for (i=1;i<argc;i++) {
    if (!strcmp(argv[i],"-N")) N=atoi(argv[++i]);
    else if (!strcmp(argv[i],"-rho")) rho=atof(argv[++i]);
    else if (!strcmp(argv[i],"-T")) T=atof(argv[++i]);
    else if (!strcmp(argv[i],"-dr")) dr=atof(argv[++i]);
    else if (!strcmp(argv[i],"-rc")) rc2=atof(argv[++i]);
    else if (!strcmp(argv[i],"-nc")) nCycles = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-ne")) nEq = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-so")) short_out=1;
    else if (!strcmp(argv[i],"+tc")) tailcorr=0;
    else if (!strcmp(argv[i],"-sh")) shift=1;
    else if (!strcmp(argv[i],"-prog")) prog = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-which-sim")) which_sim = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-hist-n")) hist_n = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-du-min")) du_min = atof(argv[++i]);
    else if (!strcmp(argv[i],"-du-max")) du_max = atof(argv[++i]);    
    else if (!strcmp(argv[i],"-s")) 
      Seed = (unsigned long)atoi(argv[++i]);
    else if (!strcmp(argv[i],"-traj_samp")) traj_samp = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-traj"))  {
      traj_fn=argv[++i];
      /* detect format */
      if (strstr(traj_fn, ".xyz") != NULL) {
        traj_out=XYZ;
      } else {
        fprintf(stdout,"File format of %s not recognized; must by xyz.\n",traj_fn);
        traj_out=NONE;
      }
    }
    else {
      fprintf(stderr,"Error.  Argument '%s' is not recognized.\n",argv[i]);
      exit(-1);
    }
  }

  /* Allocate the histogram for Delta-U */
  h = gsl_histogram_alloc(hist_n);
  if (gsl_histogram_set_ranges_uniform(h,du_min,du_max)) {
    printf("Histogram error.\n");
    exit(-1);
  }

  /* Compute the side-length */
  L = pow((V=N/rho),0.3333333);

  /* Compute the tail-corrections; assumes sigma and epsilon are both 1 */
  rr3 = 1.0/(rc2*rc2*rc2);
  ecor = 8*M_PI*rho*(rr3*rr3*rr3/9.0-rr3/3.0);
  pcor = 16.0/3.0*M_PI*rho*rho*(2./3.*rr3*rr3*rr3-rr3);
  ecut = 4*(rr3*rr3*rr3*rr3-rr3*rr3);

  /* Compute the *squared* cutoff, reusing the variable rc2 */
  rc2*=rc2;

  /* For computational efficiency, use reciprocal T */
  beta = 1.0/T;

  /* compute box volume */
  V = L*L*L;
  
  /* Output some initial information */
  fprintf(stdout,"# NVT MC Simulation of a Lennard-Jones fluid\n");
  fprintf(stdout,"# Overlapping distribution method, simulation %d\n",which_sim);
  fprintf(stdout,"# L = %.5lf; rho = %.5lf; N = %i; rc = %.5lf\n",
	  L,rho,N,sqrt(rc2));
  fprintf(stdout,"# nCycles %i, nEq %i, seed %lu, dR %.5lf\n",
	  nCycles,nEq,Seed,dr);
  
  /* Total number of cycles is number of "equilibration" cycles plus
     number of "production" cycles */
  nCycles+=nEq;

  /* Seed the random number generator */
  gsl_rng_set(r,Seed);

  /* Allocate the position arrays */
  rx = (double*)malloc((N+1)*sizeof(double));
  ry = (double*)malloc((N+1)*sizeof(double));
  rz = (double*)malloc((N+1)*sizeof(double));

  /* Generate initial positions on a cubic grid, 
     and measure initial energy */
  init(rx,ry,rz,N,L,r);

  E_old = total_e(rx,ry,rz,N,L,rc2,tailcorr,ecor,shift,ecut,&vir_old);
  if (traj_fn&&traj_out==XYZ) {
    fp=fopen(traj_fn,"w");
    write_xyz(fp,rx,ry,rz,N,L);
    fclose(fp);
  }
  nAcc = 0;
  esum = 0.0;
  nSamp = 0;
  p_sum = 0.0;
  w_sum = 0.0;
  anti_w_sum = 0.0;
  if (prog>0) {
    printf("#LABEL cycle <e>/<n> p mu_ex\n");
  }
  for (c=0;c<nCycles;c++) {
    /* Randomly select a particle */
    i=(int)gsl_rng_uniform_int(r,N);
    /* calculate displacement */
    dx = dr*(0.5-gsl_rng_uniform(r));
    dy = dr*(0.5-gsl_rng_uniform(r));
    dz = dr*(0.5-gsl_rng_uniform(r));
    ei_old=e_i(i,rx,ry,rz,N,L,rc2,tailcorr,ecor,shift,ecut,&ivir_old,0);
    if (which_sim==1&&c>nEq) {
      if (gsl_histogram_increment(h,ei_old)==GSL_EDOM) noutside++;
      anti_w_sum+=exp(beta*ei_old);
    }
    /* Save the current position of particle i */
    rxold=rx[i];
    ryold=ry[i];
    rzold=rz[i];

    /* Displace particle i */
    rx[i]+=dx;
    ry[i]+=dy;
    rz[i]+=dz;

    /* Apply periodic boundary conditions */
    if (rx[i]<0.0) rx[i]+=L;
    if (rx[i]>L)   rx[i]-=L;
    if (ry[i]<0.0) ry[i]+=L;
    if (ry[i]>L)   ry[i]-=L;
    if (rz[i]<0.0) rz[i]+=L;
    if (rz[i]>L)   rz[i]-=L;
    ei_new=e_i(i,rx,ry,rz,N,L,rc2,tailcorr,ecor,shift,ecut,&ivir_new,0);
    /* Get the new energy */
    //E_new = total_e(rx,ry,rz,N,L,rc2,tailcorr,ecor,shift,ecut,&vir);
    delta_e=ei_new-ei_old;
    delta_vir=ivir_new-ivir_old;

    /* Conditionally accept... */
    if (gsl_rng_uniform(r) < exp(-beta*delta_e)) {
      E_old+=delta_e;
      vir_old+=delta_vir;
      nAcc++;
    }
    /* ... or reject the move; reassign the old positions */
    else {
      rx[i]=rxold;
      ry[i]=ryold;
      rz[i]=rzold;
    }
    /* Sample: default frequency is once per trial move; We must
      include results of a move regardless of whether the move is
      accepted or rejected. */
    if (c>nEq) {
      esum+=E_old;
      p_sum+=vir_old/3.0/V+pcor;
      if (which_sim==0) {
        make_real(rx,ry,rz,N,L,rc2,shift,ecut,r,&we);
        if (gsl_histogram_increment(h,we)==GSL_EDOM) noutside++;
        w_sum+=exp(-beta*we);
      }
      nSamp++;
      if (prog>0&&!(c%prog)) {
        printf("% 10i % .5f % .5f % .5f % .5f\n",c,esum/nSamp/N,p_sum/nSamp+rho/beta,-T*log(w_sum/nSamp),T*log(anti_w_sum/nSamp));
        fflush(stdout);
      }
    }
    if (traj_fn&&!(c%traj_samp)) {
      if (traj_out==XYZ) {
        fprintf(stdout,"# Trajectory snapshot at %i\n",c);fflush(stdout);
        fp=fopen(traj_fn,"a");
        write_xyz(fp,rx,ry,rz,N,L);
        fclose(fp);
      }
    }
  }

  make_pdf(h,noutside);
  sprintf(pfn,"p%d.dat",which_sim);
  fp=fopen(pfn,"w");
  gsl_histogram_fprintf(fp,h,"%.5lf","%.8le");
  fclose(fp);
  fprintf(stdout,"Created %s\n",pfn);

  if (short_out)
    fprintf(stdout,"% .6lf % .5lf % .5lf % .5lf % .5lf % .5lf\n",
	    dr,((double)nAcc)/(N*nCycles),
	    esum/nSamp/N,p_sum/nSamp+rho*T,-T*log(w_sum/nSamp));
  else
    fprintf(stdout,"NVT Metropolis Monte Carlo Simulation"
	    " of the Lennard-Jones fluid.\n"
	    "---------------------------------------------\n"
	    "Number of particles:              %i\n"
	    "Number of cycles:                 %i\n"
	    "Cutoff radius:                    %8.5lf\n"
	    "Maximum displacement:             %8.5lf\n"
	    "Density:                          %8.5lf\n"
	    "Temperature:                      %8.5lf\n"
	    "Tail corrections applied?         %s\n"
	    "Shifted potential?                %s\n"
	    "Results:\n"
	    "Potential energy tail correction: %8.5lf\n"
	    "Pressure tail correction:         %8.5lf\n"
	    "Potential energy shift at cutoff  %8.5lf\n"
	    "Acceptance ratio:                 %8.5lf\n"
	    "Energy/particle:                  %8.5lf\n"
	    "Ideal gas pressure:               %8.5lf\n"
	    "Virial:                           %8.5lf\n"
	    "Total pressure:                   %8.5lf\n"
      "Widom estimate mu_ex:             %8.5lf\n"
      "Anti-Widom estimate mu_ex:        %8.5lf\n"
	    "Program ends.\n",
	    N,nCycles,sqrt(rc2),dr,rho,T,
	    tailcorr?"Yes":"No",shift?"Yes":"No",
	    ecor,pcor,ecut,
	    ((double)nAcc)/nCycles,
	    esum/nSamp/N,
	    rho/beta,p_sum/nSamp,
	    p_sum/nSamp+rho/beta,-T*log(w_sum/nSamp),T*log(anti_w_sum/nSamp));
  if (traj_fn) {
    fprintf(stdout,"Trajectory written to %s.\n",traj_fn);
  }
}

/* 
   Metropolis Monte Carlo simulation of a Lennard-Jones fluid
   in a periodic boundary and the grand canonical ensemble

   Cameron F. Abrams

   Written for the course CHE T580, Modern Molecular Simulation
   Spring 20-21

   compile using "gcc -o mclj_muvt mclj_muvt.c -lm -lgsl"

   runs as "./mclj_muvt -N <number_of_particles> -rho <density> \
                        -mu <mu-prime> -nc <numcycles(1e6)> \
                        -dr <delta-r> -s <seed(?)> \
                        -ne <#equil.cycles(100)>"

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

int mcmove ( double * rx, double * ry, double * rz, int i, int N,
	     int tailcorr, double ecor, int shift, double ecut, 
	     double L, double rc2, double dr, 
	     double T, double * E_old, double * vir_old,
	     gsl_rng * r ) {
  double dx, dy, dz;
  double rxold, ryold, rzold, E_new, vir;

  /* calculate displacement */
  dx = dr*(0.5-gsl_rng_uniform(r));
  dy = dr*(0.5-gsl_rng_uniform(r));
  dz = dr*(0.5-gsl_rng_uniform(r));

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

  E_new = total_e(rx,ry,rz,N,L,rc2,tailcorr,ecor,shift,ecut,&vir);

  /* Conditionally accept... */
  if (gsl_rng_uniform(r) < exp(-T*(E_new-(*E_old)))) {
    *E_old=E_new;
    *vir_old=vir;
    return 1;
  }
  /* ... or reject the move; reassign the old positions */
  else {
    rx[i]=rxold;
    ry[i]=ryold;
    rz[i]=rzold;
    return 0;
  }
}

enum {XYZ, NONE};
int main ( int argc, char * argv[] ) {

  double * rx, * ry, * rz;
  int N=216,c,a,p;
  int Nres=216,N0,Ntot;
  double L=0.0, E, arg, dumvir, beta;
  double rho=0.8, T=2.0, rc2 = 3.5, vir, vir_old, vir_sum, pcor, V;
  double rho_sum, p_sum, this_e;
  double mu = 0.0;
  double E_new, E_old, esum, rr3, ecor, ecut;
  double dr=0.2;
  double disp_wt=0.5, x;
  int i;
  int nCycles = 10, nSamp, nEq=1000;
  int nTransAtt, nTransAcc, nInsAtt, nInsAcc, nDelAtt, nDelAcc, acc, nAcc;
  int short_out=0;
  int shift=0;
  int tailcorr=1;
  int prog=0;
  gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937);
  unsigned long int Seed = 23410981;

  FILE * fp;
  char * traj_fn=NULL;
  int traj_out=XYZ;
  int traj_samp=100;

  /* Here we parse the command line arguments */
  for (i=1;i<argc;i++) {
    if (!strcmp(argv[i],"-N")) N=atoi(argv[++i]);
    else if (!strcmp(argv[i],"-Nres")) N=atoi(argv[++i]);
    else if (!strcmp(argv[i],"-rho")) rho=atof(argv[++i]);
    else if (!strcmp(argv[i],"-mu")) mu=atof(argv[++i]);
    else if (!strcmp(argv[i],"-disp-wt")) disp_wt=atof(argv[++i]);
    else if (!strcmp(argv[i],"-T")) T=atof(argv[++i]);
    else if (!strcmp(argv[i],"-dr")) dr=atof(argv[++i]);
    else if (!strcmp(argv[i],"-rc")) rc2=atof(argv[++i]);
    else if (!strcmp(argv[i],"-nc")) nCycles = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-ne")) nEq = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-prog")) prog = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-so")) short_out=1;
    else if (!strcmp(argv[i],"+tc")) tailcorr=0;
    else if (!strcmp(argv[i],"-sh")) shift=1;
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
  fprintf(stdout,"# muVT MC Simulation of a Lennard-Jones fluid\n");
  fprintf(stdout,"# L = %.5lf; rho0 = %.5lf; mu' = %.5lf; N0 = %i; rc = %.5lf\n",
	  L,rho,mu,N,sqrt(rc2));
  fprintf(stdout,"# nCycles %i, nEq %i, seed %lu, dR %.5lf\n",
	  nCycles,nEq,Seed,dr);
  
  /* Total number of cycles is number of "equilibration" cycles plus
     number of "production" cycles */
  nCycles+=nEq;

  /* Seed the random number generator */
  gsl_rng_set(r,Seed);

  /* Allocate the position arrays */
  Ntot = N + Nres;
  N0 = N;
  rx = (double*)malloc(Ntot*sizeof(double));
  ry = (double*)malloc(Ntot*sizeof(double));
  rz = (double*)malloc(Ntot*sizeof(double));

  /* Generate initial positions on a cubic grid, 
     and measure initial energy */
  init(rx,ry,rz,N,L,r);

  E_old = total_e(rx,ry,rz,N,L,rc2,tailcorr,ecor,shift,ecut,&vir_old);
  if (traj_fn&&traj_out==XYZ) {
    fp=fopen(traj_fn,"w");
    write_xyz(fp,rx,ry,rz,N,L);
    fclose(fp);
  }
  nTransAcc = 0;
  nInsAcc = 0;
  nDelAcc = 0;
  nTransAtt = 0;
  nInsAtt = 0;
  nDelAtt = 0;
  E = esum = E_old;
  nSamp = 1;
  vir = vir_old;
  p_sum = vir/3.0/V+pcor;
  rho_sum = N/V;
  if (prog>0) {
    printf("#LABEL cycle <e>/<n> p mu_ex\n");
  }
  for (c=0;c<nCycles;c++) {
    /* Get a random number between 0 and 1 */
    x=gsl_rng_uniform(r);
    if (x<disp_wt||c<nEq) { /* perform a displacement move */
      nTransAtt++;
      i=(int)gsl_rng_uniform_int(r,N);
      acc = mcmove(rx,ry,rz,i,N,tailcorr,ecor,shift,ecut,
		               L,rc2,dr,beta,&E,&vir,r);
      nTransAcc+=acc;
    }
    else { /* perform an insertion/deletion move */
      /* Get a random number between 0 and 1 */
      x=gsl_rng_uniform(r);
      if (x<0.5) { /* try deletion */
        nDelAtt++;
        i=(int)gsl_rng_uniform_int(r,N);
        this_e = e_i(i,rx,ry,rz,N,L,rc2,tailcorr,ecor,shift,ecut,&dumvir,0);
        arg = N/V*exp(-beta*(-this_e+mu));
        //printf("del e %.5le arg %.5le\n",this_e,arg);fflush(stdout);
        if (gsl_rng_uniform(r)<arg) {
            rx[i]=rx[N-1];
            ry[i]=ry[N-1];
            rz[i]=rz[N-1];
            N--;
            if (N==0) break;
            E-=this_e;
            vir-=dumvir;
            nDelAcc++;
        }
      } else { /* try insertion */
        nInsAtt++;
        rx[N]=gsl_rng_uniform(r)*L;
        ry[N]=gsl_rng_uniform(r)*L;
        rz[N]=gsl_rng_uniform(r)*L;
        this_e = e_i(N,rx,ry,rz,N+1,L,rc2,tailcorr,ecor,shift,ecut,&dumvir,0);
        arg = V/(N+1)*exp(-beta*(this_e-mu));
        if (gsl_rng_uniform(r)<arg) {
            N++;
            E+=this_e;
            vir+=dumvir;
            nInsAcc++;
        }
      }
    }

    /* Sample: default frequency is once per trial move; We must
       include results of a move regardless of whether the move is
       accepted or rejected. */
    if (c>nEq) {
      esum+=E;
      p_sum+=vir/3.0/V+pcor;
      rho_sum+=N/V;
      nSamp++;
      if (prog>0&&!(c%prog)) {
        printf("% 10i % .5f % .5f % .5f\n",c,esum/nSamp/N,p_sum/nSamp+rho_sum/nSamp*T,mu-T*log(rho_sum/nSamp));
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
  nAcc=nTransAcc+nInsAcc+nDelAcc;

  if (short_out)
    fprintf(stdout,"%.3lf %.3lf | %.6lf %.6lf : %.5lf | %.5lf %.5lf %.5lf %.5lf\n",
	    T,mu,dr,disp_wt,((double)nAcc)/nCycles,
	    esum/nSamp/N,p_sum/nSamp+rho_sum/nSamp*T,rho_sum/nSamp,mu-T*log(rho_sum/nSamp));
  else
    fprintf(stdout,"NPT Metropolis Monte Carlo Simulation"
	    " of the Lennard-Jones fluid in the Grand Canonical"
      " Ensemble\n"
	    "---------------------------------------------\n"
	    "Number of cycles:                 %i\n"
	    "Maximum particle displacement:    % .5lf\n"
	    "Displacement weight:              % .5lf\n"
	    "Temperature:                      % .5lf\n"
	    "Relative chemical potential:      % .5lf\n"
      "Initial number of particles:      %i\n"
	    "Tail corrections used?            %s\n"
	    "Shifted potentials used?          %s\n"
	    "Results:\n"
      "Final number of particles:        %i\n"
	    "Displacement attempts:            %i\n"
	    "Insertion attempts:               %i\n"
      "Deletion attempts:                %i\n"
	    "Acceptance ratio, ptcl displ:     % .5lf\n"
	    "Acceptance ratio, insertion:      % .5lf\n"
	    "Acceptance ratio, deletion:       % .5lf\n"
	    "Overall acceptance ratio:         % .5lf\n"
	    "Energy/particle:                  % .5lf\n"
	    "Density:                          % .5lf\n"
	    "Computed pressure:                % .5lf\n"
      "Excess chemical potential:        % .5lf\n"
	    "Program ends.\n",
	    nCycles,dr,disp_wt,T,mu,N0,
	    tailcorr?"Yes":"No",shift?"Yes":"No",N,
	    nTransAtt,nInsAtt,nDelAtt,
	    ((double)nTransAcc)/(nTransAtt?nTransAtt:1),
	    ((double)nInsAcc)/(nInsAtt?nInsAtt:1),
	    ((double)nDelAcc)/(nDelAtt?nDelAtt:1),
	    ((double)nAcc)/nCycles,
	    esum/nSamp/N,rho_sum/nSamp,
	    p_sum/nSamp+rho_sum/nSamp*T,
      mu-T*log(rho_sum/nSamp));
  if (traj_fn) {
    fprintf(stdout,"Trajectory written to %s.\n",traj_fn);
  }
}

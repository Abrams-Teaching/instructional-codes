/* 
   XYZ-trajectory analyzer: computes velocity autocorrelation

   Cameron F. Abrams

   Written for the course CHE T580, Modern Molecular Simulation
   Spring 20-21

   compile using "gcc -o vacf vacf.c -lm"

   runs as "./vacf -t <trajctory file(traj.xyz)>"

   Note, the trajectory file must contain velocities!

   Drexel University, Department of Chemical Engineering
   Philadelphia
   (c) 2004-2021
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

char * elem[] = {"Ar", "Xe", "NULL"};

/* Abstract data type representing one frame of simulation data */
typedef struct FRAME {
    double * rx, * ry, * rz;  // coordinates
    double * vx, * vy, * vz;  // velocities
    int * typ; // array of particle types 0, 1, ...
    int N; // number of particles
    double Lx, Ly, Lz; // box dimensions
    double cx, cy, cz; // center of mass positions
    double cvx, cvy, cvz; // center of mass velocities
} frametype;

/* Create and return an empty frame */
frametype * NewFrame ( int N, int hv ) {
    frametype * f = (frametype*)malloc(sizeof(frametype));
    f->N=N;
    f->rx=(double*)malloc(sizeof(double)*N);
    f->ry=(double*)malloc(sizeof(double)*N);
    f->rz=(double*)malloc(sizeof(double)*N);
    if (hv) {
        f->vx=(double*)malloc(sizeof(double)*N);
        f->vy=(double*)malloc(sizeof(double)*N);
        f->vz=(double*)malloc(sizeof(double)*N);
    } else {
        f->vz=NULL;
        f->vy=NULL;
        f->vx=NULL;
    }
    f->typ=(int*)malloc(sizeof(int)*N);
    f->cx=f->cy=f->cz=0.0;
    f->cvx=f->cvy=f->cvz=0.0;
    return f;
}

/* Compute squared scalar distance between particles i and j in 
   frame fi and fj, respectively; com_corr removes center of mass motion */
double vij2_unwrapped ( frametype * fi, int i, frametype * fj, int j, int com_corr ) {
    double dx, dy, dz;
    dx=(fi->vx[i]-(com_corr?fi->cvx:0))*(fj->vx[j]-(com_corr?fj->cvx:0));
    dy=(fi->vy[i]-(com_corr?fi->cvy:0))*(fj->vy[j]-(com_corr?fj->cvy:0));
    dz=(fi->vz[i]-(com_corr?fi->cvz:0))*(fj->vz[j]-(com_corr?fj->cvz:0));
    return dx+dy+dz;
}

/* Read an XYZ-format frame from stream fp; returns the new frame.  Note
   the non-conventional use of the comment line to read in boxsize information. */
frametype * read_xyz_frame ( FILE * fp ) {
    int N,i,j,hasvel=0;
    double x, y, z, Lx, Ly, Lz, vx, vy, vz;
    char typ[3], dummy[5];
    char ln[255];
    frametype * f = NULL;
    if (fgets(ln,255,fp)){
        sscanf(ln,"%i %i\n",&N,&hasvel);
        f = NewFrame(N,hasvel);
        fgets(ln,255,fp);
        sscanf(ln,"%s %lf %lf %lf\n",dummy,&f->Lx,&f->Ly,&f->Lz);
        for (i=0;i<N;i++) {
            fgets(ln,255,fp);
            sscanf(ln,"%s %lf %lf %lf %lf %lf %lf\n",typ,&f->rx[i],&f->ry[i],&f->rz[i],&vx,&vy,&vz);
            if (hasvel) {
                f->vx[i]=vx;
                f->vy[i]=vy;
                f->vz[i]=vz;
            }
            f->cx+=f->rx[i];
            f->cy+=f->ry[i];
            f->cz+=f->rz[i];
            if (hasvel) {
                f->cvx+=f->vx[i];
                f->cvy+=f->vy[i];
                f->cvz+=f->vz[i];
            }
            j=0;
            while(strcmp(elem[j],"NULL")&&strcmp(elem[j],typ)) j++;
            if (strcmp(elem[j],"NULL")) f->typ[i]=j;
            else f->typ[i]=-1;
        }
        f->cx/=N;
        f->cy/=N;
        f->cz/=N;
        f->cvx/=N;
        f->cvy/=N;
        f->cvz/=N;
    }
    return f;
}

#define MAXFRAMES 10000
int main (int argc, char * argv[] ) {
    frametype * Traj[MAXFRAMES];
    int M = 0;
    int i,nc;
    char * trajfile=NULL;
    FILE * fp;
    char * length_units = "sigma";
    char * time_units = "sigma*sqrt(mass/epsilon)";
    double * sd;
    int t, dt, * cnt, tmax=0;
    double md_time_step = 0.001;
    int traj_interval=1000;
    char * outfile="vacf.dat";
    int begin_frame=0, MAnalyzed=0;
    for (i=1;i<argc;i++) {
        if (!strcmp(argv[i],"-t")) trajfile=argv[++i];
        else if (!strcmp(argv[i],"-o")) outfile=argv[++i];
        else if (!strcmp(argv[i],"-begin-frame")) begin_frame=atoi(argv[++i]);
        else if (!strcmp(argv[i],"-traj-interval")) traj_interval=atoi(argv[++i]);
        else if (!strcmp(argv[i],"-md-time-step")) md_time_step=atof(argv[++i]);
        else if (!strcmp(argv[i],"-length-units")) length_units=argv[++i];
        else if (!strcmp(argv[i],"-time-units")) time_units=argv[++i];
    }
    if (!trajfile) {
        fprintf(stdout,"Error: a trajectory file must be specified with -t\n");
        exit(-1);
    }

    i=0;
    fprintf(stdout,"Reading %s\n",trajfile);fflush(stdout);
    fp=fopen(trajfile,"r");
    while (Traj[i++]=read_xyz_frame(fp));
    M=i-1;
    fclose(fp);
    if (!M) {
        fprintf(stdout,"Error: trajectory %s has no data.\n",trajfile);
        exit(-1);
    }
    fprintf(stdout,"Read %i frames from %s.\n",M,trajfile);

    /* Allocate and initialize the squared-displacement */
    sd=(double*)calloc(M-begin_frame,sizeof(double));
    cnt=(int*)calloc(M-begin_frame,sizeof(int));

    /* Compute the mean-squared displacement using
     the straightforward algorithm */
    fprintf(stdout,"# computing VACF...\n");fflush(stdout);
    for (t=begin_frame;t<M;t++) {
        for (dt=0;(t+dt)<M;dt++) {
            cnt[dt]++;  /* number of origins for interval length dt  */
            for (i=0;i<Traj[0]->N;i++) {
	            sd[dt] += vij2_unwrapped(Traj[t+dt],i,Traj[t],i,0);
            }
        }
    }
    fp=fopen(outfile,"w");
    fprintf(fp,"# VACF from %s\n",trajfile);
    fprintf(fp,"#LABEL time msd\n");
    fprintf(fp,"#UNITS %s %s^2\n",time_units,length_units);
    for (t=0;t<M-begin_frame;t++) {
        sd[t] /= cnt[t]?(Traj[0]->N*cnt[t]):1;
        fprintf(fp,"% .5lf % .8lf\n",
	        t*traj_interval*md_time_step,sd[t]);
    }
    fclose(fp);
    fprintf(stdout,"%s created.\n",outfile);
}
/* 
   XYZ-trajectory analyzer: computes velocity distribution

   Cameron F. Abrams

   Written for the course CHE T580, Modern Molecular Simulation
   Spring 20-21

   compile using "gcc -o vdf vdf.c -lm"

   runs as "./vdf -t <trajctory file(traj.xyz)>"

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

/* generic histogram allocator */
int * allocate_hist(double bs, double vmin, double vmax) {
    int * H=NULL;
    int nbins = (int)((vmax-vmin)/bs+1);
    H=calloc(nbins,sizeof(int));
    return H;
}
/* generic histogram binner */
int binit (double x, int * H, double bs, double xmin, double xmax) {
    int bin=(int)((x-xmin)/bs);
    if (x>=xmin&&x<=xmax) {
        H[bin]++;
        return 0;
    }
    return 1;
}
/* update velocity histogram that is symmetric about 0 [-vmax,vmax] */
int update_vhist ( frametype * f, int * vHist, double vbs, double vmax ) {
    int i;
    int result=0;
    for (i=0;i<f->N&&!result;i++) {
        result = binit(f->vx[i],vHist,vbs,-vmax,vmax) &&
                 binit(f->vy[i],vHist,vbs,-vmax,vmax) &&
                 binit(f->vz[i],vHist,vbs,-vmax,vmax);
    }
    return result;
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

#define MAXFRAMES 100000
int main (int argc, char * argv[] ) {
    frametype * Traj[MAXFRAMES];
    int M = 0;
    int i,nc,result;
    char * trajfile=NULL;
    FILE * fp;
    char * length_units = "sigma";
    char * time_units = "sigma*sqrt(mass/epsilon)";
    double * sd;
    int t, dt, * cnt, tmax=0;
    double md_time_step = 0.001;
    int traj_interval=1000;
    char * outfile="vdf.dat";
    int begin_frame=0, MAnalyzed=0;
    double vbs=0.01,vmax=100.0, v;
    int * vHist, histsum;
    for (i=1;i<argc;i++) {
        if (!strcmp(argv[i],"-t")) trajfile=argv[++i];
        else if (!strcmp(argv[i],"-o")) outfile=argv[++i];
        else if (!strcmp(argv[i],"-begin-frame")) begin_frame=atoi(argv[++i]);
        else if (!strcmp(argv[i],"-traj-interval")) traj_interval=atoi(argv[++i]);
        else if (!strcmp(argv[i],"-md-time-step")) md_time_step=atof(argv[++i]);
        else if (!strcmp(argv[i],"-length-units")) length_units=argv[++i];
        else if (!strcmp(argv[i],"-time-units")) time_units=argv[++i];
        else if (!strcmp(argv[i],"-vdf-bin-size")) vbs=atof(argv[++i]);
        else if (!strcmp(argv[i],"-vdf-vmax")) vmax=atof(argv[++i]);
    }
    if (!trajfile) {
        fprintf(stdout,"Error: a trajectory file must be specified with -t\n");
        exit(-1);
    }
    vHist=allocate_hist(vbs,-vmax,vmax);
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

    /* Compute velocity histogram */
    fprintf(stdout,"# computing VDF...\n");fflush(stdout);
    result=0;
    for (t=begin_frame;t<M&&!result;t++) {
        result=update_vhist(Traj[t],vHist,vbs,vmax);
    }
    if (result) {
        fprintf(stdout,"ERROR: bin violation.\n");
        exit(-1);
    }
    fp=fopen(outfile,"w");
    fprintf(fp,"# VDF from %s\n",trajfile);
    fprintf(fp,"#LABEL v freq\n");
    fprintf(fp,"#UNITS %s/%s %s/%s^(-1)\n",length_units,time_units,length_units,time_units);
    i=0;
    histsum=0;
    for (v=-vmax;v<=vmax;v+=vbs) {
        histsum+=vHist[i++];
    }
    i=0;
    for (v=-vmax;v<=vmax;v+=vbs) {
        fprintf(fp,"% .5lf % .8lf\n",v,((double)vHist[i])/(histsum*vbs));
        i++;
    }
    fclose(fp);
    fprintf(stdout,"%s created.\n",outfile);
}
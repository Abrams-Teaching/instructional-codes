/* 
   XYZ-trajectory analyzer: computes radial distribution functions

   Cameron F. Abrams

   Written for the course CHE T580, Modern Molecular Simulation
   Spring 20-21

   compile using "gcc -o rdf rdf.c -lm"

   runs as "./rdf -t <trajctory file(traj.xyz)> -dr <resolution(0.1)> -rcut <cutoff radius(3.5)>"

   Note, the trajectory file must contain the box size in the comment line for
   each frame (this is my convention, not the standard)

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
    return f;
}

/* Compute scalar distance between particles i and j in frame f;
   note the use of the minimum image convention */
double rij ( frametype * f, int i, int j ) {
    double dx, dy, dz;
    double hLx=0.5*f->Lx,hLy=0.5*f->Ly,hLz=0.5*f->Lz;
    dx=f->rx[i]-f->rx[j];
    dy=f->ry[i]-f->ry[j];
    dz=f->rz[i]-f->rz[j];
    if (dx<-hLx) dx+=f->Lx;
    if (dx> hLx) dx-=f->Lx;
    if (dy<-hLy) dy+=f->Ly;
    if (dy> hLy) dy-=f->Ly;
    if (dz<-hLz) dz+=f->Lz;
    if (dz> hLz) dz-=f->Lz;
    return sqrt(dx*dx+dy*dy+dz*dz);
}

/* An N^2 algorithm for computing interparticle separations
   and updating the radial distribution function histogram. */
void update_hist ( frametype * f, double rcut, double dr, int * H, int nbins ) {
    int i,j;
    double r;
    int bin;
    for (i=0;i<f->N-1;i++) {
        for (j=i+1;j<f->N;j++) {
            r=rij(f,i,j);
            if (r<rcut) {
                bin=(int)(r/dr);
                if (bin<0||bin>=nbins) {
                    fprintf(stderr,"Warning: bin range violation: %.3lf not on [0.0,%.3lf]\n",r,rcut);
                } else {
                    H[bin]+=2;
                }
            }
        }
    }
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
            j=0;
            while(strcmp(elem[j],"NULL")&&strcmp(elem[j],typ)) j++;
            if (strcmp(elem[j],"NULL")) f->typ[i]=j;
            else f->typ[i]=-1;
        }
    }
    return f;
}

double min(double x, double y ) {
    if (x<y) return x;
    else return y;
}

#define MAXFRAMES 10000
int main (int argc, char * argv[] ) {
    frametype * Traj[MAXFRAMES];
    int nFrames = 0;
    int i,nc;
    char * trajfile=NULL;
    FILE * fp;
    double dr=0.1,rcut=3.5, rho, r, vb, nid, L2min;
    int * H, nbins=0;
    char * outfile="rdf.dat";
    for (i=1;i<argc;i++) {
        if (!strcmp(argv[i],"-t")) trajfile=argv[++i];
        else if (!strcmp(argv[i],"-dr")) dr=atof(argv[++i]);
        else if (!strcmp(argv[i],"-rcut")) rcut=atof(argv[++i]);
        else if (!strcmp(argv[i],"-o")) outfile=argv[++i];
    }
    if (!trajfile) {
        fprintf(stdout,"Error: a trajectory file must be specified with -t\n");
        exit(-1);
    }

    i=0;
    fprintf(stdout,"Reading %s\n",trajfile);fflush(stdout);
    fp=fopen(trajfile,"r");
    while (Traj[i++]=read_xyz_frame(fp));
    nFrames=i-1;
    fclose(fp);
    if (!nFrames) {
        fprintf(stdout,"Error: trajectory %s has no data.\n",trajfile);
    }
    fprintf(stdout,"Read %i frames from %s.\n",nFrames,trajfile);

    /* Adjust cutoff and compute histogram */
    L2min=min(Traj[0]->Lx/2,min(Traj[0]->Ly/2,Traj[0]->Lz/2));
    if (rcut>L2min) rcut=L2min;
    nbins=(int)(rcut/dr)+1;
    H=(int*)malloc(sizeof(int)*nbins);
    for (i=0;i<nbins;i++) H[i]=0;
    for (i=0;i<nFrames;i++) update_hist(Traj[i],rcut,dr,H,nbins);

    /* Normalize and output g(r) to the terminal */
    /* Compute density, assuming NVT ensemble */
    fp=fopen(outfile,"w");
    fprintf(fp,"# RDF from %s\n",trajfile);
    rho=Traj[0]->N/(Traj[0]->Lx*Traj[0]->Ly*Traj[0]->Lz);
    for (i=0;i<nbins-1;i++) {
        r=dr*(i+0.5);
        vb=((i+1)*(i+1)*(i+1)-i*i*i)*dr*dr*dr;
        nid=(4./3.)*M_PI*vb*rho;
        fprintf(fp,"%.4lf %.4lf\n",i*dr,(double)(H[i])/(nFrames*Traj[0]->N*nid));
    }
    fclose(fp);
    fprintf(stdout,"%s created.\n",outfile);
}
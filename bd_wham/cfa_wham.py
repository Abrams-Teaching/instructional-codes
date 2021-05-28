# Weighted Histogram Analysis Method
# Cameron F Abrams cfa22@drexel.edu

import matplotlib.pyplot as plt
import numpy as np
import argparse as ap
import matplotlib.cm as cm
import scipy.constants as sc

kB_kcal_per_mol=sc.physical_constants['Boltzmann constant'][0]/sc.calorie*sc.Avogadro/1.e3

parser=ap.ArgumentParser()
parser.add_argument('-i',metavar='<str>',type=str,default='info.in',help='input file')
parser.add_argument('-nbins',metavar='<int>',type=int,default=1000,help='number of bins')
parser.add_argument('-zlim',metavar=('zmin','zmax'),type=float,default=[2,5],nargs=2,help='limits on histogram domain')
parser.add_argument('-T',metavar='<float>',type=float,default=310,help='temperature')
parser.add_argument('-tol',metavar='<float>',type=float,default=1.0e-6,help='WHAM tolerance')
parser.add_argument('-figsize',type=float,nargs=2,default=[6,4],help='figure size')
parser.add_argument('-skip',type=int,default=0,help='skip this many lines at beginning of raw data files')
parser.add_argument('-of',metavar='<str>',type=str,default='F.dat',help='output data file F vs z')
parser.add_argument('-o',metavar='<str>',type=str,default='plot.png',help='output plot of F vs z')
parser.add_argument('-traj-filename-format',metavar='<str>',type=str,default='frm{:d}_{:.0f}-run.colvars.traj',help='traj filename format')
parser.add_argument('-print-every',metavar='<int>',type=int,default=10,help='print output to terminal every this many WHAM iterations')
args=parser.parse_args()

# load the input file.   Input file has one line per window simulation.
# Each line reports z_0 (window potential anchor point), unique id number, and k (spring constant)
setdat=np.loadtxt(args.i).T
z0=setdat[0] # window potential centers
ids=[int(_) for _ in setdat[1]]
Ks=setdat[2] # spring constants

# Read in each *.traj file, generate histogram, keep tally 
H=[] # list of histograms
M=[] # list of histogram sums
for i,k in zip(ids,Ks):
    fn=args.traj_filename_format.format(i,k)
    print("# Reading/processing {:s}...".format(fn))
    t,z=np.loadtxt(fn,unpack=True)
    h,e=np.histogram(z[args.skip:],args.nbins,range=args.zlim)
    H.append(h)
    M.append(h.sum())
n=len(H) # number of windows is inferred
M=np.array(M)
print("# Computed {:d} histograms".format(n))

# apparent dz from user-supplied domain limits and number of bins
dz=(args.zlim[1]-args.zlim[0])/args.nbins
# generate array of bin center positions; one fewer than number of bin edges returned from np.histogram
zz=e[:-1]+0.5*dz
nz=len(zz)
print("# Z domain [{:.5f},{:.5f}] in {:d} increments of {:.5f}".format(zz[0],zz[-1],nz,dz))

plt,ax=plt.subplots(1,2,figsize=args.figsize)
ax[0].set_xlabel('$z$')
ax[0].set_ylabel('$p_i(z)$')
cmap=cm.get_cmap('inferno')
for i,(z,h) in enumerate(zip(z0,H)):
    ax[0].plot(zz,h,label='$z_0$ = {:.2f}'.format(z),color=cmap(i/len(H)))

# wham

# window potential
def W(x,x0,k):
    return 0.5*k*(x-x0)**2

F=np.zeros(n) # initial window F's

beta=1/(kB_kcal_per_mol*args.T)

# calculate sum of all histograms
sumH=np.zeros(nz)
for i in range(n):
    sumH+=H[i]

# WHAM convergence loop
iterating = True
ii=0
while iterating:
    iterating=False
    P0=sumH.copy()
    for i in range(nz):
        d=0.0
        for j,(z,k) in enumerate(zip(z0,Ks)):
            d+=M[j]*np.exp(-beta*(W(zz[i],z,k)-F[j]))
        P0[i]/=d
    sqerr=0.0
    for j,(z,k) in enumerate(zip(z0,Ks)):
        newF=-1/beta*np.log((P0*np.exp(-beta*W(zz,z,k))).sum())
        this_dF=F[j]-newF
        F[j]=newF
        sqerr+=this_dF**2
    if ii%args.print_every==0:
        print('{:d} {:.5e}'.format(ii,sqerr))
    if sqerr>args.tol:
        iterating=True
        ii+=1
        

ax[1].plot(zz,-1./beta*np.log(P0))
ax[1].set_xlabel('z')
ax[1].set_ylabel('$F(z)$ (kcal/mol)')

plt.savefig(args.o,bbox_inches='tight')
np.savetxt(args.of,np.array([zz,-1/beta*np.log(P0)]).T)
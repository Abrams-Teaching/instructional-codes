import matplotlib.pyplot as plt
import numpy as np
import argparse as ap
import matplotlib.cm as cm
import warnings

parser=ap.ArgumentParser()

parser.add_argument('-n',metavar='<int>',type=int,default=16,help='number of windows')
parser.add_argument('-pxlim',metavar=('p-x-min','p-y-min'),type=float,nargs=2,default=[-12,12],help='limits on x-axis for plots')
parser.add_argument('-zlim',metavar=('zmin','zmax'),type=float,default=[-8,8],nargs=2,help='limits on histogram domain')
parser.add_argument('-k',metavar='<float>',type=float,default=15,help='window potential spring constant')
parser.add_argument('-T',metavar='<float>',type=float,default=1.0,help='temperature')
parser.add_argument('-tol',metavar='<float>',type=float,default=1.0e-6,help='WHAM tolerance')
parser.add_argument('-fylim',metavar=('y-low','y-high'),type=float,default=[0,25],nargs=2,help='ylims for potential plot')
parser.add_argument('-pscale',metavar='<str>',type=str,default='linear',help='y-axis scaling for histograms')
parser.add_argument('-p',metavar='<str>',type=str,default='p{:d}.dat',help='filename format for window histograms')
parser.add_argument('-fpot',metavar='<str>',type=str,default='V.dat',help='file containing governing potential vs z')
parser.add_argument('-w',metavar='<str>',type=str,default='W{:d}.dat',help='files containing window potentials')
parser.add_argument('-o',metavar='<str>',type=str,default='out.png',help='output plot image file')
parser.add_argument('-of',metavar='<str>',type=str,default='F.dat',help='output data file containing F vs z')
args=parser.parse_args()

fdat=np.loadtxt(args.fpot).T

wdat=[]
pdat=[]
for i in range(args.n):
    wdat.append(np.loadtxt(args.w.format(i)).T)
    pdat.append(np.loadtxt(args.p.format(i)).T)

plt,ax=plt.subplots(1,3,figsize=(16,5))

ax[0].set_xlabel('$x$')
ax[0].set_ylabel('$V(x)$, $W_i(x)$')
ax[0].set_ylim(args.fylim)
ax[0].set_xlim(args.pxlim)
cmap=cm.get_cmap('inferno')
for i in range(args.n):
    ax[0].plot(wdat[i][0],wdat[i][1],label='',color=cmap(i/args.n),alpha=0.4)
ax[0].plot(fdat[0],fdat[1]-fdat[1].min(),label='V(x)',linewidth=2)
ax[0].legend()

ax[1].set_xlabel('$x$')
ax[1].set_ylabel('$H_i(x)$')
ax[1].set_xlim(args.pxlim)
M=[]
for i in range(args.n):
    ax[1].plot(0.5*(pdat[i][0]+pdat[i][1]),pdat[i][2],label='',color=cmap(i/args.n),alpha=0.4)
    M.append(pdat[i][2].sum())


# wham

def W(x,x0,k):
    return 0.5*k*(x-x0)**2

F=np.zeros(args.n)
M=np.array(M)
k = args.k
win_edges = np.linspace(args.zlim[0],args.zlim[1],args.n+1)
wz=win_edges[:-1]+0.5*(args.zlim[1]-args.zlim[0])/args.n
Z=fdat[0][:]

beta=1/args.T
sumH=np.zeros(len(Z))
for i in range(args.n):
    sumH+=pdat[i][2]
iterating = True
ii=0
print_every=10
while iterating:
    iterating=False
    P0=sumH.copy()
    for i in range(len(Z)):
        d=0.0
        for j,w in enumerate(wz):
            d+=M[j]*np.exp(-beta*(W(Z[i],w,k)-F[j]))
        P0[i]/=d
       # print(d)
    sqerr=0.0
    for j,w in enumerate(wz):
        newF=-args.T*np.log((P0*np.exp(-beta*W(Z,w,k))).sum())
        this_dF=F[j]-newF
        F[j]=newF
        sqerr+=this_dF**2
    if ii%print_every==0:
        print('{:d} {:.5e}'.format(ii,sqerr))
    if sqerr>args.tol:
        iterating=True
        ii+=1

ymin,ymax=ax[1].get_ylim()
pscale=ymax/P0.max()
ax[1].plot(Z,P0*pscale,label='{:.2e}$P_0$'.format(pscale),linewidth=2)
ax[1].legend()
ax[1].set_yscale(args.pscale)
ax[2].set_xlim(args.pxlim)
ax[2].set_yscale(args.pscale)
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    F=-args.T*np.log(P0)
ax[2].plot(fdat[0],fdat[1]-fdat[1].min(),label='V(x)',linewidth=2,alpha=0.2)
ax[2].plot(Z,F-F.min(),label='F(x)',linewidth=2)
ax[2].set_xlabel('x')
ax[2].set_ylabel('$F(x) = -k_BT\ln P_0$')
ax[2].set_ylim(args.fylim)
ax[2].legend()
plt.savefig(args.o,bbox_inches='tight')
np.savetxt(args.of,np.array([Z,F]).T)
import matplotlib.pyplot as plt
import numpy as np
import argparse as ap
import matplotlib.cm as cm

parser=ap.ArgumentParser()

parser.add_argument('-n',type=int,default=16,help='number of windows')
parser.add_argument('-zlim',type=float,default=[-8,8],nargs='+',help='limits on histogram domain')
parser.add_argument('-k',type=float,default=16,help='window potential spring constant')
parser.add_argument('-T',type=float,default=1.0,help='temperature')
parser.add_argument('-tol',type=float,default=1.0e-6,help='WHAM tolerance')
parser.add_argument('-fylim',type=float,default=[-15,15],nargs='+',help='ylims for potential plot')
parser.add_argument('-p',type=str,default='p{:d}.dat',help='filename format for window histograms')
parser.add_argument('-fpot',type=str,default='f-pot.dat',help='file containing governing potential vs z')
parser.add_argument('-w',type=str,default='win-pot{:d}.dat',help='files containing window potentials')
parser.add_argument('-o',type=str,default='out.png',help='output plot image file')
parser.add_argument('-of',type=str,default='F.dat',help='output data file containing F vs z')
args=parser.parse_args()

fdat=np.loadtxt(args.fpot).T

wdat=[]
pdat=[]
for i in range(args.n):
    wdat.append(np.loadtxt(args.w.format(i)).T)
    pdat.append(np.loadtxt(args.p.format(i)).T)

plt,ax=plt.subplots(1,3,figsize=(16,5))

ax[0].set_xlabel('x')
ax[0].set_ylabel('f(x),p(x)')
ax[0].set_ylim(args.fylim)
#ax[1].plot(fdat[0],fdat[1],label='f')
cmap=cm.get_cmap('inferno')

for i in range(args.n):
    ax[0].plot(wdat[i][0],wdat[i][1],label='{:d}'.format(i),color=cmap(i/args.n))
ax[0].plot(fdat[0],fdat[1],label='f(x)')

#ax[0].legend()

ax[1].set_xlabel('x')
ax[1].set_ylabel('p(x)')
M=[]
for i in range(args.n):
    ax[1].plot(0.5*(pdat[i][0]+pdat[i][1]),pdat[i][2],label='{:d}'.format(i),color=cmap(i/args.n))
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
        

ax[2].plot(Z,-args.T*np.log(P0),label='P0')
ax[2].set_xlabel('x')
ax[2].set_ylabel('$F(x)$')

plt.savefig(args.o,bbox_inches='tight')
np.savetxt(args.of,np.array([Z,-args.T*np.log(P0)]).T)
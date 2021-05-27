import matplotlib.pyplot as plt
import numpy as np
import argparse as ap
import matplotlib.cm as cm

parser=ap.ArgumentParser()
parser.add_argument('-n',type=int,default=16,help='number of windows')
parser.add_argument('-xrange',type=float,default=[-8,8],nargs='+',help='x-range')
parser.add_argument('-k',type=float,default=10,help='k')
parser.add_argument('-fylim',type=float,default=[-1,1],nargs='+',help='ylims for potential')
parser.add_argument('-w',type=str,default='win-pot{:d}.dat',help='file with window potential data')
parser.add_argument('-f',type=str,default='f-pot.dat',help='file with potential data')
parser.add_argument('-p',type=str,default='p{:d}.dat',help='file with histogram data')
#parser.add_argument('-l',type=str,default='log.dat',help='log of x')
#parser.add_argument('-hist',type=str,default='p0.dat',help='file with histogram')
#parser.add_argument('-hist-scale',type=float,default=1.0,help='scale the histogram')
parser.add_argument('-o',type=str,default='out.png',help='output file')

args=parser.parse_args()

fdat=np.loadtxt(args.f).T
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

# window potential
def W(x,x0,k):
    return 0.5*k*(x-x0)**2

F=np.zeros(args.n)
M=np.array(M)
k = args.k
wxs = np.linspace(args.xrange[0],args.xrange[1],args.n)
X=fdat[0][0][:]
T=1.0
beta=1/T
sumH=np.zeros(len(X))
for i in range(args.n):
    sumH+=pdat[i][2]
iterating = True
tolerance = 1.e-6
ii=0
print_every=10
while iterating:
    iterating=False
    P0=sumH.copy()
    for i in range(len(X)):
        d=0.0
        for j,w in enumerate(wxs):
            d+=M[j]*np.exp(-beta*(W(X[i],w,k)-F[j]))
        P0[i]/=d
       # print(d)
    sqerr=0.0
    for j,w in enumerate(wxs):
        newF=-T*np.log((P0*np.exp(-beta*W(X,w,k))).sum())
        this_dF=F[j]-newF
        F[j]=newF
        sqerr+=this_dF**2
    if ii%print_every==0:
        print('{:d} {:.5e}'.format(ii,sqerr))
    if sqerr>tolerance:
        iterating=True
        ii+=1
        
ax[2].plot(X,-T*np.log(P0),label='P0')
ax[2].set_xlabel('x')
ax[2].set_ylabel('$F(x)$')

plt.savefig(args.o,bbox_inches='tight')
np.savetxt("F.dat",np.array([X,-T*np.log(P0)]).T)
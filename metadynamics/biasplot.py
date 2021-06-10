import matplotlib.pyplot as plt
import argparse as ap
import numpy as np
import scipy.constants as sc
import matplotlib.cm as cm
import warnings

def kernel(z,z0,s,w):
    return w*np.exp(-(z-z0)**2/s)

kB_kcal_per_mol=sc.physical_constants['Boltzmann constant'][0]/sc.calorie*sc.Avogadro/1.e3

parser=ap.ArgumentParser()
parser.add_argument('-hills-traj',metavar='<str>',default='in.meta_r.hills.traj',type=str,help='hills traj file')
#parser.add_argument('-series',nargs=3,type=int,default=[0,0,0])
#parser.add_argument('-series-average',action='store_true',help='show series average in plot')
parser.add_argument('-o',metavar='<str>',default='plot.png',type=str,help='output image file')
#parser.add_argument('-ylim',default=[0,10],nargs=2,type=float,help='plot ylimits')
parser.add_argument('-xlim',default=[2.5,4.5],nargs=2,type=float,help='plot xlimits')
parser.add_argument('-label',default='Metadynamics',type=str)
parser.add_argument('-T',default=310,type=float)
#parser.add_argument('-traj-long-md',metavar='<str>',type=str,default='',help='traj file of long MD run')
#parser.add_argument('-pmf-long-md',metavar='<str>',type=str,default='',help='pmf file of long MD run')
parser.add_argument('-nbins',type=int,default=100)
parser.add_argument('-alpha',type=float,default=0.1)
parser.add_argument('-lw',type=float,default=0.5)
parser.add_argument('-figsize',type=float,nargs=2,default=[7,6])
parser.add_argument('-fontsize',type=float,default=14)
parser.add_argument('-every',type=int,default=1)
args=parser.parse_args()

T,Z,S,W=np.loadtxt(args.hills_traj,unpack=True)

plt,ax=plt.subplots(1,1,figsize=args.figsize)

X=np.linspace(args.xlim[0],args.xlim[1],args.nbins)
Vb=np.zeros(len(X))

#ax.set_ylim(args.ylim)
ax.set_xlim(args.xlim)
ax.set_xlabel('$z$ (Angstrom)',fontsize=args.fontsize)
ax.set_ylabel('$V_b(z)$ (kcal/mol)',fontsize=args.fontsize)
cmap = cm.get_cmap('plasma')
for i,(t,z,s,w) in enumerate(zip(T,Z,S,W)):
    tVb=kernel(X,z,s,w)
    Vb+=tVb
    if i%args.every==0:
        ax.plot(X,Vb-Vb.min(),color=cmap(i/len(T)),alpha=args.alpha,linewidth=args.lw)


plt.savefig(args.o,bbox_inches='tight')

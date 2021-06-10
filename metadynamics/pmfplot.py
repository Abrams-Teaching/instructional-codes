import matplotlib.pyplot as plt
import argparse as ap
import numpy as np
import scipy.constants as sc
import matplotlib.cm as cm
import warnings

kB_kcal_per_mol=sc.physical_constants['Boltzmann constant'][0]/sc.calorie*sc.Avogadro/1.e3

parser=ap.ArgumentParser()
parser.add_argument('-i',metavar='<str>',default='in.pmf',type=str,help='pmf file')
parser.add_argument('-series',nargs=3,type=int,default=[0,0,0])
parser.add_argument('-series-average',action='store_true',help='show series average in plot')
parser.add_argument('-opmf',type=str,default='')
parser.add_argument('-plot-intermediates',action='store_true')
parser.add_argument('-o',metavar='<str>',default='plot.png',type=str,help='output image file')
parser.add_argument('-ylim',default=[0,10],nargs=2,type=float,help='plot ylimits')
parser.add_argument('-xlim',default=[2.5,4.5],nargs=2,type=float,help='plot xlimits')
parser.add_argument('-label',default='Metadynamics',type=str)
parser.add_argument('-T',default=310,type=float)
parser.add_argument('-traj-long-md',metavar='<str>',type=str,default='',help='traj file of long MD run')
parser.add_argument('-pmf-long-md',metavar='<str>',type=str,default='',help='pmf file of long MD run')
parser.add_argument('-nbins',type=int,default=100)
parser.add_argument('-figsize',type=float,nargs=2,default=[7,6])
parser.add_argument('-fontsize',type=float,default=14)
parser.add_argument('-alpha',type=float,default=0.1)
args=parser.parse_args()

gotlong=False
if args.traj_long_md != '':
    gotlong=True
    fn=args.traj_long_md
    print("# Reading/processing md trajectory data from {:s}...".format(fn))
    t,z=np.loadtxt(fn,unpack=True)
    hlong,e=np.histogram(z,args.nbins,range=args.xlim)
    beta=1/(kB_kcal_per_mol*args.T)
    flong=-1./beta*np.log(hlong)
    flong=flong-flong.min()
    # for binning long md data only
    # apparent dz from user-supplied domain limits and number of bins
    dz=(args.xlim[1]-args.xlim[0])/args.nbins
    # generate array of bin center positions; one fewer than number of bin edges returned from np.histogram
    zz=e[:-1]+0.5*dz
    if args.pmf_long_md != '':
        np.savetxt(args.pmf_long_md,np.array([zz,flong]).T,header='z F(z)')
        print('Wrote {:s}.'.format(args.pmf_long_md))
elif args.pmf_long_md != '':
    gotlong=True
    zz,flong=np.loadtxt(args.pmf_long_md,unpack=True)

plt,ax=plt.subplots(1,1,figsize=args.figsize)

ax.set_ylim(args.ylim)
ax.set_xlim(args.xlim)
ax.set_xlabel('$z$ (Angstrom)',fontsize=args.fontsize)
ax.set_ylabel('$F(z)$ (kcal/mol)',fontsize=args.fontsize)
if args.series[2]!=0:
    cmap=cm.get_cmap('inferno')
    f={}
    ni=int((args.series[1]-args.series[0])/args.series[2])
    ffa=[]
    for i,s in enumerate(range(args.series[0],args.series[1],args.series[2])):
        fn=args.i.format(s)
        z,f[i]=np.loadtxt(fn,unpack=True)
        if len(ffa)==0:
            ffa=f[i].copy()
        else:
            ffa+=f[i]
        if args.plot_intermediates:
            ax.plot(z,f[i],color=cmap(i/ni),alpha=args.alpha)
    fn=args.i.format(args.series[1])
    z,ff=np.loadtxt(fn,unpack=True)
    ffa+=ff
    if args.series_average:
        ax.plot(z,ffa/ni,color='black',label='{:s}-avg'.format(args.label))
        if args.opmf!='':
            np.savetxt(args.opmf,np.array([z,ffa/ni]).T,header='z PMF')
            print('Created {:s}.'.format(args.opmf))
    else:
        ax.plot(z,ff,color='black',label=args.label)

else:
    z,f=np.loadtxt(args.i,unpack=True)
    ax.plot(z,f,label=args.label)
if gotlong:
    ax.plot(zz,flong-flong.min(),'g--',label='long MD',alpha=args.alpha)
    ax.legend()

plt.savefig(args.o,bbox_inches='tight')

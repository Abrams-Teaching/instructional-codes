import matplotlib.pyplot as plt 
import argparse as ap
import numpy as np

def get_f(p,beta=1.0):
    f=[]
    du=[]
    for l,r,pdu in zip(p[0],p[1],p[2]):
        tdu=(l+r)/2
        try:
            lpdu=np.log(pdu)
        except:
            continue
        f.append(lpdu+beta/2*tdu)
        du.append(tdu)
    return np.array(du),np.array(f)

def get_df (du0,f0,du1,f1):
    ddu=[]
    df=[]
    for d0,ff0 in zip(du0,f0):
        i=np.where(du1==d0)[0]
        if i!=-1:
            ddu.append(d0)
            df.append(f1[i]-f0[i])
    return np.array(ddu),np.array(df)


parser=ap.ArgumentParser()
parser.add_argument('-p0',type=str,default='p0.dat',help='Datafile with p0')
parser.add_argument('-p1',type=str,default='p1.dat',help='Datafile with p1')
parser.add_argument('-scale-p0',type=float,default=1.0,help='plot scale factor for p0')
parser.add_argument('-T',type=float,default=1.0,help='Temperature')
parser.add_argument('-xlim',type=float,nargs=2,default=[-30,30],help='xlim')
parser.add_argument('-fylim',type=float,nargs=2,default=[-30,10],help='fylim')
parser.add_argument('-figsize',type=float,nargs=2,default=[9,3],help='figsize')
parser.add_argument('-widom-est',type=float,default=0,help='widom estimate')
parser.add_argument('-o',type=str,default='plot.png',help='name of output graphic file')
args=parser.parse_args()

beta = 1.0/args.T

p0=np.loadtxt(args.p0).T
p1=np.loadtxt(args.p1).T

du0,f0=get_f(p0,-beta)
du1,f1=get_f(p1,beta)
ddu,df=get_df(du0,f0,du1,f1)
goodones=~np.isnan(df)
goodones*=~np.isinf(df)

fig,ax=plt.subplots(1,2,figsize=args.figsize)
plt.subplots_adjust(wspace=0.25)
ax[0].plot(0.5*(p0[0]+p0[1]),p0[2],label='$p_0$')
ax[0].plot(0.5*(p1[0]+p1[1]),p1[2],label='$p_1$')
ax[0].set_yscale('log')
ax[0].set_xlim(args.xlim)
ax[0].set_xlabel('$\Delta U$')
ax[0].set_ylabel('$p_0(\Delta U)$, $p_1(\Delta U)$')
ax[1].plot(du0,f0,label='$f_0$')
ax[1].plot(du1,f1,label='$f_1$')
ax[1].plot(ddu,df,label='$f_1-f_0$; $\mu_{{ex}}$={:.2f}'.format(df.mean(where=goodones)/beta))
ax[1].set_xlim(args.xlim)
ax[1].set_ylim(args.fylim)
ax[1].plot(args.xlim,[args.widom_est*beta]*2,label='Widom est.')
ax[1].set_xlabel('$\Delta U$')
ax[1].set_ylabel('$f_0(\Delta U)$, $f_1(\Delta U)$')

ax[0].legend()
ax[1].legend()
plt.subplots_adjust(wspace=0.25)
plt.savefig(args.o,bbox_inches='tight')

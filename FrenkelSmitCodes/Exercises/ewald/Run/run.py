import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tabulate import tabulate
ncell=8
dat={}
df={}
for pgm in ['./ewald','./ewald-l']:
    dat[pgm]={}
    for alpha in np.linspace(1.0,4.0,7):
        for kmax in [4,8,16]:
            fp=open('fort.21','w')
            fp.write('{:d} {:.2f} {:d}\n'.format(ncell,1./alpha,kmax))
            fp.close()
            os.system(pgm+' > out')
            if 'cutoff' in dat[pgm]:
                dat[pgm]['cutoff'].append(alpha)
            else:
                dat[pgm]['cutoff']=[alpha]
            with open('out') as fp:
                for l in fp:
                    tok=l.strip().split(':')
                    if len(tok)<2:
                        continue
                    k=tok[0].strip()
                    v=tok[1].strip()
                    if k in dat[pgm]:
                        dat[pgm][k].append(v)
                    else:
                        dat[pgm][k]=[v]

    df[pgm]=pd.DataFrame.from_dict(dat[pgm])
    print(tabulate(df[pgm][['cutoff','Kmax','Real Part','Fourier Part','Self Part','Total','Madelung Constant']],
          headers='keys',tablefmt='latex',showindex=False,
          floatfmt=['.2f','d','.4f','.4e','.4f','.4f','.4f']))
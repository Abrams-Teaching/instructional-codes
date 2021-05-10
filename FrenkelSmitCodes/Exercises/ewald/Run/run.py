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
    for alpha in np.linspace(0.2,1.4,7):
        for kmax in [8]:
            fp=open('fort.21','w')
            fp.write('{:d} {:.2f} {:d}\n'.format(ncell,alpha,kmax))
            fp.close()
            os.system(pgm+' > out')
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
    print(tabulate(df[pgm][['Alpha','Kmax','Real Part','Fourier Part','Self Part','Total','Madelung Constant']],headers='keys',tablefmt='latex',floatfmt='.4f',showindex=False))
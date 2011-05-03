from __future__ import division

import scipy.signal
import numpy as np

def P_factor(ux,uy,sx,sy):
    return lambda x,y: 1/(2* np.pi) * np.exp(-((x-ux)**2/(2 * sx**2))
                                       -((y-uy)**2/(2 * sx**2))
                                       )/(sx*sy)




def run(ux,uy,sx,sy,total,bkg,count,wind):
    X = np.vstack([np.arange(-wind,wind+1,1)]*(2*wind+1))
    Y = X.T
    def _run(P):
        tmp = np.random.poisson(
            np.array([[P(x,y) for x in np.arange(-wind,wind+1,1)]
                   for y in np.arange(-wind,wind+1,1)]
                  )*total + bkg)-bkg
        s = np.sum(tmp)
        
        return (np.sum(tmp*X)/s, np.sum(tmp*Y)/s)

    P = P_factor(ux,uy,sx,sy)

    dist_x,dist_y = zip(*[_run(P) for x in range(0,count)])

    dist_x,dist_y = np.array(dist_x) - ux,np.array(dist_y)-uy

    return dist_x,dist_y

#Copyright 2010 Thomas A Caswell
#tcaswell@uchicago.edu
#http://jfi.uchicago.edu/~tcaswell
#
#This program is free software; you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation; either version 3 of the License, or (at
#your option) any later version.
#
#This program is distributed in the hope that it will be useful, but
#WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
#General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program; if not, see <http://www.gnu.org/licenses>.

from util import cord_pairs
import scipy.odr as sodr
import numpy as np
import scipy.optimize
import general

def _trim_gofr(gofr,r0):
    """     Takes in gofr as a cord_pairs and r0 in
    units of 1/d, d as determined by the location
    of the first peak.

    Returns a cords_pair object with g(r>r0) 
    """
    # find the first peak
    d0 = gofr.x[np.argmax(gofr.y)]

    print "the Diameter is " + str(d0)

    # convert r0 to real units
    r0 *=d0

    
    
    # cut the data at r0
    r_arg = np.argmax([v for v in gofr.x if v<= r0])
    gofr_trim = cord_pairs(gofr.x[r_arg:],gofr.y[r_arg:])



    return gofr_trim

def fun_flipper(fun):
    def ffun(a,b):
        return fun(b,a)
    return ffun

def fun_decay_exp_inv(p,r):
    """Returns C/r exp(- r/a) cos(K(r)+phi_0) + m r + b
    evaluated at r.  p = (a,K,C,phi_0,b)"""
    return (p[2] / r) * np.exp(-r/p[0]  ) * np.cos(p[1]*r + p[3])+ p[4]

def fun_decay_exp_inv_dr(p,r):
    """ d(C/r exp(- r/a) cos(K(r)+phi_0) + m r + b)/dr
    evaluated at r.  p = (a,K,C,phi_0,b)"""
    return (np.exp(-(r/p[0]))* (-p[2]* (p[0] + r)* np.cos(p[3] + p[1]* r) - p[0] *p[2]* p[1]* r* np.sin(p[3] + p[1]* r)))/( p[0]* r**2)


def fun_decay_exp(p,r):
    """Returns C exp(- r/a) cos(K(r)+phi_0)
    evaluated at r.  p = (a,K,C,phi_0)"""
    return (p[2] ) * np.exp(-r/p[0]  ) * np.cos(p[1]*r + p[3]) + (r)*p[4] + p[5]


def fun_decay_inv(p,r):
    """Returns C/r cos(K(r)+phi_0)
    evaluated at r.  p = (K,C,phi_0)"""
    return (p[1] / r) * np.cos(p[0]*r + p[2])


def fit_gofr2(gofr,r0,func,p0=(1,7,2,0,1)):
    """
    Fits to 
    Takes in gofr as a cord_pairs and r0 in
    units of d, d as determined by the location
    of the first peak.

    Returns the tuple that is the argument for the function
    used for the fitting, the covariance matrix, and the
    sum of the residual squared
    """
    # trim the g(r) data
    gofr = _trim_gofr(gofr,r0)
    
    
    data = sodr.Data(gofr.x,gofr.y)
    model = sodr.Model(func)
    worker = sodr.ODR(data,model,p0)
    out = worker.run()
    out = worker.restart()
    return out
    


def fit_tmp_series_gofr2D(res,conn):
    '''Takes in a sample name and fits all of the 2D gofrs'''
    
    

    temps = []
    fits = []
    for r in res:
        gofr = general.get_gofr2D(r[0],conn)
        fits.append(fit_gofr2(gofr,2,fun_decay_exp,(2,7.35,1.5,0,0,1)))
        try:
            temps.append(float(r[2]))
        except (ValueError,TypeError ) :
            temps.append(25)

            
    return zip(fits,temps)

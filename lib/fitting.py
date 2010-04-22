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
import numpy as np
import scipy.optimize

def _trim_gofr(gofr,r0):
    """     Takes in gofr as a cord_pairs and r0 in
    units of 1/d, d as determined by the location
    of the first peak.

    Returns a cords_pair object with g(r>r0) 
    """
    # find the first peak
    d0 = gofr.x[np.argmax(gofr.y)]


    # convert r0 to real units
    r0 *=d0

    
    
    # cut the data at r0
    r_arg = np.argmax([v for v in gofr.x if v<= r0])
    gofr_trim = cord_pairs(gofr.x[r_arg:],gofr.y[r_arg:])



    return gofr_trim

def fun_decay_exp_inv(r,p):
    """Returns C/r exp(-a r) cos(K(r)+phi_0)
    evaluated at r.  p = (a,K,C,phi_0)"""
    return (p[2] / r) * np.exp(-p[0] * r ) * np.cos(p[1]*r + p[3])


def fun_decay_exp(r,p):
    """Returns C/r exp(-a r) cos(K(r)+phi_0)
    evaluated at r.  p = (a,K,C,phi_0)"""
    return (p[2] ) * np.exp(-p[0] * r ) * np.cos(p[1]*r + p[3])


def fun_decay_inv(r,p):
    """Returns C/r exp(-a r) cos(K(r)+phi_0)
    evaluated at r.  p = (K,C,phi_0)"""
    return (p[1] / r) * np.cos(p[0]*r + p[2])


def fit_gofr(gofr,r0,func,p0 =(1,12,.5,0)):
    """
    Fits to 
    Takes in gofr as a cord_pairs and r0 in
    units of 1/d, d as determined by the location
    of the first peak.

    Returns the tuple that is the argument for the function
    used for the fitting, the covariance matrix, and the
    sum of the residual squared
    """
    # trim the g(r) data
    gofr = _trim_gofr(gofr,r0)
    

    def local_fun(p):
        return func(gofr.x,p) - gofr.y

    (p_out,p_cov,id,m,flg) = scipy.optimize.leastsq(local_fun,p0,full_output=1)
    print p_cov
    return p_out,p_cov,sum(local_fun(p_out))**2
    # shove into the numpy fitting code

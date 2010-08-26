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
    d0 = gofr.x[np.argmax(gofr.y[15:])+15]

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

def fit_quad_to_peak(x,y):
    """Fits a quadratic to the data points handed in """
    def quad(B,x):
        return B[0] *(x -B[1]) ** 2 + B[2]

    beta = (0,np.mean(x),y[val_to_indx(x,np.mean(x))])

    data = sodr.Data(x,y)
    model = sodr.Model(quad)
    worker = sodr.ODR(data,model,beta)
    out = worker.run()


    
    ## plts.figure()
    ## plts.plot(x,y)
    ## plts.plot(x,quad(out.beta,x))
    ## plts.title(out.beta[1])
    return out

def try_fits(dset_key,conn):
    """Try's a couple of fits and plots the results """

    # get out the computation number
    res = conn.execute("select comp_key from comps where function = 'gofr3D' and dset_key = ?",(dset_key,)).fetchall()
    if not len(res) == 1:
        raise "die"

    # get gofr
    gofr = gen.get_gofr3D(res[0][0],conn)
    gofr = fitting._trim_gofr(gofr,.2)
    
    # fits
   
    (p_out1_2,cov1_2,err1_2) = fitting.fit_gofr(gofr,2,fitting.fun_decay_exp_inv,(2,7.35,1.5,0,0,0))
    (p_out2_2,cov2_2,err2_2) = fitting.fit_gofr(gofr,2,fitting.fun_decay_exp,(1.5,7.35,1.5,0,0,0))

    
    # plots
    

    
    # check interactive plotting and turn it off
    istatus = plt.isinteractive();
    print istatus
    if istatus:plt.ioff()

    leg_hands = []
    leg_str = []

    fig = plt.figure()
    ax = fig.add_axes([.1,.1,.8,.8])
    ax.hold(True)
    ax.grid(True)
    #ax.set_aspect('equal')
    leg_hands.append(ax.step(gofr.x,gofr.y-1))
    leg_str.append("g(r)")


    leg_hands.append(ax.step(gofr.x,fitting.fun_decay_exp_inv(p_out1_2,gofr.x)))
    leg_str.append("exp inv 2")


    leg_hands.append(ax.step(gofr.x,fitting.fun_decay_exp(p_out2_2,gofr.x)))
    leg_str.append("exp 2")



    print p_out1_2
    print "exp inv 2 err: " + str(err1_2)
    print p_out2_2
    print "exp 2 err: " + str(err2_2)



    ax.legend(leg_hands,leg_str)
    ax.set_title('g(r) fitting')
    ax.set_xlabel(r' r [$\mu$m]')
    ax.set_ylabel('g(r)')
    
            
    if istatus:
        print "displaying figure"
        plt.ion()
        plt.show()
    else:
        print "closing figure"
        plt.close(fig)

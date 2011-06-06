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
import plots

def fun_flipper(fun):
    def ffun(a,b):
        return fun(b,a)
    return ffun

###################
# functional forms#
###################

# for g(r)
def fun_decay_exp_inv_gen(r0):
    """Returns C/r exp(- (r-r0)/a) cos(K(r)+phi_0) + 1
    evaluated at r.  p = (a,K,C,phi_0)"""
    def fit_fun(p,r):
        return (p[2] / (r-p[3])) * np.exp(-(r-p[3])/p[0]  ) * np.sin(p[1]*(r- p[3])) +1
    return fit_fun

def fun_decay_exp_inv_dr_gen(r0):
    """ d(C/r exp(- (r-r_0)/a) cos(K(r)+phi_0) + b)/dr
    evaluated at r.  p = (a,K,C,phi_0)"""
    def ret_fun(p,r):
        return (p[2] /( r-p[3])) * np.exp(-(r-r0)/p[0]) * np.sin(p[1]*(r - p[3]))*\
               -1*(1/p[0] + 1/r + p[1]/ np.tan(p[1] * (r - p[3])))

    return ret_fun


def fun_decay_exp_gen(r0):
    """Returns C exp(- r/a) cos(K(r)+phi_0)
    evaluated at r.  p = (a,K,C,phi_0)"""
    def fit_fun(p,r):
        return (p[2] ) * np.exp(-(r-r0)/p[0]  ) * np.cos(p[1]*r + p[3]) + (r)*p[4] + p[5]
    return fit_fun

def fun_decay_inv(p,r):
    """Returns C/r cos(K(r)+phi_0)
    evaluated at r.  p = (K,C,phi_0)"""
    return (p[1]*r0 / r) * np.cos(p[0]*r + p[2])

# for van Hove

def fun_lorentzian(p,r):
    """Returns C/((r/a)^2 +1) evaluated at r, p = (a,C)"""
    return p[1] / ((r/p[0])**2 + 1)

def fun_lorentzian_p_gauss(p,r):
    """Returns C/((r/a)^2 +1) + C_2 exp(-(r/a_2)^2) evaluated at r, p = (a,C,C2)"""
    return p[1] / ((r/p[0])**2 + 1) + p[2] * np.exp(-((r/p[0])**2)/2)

def fun_lorentzian_t_gauss(p,r):
    """Returns C/((r/a)^2 +1) + C_2 exp(-(r/a_2)^2) evaluated at r, p = (a,C,C2)"""
    return p[2] / ((r/p[0])**2 + 1) * np.exp(-((r/p[1])**2)/2)

def fun_gauss(p,r):
    """Returns C/((r/a)^2 +1) + C_2 exp(-(r/a_2)^2) evaluated at r, p = (a,C,a2,C2)"""
    return  p[1] * np.exp(-((r/p[0])**2))

def fun_exp_p_gauss(p,r):
    """Returns C/((r/a)^2 +1) + C_2 exp(-(r/a_2)^2) evaluated at r, p = (a,C,a2,C2)"""
    return  p[1] * np.exp(-((r**2/p[0]))) +  p[3] * np.exp(-((np.abs(r)/p[2])))


def fun_exp_p_exp(p,r):
    """Returns C/((r/a)^2 +1) + C_2 exp(-(r/a_2)^2) evaluated at r, p = (a,C,a2,C2)"""
    return  p[1] * np.exp(-((np.abs(r)/p[0]))) +  p[3] * np.exp(-((np.abs(r)/p[2])))


def fun_exp_t_gauss(p,r):
    """Returns C/((r/a)^2 +1) + C_2 exp(-(r/a_2)^2) evaluated at r, p = (a,C,a2,C2)"""
    return  p[2] * np.exp(-((r**2/p[0]))-((np.abs(r)/p[1]))) 

def fun_gauss_gauss(p,r):
    """Returns C/((r/a)^2 +1) + C_2 exp(-(r/a_2)^2) evaluated at r, p = (a,C,a2,C2)"""
    return  p[1] * np.exp(-((r/p[0])**2)) + p[3] * np.exp(-((r/p[2])**2))


####################
# fitting functions#
####################
def fit_curve(x,y,p0,func):
    """Fits y = func(x|p) with the initial parameters of p0.  func
    must be of the form y = func(p,x).  uses scipy odr code  """
    
    data = sodr.Data(x,y)
    model = sodr.Model(func)
    worker = sodr.ODR(data,model,p0)
    out = worker.run()
    out = worker.restart()
    return out
    

def display_fit(x,y,p,func,fig=None):
    """Displays the raw data and the fit. fig is an plots.tac_figure
    object"""
    if fig is None:
        fig = plots.tac_figure('x','y','fitting')
        fig.plot(x,np.log(y),label='data')
        
    
    fig.plot(x,np.log(func(p,x)),'--x',label=func.__name__ +  '('+
             ','.join(['%.1e'%k for k in p])+ ')')
    
    return fig

#########################
# functions to be killed#
#########################

    


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

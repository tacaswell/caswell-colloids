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

import plots as plots
import h5py
import numpy as np
import pdb
import matplotlib
_C_to_K = 273.15
def _ff(n):
    """Formats frame names """
    return "frame%(#)06d"%{"#":n}

def _fd(str_,n):
    """ formats dset names"""
    return str_ + "_%(#)07d"%{"#":n}


def write_delta_T(comp,conn):
    """Goes back and adds the delta T md to msd results"""
    (fin,fout) = conn.execute("select fin,fout from comps where comp_key = ?",(comp,)).fetchone()
    (iden_comp,) = conn.execute("select comp_key from comps where fout = ? and function = 'Iden'",(fin,)).fetchone()
    Fin = h5py.File(fin,'r')
    # get delta T
    dt = Fin['frame000006'].attrs['dtime']
    Fin.close()
    del Fin
    # set it
    Fout = h5py.File(fout,'r+')
    Fout[_fd('msd_squared',comp)].attrs['dtime'] = dt
    Fout[_fd('mean_squared_disp',comp)].attrs['dtime'] = dt
    Fout[_fd('mean_disp',comp)].attrs['dtime'] = dt
    Fout.close()
    del Fout
    
def plot_msd_old(comp_key,conn,fig=None):
    if fig is None:
        fig = plots.tac_figure('t[s]','msd','msd')
        pass
    
    (fin,dset) = conn.execute("select fout,dset_key from comps where comp_key = ?",(comp_key,)).fetchone()
    (temp,) = conn.execute("select temp from dsets where key = ?",(dset,)).fetchone()
    
    
    
    Fin = h5py.File(fin,'r')
    g_name = _fd('mean_squared_disp',comp_key)
    msd = Fin[g_name]
    msd = Fin[g_name]['data'][:]
    dt = Fin[g_name].attrs['dtime']
    print 
    print 'the delta is ',  dt, 'for comp ' ,comp_key
    t = (np.arange(len(msd))+1)*dt

    cm = plots.color_mapper(27,33)
    
    fig.plot(t,msd,label=str(temp),color=cm.get_color(temp))
    Fin.close()
    del Fin
    return fig


def plot_msd(comp_key,conn,fig=None):
    if fig is None:
        fig = plots.tac_figure('t[s]','msd','msd')
        pass
    
    (fin,) = conn.execute("select fout from msd where comp_key = ?",(comp_key,)).fetchone()

    
    
    
    Fin = h5py.File(fin,'r')
    g_name = _fd('mean_squared_disp',comp_key)
    
    msd = Fin[g_name]['msd'][:]
    dt = Fin[g_name].attrs['dtime']
    temp = Fin[g_name].attrs['temperature']
    mtl = Fin[g_name].attrs['min_track_length']
    print 
    print 'the delta is ',  dt, 'for comp ' ,comp_key
    t = (np.arange(len(msd))+1)*dt

    cm = plots.color_mapper(27,33)
    
    fig.plot(t,msd,label='%(#).2fC, %(!)d'%{"!":mtl,"#":temp},color=cm.get_color(temp))
    Fin.close()
    del Fin
    return fig


def plot_msd_series(comp_key_lst,conn,sname=None):
    """ Takes in a list of msd computations and makes a nice plot"""

    tltstr = 'mean squared displacement'
    if not sname is None:
        tltstr = tltstr + " " + sname
        

        
    
    fig = plots.tac_figure('t[ms]',r'$\langle \Delta \rangle ^2$',tltstr,count=len(comp_key_lst),func=matplotlib.axes.Axes.plot)
    for c in comp_key_lst:
        plot_msd(c[0],conn,fig)

def series_fit(comp_key_lst,conn):
    """Takes a list of computation key  """

    
    
    
    
    

def fit_msd(comp_key,conn,make_plot = False):
    (fin,) = conn.execute("select fout from msd where comp_key = ?",(comp_key,)).fetchone()
    Fin = h5py.File(fin,'r')
    g_name = _fd('mean_squared_disp',comp_key)
    
    msd = Fin[g_name]['msd'][:]
    dt = Fin[g_name].attrs['dtime']
    temp = Fin[g_name].attrs['temperature']
    Fin.close()
    del Fin
    
    t = (np.arange(len(msd))+1)*dt
     
    (x,r,rn,s) = np.linalg.lstsq(np.transpose(np.array([t,np.ones(len(msd))])),msd)

    if(make_plot):
        fig = plots.tac_figure('t [ms]','msd [px^2]','msd and fit')
        fig.draw_line(t,msd,'x',label= 'msd')
        print x[1]
        fig.draw_line(t,t*x[0] + x[1],label= 'fit')

    scale = 6.45/60                     # u[m]/[pix]

    # lamda_x = sqrt(2Dt) (Einstein 1905)
    # msd = 2*lambda_x^2 = 4 D t (Einstein 1905)
    # the extra 2 is because I am looking at total displacement in 2D
    D = (x[0]/4)*(scale**2)*1000      # u[m]^2/sec
    r = _dilute_D_to_rad(D,temp)         # u[m]
    
    return (x[0],temp,r)

def _dilute_D_to_rad(D,T):
    """Does the computation to get the radius of a particle out of a
    diffusion constant assuming a dilute system in water.
    Assumes that D has units of (u[m])^2/[s].

    returns the radius in u[m]
    """
    nu = _nu(T)                         # m[Pa]*[s] 
    
    nu = nu/1000                        # [Pa] *[s]
    
    k = 1.3806504*(10**-23)             # [J]/[K]  = Pa * [m]^3 / [K]
    
    T = T + _C_to_K                     # [K]
    #D = (k*T )/6pi nu r
    print (T,nu*1000)
    D = D*(10**-12)                     # ([m])^2/[s]
    
    r = (k * T)/(6 * np.pi * nu * D)

    return r*(10**6)

def _nu(T):
    """ computes the viscosity of water.  Taken from:
    http://www.ddbst.com/en/online/Online_Calc_visc_Form.php"""
    A = -3.7188
    B = 578.919
    C = -137.546
    K_min = 273
    K_max = 373
    

    T = T + _C_to_K

    
    
    if T> K_max or T< K_min:
        raise Exception("out side of validity range")

    return np.exp(A + B/(C + (T)))

    
    

def delete_msd(comp_key,conn):
    """Deletes all record of a computation from both the data files and the data base """
    # get data file name
    (fname,) = conn.execute("select fout from comps where comp_key = ?",(comp_key,)).fetchone()
    # open file
    F = h5py.File(fname)
    # delete data sets
    del F[_fd('mean_squared_disp',comp_key)]
    del F[_fd('mean_disp',comp_key)]
    del F[_fd('msd_squared',comp_key)]
    # make sure changes take and clean up file
    F.flush()
    F.close()
    del F
    # delete from database
    conn.execute("delete from comps where comp_key = ?",(comp_key,))
    conn.execute("delete from msd_prams where comp_key = ?",(comp_key,))
    conn.commit()
                 

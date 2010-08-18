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


import sqlite3
import h5py
import pov 
import plots 
import util 
import numpy as np
from util import cord_pairs
import scipy.optimize as sopt
import scipy.odr as sodr
import bisect
import itertools
import matplotlib.pyplot as plts
import fitting



def make_gofr_2Dv3D(conn):
    keys = conn.execute("select key from dsets where dtype = 'z'").fetchall()
    for k in keys:
        plots.make_2dv3d_plot(k[0],conn,'figures/2Dv3D/%(#)02d.png'%{"#":k[0]})


def get_gofr3D(comp_num,conn):
    '''Takes in computation number and database connection and extracts the given g(r) and returns it as a
    numpy array'''

    gname='gofr3D'

    dset_names = ['bin_count', 'bin_edges']
    
    res = conn.execute("select fout from comps where comp_key = ?",(comp_num,)).fetchall()
    if not len(res) == 1:
        print len(res)
        raise util.dbase_error("error looking up computation")

    
    g = get_gofr_group(res[0][0],gname,comp_num)
    gofr = np.array(g[dset_names[0]])
    bins = np.array(g[dset_names[1]])
    return cord_pairs(bins,gofr)
    

def get_gofr2D(comp_num,conn):
    '''Takes in computation number and database connection and extracts the given g(r) and returns it as a
    numpy array'''

    dset_names = ['bin_count', 'bin_edges']
    
    res = conn.execute("select fout from comps where comp_key = ?",(comp_num,)).fetchall()
    if not len(res) == 1:
        print len(res)
        raise util.dbase_error("error looking up computation")
    
    print comp_num
    F = h5py.File(res[0][0],'r')
    print "gofr_%(#)07d"%{"#":comp_num}
    print "gofr_%(#)07d"%{"#":comp_num} in F.keys()
    g = F["gofr_%(#)07d"%{"#":comp_num}]
    gofr = np.array(g[dset_names[0]])
    bins = np.array(g[dset_names[1]])*6.45/60
    F.close()
    return cord_pairs(bins,gofr)

    
def find_peaks(gofr_pairs,expt_spc = None):
    '''Takes in a pairs object and returns a list of tuples with (x[indx],y[indx],indx)

    expt_spc is the expected spacing of the peaks

    Assumes that the largest peak comes first
    '''

    # find the largest (assumed to be first) peak
    
    
        
    if expt_spc is  None:
        # if no guess given, make a guess based on the location of the first peak
        pass


def find_peaks_fit(gp,dfun,p):
    """Looks for peaks based on hints from the derivative of the expected function and
    the fit parameters"""

    dfun = fitting.fun_flipper(dfun)
    wind = 20;

    lmax = []
    lmin = []
    diffs = []
    sz = []
    # find the first max
    indx = np.argmax(gp.y[15:])+15
    pfit = fit_peak(gp.x[indx-wind:indx+wind],gp.y[indx-wind:indx+wind])

    lmax.append((pfit.beta[1],pfit.beta[2]))
    diffs.append(gp.x[indx]-pfit.beta[1])
    # get expected spacing from
    e_spc = (2*np.pi)/p[1]

    print e_spc

    cur_state = np.sign(pfit.beta[0])
    sz.append(pfit.beta[0])
    #start at first peak + 1/4 e_spc
    cur_pos = pfit.beta[1] + e_spc/4
    # step over x range in e_spc/2 steps
    max_pos = np.max(gp.x)
    print max_pos
    print cur_pos
    while cur_pos < max_pos:
        # get zero crossing from dfun
        try:
            crit_p = sopt.brentq(dfun,cur_pos,cur_pos + e_spc/2,args=p)
            print 'diff from expected ' + str(crit_p - (cur_pos + e_spc/4))
        except ValueError:
            print "no zero found"
            break
        # convert value into indx
        indx = val_to_indx(gp.x,crit_p)
        # pick out window around that box
        pfit = fit_peak(gp.x[indx-wind:indx+wind],gp.y[indx-wind:indx+wind])
        

        # determine if max or min
        # add center/value to max or min output
        if np.abs(pfit.beta[0])<.05:
            print "peak too small"
            break

        print 'diff from crossing and min/max ' + str( crit_p - pfit.beta[1])
        diffs.append(crit_p - pfit.beta[1])
        sz.append(pfit.beta[0])
        
        if pfit.beta[0] >0:
            if cur_state != -1:
                print "found two peaks in a row"
                break
            lmin.append((pfit.beta[1],pfit.beta[2]))
            cur_state = np.sign(pfit.beta[0])
        elif pfit.beta[0]<0:
            if cur_state != 1:
                print "found two troughs in a row"
                break
            lmax.append((pfit.beta[1],pfit.beta[2]))
            cur_state = np.sign(pfit.beta[0])

            
        # increment cur_pos
        cur_pos = crit_p +  e_spc/4
        wind +=1
        pass
    
    return lmax,lmin,diffs,sz



def fit_peak(x,y):
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

def get_list_gofr (sname,conn,date = None,gtype = 'gofr'):
    """Takes in a  sample name and an optional string specifying
    the date of the computations to use"""

    print date
    print gtype
    if date is None:
        return conn.execute("select comps.comp_key,comps.fout,dsets.temp,comps.fin \
        from comps,dsets where comps.dset_key = dsets.key and comps.function=? and\
        dsets.sname = ?",(gtype,sname,)).fetchall()
    else:
        return conn.execute("select comps.comp_key,comps.fout,dsets.temp,comps.fin \
        from comps,dsets where comps.dset_key = dsets.key and comps.function=? and\
        dsets.sname = ? and dsets.dtype = 't' and comps.date = ?",(gtype,sname,date)).fetchall()

    

def get_gofr_by_plane_cps(comp_num,conn):
    fname = conn.execute("select fout from comps where comp_key = ?",(comp_num,)).fetchall()
    fname = fname[-1][0]
    F =  h5py.File(fname)
    g = F['gofr_by_plane_%(#)07d'%{'#':comp_num}]

    g_l = [cord_pairs(g[c]['bin_edges'][:],g[c]['bin_count'][:]) for c in g]

    
    del g
    F.close()
    return g_l

def get_gofr_by_plane_tmp(comp_num,conn):
    (fname,) = conn.execute("select fout from comps where comp_key = ?",(comp_num,)).fetchone()
    
    F =  h5py.File(fname)
    g = F['gofr_by_plane_%(#)07d'%{'#':comp_num}]

    temps = [g[c].attrs['temperature'] for c in g]

    
    del g
    F.close()
    return temps

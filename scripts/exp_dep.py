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
import lib.pov 
import lib.plots as lplts
import lib.util 
import lib.general as gen
import trackpy.cpp_wrapper
import numpy as np

def shift_mod(F,fr,comp_num,dr,bin_edges):
    tmp_data = F[fr][dr + "_shift_%(#)07d"%{'#':comp_num}][:]
    (tmp_hist, junk) = np.histogram(np.mod(tmp_data,1),bin_edges,new=True)
    return tmp_hist

def x_shift_mod(F,fr,comp_num,bin_edges):
    return shift_mod(F,fr,comp_num,'x',bin_edges)

def y_shift_mod(F,fr,comp_num,bin_edges):
    return shift_mod(F,fr,comp_num,'y',bin_edges)

def xy_shift_mod(F,fr,comp_num,bin_edges):
    return x_shift_mod(F,fr,comp_num,bin_edges) + x_shift_mod(F,fr,comp_num,bin_edges)

def _shift(F,fr,comp_num,dr,bin_edges):
    tmp_data = F[fr][dr + "_shift_%(#)07d"%{'#':comp_num}][:]
    (tmp_hist, junk) = np.histogram(tmp_data,bin_edges,new=True)
    return tmp_hist

def x_shift(F,fr,comp_num,bin_edges):
    return _shift(F,fr,comp_num,'x',bin_edges)

def y_shift(F,fr,comp_num,bin_edges):
    return _shift(F,fr,comp_num,'y',bin_edges)

def xy_shift(F,fr,comp_num,bin_edges):
    return _shift(F,fr,comp_num,'x',bin_edges) + _shift(F,fr,comp_num,'y',bin_edges)

def shift_mag(F,fr,comp_num,bin_edges):
    tmp_data = np.sqrt(F[fr]["x_shift_%(#)07d"%{'#':comp_num}][:]**2 +
                       F[fr]["y_shift_%(#)07d"%{'#':comp_num}][:]**2)
    (tmp_hist, junk) = np.histogram(tmp_data,bin_edges,new=True)
    return tmp_hist

def shift_mag_mod(F,fr,comp_num,bin_edges):
    tmp_data = np.sqrt(F[fr]["x_shift_%(#)07d"%{'#':comp_num}][:]**2 +
                       F[fr]["y_shift_%(#)07d"%{'#':comp_num}][:]**2)
    (tmp_hist, junk) = np.histogram(np.mod(tmp_data,1),bin_edges,new=True)
    return tmp_hist

def hist_shifts(key,conn,fun,range_bn = None, fig = None):
    """makes histograms of the shift in the x and y directions """

    # get file name/comp_num
    (comp_num,fname) = conn.execute("select comp_key,fout from comps\
    where dset_key = ? and function = 'Iden'",(key,)).fetchone()

    # open file/group

    F = h5py.File(fname,'r')

    nbins = 100

    if range_bn is None:
        bin_edges = np.linspace(-2,2,nbins + 1)
    else:
        bin_edges = np.linspace(*(range_bn +  (nbins + 1,)))
    bin_counts = np.zeros(nbins)
    # extract the relevant data
    for fr in F:
        if fr == 'parameters':
            continue
        bin_counts += fun(F,fr,comp_num,bin_edges)
        
    # plot
    istatus = lplts.non_i_plot_start()
    if fig is None:
        (fig,ax) = lplts.set_up_plot()
    else:
        ax = fig.get_axes()[0]
        
    sh = ax.step(bin_edges[:-1],bin_counts/np.sum(bin_counts))

    if ax.get_legend() is None:
        print 'attempt to set leg'
        ax.legend([sh],[F.attrs['Exposure']],loc = 3)
    else:
        #leg =aff  ax.get_legend()
        pass
    lplts.non_i_plot_stop(istatus)
    # clean up
    F.close()
    del F

    return fig

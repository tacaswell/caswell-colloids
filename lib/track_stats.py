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
from __future__ import division

import h5py

import matplotlib

import plots as plt


import plots
import numpy as np
from general import fd
Fig = plots.Figure
color_mapper = plots.color_mapper



def _fd(str_,n):
    """ formats dset names"""
    return str_ + "_%(#)07d"%{"#":n}


def make_disp_sq_plots(comp_key,conn):
    """Takes in a comp_key for a track_stat computation and plots the
    square displacement histograms"""

    (fin,) = conn.execute("select fout from comps where function = 'track_stats' and comp_key = ?",
                          (comp_key,)).fetchone()

    F = h5py.File(fin,'r')
    
    g = F[_fd('disp_sq_hist',comp_key)]
    cmap = color_mapper(0,len(g))

    (fig,ax) = plots.set_up_plot()
    istatus = plots.non_i_plot_start();
    for s in g:
        step = int(s[-7:])
        val = g[s]['bin_value'][:]
        
        ax.semilogy(g[s]['bin_edges'],val,color = cmap.get_color(step))
    F.close()

    (iden_fun,dset_key ) = conn.execute("select function,dset_key from comps where " +
                                "comp_key in " +
                                "(select iden_key from trk_stat_prams where "+
                                "comp_key = ?)",(comp_key,)).fetchone()

    ax.set_title("dset: " + str(dset_key) + " " + iden_fun)

    plots.non_i_plot_stop(istatus)

def track_len_hists(comp_lst,conn):
    """makes histograms of the track lengths for the tracking computations in comp_lst """

    res = [_extract_track_lengths(c,conn) for c in comp_lst]
    hist_res = [np.histogram(r[0],bins=np.max(r[0])) + (r[1],) for r in res]

    tmps = [d[1] for d in res]
    cm = plt.color_mapper(np.min(tmps),np.max(tmps))
    
    
    fig = Fig('length','counts','Track length histogram',func = matplotlib.axes.Axes.loglog)
    for hr in hist_res:
        temp = hr[2]
        fig.plot(hr[1],hr[0],label='%(#)0.1f C'%{'#':temp},color=cm.get_color(temp))

def _extract_track_lengths(track_key,conn):
    """Extracts the array of track lengths
    
    Assumes track_key is a length 1 tuple

    Assumes that the tracking and iden data are in the same file
    """
    print track_key
    
    (fname,iden_key,track_key) = conn.execute("select fout,iden_key,comp_key from tracking where comp_key = ?",
                                    track_key).fetchone()
    
    F = h5py.File(fname,'r')
    len_vec = F[fd('tracking',track_key)]['length'][:]
    
    temp = 0
    fr_count = 0
    for g in F.keys():
        if g[0:5] == 'frame':
            temp += F[g].attrs['temperature']
            fr_count += 1

    
    F.close()
    del F
    return len_vec, temp/fr_count

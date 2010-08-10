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

import matplotlib.pyplot as plt

import plots
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

    
def make_trk_len_plots(comp_key,conn,fig,scale=1):
    """ takes in an comp_key for a track_stat computation and plot the
    track length histogram"""

    
    (fin,) = conn.execute("select fout from comps where function = 'track_stats' and comp_key = ?",
                          (comp_key,)).fetchone()

    
    (iden_fun,dset_key ) = conn.execute("select function,dset_key from comps where " +
                                "comp_key in " +
                                "(select iden_key from trk_stat_prams where "+
                                "comp_key = ?)",(comp_key,)).fetchone()

    

    F = h5py.File(fin,'r')
    
    g = F[_fd('track_length',comp_key)]
    

    

    val = g['bin_value'][:]
    fig.plot(g['bin_edges'][:]*scale,val,label=iden_fun)
    print 'over: ' + str(g.attrs['over_count'])
    F.close()



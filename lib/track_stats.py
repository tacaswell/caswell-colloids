#Copyright 2010,2011 Thomas A Caswell
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
import matplotlib.pyplot as plt
import numpy as np

import plots
from general import fd,ff


Fig = plots.tac_figure
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



def track_len_hists(comp_lst,conn,ax=None,*args,**kwargs):
    """makes histograms of the track lengths for the tracking computations in comp_lst """

    res = [_extract_track_lengths(c,conn) for c in comp_lst]
    hist_res = [np.histogram(r[0],bins=np.max(r[0])) + (r[1],r[2]) for r in res]

    tmps = [d[1] for d in res]
    cm = plots.color_mapper(np.min(tmps),np.max(tmps))
    
    if ax is None:
        ax = plots.set_up_axis('n*dtime [s]','cumsum $n P(n)$','',func = matplotlib.axes.Axes.plot)
        
    for hr in hist_res:
        temp = hr[2]
        ax.plot(np.cumsum((np.diff(hr[1]))+hr[1][0])*hr[3]/1000,np.cumsum(hr[0]*(hr[1][:-1]))/np.sum(hr[0]*(hr[1][:-1]))
                      ,label='%(#)0.1f C'%{'#':temp},color=cm.get_color(temp),*args,**kwargs)

class Hist2D_accumlator (object):
    def __init__(self,data_dims,hist_dims):
        self.data_dims = data_dims
        self.hist_dims = hist_dims
        self.data = np.zeros(self.hist_dims)
    def add_point(self,point):
        self.data[tuple(np.floor((point/self.data_dims)*self.hist_dims).tolist())]+=1
        
    pass


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
    dtime = 0
    fr_count = 0
    for g in F.keys():
        if g[0:5] == 'frame':
            temp += F[g].attrs['temperature']
            dtime += F[g].attrs['dtime']
            fr_count += 1

    
    F.close()
    del F
    return len_vec, temp/fr_count, dtime/fr_count


def new_track_density(track_key,hist_dims,conn):
    """This makes a 2D histogram of where new tracks are found in the
    sample.  This is useful to check if there are regions of the
    sample that are tracking badly"""

    # extract all of the md needed to look up the data

    (fname,iden_key,track_key,dset_key) = conn.execute("select fout,iden_key,comp_key,dset_key from tracking where comp_key = ?",
                                    track_key).fetchone()
    print fname
    F = h5py.File(fname,'r')
    print F.keys()[:5]
    try:
        start_plane = F[fd('tracking',track_key)]['start_plane'][:]
        start_part = F[fd('tracking',track_key)]['start_particle'][:]

        print len(start_plane)
                
        # figure out the right size to make the array
        dims = F.attrs['dims']
        print dims
        # make data collection object
        hist2D_ac = Hist2D_accumlator(dims,hist_dims)
        # loop over the heads of track index and hash result
        cur_plane = None
        cur_x = None
        cur_y = None
        temp = 0
        fr_count = 0
        for plane,part in zip(start_plane,start_part):
            if not plane == cur_plane:
                cur_plane = plane
                cp =  F[ff(cur_plane)]
                cur_x = cp[fd('x',iden_key)]
                cur_y = cp[fd('y',iden_key)]
                temp += cp.attrs['temperature']
                fr_count += 1

            hist2D_ac.add_point(
                (cur_x[part],
                 cur_y[part])
                )
            pass
    except ValueError,er:
        print ff(cur_plane)
        
        
    finally:
        F.close()
        del F

    f = plt.figure()
    ax = f.add_axes([.1,.1,.8,.8])
    c = ax.imshow(np.flipud(hist2D_ac.data.T),interpolation='nearest')
    plt.colorbar(c)
    ax.set_title('%.2f C '%(temp/fr_count) + str(dset_key))
    return hist2D_ac.data

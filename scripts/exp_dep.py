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

import sqlite3
import h5py
import lib.pov 
import lib.plots as lplts
import lib.util 
import lib.general as gen
import appwrappy.cpp_wrapper
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

def shift_mod(F,fr,comp_num,dr,bin_edges):
    tmp_data = F[fr][dr + "_%(#)07d"%{'#':comp_num}][:]

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

def hist_shifts_frame(key,frame,conn,fun,range_bn = None, fig = None):
    
    # get file name/comp_num
    (comp_num,fname) = conn.execute("select comp_key,fout from comps\
    where dset_key = ? and function = 'Iden'",(key,)).fetchone()

    # open file/group

    F = h5py.File(fname,'r')

    nbins = 1000
    
    if range_bn is None:
        bin_edges = np.linspace(-2,2,nbins + 1)
    else:
        bin_edges = np.linspace(*(range_bn +  (nbins + 1,)))


    fr = 'frame%(#)06d'%{'#':frame}
    bin_counts = fun(F,fr,comp_num,bin_edges)


        
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


def hist_shifts(key,conn,fun,range_bn = None, fig = None):
    """makes histograms of the shift in the x and y directions """

    # get file name/comp_num
    (comp_num,fname) = conn.execute("select comp_key,fout from comps\
    where comp_key = ? and function = 'Iden'",(key,)).fetchone()

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


def I_hist(F,fr,comp_num,bins):
    tmp_data = F[fr]["intensity_%(#)07d"%{'#':comp_num}][:]
    (tmp_hist, tmp_bins) = np.histogram(tmp_data,bins,new=True)
    return tmp_hist,tmp_bins


def e_hist(F,fr,comp_num,bins):
    
    tmp_data = F[fr]["eccentricity_%(#)07d"%{'#':comp_num}][:]
    (tmp_hist, tmp_bins) = np.histogram(tmp_data,bins,new=True)
    return tmp_hist,tmp_bins


def comp_hists(key,frames,fun,conn):
    """Makes histograms plots comparing the values from various frames
    in the same data set"""

    
    # get file name/comp_num
    (comp_num,fname) = conn.execute("select comp_key,fout from comps\
    where dset_key = ? and function = 'Iden'",(key,)).fetchone()

    # open file/group

    F = h5py.File(fname,'r')

    res = [fun(F,'frame%(#)06d'%{'#':fr},comp_num,100) for fr in frames]
    
    
    istatus = lplts.non_i_plot_start()
    (fig,ax) = lplts.set_up_plot()
    hands = [ax.step(r[1][:-1],r[0]) for r in res]
    ax.legend(hands,[str(f) for f in frames])
    lplts.non_i_plot_stop(istatus)

    F.close()
    del F

    return res



def comp_hists_avg(key,frames_a,frames_b,fun,conn):
    """Makes histograms plots comparing the values from various frames
    in the same data set

    not written

    """

    
    # get file name/comp_num
    (comp_num,fname) = conn.execute("select comp_key,fout from comps\
    where dset_key = ? and function = 'Iden'",(key,)).fetchone()

    # open file/group

    F = h5py.File(fname,'r')
    hists_a = np.zeros()
    for fr in frames_a:
        tmp_h,tmp_b = fun(F,'frame%(#)06d'%{'#':fr},comp_num)
        hists_a += tmp_h
    
    
    istatus = lplts.non_i_plot_start()
    (fig,ax) = lplts.set_up_plot()
    hands = [ax.step(r[1][:-1],r[0]) for r in res]
    ax.legend(hands,[str(f) for f in frames])
    lplts.non_i_plot_stop(istatus)

    F.close()
    del F

    return res

def I_v_shift_scatter(comp_key,conn):
    """makes a scatter plot of intensity vs shift (which one determined by func """
        
    # get file name/comp_num
    (comp_num,fname) = conn.execute("select comp_key,fout from comps\
    where comp_key = ? and function = 'Iden'",(comp_key,)).fetchone()
    (hwhm,p_rad) = conn.execute("select hwhm,p_rad from Iden_prams \
    where comp_key = ?",(comp_num,)).fetchone()
    frame = 15
    fr = 'frame%(#)06d'%{'#':frame}
    # open file/group

    F = h5py.File(fname,'r')


    shift_data = np.sqrt(F[fr]["x_shift_%(#)07d"%{'#':comp_num}][:]**2 +
                       F[fr]["y_shift_%(#)07d"%{'#':comp_num}][:]**2)
    
    I_data = F[fr]["eccentricity_%(#)07d"%{'#':comp_num}][:]
    F.close()

    fig = lplts.Figure('E','shift mag','I v shift_mag',func = matplotlib.axes.Axes.scatter)
    fig.plot(I_data,shift_data,lab = '(' + str(hwhm) + ',' + str(p_rad)+ ')')

    ## stolen from numpy manual
    
    H, (edges) = np.histogramdd( (I_data/np.max(I_data),shift_data), bins=(100, 200))
    print H.shape

    # We can now use the Matplotlib to visualize this 2-dimensional histogram:

    
    import matplotlib.pyplot as plt
    plt.figure()
    plt.imshow(np.flipud(H),extent=[0,2,0,1])
    # <matplotlib.image.AxesImage object at ...>
    plt.show()
    plt.draw()
    ##

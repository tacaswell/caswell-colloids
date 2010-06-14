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

def hist_shifts(key,conn,fig = None):
    """makes histograms of the shift in the x and y directions """

    # get file name/comp_num
    (comp_num,fname) = conn.execute("select comp_key,fout from comps\
    where dset_key = ? and function = 'Iden'",(key,)).fetchone()

    # open file/group

    F = h5py.File(fname,'r')

    nbins = 100
    
    bin_edges = np.linspace(-2,2,nbins + 1)
    bin_counts = np.zeros(nbins)
    # extract the relevant data
    for fr in F:
        if fr == 'parameters':
            continue
        tmp_data = F[fr]["x_shift_%(#)07d"%{'#':comp_num}][:]
        (tmp_hist, junk) = np.histogram(np.mod(tmp_data,1),bin_edges,new=True)
        bin_counts += tmp_hist
        
    # plot
    istatus = lplts.non_i_plot_start()
    if fig is None:
        (fig,ax) = lplts.set_up_plot()
    else:
        ax = fig.get_axes()[0]
        
    sh = ax.step(bin_edges[:-1],bin_counts)

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

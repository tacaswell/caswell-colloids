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

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm

import numpy as np
def figure_out_grid(tot):
    """Computes the 'ideal' grid dimensions for the total number of
    graphs handed in"""
    guess = np.round(np.sqrt(tot))

    dims = [guess,guess]

    flag = True

    while(flag):
        if dims[0]*dims[1] < tot:
            dims[0] += 1
            flag = True
        elif dims[0]*dims[1] >tot and dims[0]*dims[1] <(tot-dims[1]):
            dims[0] -= 1
            flag = True
        else:
            flag = False
    return tuple(dims)
        

def disp_all_cmap():

    an = cm.cmapnames
    fig = plt.figure()
    fig.suptitle('all color maps')

    dims = figure_out_grid(len(an))
    plt_count = 1
    for name in an:
        sp_arg = dims +(plt_count,)
        ax = fig.add_subplot(*sp_arg)

        ax.imshow([range(0,30)]*6,interpolation='nearest',cmap=cm.get_cmap(name))
        ax.set_title(name)
        
        plt_count += 1
    plt.show()
    plt.draw()

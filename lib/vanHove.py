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

import plots as plt
import h5py
import numpy as np
import pdb
import matplotlib
import matplotlib.pyplot as mplt
from general import fd

def read_vanHove(comp_key,conn):
    '''Takes in a comp_key to a van Hove computation and returns a
    sensible data structure'''

    (fin,) = conn.execute("select fout from comps where comp_key = ?",(comp_key,)).fetchone()

    Fin = h5py.File(fin,'r')

    g = Fin[fd('vanHove',comp_key)]
    vanHove = [(g[fd('step',j+1)]['x']['disp_count'][:],g[fd('step',j+1)]['x']['disp_edges'][:])   for j in range(0,g.attrs['max_step'])]
    del g
    Fin.close()
    del Fin
    return vanHove
               
               
def plot_vanHove(comp_lst,time_step,conn):
    ''' '''

    cm = plt.color_mapper(27,33)
    fig =  plt.Figure(r'$\Delta$','log(N)','van Hove',func=matplotlib.axes.Axes.step)
  
    
    for c in comp_lst:
        (fin,) = conn.execute("select fout from comps where comp_key = ?",c).fetchone()
        Fin = h5py.File(fin,'r')
        g = Fin[fd('vanHove',c[0])]
        
        temp = g.attrs['temperature']
        dtime = g.attrs['dtime']
        j = time_step//dtime

        count = g[fd('step',j)]['x']['disp_count'][:]
        edges = g[fd('step',j)]['x']['disp_edges'][:]
        fig.plot(edges,numpy.log(count+1),label='%.2f'%temp,color=cm.get_color(temp))
        
        del g
        Fin.close()
        del Fin

    

        
               
def plot_vanHove_pn(comp_lst,time_step,conn):
    ''' '''

    cm = plt.color_mapper(27,33)
    fig =  plt.Figure(r'$\Delta$','P(N)','van Hove',func=matplotlib.axes.Axes.step)
  
    
    for c in comp_lst:
        (fin,) = conn.execute("select fout from comps where comp_key = ?",c).fetchone()
        Fin = h5py.File(fin,'r')
        g = Fin[fd('vanHove',c[0])]
        
        temp = g.attrs['temperature']
        dtime = g.attrs['dtime']
        j = int(time_step/dtime)

        count = g[fd('step',j)]['x']['disp_count'][:]
        edges = g[fd('step',j)]['x']['disp_edges'][:]
        fig.plot(edges,count/numpy.sum(count),label='%.2f'%temp,color=cm.get_color(temp))
        
        del g
        Fin.close()
        del Fin

    

        
             
def plot_vanHove_sp(comp_lst,time_step,conn):
    ''' '''

    cm = plt.color_mapper(27,33)

    fig = mplt.figure()
    fig.suptitle(r'van Hove dist $\tau = %.2f \rm{s}$'% (time_step/1000))
    dims = figure_out_grid(len(comp_lst))
    
    plt_count = 1
    
    for c in comp_lst:
        (fin,) = conn.execute("select fout from comps where comp_key = ?",c).fetchone()
        Fin = h5py.File(fin,'r')
        g = Fin[fd('vanHove',c[0])]
        
        temp = g.attrs['temperature']
        dtime = g.attrs['dtime']
        j = int(time_step/dtime)
        print j
        count = g[fd('step',j)]['y']['disp_count'][:]
        edges = g[fd('step',j)]['y']['disp_edges'][:]
        sp_arg = dims +(plt_count,)
        ax = fig.add_subplot(*sp_arg)
        ax.grid(True)
        
        ax.set_ylabel(r'$\log{N}$')
        ax.step(edges,np.log((count+1)/np.sum(count)),color=cm.get_color(temp))
        ax.set_title('%.2f'%temp)
        del g
        Fin.close()
        del Fin
        plt_count += 1

    mplt.draw()

        
def figure_out_grid(tot):

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
        
    

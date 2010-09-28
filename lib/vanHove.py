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
import fitting
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

def extract_vanHove(c,conn,time_step,count_cut,wind=1):
    """
    Function to encapsulate extracting a vanHove at a given time step.
    The actual time difference will always be less than requested due
    to rounding and finite sampling rates.

    This function takes care of looking up the correct files and dealing with
    the hdf layer.
    """
    (fin,) = conn.execute("select fout from comps where comp_key = ?",c).fetchone()
    Fin = h5py.File(fin,'r')
    g = Fin[fd('vanHove',c[0])]

    temp = g.attrs['temperature']
    dtime = g.attrs['dtime']
    j = int(time_step/dtime)
    print j

    (edges,count,x_lim) = _extract_vanHove(g,j,count_cut,wind)
    
    del g
    Fin.close()
    del Fin
        
    

    return edges,count,temp,dtime,x_lim

def _extract_vanHove(g,j,count_cut,wind):
    """
    Does the actual work of extracting the proper data set from an open
    hdf group object.

    The x and y data are merged, further binned, and trimmed

    g hdf group object
    j the time step to extract
    count_cut the minimum number of counts to have for the bin to be kept
    wind how many bins to sum togethur
    """
    count = g[fd('step',j)]['y']['disp_count'][:]
    # get better statistics
    count += g[fd('step',j)]['x']['disp_count'][:]
    #        count += count[::-1]

    edges = g[fd('step',j)]['y']['disp_edges'][:]
    
    edges = np.array([np.mean(edges[j:j+wind]) for j in range(1,len(count)-wind,wind)])
    count = np.array([np.sum(count[j:j+wind]) for j in range(1,len(count)-wind,wind)])

    x_lim = (np.min(edges),np.max(edges))

    edges = edges[count>count_cut]
    count = count[count>count_cut]

    return edges,count,x_lim

def _alpha2(edges,count):
    """
    Computes the \alpha_2 value for the distribution
    """
    # normalize count
    tmp_count = count/np.sum(count)
    return (
        np.sum((edges**4) *tmp_count )/
            (
            (5/3)*
                np.sum((edges**2)*tmp_count)
                **2)
        ) -1


def compute_alpha(comp,conn,wind ,min_count ):
    """
    Takes in a computation number, opens the hdf file, extracts all the vanHove distributions
    an computes alpha2 for them.

    Returns a list of tuples(wait_time,alpha2) and the temperature
    """
    
    (fin,) = conn.execute("select fout from comps where comp_key = ?",comp).fetchone()
    (max_step,) = conn.execute("select max_step from vanHove_prams where comp_key = ?",comp).fetchone()
    Fin = h5py.File(fin,'r')
    g = Fin[fd('vanHove',comp[0])]

    temp = g.attrs['temperature']
    dtime = g.attrs['dtime']
    dtime = 1
    a = []
    for j in range(0,max_step):
        (edges,count,junk) = _extract_vanHove(g,j+1,min_count,wind)
        a_tmp = _alpha2(edges,count)
        if not np.isnan(a_tmp):
            a.append((dtime*(j+1),a_tmp))
    
    del g
    Fin.close()
    del Fin

    return a,temp
def plot_alpha2(comp_list,conn,wind,min_count):
    """
    Takes in a comp_list and plots the alpha2().

    returns a list of the compute_alpha results
    """
    a_lst = [compute_alpha(c,conn,wind,min_count) for c in comp_list]
    _plot_alpha2(a_lst)

    return a_lst

def _plot_alpha2(a_list):
    """
    plots alpha2 from the list of outputs of compute_alpha handed in
    """
    cm = plt.color_mapper(27,33)
    fig =  plt.Figure(r'$\Delta$','P(N)','van Hove',func=matplotlib.axes.Axes.step)

    for a,temp in a_list:
        dt,a2 = zip(*a)
        fig.plot(dt,a2,
                 label='%.2f'%temp,
                 color=cm.get_color(temp)
                 )
    


def plot_vanHove_dt(comp,conn,start,step_size,steps):
    """
    Plots a grid array of the Von Hove distributions of a single computation
    at different times in steps even steps from the first to last.
    """
        
    (fin,) = conn.execute("select fout from comps where comp_key = ?",comp).fetchone()
    (max_step,) = conn.execute("select max_step from vanHove_prams where comp_key = ?",comp).fetchone()
    Fin = h5py.File(fin,'r')
    g = Fin[fd('vanHove',comp[0])]

    temp = g.attrs['temperature']
    dtime = g.attrs['dtime']


    #  istatus = plt.non_i_plot_start()
    
    fig = mplt.figure()
    fig.suptitle(r'van Hove dist temp: %.2f dtime: %d'% (temp,dtime))
    dims = figure_out_grid(steps)
    
    plt_count = 1
    outs = []
    tmps = []
    for j in range(start,start+step_size*steps, step_size):
        (edges,count,x_lim) = _extract_vanHove(g,j+1,1,5)
        if len(count) < 50:
            plt_count += 1
            continue
        #count = count/np.sum(count)
        
        sp_arg = dims +(plt_count,)
        ax = fig.add_subplot(*sp_arg)
        ax.grid(True)

        
        alpha = _alpha2(edges,count)
        
        ax.set_ylabel(r'$\log{P(N)}$')
        ax.step(edges,np.log((count/np.sum(count))),lw=2)
        ax.set_title(r'$\alpha_2 = %.2f$'%alpha + ' j:%d '%j  )
        ax.set_xlim(x_lim)
        plt_count += 1

    mplt.draw()

    # plt.non_i_plot_start(istatus)

    del g
    Fin.close()
    del Fin

def plot_vanHove_sp(comp_lst,time_step,conn,wind =1,func = fitting.fun_exp_p_gauss):
    '''Plots a grid array of the Von Hove functions at the time step
    given for all of the comps in the comp_lst'''

    cm = plt.color_mapper(27,33)
    
    fig = mplt.figure()
    fig.suptitle(r'van Hove dist $\tau = %.2f \rm{s}$'% (time_step/1000))
    dims = figure_out_grid(len(comp_lst))
    
    plt_count = 1
    outs = []
    tmps = []
    for c in comp_lst:
        (edges,count,temp,dtime,x_lim) = extract_vanHove(c,conn,time_step,50,wind)
        if len(count) < 50:
            break
        #count = count/np.sum(count)
        
        sp_arg = dims +(plt_count,)
        ax = fig.add_subplot(*sp_arg)
        ax.grid(True)


        
        
        p0 = (2,np.max(count)*.99,1,np.max(count)*.01)
    
    
    
        out = fitting.fit_curve(edges,count,p0,func)
        out = fitting.fit_curve(edges,count,out.beta,func)
        outs.append(out)
        tmps.append(temp)
        ax.plot(edges,np.log(func(out.beta,edges)),'--k')
        ax.set_ylabel(r'$\log{N}$')
        ax.step(edges,np.log((count)),color=cm.get_color(temp))
        ax.set_title('%.2f'%temp + ' dtime:%d ms'%dtime)
        ax.set_xlim(x_lim)
        plt_count += 1

    mplt.draw()

    return outs,tmps

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
        
    
def fit_lorentzian(comp_key,p0,time_step,conn,func = fitting.fun_lorentzian,fig=None,wind=3):
    """Fits a Lorentzian to the van Hove function handed in """

    (fin,) = conn.execute("select fout from comps where comp_key = ?",(comp_key,)).fetchone()
    Fin = h5py.File(fin,'r')
    g = Fin[fd('vanHove',comp_key)]

    temp = g.attrs['temperature']
    dtime = g.attrs['dtime']
    j = int(time_step/dtime)
    count = g[fd('step',j)]['y']['disp_count'][:]
    count += g[fd('step',j)]['x']['disp_count'][:]
    edges = g[fd('step',j)]['y']['disp_edges'][:]

    
    edges = np.array([np.mean(edges[j:j+wind]) for j in range(1,len(count)-wind,wind)])
    count = np.array([np.sum(count[j:j+wind]) for j in range(1,len(count)-wind,wind)])
    
        
    edges = edges[count>30]
    count = count[count>30]
    

    out = fitting.fit_curve(edges,count,p0,func)
    fig = fitting.display_fit(edges,count,out.beta,func,fig)
    print out.beta
    return out,fig

def find_knee(x,y):
    """Takes in a curve. fits a line to the top and bottom parts and
    tries to find the inflection point"""

    # find ranges
    if len(x) != len(y):
        raise Exception("bad data")
    tot_len = len(x)
    
    
    
    # fit strait lines to both

    # find intercept
    knee_r = (f_top.beta[1] - f_bottom.beta[1])/(-f_top.beta[0] + f_bottom.beta[0])
    

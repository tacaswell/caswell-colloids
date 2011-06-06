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
import matplotlib.pyplot as mplt
import numpy as np

import fitting
import general as gen
import plots as plots
from general import fd

import T_convert as ltc


def read_vanHove(comp_key,conn):
    '''Takes in a comp_key to a van Hove computation and returns a
    sensible data structure'''

    (fin,) = conn.execute("select fout from vanHove where comp_key = ?",(comp_key,)).fetchone()

    Fin = h5py.File(fin,'r')

    g = Fin[fd('vanHove',comp_key)]
    vanHove = [(g[fd('step',j+1)]['x']['disp_count'][:],g[fd('step',j+1)]['x']['disp_edges'][:])   for j in range(0,g.attrs['max_step'])]
    del g
    Fin.close()
    del Fin
    return vanHove
               
               
def plot_vanHove(comp_lst,time_step,conn):
    ''' '''

    cm = plots.color_mapper(27,33)
    fig =  plots.tac_figure(r'$\Delta$','log(N)','van Hove',func=matplotlib.axes.Axes.step)
  
    
    for c in comp_lst:
        (fin,) = conn.execute("select fout from vanHove where comp_key = ?",c).fetchone()
        Fin = h5py.File(fin,'r')
        g = Fin[fd('vanHove',c[0])]
        
        temp = g.attrs['temperature']
        dtime = g.attrs['dtime']
        j = time_step//dtime

        count = g[fd('step',j)]['x']['disp_count'][:]
        edges = g[fd('step',j)]['x']['disp_edges'][:]
        fig.draw_line(edges,numpy.log(count+1),label='%.2f'%temp,color=cm.get_color(temp))
        
        del g
        Fin.close()
        del Fin

    

        
               
def plot_vanHove_pn(comp_lst,time_step,conn):
    ''' '''
    
    cm = plots.color_mapper(27,33)
    fig =  plots.tac_figure(r'$\Delta$','P(N)','van Hove',func=matplotlib.axes.Axes.step)
  
    
    for c in comp_lst:
        (fin,) = conn.execute("select fout from comps where comp_key = ?",c).fetchone()
        Fin = h5py.File(fin,'r')
        g = Fin[fd('vanHove',c[0])]
        
        temp = g.attrs['temperature']
        dtime = g.attrs['dtime']
        if dtime != 0:
            j = int(time_step/dtime)

            count = g[fd('step',j)]['x']['disp_count'][:]
            edges = g[fd('step',j)]['x']['disp_edges'][:]
            fig.draw_line(edges,count/np.sum(count),label='%.2f'%temp,color=cm.get_color(temp))
        
        del g
        Fin.close()
        del Fin

def extract_vanHove(c,conn,time_step,count_cut,wind=1,norm=False):
    """
    Function to encapsulate extracting a vanHove at a given time step.
    The actual time difference will always be less than requested due
    to rounding and finite sampling rates.

    This function takes care of looking up the correct files and dealing with
    the hdf layer.
    """
    (fin,) = conn.execute("select fout from vanHove where comp_key = ?",c).fetchone()
    Fin = h5py.File(fin,'r')
    g = Fin[fd('vanHove',c[0])]

    temp = g.attrs['temperature']
    dtime = g.attrs['dtime']
    max_step = g.attrs['max_step']
    j = int(time_step/dtime)
    
    if j > max_step or j <= 0:
        return np.array([]),np.array([]),0,0,[]
    
    (edges,count,x_lim) = _extract_vanHove(g,j,count_cut,wind)
    
    del g
    Fin.close()
    del Fin
        
    if norm:
        count = count/np.sum(count)

    return edges,count,temp,dtime,x_lim

def extract_vanHove_all(c,conn,wind=1):
    """
    Alternate extraction function that returns all of the 
    
    """

    (fin,min_track_length) = conn.execute("select fout,min_track_length from vanHove where comp_key = ?",c).fetchone()
    Fin = h5py.File(fin,'r')
    g = Fin[fd('vanHove',c[0])]

    temp = g.attrs['temperature']
    dtime = g.attrs['dtime']
    
    vanHove = [_extract_vanHove(g,step,0,wind) for step in range(1,min_track_length)]
    
    del g
    Fin.close()
    del Fin
    
    return vanHove,temp,dtime


    

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
    
    edges = np.array([np.mean(edges[j:j+wind]) for j in range(0,len(count)-wind,wind)])
    count = np.array([np.sum(count[j:j+wind]) for j in range(0,len(count)-wind,wind)])

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
        np.sum(
            (edges**4) *tmp_count )/
        (3*(np.sum((edges**2)*tmp_count)
                **2))
        ) -1


def compute_alpha(comp,conn,wind ,min_count ):
    """
    Takes in a computation number, opens the hdf file, extracts all the vanHove distributions
    an computes alpha2 for them.

    Returns a list of tuples(wait_time,alpha2) and the temperature
    """
    
    (fin,max_step) = conn.execute("select fout,max_step from vanHove where comp_key = ?",comp).fetchone()
    
    Fin = h5py.File(fin,'r')
    g = Fin[fd('vanHove',comp[0])]

    temp = g.attrs['temperature']
    dtime = g.attrs['dtime']
    
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
def plot_alpha2(comp_list,conn,wind,min_count,Tc=0,ax=None):
    """
    Takes in a comp_list and plots the alpha2().

    returns a list of the compute_alpha results
    """
    a_lst = [compute_alpha(c,conn,wind,min_count) for c in comp_list]
    _plot_alpha2(a_lst,ax)

    return a_lst

def _plot_alpha2(a_list,ax):
    """
    plots alpha2 from the list of outputs of compute_alpha handed in
    """
    cm = plots.color_mapper(27,33)
    if ax is None:
        ax =  plots.set_up_axis(r'$\Delta \tau$ [ms]',r'$\alpha_2$','')

    for a,temp in a_list:
        dt,a2 = zip(*a)
        ax.step(dt,a2,
                label='%.2f'%temp,
                color=cm.get_color(temp),
                where='post')
    
def plot_alpha2_grid(lsts,conn,title):
    ''' plots a grid of alpha2 plots.  Each entry in lsts is a pair
    the first entry is a list of 3-tuples that are (vanHove key, wind,min_count) and
    the second entry is a string for the title'''
    cm = plots.color_mapper(27,33)
    fig = mplt.figure()
    fig.suptitle(title)
    dims = figure_out_grid(len(lsts))
    
    plt_count = 1
    for lst,t in lsts:
        sp_arg = dims +(plt_count,)
        ax = fig.add_subplot(*sp_arg)
        ax.grid(True)
        a_list = [compute_alpha(c,conn,wind,min_count) for (c,wind,min_count) in lst]


        for a,temp in a_list:
            dt,a2 = zip(*a)
            ax.step(dt,a2,color=cm.get_color(temp))
            
        plt_count +=1
        
    mplt.draw()
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


    #  istatus = plots.non_i_plot_start()
    
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

    # plots.non_i_plot_start(istatus)

    del g
    Fin.close()
    del Fin

def plot_vanHove_sp(comp_lst,time_step,conn,wind =1,func = fitting.fun_exp_p_gauss,norm=False):
    '''Plots a grid array  of the Von Hove functions at the time step
    given for all of the comps in the comp_lst'''

    
    
    fig = mplt.figure()
    fig.suptitle(r'van Hove dist $\tau = %.2f \rm{s}$'% (time_step/1000))
    dims = figure_out_grid(len(comp_lst))
    
    plt_count = 1
    outs = []
    
    data = [extract_vanHove(c,conn,time_step,1,wind,norm=norm) for c in comp_lst]
    tmps = [d[2] for d in data]
    cm = plots.color_mapper(np.min(tmps),np.max(tmps))
    data.sort(key=lambda x: x[2])
    
    extream = [(np.min(d[1]), np.max(d[1])) for d in data]
    print extream
    y_lim = [np.min([d[0] for d in extream]), np.max([d[1] for d in extream])]
    print y_lim
    y_lim = [np.log(y_lim[0]*.9), np.log(y_lim[1]*1.1)]
    print y_lim
    for (edges,count,temp,dtime,x_lim) in data:
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
        
        ax.plot(edges,np.log(func(out.beta,edges)),'--k')
        if norm:
            ax.set_ylabel(r'$\ln{P(\Delta)}$')
        else:
            ax.set_ylabel(r'$\ln{N}$')
        ax.step(edges,np.log((count)),color=cm.get_color(temp))
        ax.set_title('%.2f'%temp + ' dtime:%d ms'%dtime)
        ax.set_xlim(x_lim)
        ax.set_ylim(y_lim)
        plt_count += 1

    mplt.draw()

    return outs,tmps


def plot_vanHove_single_axis(comp_lst,time_step,conn,title=None,ax=None,wind =1,norm=False,**kwargs):
    ''''''
    r_scale = 6.45/60
    if 'Tc' in kwargs:
        T_conv_fun = ltc.T_to_phi_factory(kwargs['Tc'],ltc.linear_T_to_r_factory(-.011,0.848))
        del kwargs['Tc']
        label_str = '%(#)0.2f'
    else:
        T_conv_fun = lambda x:x
        label_str = '%(#)0.1f C'

    bin_min = 0
    if 'bin_min' in kwargs:
        bin_min = kwargs['bin_min']
        del kwargs['bin_min']

    norm_fun = np.max
    if 'norm_fun' in kwargs:
        norm_fun = kwargs['norm_fun']
        del kwargs['norm_fun']
    
    
    data = [extract_vanHove(c,conn,time_step,bin_min,wind,norm=norm) for c in comp_lst]
    tmps = [d[2] for d in data]
    cm = plots.color_mapper(np.min(tmps),np.max(tmps))
    data.sort(key=lambda x: x[2])
    if title is None:
        title = r'van Hove $\tau$: %g s'%(time_step/1000)
    cm = plots.color_mapper(np.min(tmps),np.max(tmps))
    print tmps
    if ax is None:
        ax = plots.set_up_axis(r'$\Delta$ [$\mu$m]',r'$N/N_{max}$',title)
        
        
    [ax.step(edges*r_scale,count/norm_fun(count),
                 label=label_str%{'#':T_conv_fun(temp)},color=cm.get_color(temp),where='post',**kwargs)
     for (edges,count,temp,dtime,x_lim) in data if len(count)>0]
    ax.set_yscale('log')
    ax.set_ylim(1e-3,1)



def plot_vanHove_diff(vh_a,vh_b,time_step,conn,ax,wind =1,**kwargs):
    ''''''

    
    
    
    
    data = [extract_vanHove(c,conn,time_step,1,wind,norm=norm) for c in [vh_a,vh_b]]
        
    data.sort(key=lambda x: x[2])
    if title is None:
        title = r'van Hove $\tau$: %g s'%(time_step/1000)
    cm = plots.color_mapper(np.min(tmps),np.max(tmps))
    print tmps
    if ax is None:
        ax = plots.set_up_axis(r'$\Delta$ [px]',r'$N/N_{max}$',title)
        
    count = data[0][1] - data[1][1]

    edges = (data[0][0] + data[1][0])/2
        
    ax.plot(edges,count/np.max(count),**kwargs)
    
    


def plot_vanHove(c,conn,ax,tau,*args,**kwargs):
    '''Extracts and plots a van Hove plot on to the given axis'''
    
    wind = 1
    if 'wind' in kwargs:
        wind = kwargs['wind']
        del kwargs['wind']
        
    (edges,count,temp,dtime,x_lim) = extract_vanHove(c,conn,tau,1,wind,False)
    ax.plot(edges,count/np.max(count),*args,**kwargs)

def plot_traking_variation(dset_key,time_step,conn,wind=1,title=None,norm=None):
    comp_lst = conn.execute("select comp_key from vanHove where dset_key = ? and min_track_length = 45",(dset_key,)).fetchall()
    print comp_lst
    data = [extract_vanHove(c,conn,time_step,-1,wind,norm=norm) for c in comp_lst]
    if title is None:
        title = 'van Hove %(#)0.1f C'%{'#':data[0][2]}

        

    fig = plots.tac_figure(r'$\Delta$ [px]',r'$N/N_{max}$',title,func=matplotlib.axes.Axes.semilogy)
    for (a,b) in zip(data,comp_lst):
        (edges,count,temp,dtime,x_lim,comp_key) = a+b
        fig.draw_line(edges,count/(np.max(count) ),label=str(comp_key))
    
    

def plot_hwhm(comp_lst,time_step,conn,title=None,wind =1,norm=True,ax= None,**kwargs):
    '''tramp function to keep from needing to re-write a lot of code
    due to the below renaming'''
    plot_vanHove_reduction(comp_lst,time_step,conn,fun=_vh_hwhm,title=title,wind=wind,norm=norm,ax=ax,**kwargs)
    
def plot_vanHove_reduction(comp_lst,time_step,conn,fun=None,title=None,wind =1,norm=True,ax= None,**kwargs):
    '''the half-width half max of the van Hove distributions vs
    temperature at a fixed time step'''

    if fun is None:
        fun = _vh_msd
    
    r_scale = 6.45/60
    if 'Tc' in kwargs:
        T_conv_fun = ltc.T_to_phi_factory(kwargs['Tc'],ltc.linear_T_to_r_factory(-.011,0.848))
        del kwargs['Tc']
        x_lab = r'$\phi/\phi^*$'
        if title is None:
            title = r'Half width half max vs $\phi/\phi^*$'
    
    else:
        T_conv_fun = lambda x:x
        x_lab = 'T [C]'
        if title is None:
            title = 'Half width half max vs Temperature'

        
    data = [extract_vanHove(c,conn,time_step,1,wind,norm=norm) for c in comp_lst]
    tmps = [d[2] for d in data]
    cm = plots.color_mapper(np.min(tmps),np.max(tmps))
    data.sort(key=lambda x: x[2])
    
    T = np.array([d[2] for d in data])
    hwhm = np.array([fun(d[0]*r_scale,d[1]) for d in data])
    
    print T
    print hwhm
    if ax is None  and fun ==_vh_hwhm:
        ax = plots.set_up_axis(x_lab,'hwhm [$\mu$m]',title)
    ax.plot(T_conv_fun(T),hwhm,'x-',label = '%.1f S'%(time_step/1000),**kwargs)
        
    return ax


def plot_hwhm_v_tau(comp,conn,steps,ax,wind =1,*args,**kwargs):
    '''the half-width half max of the van Hove distributions vs tau at
    a fixed temperature'''
    

        
        
    (vanHove,temp,dtime) = extract_vanHove_all(comp,conn,wind)
    
    
    if not 'label' in kwargs:
        kwargs['label'] = '%.2f'%temp
    if not 'color' in kwargs:
        cm = plots.color_mapper(27,32)
        kwargs['color'] = cm.get_color(temp)
    
    hwhm = [_vh_hwhm(v[0],v[1]) for v in vanHove]
      
        
        
    ax.plot(dtime*(np.arange(0,len(hwhm))+1),np.array(hwhm),*args,**kwargs)
        
    tmp_fig = mplt.gcf()
    mplt.figure(ax.get_figure().number)
    mplt.draw()
    mplt.figure(tmp_fig.number)

def set_up_hwhm_tau_axis():
    """
    Sets up the axis
    """
    fig = mplt.figure()
    ax = mplt.gca()
    ax.set_xlabel(r'$\tau$ [ms]')
    ax.set_ylabel('hwhm [pix]')
    ax.hold(True)
    ax.grid(True)
    ax.set_yscale('log')
    ax.set_xscale('log')
    return ax

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
    

def fix_vanHove_dtime(vh_key,conn):
    """Fixes vanHove data sets that have the dtime messed up"""

    # look up trk_key
    (trk_key,) = conn.execute("select track_key from vanHove_prams where comp_key = ?"
                              ,(vh_key,)).fetchone()
    # get avgerage dtime
    dtime = gen.avg_dtime(trk_key,conn)
    # set dtime attribute
    (fname,) = conn.execute("select fout from comps where comp_key = ? and function = 'vanHove'"
                              ,(vh_key,)).fetchone()
    Fin = h5py.File(fname,'r+')
    g = Fin[fd('vanHove',vh_key)]
    
    g.attrs['dtime'] = dtime

    Fin.close()
    del Fin


def _vh_hwhm(edges,count):
    """Finds the half width half max of the van Hove """

    # find the maximum value

    c_max = np.max(count)
    indx = np.nonzero(count >= (c_max/2))[0]
    return (edges[indx[-1]] - edges[indx[0]])/2

def _vh_msd(edges,count):
    """Computes the MSD  """
    
    # get bin centers
    cents = gen.get_bin_centers(edges)
    
    return np.sum(count*(cents**2))/np.sum(count)

def F_qt_plot(comp_lst,conn,q,wind=5,**kwargs):
    """Plots F(q,\tau)"""
    if 'Tc' in kwargs:
        T_conv_fun = ltc.T_to_phi_factory(kwargs['Tc'],ltc.linear_T_to_r_factory(-.011,0.848))
        del kwargs['Tc']
    
    else:
        T_conv_fun = lambda x:x
    

    x_lab = r'$\tau$ [ms]'
    y_lab = r'$F(q,\tau)$'

    fig,ax = plots.set_up_plot()
    plots.add_labels(ax,'',x_lab,y_lab)
        
        
    
    cm = plots.color_mapper(27,32)
    
    
    
    for vh_comp in comp_lst:
        Fqt,temp,dtime = _F_qt(vh_comp,conn,q,wind)
        Fqt = np.abs(Fqt)
        time = np.cumsum(dtime * (np.arange(0,len(Fqt))+1))
        ax.plot(time,Fqt,label='%.2f'%T_conv_fun(temp),color=cm.get_color(temp))

    ax.legend(loc=0)

    ax.set_ylim(0,1)
    
def _F_qt(vh_comp,conn,q,wind=5):
    """implements F(q,\tau) as described in the zexin nature paper """

    # get the van Hove data
    (vanHove,temp,dtime) = extract_vanHove_all(vh_comp,conn,wind)
    
    Fqt = [_Fqt_comp(vh,q) for vh in vanHove]

    return Fqt,temp,dtime
    
    
    

def _Fqt_comp(vh,q):
    """worker function to do the F(q,t) computation """
    r_scale = 6.45/60
    edges,count,x_lim = vh
    # make sure that vh is normalized
    count = count/np.sum(count)

    return np.sum(count * np.exp(1j*q*edges*r_scale))
    
def remove_vanHove_computation(comp_key,conn):
    '''Completely removes the computation from the results and the
    database'''
    (fin,) = conn.execute("select fout from vanHove where comp_key = ?",
                          comp_key).fetchone()
    print fin
    conn.execute("delete from vanHove where comp_key = ?",comp_key)
    conn.execute("delete from comps where comp_key = ?",comp_key)
    conn.commit()
    
    Fin = h5py.File(fin,'r+')
    del Fin[fd('vanHove',comp_key[0])]
    Fin.close()
    del Fin
    

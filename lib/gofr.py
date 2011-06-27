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

import itertools
import os.path
import re

import h5py
import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plts
import numpy as np
import scipy.odr as sodr
import scipy.optimize as sopt



import T_convert as ltc
import fitting
import general as gen
import plots 
import util 
from util import cord_pairs


##################
#helper functions#
##################
def find_peaks_fit(gp,dfun,p):
    """Looks for peaks based on hints from the derivative of the expected function and
    the fit parameters
    returns 4 lists, one entry per peak found
    lmax [(location, height),]
    lmin [(location, depth),]
    diffs [expected-actual, ] (min and max interleaved)
    sz [coeff in quad,]
    """

    dfun = fitting.fun_flipper(dfun)
    wind = 20;

    lmax = []
    lmin = []
    diffs = []
    sz = []
    # find the first max
    indx = np.argmax(gp.y[15:])+15
    pfit = fit_quad_to_peak(gp.x[indx-wind:indx+wind],gp.y[indx-wind:indx+wind])

    lmax.append((pfit.beta[1],pfit.beta[2]))
    diffs.append(gp.x[indx]-pfit.beta[1])
    # get expected spacing from
    e_spc = (2*np.pi)/p[1]

    print e_spc

    cur_state = np.sign(pfit.beta[0])
    sz.append(pfit.beta[0])
    #start at first peak + 1/4 e_spc
    cur_pos = pfit.beta[1] + e_spc/4
    # step over x range in e_spc/2 steps
    max_pos = np.max(gp.x)
    print max_pos
    print cur_pos
    while cur_pos < max_pos:
        # get zero crossing from dfun
        try:
            crit_p = sopt.brentq(dfun,cur_pos,cur_pos + e_spc/2,args=p)
            #print 'diff from expected ' + str(crit_p - (cur_pos + e_spc/4))
        except ValueError:
            print "no zero found"
            break
        # convert value into indx
        indx = gen.val_to_indx(gp.x,crit_p)
        # pick out window around that box
        pfit = fit_quad_to_peak(gp.x[indx-wind:indx+wind],gp.y[indx-wind:indx+wind])
        

        # determine if max or min
        # add center/value to max or min output
        if np.abs(pfit.beta[0])<.05:
            print "peak too small"
            break

        print 'diff from crossing and min/max ' + str( crit_p - pfit.beta[1])
        diffs.append(crit_p - pfit.beta[1])
        sz.append(pfit.beta[0])
        
        if pfit.beta[0] >0:
            if cur_state != -1:
                print "found two peaks in a row"
                break
            lmin.append((pfit.beta[1],pfit.beta[2]))
            cur_state = np.sign(pfit.beta[0])
        elif pfit.beta[0]<0:
            if cur_state != 1:
                print "found two troughs in a row"
                break
            lmax.append((pfit.beta[1],pfit.beta[2]))
            cur_state = np.sign(pfit.beta[0])

            
        # increment cur_pos
        cur_pos = crit_p +  e_spc/4
        wind +=1
        pass
    
    return lmax,lmin,diffs,sz

def fit_quad_to_peak(x,y):
    """Fits a quadratic to the data points handed in """
    def quad(B,x):
        return B[0] *(x -B[1]) ** 2 + B[2]

    beta = (0,np.mean(x),y[gen.val_to_indx(x,np.mean(x))])

    data = sodr.Data(x,y)
    model = sodr.Model(quad)
    worker = sodr.ODR(data,model,beta)
    out = worker.run()
    return out

def _draw_gofr_hex_lines(ax,r0):
    '''This function will draw on a graph vertical linse where peaks in g(r) are expected
    for a hex packing'''

    lin_scale = .85#((4+2*math.sqrt(3))/2 -2)/2
    irr_pos = [2*math.sqrt(3),  math.sqrt(28)]
    irr_pos_txt = [r'$2\, \sqrt{3}$',r'$\sqrt{28}$']
    for s in range(2,9):
        ax.plot(s*r0*np.ones(2),[0 ,3],'r')
        ax.annotate(str(s),xy=(s*r0,2.5),xycoords='data',
                    xytext=(-1,0),textcoords='offset points')
    for s,t in zip(irr_pos,irr_pos_txt):
        ax.annotate(t,xy=(s*r0,2.75),xycoords='data',
                    xytext=(-1,0),textcoords='offset points')
        ax.plot(s*r0*np.ones(2),[0 ,3],'k')
    for s in range(1,6):
        ax.plot((1+ s*lin_scale)*2*r0*np.ones(2),[0 ,3],'m')
        ax.annotate(str(s),xy=(2*r0*(1+s*lin_scale),2.25),xycoords='data',
                    xytext=(-1,0),textcoords='offset points')
######################
#extraction functions#
######################
def get_gofr3D(comp_num,conn):
    '''Takes in computation number and database connection and extracts the given g(r) and returns it as a
    cord_pairs object (left bin edge,value)'''

    gname='gofr3D'

    dset_names = ['bin_count', 'bin_edges']
    
    res = conn.execute("select fout from comps where comp_key = ?",(comp_num,)).fetchall()
    if not len(res) == 1:
        print len(res)
        raise util.dbase_error("error looking up computation")

    
    g = get_gofr_group(res[0][0],gname,comp_num)
    gofr = np.array(g[dset_names[0]])
    bins = np.array(g[dset_names[1]])
    return cord_pairs(bins,gofr)

def get_gofr2D(comp_num,conn):
    '''Takes in computation number and database connection and extracts the given g(r) and returns it as a
    cord_pairs object (left bin edge,value)'''
    scale = 6.45/60
    dset_names = ['bin_count', 'bin_edges']
    
    res = conn.execute("select fout from gofr where comp_key = ?",(comp_num,)).fetchall()
    if not len(res) == 1:
        print len(res)
        raise util.dbase_error("error looking up computation")
    
    print comp_num
    F = h5py.File(res[0][0],'r')
    print "gofr_%(#)07d"%{"#":comp_num}
    print "gofr_%(#)07d"%{"#":comp_num} in F.keys()
    g = F["gofr_%(#)07d"%{"#":comp_num}]
    gofr = np.array(g[dset_names[0]])
    bins = np.array(g[dset_names[1]])*scale
    F.close()
    return cord_pairs(bins,gofr)

def get_gofr2D_T(comp_num,conn):
    '''Takes in computation number and database connection and extracts the given g(r) and returns it as a
    cord_pairs object (left bin edge,value)'''

    dset_names = ['bin_count', 'bin_edges']
    
    res = conn.execute("select fout from gofr where comp_key = ?",(comp_num,)).fetchall()
    if not len(res) == 1:
        print len(res)
        raise util.dbase_error("error looking up computation")
    
    print comp_num
    F = h5py.File(res[0][0],'r')
    
    

    # data
    g = F["gofr_%(#)07d"%{"#":comp_num}]
    gofr = np.array(g[dset_names[0]])
    bins = np.array(g[dset_names[1]])*(6.45/60)

    # temperature
    if 'temperature' in  F["gofr_%(#)07d"%{"#":comp_num}].attrs:
        t = F["gofr_%(#)07d"%{"#":comp_num}].attrs['temperature']
        print 't from file',t
    else:
        # hack to get the file format from 2010-05-24 to work
        (data_fname,) = conn.execute("select fname from dsets where dset_key in " +
                                     "(select dset_key from gofr where comp_key = ?)",
                                     (comp_num,)).fetchone()

        pos_temps = re.findall('[23][0-9]-[0-9]',data_fname)
        if len(pos_temps) > 0:
            t = float(pos_temps[0].replace('-','.'))
            print 't from filename'
        elif data_fname.find('room')>0:
            t = 27.0
        
        else:
            pos_temps = re.findall('[23][0-9]\.[0-9]',data_fname)
            if len(pos_temps) > 0:
                t = float(pos_temps[0].replace('-','.'))
                print 't from filename'
            else:
                print data_fname
                raise Exception("can't figure out temperature")
    
    
    F.close()
    # kludge to deal with iden errors
    gofr[0:3] = 0
    return cord_pairs(bins,gofr),t

def get_gofr2D_rho(comp_num,conn):
    '''Takes in computation number and database connection and extracts the given g(r) and returns it as a
    cord_pairs object (left bin edge,value)'''

    dset_names = ['bin_count', 'bin_edges']
    
    res = conn.execute("select fout from gofr where comp_key = ?",(comp_num,)).fetchall()
    if not len(res) == 1:
        print len(res)
        raise util.dbase_error("error looking up computation")
    
    
    try:
        scale = 6.45/60
        F = h5py.File(res[0][0],'r')
        g = F["gofr_%(#)07d"%{"#":comp_num}]
        rho = g.attrs['rho']/(scale**2)
        
        gofr = np.array(g[dset_names[0]])
        bins = np.array(g[dset_names[1]])*scale
 
    finally:
        F.close()

    return cord_pairs(bins,gofr),rho



def get_gofr2D_dict(comp_num,conn):
    '''Takes in computation number and database connection and
    extracts the given g(r) and returns it as a cord_pairs object
    (left bin edge,value)
    Also return a dictionary of all the meta data in the file.
    '''

    dset_names = ['bin_count', 'bin_edges']
    
    
    
    res = conn.execute("select fout from gofr where comp_key = ?",comp_num).fetchall()
    
    if not len(res) == 1:
        print len(res)
        raise util.dbase_error("error looking up computation")
    
    
    try:
        
        F = h5py.File(res[0][0],'r')
        g = F["gofr_%(#)07d"%{"#":comp_num[0]}]
        attr = dict(g.attrs)
        # see if scale is in the meta data, if not add it
        if 'scale' not in attr:
            attr['scale'] = 6.45/60

        # convert the units on rho
        attr['rho'] = attr['rho']/(attr['scale']**2)
        
        # deal with temperature
        if 'temperature' not in  attr:
            # hack to get the file format from 2010-05-24 to work
            (data_fname,) = conn.execute("select fname from dsets where dset_key in " +
                                         "(select dset_key from gofr where comp_key = ?)",
                                         (comp_num,)).fetchone()

            pos_temps = re.findall('[23][0-9]-[0-9]',data_fname)
            if len(pos_temps) > 0:
                t = float(pos_temps[0].replace('-','.'))
            elif data_fname.find('room')>0:
                t = 27.0
            else:
                pos_temps = re.findall('[23][0-9]\.[0-9]',data_fname)
                if len(pos_temps) > 0:
                    t = float(pos_temps[0].replace('-','.'))
                else:
                    print data_fname
                    raise Exception("can't figure out temperature")
            attr['temperature'] = t


        
        gofr = np.array(g[dset_names[0]])
        bins = np.array(g[dset_names[1]])*attr['scale']

        
 
    finally:
        F.close()

    return cord_pairs(bins,gofr),attr

def get_gofr_tmp(fname,comp_num,conn):
    F = h5py.File(fname,'r')
    print fname
    
    if 'temperature' in  F["gofr_%(#)07d"%{"#":comp_num}].attrs:
        t = F["gofr_%(#)07d"%{"#":comp_num}].attrs['temperature']
        print 't from file',t
    else:
        # hack to get the file format from 2010-05-24 to work
        (data_fname,) = conn.execute("select fname from dsets where dset_key in " +
                                     "(select dset_key from gofr where comp_key = ?)",
                                     (comp_num,)).fetchone()

        pos_temps = re.findall('[23][0-9]-[0-9]',data_fname)
        if len(pos_temps) > 0:
            t = float(pos_temps[0].replace('-','.'))
            print 't from filename'
        elif data_fname.find('room')>0:
            t = 27.0
        
        else:
            pos_temps = re.findall('[23][0-9]\.[0-9]',data_fname)
            if len(pos_temps) > 0:
                t = float(pos_temps[0].replace('-','.'))
                print 't from filename'
            else:
                print data_fname
                raise Exception("can't figure out temperature")
    F.close()
    return t


def get_rhoT(comp,conn):
    '''Returns the tuple of (rho,temperature) for the given computation'''

    
    (dset_key,fname) = conn.execute("select dset_key,fout from gofr where" +
                                    " comp_key = ?",comp).fetchone()

    
    cmap = plots.color_mapper(26,31)


    
    F = h5py.File(fname,'r')
    g = F[gen.fd('gofr',comp[0])]
    rho = g.attrs['rho'] 
    temp = g.attrs['temperature'] 
    
    del g
    F.close()
    del F
    
    return (rho,temp)

    
def rebin_gofr(comp_key,conn,bin_count):
    ''' Takes a g(r) computation and re bins the data to wider bins
    (integer multiple of existing bins

    returns the same data is the same format as the other g(r) extraction functions.

    '''

    # extract the raw g(r)
    data,attr = get_gofr2D_dict(comp_key,conn)
    # convert the bin values to un-normalized counts

    
    bins = np.append(data.x ,[attr['max_range'] * attr['scale']])
    vals = data.y * np.diff(bins**2)    # we don't need pi because it will be divided out below
    
    # calculate new bin edges
    n_bins = data.x[::bin_count]
    if len(data.x)%bin_count == 0:
        n_bins = np.append(n_bins, bins[(len(data.x)//bin_count)*bin_count])
    # sum bin counts
    n_vals = np.array([np.sum(vals[j:j+bin_count]) for j in range(0,len(vals)-bin_count+1,bin_count)])
    # re-normalize bin values

    n_vals = n_vals/np.diff(n_bins**2)  # we don't need pi because it is omitted above as well

    return cord_pairs(n_bins[:-1],n_vals),attr


####################
#plotting functions#
####################
def make_2d_gofr_plot(comp_key,conn,fname = None):
    # add error handling on all of these calls
    
    # get comp_number of gofr
    (key,g_fname) = conn.execute("select dset_key,fout from comps where comp_key == ? and function == 'gofr'",(comp_key,)).fetchone()
    
        
    # get dset name
    (sname, temp) = conn.execute("select sname,temp from dsets where key == ? ",(key,)).fetchone()

    print sname + " " + str(temp)
    group = gen.get_gofr_group(g_fname,'gofr',comp_key)



    # make plot
    istatus = plts.isinteractive();
    if istatus:plts.ioff()
    
    dset_names = ['bin_count', 'bin_edges']
    fig = plts.figure()
    
    ax = fig.add_axes([.1,.1,.8,.8])
    ax.hold(True)
    ax.plot(group[dset_names[1]][:],group[dset_names[0]])


    # finds the location of the maximum, assume to be the first peak
    d0 = group[dset_names[1]][np.argmax(group[dset_names[0]])]
    print np.argmax(group[dset_names[0]])
    print d0
    #_draw_gofr_hex_lines(ax,d0/2)
    ax.set_title(sname + " temp: " + str(temp) + ' d0' + str(d0))
    ax.set_ylim(0,3)

    ax.grid(True)

    # save figure
    if not fname == None:
        f_path = '/home/tcaswell/python/figures/' + sname + '/'
        if not os.path.isdir(f_path):
            os.makedirs(f_path)
        
        fname = f_path + str(key) + '_gofr2D.png'
        fig.savefig(fname)


     
        
    if istatus:
        plts.ion()
        plts.show()
    else:
        plts.close(fig)


def plot_with_fitting(g_key,conn):
    """Plots a single g(r) with the fitting function on top
    
    
    """

    (comp_key,fname,dset_key) = conn.execute("select comp_key,fout,dset_key\
    from gofr where comp_key = ?",(g_key,)).fetchone()
    

    
    
    temps = get_gofr_tmp(fname,comp_key,conn) 
    print temps

    g = get_gofr2D(comp_key,conn) 
    (fits,r0um) = fit_gofr3(g,1,fitting.fun_decay_exp_inv_gen) 

    
    
    istatus = plots.non_i_plot_start()

    fig = plots.tac_figure('r[$\mu m$]','g(r)','g(r) + fitting temp: %.2f, '%temps + str(dset_key) ,
                       func = matplotlib.axes.Axes.step) 
    print fits.beta
    
    fig.draw_line(g.x,g.y,label='data')
    fig.draw_line(g.x[25:],fitting.fun_decay_exp_inv_gen(r0um)(fits.beta,g.x[25:]),label='fitting')
    fig.axis((0,12),(0,3))
    plots.non_i_plot_stop(istatus)



def peak_locations(comp,conn,ax,**kwargs):
    """
    Plots the peak locations for the given computation
    """

 

    
    res = conn.execute("select comp_key,fout " + 
                       "from gofr where comp_key = ?",comp).fetchone()
          

    
    
    
    
    gofr = get_gofr2D(res[0],conn) 
    fit,r0= fit_gofr3(gofr,2.1,fitting.fun_decay_exp_inv_gen)
    
    peaks = find_peaks_fit(gofr,fitting.fun_decay_exp_inv_dr_gen(r0),fit.beta)
        


    
    istatus = plots.non_i_plot_start()
    
    # make plot for trough locations
    
    leg_strs = []
    leg_hands = []
    print peaks
    pk_lst,pk_ht = zip(*peaks[0])
    print pk_lst
    tr_lst,tr_ht = zip(*peaks[1])
    ax.plot(np.arange(len(pk_lst))+1,pk_lst,linestyle='',**kwargs)
    ax.plot(np.arange(0,len(tr_lst))+1.5,tr_lst,linestyle='',**kwargs)
    
    plots.non_i_plot_stop(istatus)
    

###########################
#function of tmp functions#
###########################
def tmp_series_gn2D(comp_list,conn,**kwargs):
    """Makes plots for g_1 to g_n
    
    comp_list : list of 1 element tuple that contain the computation
    number to be plotted
    """


    r_scale = 6.45/60
    if 'Tc' in kwargs:
        T_conv_fun = ltc.T_to_phi_factory(kwargs['Tc'],ltc.linear_T_to_r_factory(-.011,0.848))
        del kwargs['Tc']
        x_lab = r'$\phi/\phi^*$'
        
    else:
        T_conv_fun = lambda x:x
        x_lab = 'T [C]'
        

    
    res = [conn.execute("select comp_key,fout\
    from gofr where comp_key = ?",c).fetchone()
           for c in comp_list]

    
    
    temps = [get_gofr_tmp(r[1],r[0],conn) for r in res]
    print temps
    
    gofrs = [get_gofr2D(r[0],conn) for r in res]
    fits_r0 = [fit_gofr3(g,2.1,fitting.fun_decay_exp_inv_gen) for g in gofrs]
    fits,r0 = zip(*fits_r0)
    peaks = [find_peaks_fit(g,fitting.fun_decay_exp_inv_dr_gen(r),p.beta)
             for g,p,r in zip(gofrs,fits,r0)]


    
    istatus = plots.non_i_plot_start()
    
    # make g_n plots for peaks
    fig,ax = plots.set_up_plot()
    
    for j in range(max([len(p[0]) for p in peaks])):
        pairs = [(t,p[0][j][1]-1) for (t,p) in itertools.izip(temps,peaks) if len(p[0]) > j]
        t,v = zip(*pairs)
        ax.plot(T_conv_fun(np.array(t)),v,'x-',label='$g_{'+str(j)+'}$' )
        

    l = ax.legend()
    l.set_bbox_to_anchor((1.125,1))
    plots.add_labels(ax,'$g_n$ maximums',r'T [C]','g(peak) -1')


    
    # make g_n plots for troughs
    fig,ax = plots.set_up_plot()
    leg_strs = []
    leg_hands = []
    for j in range(max([len(p[1]) for p in peaks])):
        pairs = [(t,p[1][j][1]-1) for (t,p) in itertools.izip(temps,peaks) if len(p[1]) > j]
        t,v = zip(*pairs)
        ax.plot(T_conv_fun(np.array(t)),v,'x-',label='$g_{'+str(j)+'}$' )

    ax.legend()
    plots.add_labels(ax,'$g_n$ minimums',r'T [C]','g(peak) -1')
    

    # make plot for trough locations
    fig,ax = plots.set_up_plot()
    leg_strs = []
    leg_hands = []
    for j in range(max([len(p[1]) for p in peaks])):
        pairs = [(t,p[1][j][0]) for (t,p) in itertools.izip(temps,peaks) if len(p[1]) > j]
        t,v = zip(*pairs)
        ax.plot(T_conv_fun(np.array(t)),v,'x-',label='$g_{'+str(j)+'}$' )

    ax.legend()
    plots.add_labels(ax,'minimum locations',r'T [C]','r [$\mu m$]')


    # make plot for peak locations
    fig,ax = plots.set_up_plot()
    
    for j in range(max([len(p[0]) for p in peaks])):
        pairs = [(t,p[0][j][0]) for (t,p) in itertools.izip(temps,peaks) if len(p[0]) > j]
        t,v = zip(*pairs)
        ax.plot(T_conv_fun(np.array(t)),v,'x-',label='$g_{'+str(j)+'}$' )

    ax.legend()
    plots.add_labels(ax,'peak locations',r'T [C]','r [$\mu m$]')

    
    plots.non_i_plot_stop(istatus)
    
    return peaks


def plot_gofr_inset(comp_key,conn,main_lim_max = None, inset_lim = None,*args,**kwargs):
    '''Plots a single g(r) plot with an inset of the higher peaks '''
    if 'Tc' in kwargs:
        T_conv_fun = ltc.T_to_phi_factory(kwargs['Tc'],ltc.linear_T_to_r_factory(-.011,0.848))
        del kwargs['Tc']
        x_lab = r'$\phi/\phi^*$'
        
    else:
        T_conv_fun = lambda x:x
        x_lab = 'T [C]'
        

    if 'r0' in kwargs:
        r_0 = kwargs['r0']
        del kwargs['r0']
    else:
        r_0 = 1.8

    res = conn.execute("select comp_key,dset_key from gofr where comp_key = ?",comp_key).fetchone()
    (d_fname,) = conn.execute("select fname from dsets where dset_key = ?",res[1:]).fetchone()
    (g,temp) = get_gofr2D_T(res[0],conn)

    istatus = plots.non_i_plot_start()
    
    f = plts.figure()
    #f.set_size_inches(4,4,forward=True)
    
    a1 = f.add_axes([.1,.1,.85,.8])
    a2 = f.add_axes([.505+.055,.575+.01,.38,.28])

    a1.step(g.x,g.y-1)
    a1.grid(True)
    a1.set_xlabel(r'r [$\mu m$]')
    a1.set_ylabel(r'$g(r) - 1$')
    

    a2.step(g.x[1000:],g.y[1000:]-1)
    ## a2.set_xlabel(r'r [$\mu m$]')
    ## a2.set_ylabel(r'G(r) - 1')
    a2.grid(True)
    a2.get_yaxis().get_major_formatter().set_powerlimits((2,2))
    a2.get_yaxis().set_major_locator(matplotlib.ticker.LinearLocator(5))

    a1.set_title((str(res[1]) + ', '
                 + os.path.basename(d_fname)
                 + ', %.2f, '%T_conv_fun(temp)
                 + gen.get_acq_time(res[1],conn)).replace('_','\_'))

    
    fit,r0 = fit_gofr3(g,r_0,fitting.fun_decay_exp_inv_gen)
    fit_fun = fitting.fun_decay_exp_inv_gen(r0)
    
    x_range = np.arange(r0,np.max(g.x),.01)
    a1.plot(x_range,fit_fun(fit.beta,x_range)-1,'k-')
    r_ind = np.argmax(g.x>r0)
    a1.plot(g.x[r_ind:],fit_fun(fit.beta,g.x[r_ind:])-g.y[r_ind:],'m-')
    a1.plot(g.x,fit_fun(fit.beta,g.x)-1,'g--')
    x_range = np.arange(g.x[1000],np.max(g.x),.01)
    a2.plot(x_range,fit_fun(fit.beta,x_range)-1,'k-')
    a2.plot(g.x[1000:],fit_fun(fit.beta,g.x[1000:])-g.y[1000:],'m-')

    print 'xi: %.4f'%fit.beta[0]
    print 'K: %.4f'%(2*np.pi/fit.beta[1])
    print 'C: %.4f'%fit.beta[2]
    print 'r*: %.4f'%fit.beta[3]

    if inset_lim is None:
        y_lim = np.max(np.abs(a2.get_ylim()))
        a2.set_ylim(-y_lim,y_lim)
    else:
        a2.set_ylim(-inset_lim,inset_lim)
    
    if main_lim_max is not None:
        a1.set_ylim(-1,main_lim_max)
    
    plots.non_i_plot_stop(istatus)

    return f,a1,a2

    
def plot_gofr_diff(comp_key_1,comp_key_2,conn,main_lim_max = None, inset_lim = None,*args,**kwargs):
    '''Plots a single g(r) plot with an inset of the higher peaks '''
    if 'Tc' in kwargs:
        T_conv_fun = ltc.T_to_phi_factory(kwargs['Tc'],ltc.linear_T_to_r_factory(-.011,0.848))
        del kwargs['Tc']
        x_lab = r'$\phi/\phi^*$'
        
    else:
        T_conv_fun = lambda x:x
        x_lab = 'T [C]'
        

    if 'r0' in kwargs:
        r_0 = kwargs['r0']
        del kwargs['r0']
    else:
        r_0 = 1.8

    def helper(c_k):
        res = conn.execute("select comp_key,dset_key from gofr where comp_key = ?",c_k).fetchone()
        return   get_gofr2D_T(res[0],conn)
    
    g_1,t_1 = helper(comp_key_1)
    g_2,t_2 = helper(comp_key_2)

    if len(g_1.x) != len(g_2.x):
        die
        return

    istatus = plots.non_i_plot_start()
    
    f = plts.figure()
    #f.set_size_inches(4,4,forward=True)
    
    a1 = f.add_axes([.1,.1,.85,.8])
    a2 = f.add_axes([.505+.055,.575+.01,.38,.28])

    a1.step(g_1.x,g_1.y-g_2.y,label='diff')
    a1.step(g_1.x,g_1.y-1,'-',label='%.2f'%T_conv_fun(t_1))
    a1.step(g_2.x,g_2.y-1,'--',label='%.2f'%T_conv_fun(t_2))
    a1.grid(True)
    a1.set_xlabel(r'r [$\mu m$]')
    a1.set_ylabel(r'')
    

    a2.step(g_1.x[1000:],g_1.y[1000:]-g_2.y[1000:])
    a2.step(g_1.x[1000:],g_1.y[1000:]-1)
    a2.step(g_2.x[1000:],g_2.y[1000:]-1)
    ## ## a2.set_xlabel(r'r [$\mu m$]')
    ## ## a2.set_ylabel(r'G(r) - 1')
    a2.grid(True)
    a2.get_yaxis().get_major_formatter().set_powerlimits((2,2))
    a2.get_yaxis().set_major_locator(matplotlib.ticker.LinearLocator(5))

    ## a1.set_title((str(res[1]) + ', '
    ##              + os.path.basename(d_fname)
    ##              + ', %.2f, '%temp
    ##              + gen.get_acq_time(res[1],conn)).replace('_','\_'))

    
    if inset_lim is None:
        y_lim = np.max(np.abs(a2.get_ylim()))
        a2.set_ylim(-y_lim,y_lim)
    else:
        a2.set_ylim(-inset_lim,inset_lim)
    
    if main_lim_max is not None:
        a1.set_ylim(-1,main_lim_max)
    a1.legend(loc=4)
    plots.non_i_plot_stop(istatus)
    
    ## return f,a1,a2

def tmp_series_fit_plots(comp_list,conn,**kwargs):
    """Makes plots for g_1 to g_n
    
    comp_list : list of 1 element tuple that contain the computation
    number to be plotted
    """
    
    r_scale = 6.45/60
    if 'Tc' in kwargs:
        T_conv_fun = ltc.T_to_phi_factory(kwargs['Tc'],ltc.linear_T_to_r_factory(-.011,0.848))
        del kwargs['Tc']
        x_lab = r'$\phi/\phi^*$'
    else:
        T_conv_fun = lambda x:x
        x_lab = 'T [C]'
    if 'r0' in kwargs:
        r_0 = kwargs['r0']
        del kwargs['r0']
    else:
        r_0 = 1.8

    
    res = [conn.execute("select comp_key,fout " +
                        "from gofr where comp_key = ?",c).fetchone()
           for c in comp_list]


    
    temps = np.array([get_gofr_tmp(r[1],r[0],conn) for r in res])
    zip_lst = zip(temps,res)
    zip_lst.sort(key=lambda x:x[0])
    temps,res = zip(*zip_lst)
    temps = np.array(temps)
    print temps
    
    gofrs = [get_gofr2D(r[0],conn) for r in res]
    fits = [fit_gofr3(g,r_0,fitting.fun_decay_exp_inv_gen) for g in gofrs]
    fits,r0s = zip(*fits)
    
    istatus = plots.non_i_plot_start()

    fig_xi = plots.tac_figure(x_lab,r'fitting parameters [$\mu m$]','Fitting parameters') 

    print [p.beta[0] for p in fits ]
     
    fig_xi.draw_line(T_conv_fun(temps),[p.beta[0] for p in fits ],
                    marker='x',label=r'$\xi$')
    fig_xi.draw_line(T_conv_fun(temps),[(np.pi*2)/p.beta[1] for p in fits ],
                    marker='x',label=r'$K$')
    fig_xi.draw_line(T_conv_fun(temps),[p.beta[2] for p in fits ],
                    marker='x',label=r'$C$')
    fig_xi.draw_line(T_conv_fun(temps),[p.beta[3] for p in fits ],
                    marker='x',label=r'$r^*$')
    fig_xi.draw_line(T_conv_fun(temps),r0s,
                    marker='x',label=r'$r_0$')
    
    
    fig_err = plots.tac_figure('T [C]','errors [diff units]','fitting error')
    fig_err.draw_line(T_conv_fun(temps),[p.sum_square for p in fits ],'-x',label=r'sum square')
    fig_err.draw_line(T_conv_fun(temps),[p.sum_square_delta for p in fits ],'-x',label=r'sum square delta')
    fig_err.draw_line(T_conv_fun(temps),[p.res_var for p in fits ],'-x',label=r'res var')
    
    
    plots.non_i_plot_stop(istatus)


    
    return fits
def plot_residue(comp_list,conn):
    '''Plots the residue of the fitting and data '''
    
    res = [conn.execute("select comp_key,fout\
    from gofr where comp_key = ?",c).fetchone()
           for c in comp_list]


    
    temps = [get_gofr_tmp(r[1],r[0],conn) for r in res]
    zip_lst = zip(temps,res)
    zip_lst.sort(key=lambda x:x[0])
    temps,res = zip(*zip_lst)

    tmax = max(temps)
    tmin = min(temps)
    cmap = plots.color_mapper(tmin,tmax)
    tmax = max(temps)
    tmin = min(temps)
    
    cmap = plots.color_mapper(tmin,tmax)


    r0 = 2.1
    gofrs = [get_gofr2D(r[0],conn) for r in res]
    fit_r0 = [fit_gofr3(g,r0,fitting.fun_decay_exp_inv_gen) for g in gofrs]
    fits,r0_ums = zip(*fit_r0)
    
    
    istatus = plots.non_i_plot_start()
    fig = plots.tac_figure('r [um]',r'$\frac{g(r) - fit(r)}{g(r)}$','Fitting residue') 
    for f,g,t,r0um in zip(fits,gofrs,temps,r0_ums):
        fun = fitting.fun_decay_exp_inv_gen(r0um)
        g_tmp,r = _trim_gofr(g,r0)
        fig.draw_line(g_tmp.x,(g_tmp.y - fun(f.beta,g_tmp.x))/g_tmp.y,
                 label='%.2f'%t,
                 color = cmap.get_color(t) )
    
    plots.non_i_plot_stop(istatus)

def set_up_gn_plots(sname,T=True):
    gn_fig = plts.figure()
    gn_ax = gn_fig.add_axes([.1,.1,.8,.8])
    fig = plts.figure()
    ax = fig.add_axes([.1,.1,.8,.8])
    ax.hold(True)
    ax.grid(True)
    ax.set_title(sname)
    ax.set_xlabel(r'r [$\mu m$]')
    ax.set_ylabel(r'$G(r)-1$')
    if T:
        gn_ax.set_title(sname + r' $g_1(T)$')
        gn_ax.set_xlabel('T')
    else:
        gn_ax.set_title(sname + r' $g_1(\phi/\phi^*)$')
        gn_ax.set_xlabel('$\phi/\phi^*$')
    gn_ax.set_ylabel('$g_1-1$')
    
    
    return (ax,gn_ax)

def make_gofr_tmp_series(comp_list,conn,ax,gn_ax,T_correction = 0,*args,**kwargs):
    '''Takes in a sample name and plots all of the g(r) for it '''
    r_scale = 6.45/60
    if 'Tc' in kwargs:
        T_conv_fun = ltc.T_to_phi_factory(kwargs['Tc'],ltc.linear_T_to_r_factory(-.011,0.848))
        del kwargs['Tc']
    else:
        T_conv_fun = lambda x:x
    
    dset_names = ['bin_count', 'bin_edges']

    
    res = [conn.execute("select comp_key,fout,fin,dset_key from gofr where comp_key = ? ",c).fetchone()
           for c in comp_list]

    #res.sort(key=lambda x:x[2])
    print res
    
    #(sname,) = conn.execute("select sname from dsets where dset_key = ?",(res[0][3],)).fetchone()

    print "there are " + str(len(res)) + " entries found"
    # check interactive plotting and turn it off
    istatus = plts.isinteractive();
    print istatus
    if istatus:plts.ioff()

    
    fig = ax.get_figure()
    
    gn_g = []
    gn_t = []
    gn_p = []

    
    gn_fig = gn_ax.get_figure()
        
        
    temps = [get_gofr_tmp(r[1],r[0],conn)+T_correction for r in res]
    tmax = max(temps)
    tmin = min(temps)
    
    cmap = plots.color_mapper(tmin,tmax)
    for r,temperature in zip(res,temps):

        g = get_gofr2D(r[0],conn)
        
        
        ax.step(g.x,g.y-1,
                color = cmap.get_color(temperature)
                ,label = '%.2f'%T_conv_fun(temperature) + ', %d'%r[3])
        
            
        

        wind = 30
        m = np.argmax(g.y[5:])+5
        betas = fit_quad_to_peak(g.x[m-wind:m+wind],g.y[m-wind:m+wind])
        
        g1 = betas.beta[2]
        print g1,np.max(g.y[5:])
        if not np.isnan(g1):
            try:
                gn_p.append((temperature,g1))
                # gn_t.append((float(r[2]))
                # gn_g.append(np.max(g[dset_names[0]]))
            except (ValueError,TypeError ) :
                gn_p.append((25,np.max(g[dset_names[0]])))
                # gn_t.append(25)
                # gn_g.append(np.max(g[dset_names[0]]))
                pass
            
            

    gn_p.sort(lambda x,y: int(np.sign(x[0]-y[0])))
    for p in gn_p:
        gn_t.append(p[0])
        gn_g.append(p[1])
    print gn_t
    print gn_g
    ax.legend()
    
    
    gn_ax.plot(T_conv_fun(np.array(gn_t)),np.array(gn_g)-1,'x-',**kwargs)
    gn_ylim = list(gn_ax.get_ylim())
    gn_ylim[0] = 0
    gn_ax.set_ylim(gn_ylim)
        
    print 
    gn_ax.grid(True)

    plts.figure(fig.number)
    plts.draw()
    plts.figure(gn_fig.number)
    plts.draw()
        
    if istatus:
        print "displaying figure"
        plts.ion()
        plts.show()
    else:
        print "closing figure"
        plts.close(fig)
        plts.close(gn_fig)

    return zip(gn_t,gn_g)






    
def make_gn_plots(comp_list,conn,axs,*args,**kwargs):
    '''Takes in a computation list.  Plots the g_n line on to the nth
    axes in axs.  args and kwargs are passed to plot'''
    

    
    res = [conn.execute("select comp_key,fout\
    from gofr where comp_key = ?",c).fetchone()
           for c in comp_list]

    
    
    temps = [get_gofr_tmp(r[1],r[0],conn) for r in res]
    
    gofrs = [get_gofr2D(r[0],conn) for r in res]
    fits_r0 = [fit_gofr3(g,2.1,fitting.fun_decay_exp_inv_gen) for g in gofrs]
    fits,r0 = zip(*fits_r0)
    peaks = [find_peaks_fit(g,fitting.fun_decay_exp_inv_dr_gen(r),p.beta)
             for g,p,r in zip(gofrs,fits,r0)]
    j = 0
    # loop over axes
    for ax in axs:
        (t,val) = zip(*[(t,p[0][j][1]-1) for (t,p) in itertools.izip(temps,peaks) if len(p[0]) > j])
        
        ax.plot(t,val,*args,**kwargs)
        j+=1
        plts.figure(ax.get_figure().number)
        plts.draw()
        

    # check interactive plotting and turn it off
    istatus = plts.isinteractive();
    print istatus
    #if istatus:plts.ioff()
    
    
    
    
    #plts.draw()

def set_up_gn_plot(t_range,g_range,title,n=1):
    '''helper function to set up g_n axes '''

    fig = plts.figure()
    ax = fig.add_axes([.1,.1,.8,.8])
    ax.hold(True)
    ax.grid(True)
    ax.grid(True)
    ax.set_title(title)
    ax.set_ylim(g_range)
    ax.set_xlim(t_range)
    ax.set_xlabel('T')
    ax.set_ylabel('$g_%d-1$'%n)
    return ax
    
def tmp_series_tables(comp_list,conn):
    
    res = [conn.execute("select comp_key,fout\
    from comps where comps.comp_key = ?",c).fetchone()
           for c in comp_list]

    
    gofrs = [get_gofr2D(r[0],conn) for r in res]
        
    temps = [get_gofr_tmp(r[1],r[0],conn) for r in res]

    fits = [fit_gofr2(g,2.1,fitting.fun_decay_exp_inv) for g in gofrs]
    peaks = [find_peaks_fit(g,fitting.fun_decay_exp_inv_dr,p.beta)
             for g,p in zip(gofrs,fits)]

    # make max location table
    max_table = ['|%.2f|'%t + '%.3f|'%k[0][0][0] + '|'.join(['%.3f'%i for i in diff([p[0]/k[0][0][0] for p in k[0]])]) + '|' for k,t in zip(peaks,temps)]
    min_table = ['|%.2f|'%t + '%.3f|'%k[0][0][0]+'%.3f|'%(k[1][0][0]/k[0][0][0]) + '|'.join(['%.3f'%i for i in diff([p[0]/k[0][0][0] for p in k[1]])]) + '|' for k,t in zip(peaks,temps)]
    
####################    
#by plane functions#
####################    
def first_comp_g1_plot(comp_lst,conn,n=1,g1_fig=None,lab=None):
    '''
    Takes in a list of gofr_by_plane computations.  Makes the g_1 plot
    by averaging the first n computations.
    '''
    wind = 30
    points = []
    points2 = []
    for comp in comp_lst:
        # get out the g(r) values of interest
        gc = get_gofr_by_plane_cps(comp[0],conn)
        temps = get_gofr_by_plane_tmp(comp[0],conn)

        # pull out the frames of interest
        gc = gc[:n]
        temp = np.mean(temps[:n])
        
        # find the location of the first peak
        max_indx = [np.argmax(g.y[5:])+5 for g in gc]
        betas = [fit_quad_to_peak(g.x[m-wind:m+wind],g.y[m-wind:m+wind]) for m,g in zip(max_indx,gc)]

        max_vals = [np.max(g.y[5:]) for g in gc]
        #points.append((np.mean(np.array([b.beta[2] for b in betas]))-1,temp))
        for b,t in zip(betas,temps):
            points.append((b.beta[2]-1,t))
        #points2.append((np.mean(max_vals)-1,temp))
        for m,t in zip(max_vals,temps):
            points2.append((m-1,t))
        
    points.sort(key=lambda x:x[1])
    points2.sort(key=lambda x:x[1])
    p,t = zip(*points)
    p2,t2 = zip(*points2)
    
    # make plot
    if g1_fig is  None:
        g1_fig = plots.tac_figure('Temperature [C]','$g_1 -1$','')
        
        
    
    g1_fig.plot(t2,p2,'-o',label=lab)
    g1_fig.plot(t,p,'-x',label=lab)
    
    return g1_fig
            


def make_gofr_by_plane_plots(comp,conn,ax=None,*args,**kwargs):
    istatus = plots.non_i_plot_start()
    if ax is None:
        (fig,ax) = plots.set_up_plot()
    else:
        fig = ax.get_figure()
    if 'Tc' in kwargs:
        T_conv_fun = ltc.T_to_phi_factory(kwargs['Tc'],ltc.linear_T_to_r_factory(-.011,0.848))
        del kwargs['Tc']
        xlab = r'$\phi^\prime$'
    else:
        T_conv_fun = lambda x:x
        xlab = 'T[C]'
    gc = get_gofr_by_plane_cps(comp,conn)
    temps = get_gofr_by_plane_tmp(comp,conn)
    dset = conn.execute("select dset_key from gofr_by_plane where comp_key = ?",(comp,)).fetchone()[0]
    (dtype,fname) = conn.execute("select dtype,fname from dsets where dset_key = ?",(dset,)).fetchone()

    wind = 30
    max_indx = [np.argmax(g.y[5:])+5 for g in gc]
    betas = [fit_quad_to_peak(g.x[m-wind:m+wind],g.y[m-wind:m+wind]) for m,g in zip(max_indx,gc)]
    ax.plot(T_conv_fun(np.array(temps)),np.array([b.beta[2] for b in betas])-1,'x',*args,**kwargs)
    
    ax.set_xlabel(xlab)
    ax.set_ylabel(r'$g_1-1$')
    ax.set_ylim([0, 3])
    #    ax.set_title("dset: " + str(dset) + " temp: " + str(temp) + "C  dtype:" + dtype)
    ax.set_title("dset: " + str(dset) + " " + fname.split('/')[-1].replace('_','\_'))
    

    (fig,ax) = plots.set_up_plot()
    plane_num = range(0,len(betas))
    ax.plot(plane_num,temps)
    ax.set_xlabel('comp index')
    ax.set_ylabel(r'temperature')
    ax.set_title("dset: " + str(dset) + " " + fname.split('/')[-1].replace('_','\_'))

    # (fig,ax) = plots.set_up_plot()
    # [ax.plot(g.x,g.y) for g in gc]
    
    plots.non_i_plot_stop(istatus)


def make_gn_by_plane_plots(comp_lst,conn):
    '''
    no idea why I wrote this
    '''
    istatus = plots.non_i_plot_start()
    (fig,ax) = plots.set_up_plot()
    for comp in comp_lst:
        gc = get_gofr_by_plane_cps(comp[0],conn)
        temps = get_gofr_by_plane_tmp(comp[0],conn)

        dset = conn.execute("select dset_key from gofr_by_plane where comp_key = ?",
                            comp).fetchone()
        (dtype,fname) = conn.execute("select dtype,fname from dsets where dset_key = ?",
                                     dset).fetchone()

        wind = 30
        max_indx = [np.argmax(g.y[5:])+5 for g in gc]
        betas = [fit_quad_to_peak(g.x[m-wind:m+wind],g.y[m-wind:m+wind]) for m,g in zip(max_indx,gc)]
        temps = range(0,len(betas))
        ax.step(temps,[b.beta[2] for b in betas],label='fit')
        ax.step(temps,[np.max(g.y) for g in gc],label='max')
    #    ax.set_title("dset: " + str(dset) + " temp: " + str(temp) + "C  dtype:" + dtype)
    ax.set_title("dset: " + str(dset) + " " + fname.split('/')[-1])
    
    plots.non_i_plot_stop(istatus)

def get_gofr_by_plane_tmp(comp_num,conn):
    (fname,) = conn.execute("select fout from gofr_by_plane where comp_key = ?",(comp_num,)).fetchone()
    
    F =  h5py.File(fname)
    g = F['gofr_by_plane_%(#)07d'%{'#':comp_num}]

    temps = [g[c].attrs['temperature'] for c in g]

    
    del g
    F.close()
    return temps
def get_gofr_by_plane_cps(comp_num,conn):
    fname = conn.execute("select fout from gofr_by_plane where comp_key = ?",(comp_num,)).fetchall()
    fname = fname[-1][0]
    F =  h5py.File(fname)
    g = F['gofr_by_plane_%(#)07d'%{'#':comp_num}]

    g_l = [cord_pairs(g[c]['bin_edges'][:],g[c]['bin_count'][:]) for c in g]

    
    del g
    F.close()
    return g_l


def plot_bp_rho(comp,ax,conn):
    ''' Extracts the rho value from gofr_by_plane computations and plots the results '''
    
    (dset_key,fname) = conn.execute("select dset_key,fout from gofr_by_plane where" +
                                    " comp_key = ?",comp).fetchone()

    
    cmap = plots.color_mapper(26,31)


    
    F = h5py.File(fname,'r')
    g = F[gen.fd('gofr_by_plane',comp[0])]
    rho = [g[c].attrs['rho'] for c in g]
    temp = np.mean(np.array([g[c].attrs['temperature'] for c in g]))
    
    ax.plot(rho,label=str(dset_key),color=cmap.get_color(temp))

    del g
    F.close()
    del F
    
    return (np.mean(rho),temp)

def plot_rho_lst(comp_lst,ax,conn,*args,**kwargs):
    '''
    Plots rho(T) from a list of g(r) plots
    '''

    (rho,temp) = zip(*[get_rhoT(c,conn) for c in comp_lst])
    print kwargs
    ax.plot(temp,rho,**kwargs)

    plts.draw()

###################
#2D v 3D functions#
###################
                   
def make_gofr_2Dv3D(conn):
    keys = conn.execute("select key from dsets where dtype = 'z'").fetchall()
    for k in keys:
        plots.make_2dv3d_plot(k[0],conn,'figures/2Dv3D/%(#)02d.png'%{"#":k[0]})

    






    
def make_2dv3d_plot(key,conn,fname = None):
    # add error handling on all of these calls
    
    # get comp_number of gofr
    res = conn.execute("select comp_key,fout from comps where dset_key == ?\
    and function == 'gofr'",(key,)).fetchall()
    (g_ck,g_fname) = res[-1]
        
        
    # get comp_number of 3D gofr
    res = conn.execute("select comp_key,fout from comps where dset_key == ?\
    and function == 'gofr3D'",(key,)).fetchall()
    (g3D_ck,g3D_fname) = res[-1]

    # get dset name
    (sample_name, temp) = conn.execute("select sname,temp from dsets where key == ? ",(key,)).fetchone()

    print sample_name + " " + str(temp)
    group = gen.get_gofr_group(g_fname,'gofr',g_ck)
    group3D = gen.get_gofr_group(g3D_fname,'gofr3D',g3D_ck)


    # make plot
    istatus = plts.isinteractive();
    if istatus:plts.ioff()
    
    dset_names = ['bin_count', 'bin_edges']
    fig = plts.figure()
    ax = fig.add_axes([.1,.1,.8,.8])
    ax.hold(True)
    ax.plot(group[dset_names[1]][:]*6.45/60,group[dset_names[0]])
    ax.plot(group3D[dset_names[1]],group3D[dset_names[0]])

    # finds the location of the maximum, assume to be the first peak
    d0 = group3D[dset_names[1]][np.argmax(group3D[dset_names[0]])]
    print np.argmax(group3D[dset_names[0]])
    print d0
    _draw_gofr_hex_lines(ax,d0/2)
    ax.set_title(sample_name + " temp: " + str(temp))
    ax.set_ylim(0,3)
    ax.legend(['2D','3D'])
    ax.grid(True)

    # save figure
    if not fname == None:
        fig.savefig(fname)
     
        
    if istatus:
        plts.ion()
        plts.show()
    else:
        plts.close(fig)
        
        

    


######
#s(q)#
######

def plot_sofq(comp_key,conn,ax=None,cmap=None):
    
    r = conn.execute("select comp_key,fout\
    from gofr where gofr.comp_key = ?",(comp_key,)).fetchone()
    
    
    
    
    temp = get_gofr_tmp(r[1],r[0],conn)
    

    (gofr,rho) = get_gofr2D_rho(r[0],conn)
    wind = 30
    indx = np.argmax(gofr.y[15:])+15
    pfit = fit_quad_to_peak(gofr.x[indx-wind:indx+wind],gofr.y[indx-wind:indx+wind])
    
    print pfit.beta[1]
    print rho
    
    q_vec = 2*np.pi *np.linspace(.2,5 ,500)

    # kludge to deal with erroneous particles 
    gofr.y[0:5] = 0
    

    S = compute_sofq(gofr,rho,q_vec)

    if ax is None:
        fig,ax = plots.set_up_plot()
        plots.add_labels(ax,'',r'$k\sigma/2\pi$','$S(k)$')
        
    
    ax.plot(pfit.beta[1]*q_vec/(2*np.pi),S,label='%.2f'%temp)
    
    

    return ax

def plot_sofq_series(comp_list,conn,cmap=None):


    fig,ax = plots.set_up_plot()
    plots.add_labels(ax,'',r'$k\sigma/2\pi$','$S(k)$')

    count = len(comp_list)
    cmap = cm.get_cmap('winter')
    ax.set_color_cycle([cmap(j/(count-1)) for j in range(count)] )
        
    for c in comp_list:
        plot_sofq(c[0],conn,ax)

    ax.legend(loc=0)
    
        

def get_max_sofq_val(comp_list,conn):
    """Returns the maximum value of s(q)"""
    res = [conn.execute("select comp_key,fout\
    from gofr where comp_key = ?",comp_key).fetchone()
           for comp_key in comp_list]
       
    
    temps = [get_gofr_tmp(r[1],r[0],conn) for r in res]
    
    
    g_r = [get_gofr2D_rho(r[0],conn) for r in res]
    
        
    q_vec = 2*np.pi *np.linspace(.2,5 ,500)

    S = [compute_sofq(gofr,rho,q_vec) for gofr,rho in g_r]
        
    return [np.max(s)for s in S],temps


def get_max_sofq_q(comp_list,conn):
    """Returns the maximum value of s(q)"""
    res = [conn.execute("select comp_key,fout\
    from gofr where comp_key = ?",comp_key).fetchone()
           for comp_key in comp_list]
       
    
    temps = [get_gofr_tmp(r[1],r[0],conn) for r in res]
    
    
    g_r = [get_gofr2D_rho(r[0],conn) for r in res]
    
        
    q_vec = 2*np.pi *np.linspace(.2,5 ,500)

    S = [compute_sofq(gofr,rho,q_vec) for gofr,rho in g_r]
        
    return [q_vec[np.argmax(s)] for s in S],temps


    
def compute_sofq(gofr,rho,q_vec):
    '''Computes the structure factor from the 
    see PRUE 81 041305

    s(q) = 1 + \rho \tilde{h}(q)
    \tilde{h}(q) = (2 \pi)^{d/2} \int_0^\infty r^{d-1} h(r) \frac{J_{d/2}}{(qr)^{d/2}} dr

    which for a 3D sample,

    s(q) = 1 + \rho 4\pi \int_0^\infty \frac{r}{k} \sin(kr) h(r)
    
    gofr : a cord_pair object
    q_vec : an iterable of q values to compute the transform at
    '''

    
    r = gofr.x
    h = gofr.y-1
    dr = np.diff(r)
    dr = np.append(dr,np.mean(dr))
    # kludge to fix issue with doubled particles
    #h[0:5] = 0


    
    S = [ 1 + (rho / q) *4* np.pi * np.sum(r * np.sin(q*r) * h *dr) for q in q_vec]

    return S

def plot_s1_series(c_lst,conn,ax=None,label=None,**kwargs):
    '''Plots the s_1 series on to the given axis '''
    if 'Tc' in kwargs:
        T_conv_fun = ltc.T_to_phi_factory(kwargs['Tc'],ltc.linear_T_to_r_factory(-.011,0.848))
        del kwargs['Tc']
        xlab = r'$\phi^\prime$'
    else:
        T_conv_fun = lambda x:x
        xlab = 'T[C]'
    ylab = r'$s_1-1$'
    
    (s1,T) = get_max_sofq_val(c_lst,conn)

    if ax is None:
        fig,ax = plots.set_up_plot()
        plots.add_labels(ax,'',xlab,ylab)
        
    ax.plot(T_conv_fun(np.array(T)),np.array(s1)-1,'x-',label=label)
    ax.set_ylim(bottom=0)
    if label is not None:
        ax.legend(loc=0)
    
    plts.draw()

def plot_compressibility(comp_list,conn,ax=None,label=None,**kwargs):
    if 'Tc' in kwargs:
        T_conv_fun = ltc.T_to_phi_factory(kwargs['Tc'],ltc.linear_T_to_r_factory(-.011,0.848))
        del kwargs['Tc']
        xlab = r'$\phi^\prime$'
    else:
        T_conv_fun = lambda x:x
        xlab = 'T[C]'
    ylab = r'$KT\kappa$'

    res = [conn.execute("select comp_key,fout\
    from gofr where comp_key = ?",comp_key).fetchone()
           for comp_key in comp_list]
       
    
    T = [get_gofr_tmp(r[1],r[0],conn) for r in res]
    
    
    g_r = [get_gofr2D_rho(r[0],conn) for r in res]
    
      
    

    C = [compute_compressibility(gofr,rho) for gofr,rho in g_r]
       
        

    if ax is None:
        fig,ax = plots.set_up_plot()
        plots.add_labels(ax,'',xlab,ylab)
        
    ax.plot(T_conv_fun(np.array(T)),C,'x-',label=label)
    
    if label is not None:
        ax.legend(loc=0)
    
    plts.draw()

    
def compute_compressibility(gofr,rho):
    """
    see PRE 81 041305

    s(q) = 1 + \rho \tilde{h}(q)
    \tilde{h}(q) = (2 \pi)^{d/2} \int_0^\infty r^{d-1} h(r) \frac{J_{d/2}}{(qr)^{d/2}} dr

    which for a 3D sample,

    s(q) = 1 + \rho 4\pi \int_0^\infty \frac{r}{k} \sin(kr) h(r)

    s(0) = 1 ) \rho 4\pi \int_0^\infty r^2 h(r)

    taking the limit as k->0
    
    gofr : a cord_pair object
    rho : density
    """
    
    r = gofr.x
    h = gofr.y-1
    dr = np.diff(r)
    dr = np.append(dr,np.mean(dr))
    # kludge to fix issue with doubled particles
    h[0:5] = 0


    
    c = 1 + (rho ) * np.sum(r * r * h *dr)*4*np.pi

    
    ## plts.figure()
    ## plts.plot( (rho ) * r * r * h *dr *4*np.pi)
    return c

    
def make_sofq_3D_plot(key,conn,Q):
    '''From the key plots s(q) as computed from the 3D g(r)'''
    res = conn.execute("select comp_key from comps where dset_key=? and function='gofr3D'",(key,)).fetchall()
    if not len(res)==1:
        raise util.dbase_error("can't find 3D gofr")

    plts.figure()
    g = gen.get_gofr3D(res[0][0],conn)

    S = sofQ(g,Q)

    res2 = conn.execute("select comp_key from comps where dset_key=? and function='gofr'",(key,)).fetchall()
    if not len(res2)==1:
        raise util.dbase_error("can't find gofr")
    
    g2 = gen.get_gofr2D(res2[0][0],conn)
    S2 = gen.sofQ(g2,Q)


    istatus = plts.isinteractive();
    if istatus:plts.ioff()

    # plot s(q)
    leg_hands = []
    leg_strs = []
    
    fig = plts.figure()
    ax = fig.add_axes([.1,.1,.8,.8])
    ax.hold(True)
    leg_hands.append(ax.plot(Q,S))
    leg_strs.append('3D based')
    leg_hands.append(ax.plot(Q,S2))
    leg_strs.append('2D based')
    ax.legend(leg_hands,leg_strs)
    ax.set_title('S(q)')
    ax.set_xlabel(r' q [$1/\mu$m]')
    #ax.set_xlabel(r' qR')
    ax.set_ylabel('S(q) [arb units]')
    
    if istatus:
        print "displaying figure"
        plts.show()
        plts.ion()
    
    else:
        print "closing figure"
        plts.close(fig)
        plts.close(gn_fig)






###########
# fitting # 
###########

def _trim_gofr(gofr,r0):
    """     Takes in gofr as a cord_pairs and r0 in
    units of 1/d, d as determined by the location
    of the first peak.

    Returns a cords_pair object with g(r>r0) 
    """
    # find the first peak
    d0 = gofr.x[np.argmax(gofr.y[15:])+15]

    print "the Diameter is " + str(d0)

    # convert r0 to real units
    r0 *=d0

    
    
    # cut the data at r0
    r_arg = np.argmax([v for v in gofr.x if v<= r0])
    gofr_trim = cord_pairs(gofr.x[r_arg:],gofr.y[r_arg:])



    return gofr_trim,r0

def fit_gofr2(gofr,r0,func,p0=(1,7,2,0,1)):
    """
    Fits to 
    Takes in gofr as a cord_pairs and r0 in
    units of d, d as determined by the location
    of the first peak.

    Returns the tuple that is the argument for the function
    used for the fitting, the covariance matrix, and the
    sum of the residual squared
    """
    # trim the g(r) data
    (gofr,r0) = _trim_gofr(gofr,r0)
    
    return fitting.fit_curve(gofr.x,gofr.y,p0,func)

def fit_gofr3(gofr,r0,gen_func,p0=(2,7,1,.7)):
    (gofr,r0) = _trim_gofr(gofr,r0)

    return fitting.fit_curve(gofr.x,gofr.y,p0,gen_func(r0)),r0




###################
#format correction#
###################
def remove_gofr_computation(comp_number,conn):
    (f_gofr,) = conn.execute("select fout from gofr where comp_key = ?"
                                   ,comp_number).fetchone()

    # the order is important to keep the foreign constraints happy
    # kill gofr_prams entry
    conn.execute("delete from gofr where comp_key = ?",comp_number)
    # kill comps entry
    conn.execute("delete from comps where comp_key = ?",comp_number)
    # commit to db, commit before deleting the data as unmarked data is less irritating
    # than non-existing data
    conn.commit()
    

    
    # remove group from hdf file
    F = h5py.File(f_gofr,'r+')
    del F["gofr_%(#)07d"%{"#":comp_number[0]}]
    F.close()
    del F
    
    pass

def remove_gofr_bp_computation(comp_number,conn):
    (f_gofr,) = conn.execute("select fout from gofr_by_plane where comp_key = ?"
                                   ,comp_number).fetchone()

    # the order is important to keep the foreign constraints happy
    # kill gofr_prams entry
    conn.execute("delete from gofr_by_plane where comp_key = ?",comp_number)
    # kill comps entry
    conn.execute("delete from comps where comp_key = ?",comp_number)
    # commit to db, commit before deleting the data as unmarked data is less irritating
    # than non-existing data
    conn.commit()
    

    
    # remove group from hdf file
    F = h5py.File(f_gofr,'r+')
    del F["gofr_by_plane_%(#)07d"%{"#":comp_number[0]}]
    F.close()
    del F
    
    pass


def fix_temperature(comp_key,conn):
    '''Fixes the temperature meta data (for exmaple, if you forget to
    set it in the Iden* files)'''

    (iden_key,g_fname) = conn.execute("select iden_key,fout from gofr where comp_key = ?",
                                      (comp_key,)).fetchone()
    (i_fname,) = conn.execute("select fout from iden where comp_key = ?",(iden_key,)).fetchone()
    print iden_key,i_fname
    
    F_iden = h5py.File(i_fname,'r')
    temp_sum = 0
    frame_count = 0
    for fr in F_iden:
        if fr[0:5] == 'frame':
            temp_sum += F_iden[fr].attrs['temperature']
            frame_count += 1




    F_iden.close()
    del F_iden

    t = temp_sum/frame_count;
    print t
    F = h5py.File(g_fname,'r+')
    grp = F["gofr_%(#)07d"%{"#":comp_key}]
    if 'temperature' in  grp.attrs:
        del grp.attrs['temperature']
    
    grp.attrs.create('temperature',t,None,'float32')
    F.close()
    return t



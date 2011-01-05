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
import pov 
import plots 
import util 
import numpy as np
from util import cord_pairs
import scipy.optimize as sopt
import scipy.odr as sodr
import bisect
import itertools
import matplotlib.pyplot as plts
import matplotlib.cm as cm
import fitting
import general as gen
import matplotlib
##################
#helper functions#
##################
def find_peaks_fit(gp,dfun,p):
    """Looks for peaks based on hints from the derivative of the expected function and
    the fit parameters"""

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
            print 'diff from expected ' + str(crit_p - (cur_pos + e_spc/4))
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
    bins = np.array(g[dset_names[1]])*6.45/60
    F.close()
    return cord_pairs(bins,gofr)

def get_gofr2D_rho(comp_num,conn):
    '''Takes in computation number and database connection and extracts the given g(r) and returns it as a
    cord_pairs object (left bin edge,value)'''

    dset_names = ['bin_count', 'bin_edges']
    
    res = conn.execute("select fout from comps where comp_key = ?",(comp_num,)).fetchall()
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

def get_gofr_tmp(fname,comp_num,conn):
    F = h5py.File(fname,'r')
    if 'temperature' in  F["gofr_%(#)07d"%{"#":comp_num}].attrs:
        t = F["gofr_%(#)07d"%{"#":comp_num}].attrs['temperature']
        print 't from file',t
    else:
        (t,) = conn.execute("select temp from dsets " +
                            "where key in (select dset_key from comps where comp_key = ? )",
                            (comp_num,)).fetchone()
        print ' t from db',t
    F.close()
    return t




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

    res = conn.execute("select comp_key,fout\
    from comps where comps.comp_key = ?",(g_key,)).fetchone()
    

    
    
    temps = get_gofr_tmp(res[1],res[0],conn) 
    print temps

    g = get_gofr2D(res[0],conn) 
    fits = fit_gofr2(g,2.1,fitting.fun_decay_exp_inv) 
    
    
    istatus = plots.non_i_plot_start()

    fig = plots.Figure('r[$\mu m$]','g(r)','g(r) + fitting temp: %.2f'%temps,func = matplotlib.axes.Axes.step) 
    
    
    fig.plot(g.x,g.y,label='data')
    fig.plot(g.x[25:],fitting.fun_decay_exp_inv(fits.beta,g.x[25:]),label='fitting')
    fig.axis((0,12),(0,3))
    plots.non_i_plot_stop(istatus)


###########################
#function of tmp functions#
###########################
def tmp_series_gn2D(comp_list,conn):
    """Makes plots for g_1 to g_n
    
    comp_list : list of 1 element tuple that contain the computation
    number to be plotted
    """
    
    res = [conn.execute("select comp_key,fout\
    from comps where comps.comp_key = ?",c).fetchone()
           for c in comp_list]

    
    
    temps = [get_gofr_tmp(r[1],r[0],conn) for r in res]
    print temps
    
    gofrs = [get_gofr2D(r[0],conn) for r in res]
    fits = [fit_gofr2(g,2.1,fitting.fun_decay_exp_inv) for g in gofrs]
    peaks = [find_peaks_fit(g,fitting.fun_decay_exp_inv_dr,p.beta)
             for g,p in zip(gofrs,fits)]


    
    istatus = plots.non_i_plot_start()
    
    # make g_n plots for peaks
    fig,ax = plots.set_up_plot()
    leg_strs = []
    leg_hands = []
    for j in range(max([len(p[0]) for p in peaks])):
        pairs = [(t,p[0][j][1]) for (t,p) in itertools.izip(temps,peaks) if len(p[0]) > j]
        leg_hands.append(ax.plot([p[0] for p in pairs],[p[1]-1 for p in pairs],'x-'))
        leg_strs.append('$g_'+str(j)+'$' )

    #ax.legend(leg_hands,leg_strs,loc=0)
    plots.add_labels(ax,'$g_n$ maximums',r'T [C]','g(peak) -1')


    
    # make g_n plots for troughs
    fig,ax = plots.set_up_plot()
    leg_strs = []
    leg_hands = []
    for j in range(max([len(p[1]) for p in peaks])):
        pairs = [(t,p[1][j][1]) for (t,p) in itertools.izip(temps,peaks) if len(p[1]) > j]
        leg_hands.append(ax.plot([p[0] for p in pairs],[p[1]-1 for p in pairs],'x-'))
        leg_strs.append('$g_'+str(j)+'$' )

    #ax.legend(leg_hands,leg_strs,loc=0)
    plots.add_labels(ax,'$g_n$ minimums',r'T [C]','g(peak) -1')
    

    # make plot for trough locations
    fig,ax = plots.set_up_plot()
    leg_strs = []
    leg_hands = []
    for j in range(max([len(p[1]) for p in peaks])):
        pairs = [(t,p[1][j][0]) for (t,p) in itertools.izip(temps,peaks) if len(p[1]) > j]
        leg_hands.append(ax.plot([p[0] for p in pairs],[p[1] for p in pairs],'x-'))
        leg_strs.append('$g_'+str(j)+'$' )

    #ax.legend(leg_hands,leg_strs,loc=0)
    plots.add_labels(ax,'minimum locations',r'T [C]','r [$\mu m$]')


    # make plot for peak locations
    fig,ax = plots.set_up_plot()
    leg_strs = []
    leg_hands = []
    for j in range(max([len(p[0]) for p in peaks])):
        pairs = [(t,p[0][j][0]) for (t,p) in itertools.izip(temps,peaks) if len(p[0]) > j]
        leg_hands.append(ax.plot([p[0] for p in pairs],[p[1] for p in pairs],'x-'))
        leg_strs.append('$g_'+str(j)+'$' )

    #ax.legend(leg_hands,leg_strs,loc=0)
    plots.add_labels(ax,'peak locations',r'T [C]','r [$\mu m$]')

    
    plots.non_i_plot_stop(istatus)
    
    return peaks


def tmp_series_fit_plots(comp_list,conn):
    """Makes plots for g_1 to g_n
    
    comp_list : list of 1 element tuple that contain the computation
    number to be plotted
    """

    res = [conn.execute("select comp_key,fout\
    from gofr where comp_key = ?",c).fetchone()
           for c in comp_list]


    
    temps = [get_gofr_tmp(r[1],r[0],conn) for r in res]
    zip_lst = zip(temps,res)
    zip_lst.sort(key=lambda x:x[0])
    temps,res = zip(*zip_lst)
    
    print temps

    gofrs = [get_gofr2D(r[0],conn) for r in res]
    fits = [fit_gofr3(g,2.1,fitting.fun_decay_exp_inv_gen) for g in gofrs]
    fits,r0s = zip(*fits)
    
    istatus = plots.non_i_plot_start()

    fig_xi = plots.Figure('T [C]',r'fitting parameters [$\mu m$]','Fitting parameters') 

    print [2*np.pi/p.beta[1] for p in fits ]
    
    fig_xi.plot(temps,[p.beta[0] for p in fits ],
                    marker='x',label=r'$\xi$')
    fig_xi.plot(temps,[(np.pi*2)/p.beta[1] for p in fits ],
                    marker='x',label=r'$K$')
    fig_xi.plot(temps,[p.beta[2] for p in fits ],
                    marker='x',label=r'$C [NONE]$')
    fig_xi.plot(temps,r0s,
                    marker='x',label=r'$r_0$')
    
    
    fig_err = plots.Figure('T [C]','errors [diff units]','fitting error')
    fig_err.plot(temps,[p.sum_square for p in fits ],'-x',label=r'sum_square')
    fig_err.plot(temps,[p.sum_square_delta for p in fits ],'-x',label=r'sum_square_delta')
    fig_err.plot(temps,[p.res_var for p in fits ],'-x',label=r'res_var')

    fig_b = plots.Figure('T [C]',r'$b-1$','$g(r)$ shift') 
    fig_b.plot(temps,[p.beta[4]-1 for p in fits ],
                    marker='x',label=r'$b-1$')


    
    plots.non_i_plot_stop(istatus)


    print [(np.pi*2)/p.beta[1] for p in fits ]

def make_gofr_tmp_series(comp_list,conn,fnameg=None,fnamegn=None,gtype='gofr',date = None,g1_fig = None,T_correction = 0):
    '''Takes in a sample name and plots all of the g(r) for it '''
    r_scale = 6.45/60
    dset_names = ['bin_count', 'bin_edges']

    
    res = [conn.execute("select comp_key,fout,fin,dset_key from gofr where comp_key = ? ",c).fetchone()
           for c in comp_list]

    #res.sort(key=lambda x:x[2])

    
    (sname,) = conn.execute("select sname from dsets where dset_key = ?",(res[0][3],)).fetchone()

    print "there are " + str(len(res)) + " entries found"
    # check interactive plotting and turn it off
    istatus = plts.isinteractive();
    print istatus
    if istatus:plts.ioff()

    leg_hands = []
    leg_strs = []

    fig = plts.figure()
    ax = fig.add_axes([.1,.1,.8,.8])
    ax.hold(True)
    ax.grid(True)
    
    gn_g = []
    gn_t = []
    gn_p = []

    if g1_fig is None:
        gn_fig = plts.figure()
        gn_ax = gn_fig.add_axes([.1,.1,.8,.8])
    else:
        gn_fig = g1_fig
        gn_ax = g1_fig.gca()
        
    temps = [get_gofr_tmp(r[1],r[0],conn)+T_correction for r in res]
    tmax = max(temps)
    tmin = min(temps)
    ax.set_color_cycle([cm.jet((t-tmin)/(tmax-tmin)) for t in temps])
    for r,temperature in zip(res,temps):


        
        F = h5py.File(r[1],'r')
        print gtype + "_%(#)07d"%{"#":r[0]}
        print gtype + "_%(#)07d"%{"#":r[0]} in F
        
        g = F[gtype + "_%(#)07d"%{"#":r[0]}]
        
        
        leg_hands.append(ax.plot(g[dset_names[1]][:]*r_scale,g[dset_names[0]][:]-1))
        leg_strs.append(str(np.round(temperature,2)))
        g1 = np.max(g[dset_names[0]])
        
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
            
        del g
        F.close()
        del F
    gn_p.sort(lambda x,y: int(np.sign(x[0]-y[0])))
    for p in gn_p:
        gn_t.append(p[0])
        gn_g.append(p[1])
    print gn_t
    print gn_g
    fig.legend(leg_hands,leg_strs,loc=0)
    ax.set_title(sname)
    ax.set_xlabel(r'r [$\mu m$]')
    ax.set_ylabel(r'$G(r)-1$')
    
    
    gn_ax.plot(gn_t,np.array(gn_g)-1,'x-')
    gn_ylim = gn_ax.get_ylim()
    gn_ylim[0] = 0
    gn_ax.set_ylim(gn_ylim)
        
    print 
    gn_ax.grid(True)
    gn_ax.set_title(sname + r' $g_1(T)$')
    gn_ax.set_xlabel('T')
    gn_ax.set_ylabel('$g_1-1$')
    if g1_fig is not None:
        plts.draw()

    if fnameg is not None:
        fig.savefig(fnameg)

    if  fnamegn is not None:
        gn_fig.savefig(fnamegn)

        
    if istatus:
        print "displaying figure"
        plts.ion()
        plts.show()
    else:
        print "closing figure"
        plts.close(fig)
        plts.close(gn_fig)

    return zip(gn_t,gn_g)

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

def make_gofr_by_plane(d_lst,conn):
    istatus = plots.non_i_plot_start()
    (fig,ax) = plots.set_up_plot()
    leg_strs = []
    leg_hand = []
    cmap = color_mapper(np.min([d[1] for d in d_lst]),np.max([d[1] for d in d_lst]))
    for d in d_lst:
        dkey = d[0]
        ckeys = conn.execute("select comp_key from comps where dset_key = ? and function = 'gofr_by_plane'",(dkey,)).fetchall()
        for c in ckeys:
            print c
            gc = gen.get_gofr_by_plane_cps(c[0],conn)
            leg_hand.append(ax.plot([np.max(g.y) for g in gc],color = cmap.get_color(d[1])))
            leg_strs.append(str(d[1]))
            
    ax.legend(leg_hand,leg_strs)
    plots.non_i_plot_stop(istatus)

def make_gofr_by_plane_plots(comp,conn):
    istatus = plots.non_i_plot_start()
    (fig,ax) = plots.set_up_plot()
    
    gc = gen.get_gofr_by_plane_cps(comp,conn)
    temps = gen.get_gofr_by_plane_tmp(comp,conn)
    dset = conn.execute("select dset_key from comps where comp_key = ?",(comp,)).fetchone()[0]
    (temp,dtype,fname) = conn.execute("select temp,dtype,fname from dsets where key = ?",(dset,)).fetchone()

    wind = 30
    max_indx = [np.argmax(g.y[5:])+5 for g in gc]
    betas = [gen.fit_quad_to_peak(g.x[m-wind:m+wind],g.y[m-wind:m+wind]) for m,g in zip(max_indx,gc)]
    ax.step(temps,[b.beta[2] for b in betas])
    ax.step(temps,[np.max(g.y) for g in gc])
    #    ax.set_title("dset: " + str(dset) + " temp: " + str(temp) + "C  dtype:" + dtype)
    ax.set_title("dset: " + str(dset) + " " + fname.split('/')[-1])
    
    plots.non_i_plot_stop(istatus)

def get_gofr_by_plane_tmp(comp_num,conn):
    (fname,) = conn.execute("select fout from comps where comp_key = ?",(comp_num,)).fetchone()
    
    F =  h5py.File(fname)
    g = F['gofr_by_plane_%(#)07d'%{'#':comp_num}]

    temps = [g[c].attrs['temperature'] for c in g]

    
    del g
    F.close()
    return temps
def get_gofr_by_plane_cps(comp_num,conn):
    fname = conn.execute("select fout from comps where comp_key = ?",(comp_num,)).fetchall()
    fname = fname[-1][0]
    F =  h5py.File(fname)
    g = F['gofr_by_plane_%(#)07d'%{'#':comp_num}]

    g_l = [cord_pairs(g[c]['bin_edges'][:],g[c]['bin_count'][:]) for c in g]

    
    del g
    F.close()
    return g_l
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

def plot_sofq(comp_key,conn,fig_sofq=None,length=None,cmap=None):
    
    r = conn.execute("select comp_key,fout\
    from comps where comps.comp_key = ?",(comp_key,)).fetchone()
    
    
    
    
    temp = get_gofr_tmp(r[1],r[0],conn)
    

    (gofr,rho) = get_gofr2D_rho(r[0],conn)
    wind = 30
    indx = np.argmax(gofr.y[15:])+15
    pfit = fit_quad_to_peak(gofr.x[indx-wind:indx+wind],gofr.y[indx-wind:indx+wind])
    
    print pfit.beta[1]
    print rho
    
    q_vec = 2*np.pi *np.linspace(.2,5 ,500)

    S = compute_sofq(gofr,rho,q_vec)

    if fig_sofq is None:
        if length is None:
            fig_sofq = plots.Figure(r'$k\sigma/2\pi$',r'$S(k)$','test',cmap=cmap)
        else:
            fig_sofq = plots.Figure(r'$k\sigma/2\pi$',r'$S(k)$','test',count=length,cmap=cmap)
    
    fig_sofq.plot(pfit.beta[1]*q_vec/(2*np.pi),S,label='%.2f'%temp)
    
    

    return fig_sofq

def plot_sofq_series(comp_list,conn,cmap=None):


    c0 = comp_list.pop(0)
    fig = plot_sofq(c0[0],conn,None,len(comp_list)+1,cmap=cmap)
    for c in comp_list:
        plot_sofq(c[0],conn,fig)

        
    comp_list.reverse()
    comp_list.append(c0)
    comp_list.reverse()

def get_max_sofq_loc(comp_list,conn):
    
    res = [conn.execute("select comp_key,fout\
    from comps where comps.comp_key = ?",comp_key).fetchone()
           for comp_key in comp_list]
    
    
    
    
    temps = [get_gofr_tmp(r[1],r[0],conn) for r in res]
    

    g_r = [get_gofr2D_rho(r[0],conn) for r in res]
    
        
    q_vec = 2*np.pi *np.linspace(.2,5 ,500)

    S = [compute_sofq(gofr,rho*.25,q_vec) for gofr,rho in g_r]
    

    
    #    return [q_vec[np.argmax(s)]/(2*np.pi) for s in S],temps
    return [np.max(s)for s in S],temps
    
def compute_sofq(gofr,rho,q_vec):
    '''Computes the structure factor from the 
    see PRE 81 041305

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


    
    S = [ 1 + (rho / q) * np.sum(r * np.sin(q*r) * h *dr) for q in q_vec]

    return S
    
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

def fit_gofr3(gofr,r0,gen_func,p0=(1,7,2,0,1)):
    (gofr,r0) = _trim_gofr(gofr,r0)

    return fitting.fit_curve(gofr.x,gofr.y,p0,gen_func(r0)),r0
# kill?

def gn_type_plots(sname,conn):

    istatus = plots.non_i_plot_start()
    
    # make g_n plots for peaks
    fig,ax = plots.set_up_plot()
    leg_strs = []
    leg_hands = []


    for t in ['t','z']:
        res = conn.execute("select comps.comp_key,dsets.temp\
        from comps,dsets where comps.dset_key = dsets.key and comps.function=? and\
        dsets.sname = ? and dtype = ?",('gofr',sname,t,)
                       ).fetchall()


        temps = [r[1] for r in res]
        gofrs = [gen.get_gofr2D(r[0],conn) for r in res]
        fits = [fit_gofr2(g,2.1,fitting.fun_decay_exp_inv) for g in gofrs]
        peaks = [gen.find_peaks_fit(g,fitting.fun_decay_exp_inv_dr,p.beta)
                 for g,p in zip(gofrs,fits)]


    

        for j in range(min([len(p[0]) for p in peaks])):
            leg_hands.append(ax.plot(temps,np.array([p[0][j][1] for p in peaks])-1,'x-'))
            leg_strs.append('$g_'+str(j)+'$ ' + t )



    ax.legend(leg_hands,leg_strs)
    plots.add_labels(ax,'$g_n$ maximums',r'T [C]','g(peak) -1')


###################
#format correction#
###################
def remove_gofr_computation(comp_number,conn):
    (f_gofr,) = conn.execute("select fout from comps where comp_key = ?"
                                   ,comp_number).fetchone()

    # the order is important to keep the foreign constraints happy
    # kill gofr_prams entry
    conn.execute("delete from gofr_prams where comp_key = ?",comp_number)
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

def fix_temperature(comp_key,conn):
    '''Fixes the temperature meta data (for exmaple, if you forget to
    set it in the Iden* files)'''

    (iden_key,g_fname) = conn.execute("select iden_key,fout from gofr where comp_key = ?",
                                      (comp_key,)).fetchone()
    (i_fname,) = conn.execute("select fout from iden where comp_key = ?",(iden_key,)).fetchone()

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


    

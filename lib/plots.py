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
import trackpy.cpp_wrapper as cpp_wrapper
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import random
import itertools
import h5py
import numpy as np
import math
import util
import itertools
import os
import os.path
import general as gen
import fitting



class cord_pairs:
    def __init__(self,x,y):
        self.x = x
        self.y = y
    def __iter__(self):
        return itertool.izip(x,y)

# change to take 
def _plot_file_frame_phi6(key,conn,fr_num,fnameg=None):
    '''

    '''
    (fname,p_comp) = conn.execute("select fout,comp_key from comps where function = 'phi6' and dset_key = ?;",(key,)).fetchone()
    (comp_number,) = conn.execute("select comp_key from comps where function = 'Iden' and fout = ?",(fname,)).fetchone()

    (sname,stype,temp) = conn.execute("select sname,dtype,temp from dsets where key = ?",(key,)).fetchone()

    
    f = h5py.File(fname,'r')

    fc = f.attrs['number-of-planes']
    if fr_num <fc:
        fr_num = fc-1
    
    x = f["/frame%(#)06d"%{"#":fr_num}+"/x_%(#)07d"%{"#":comp_number}]
    y = f["/frame%(#)06d"%{"#":fr_num}+"/y_%(#)07d"%{"#":comp_number}]
    phi = f["/frame%(#)06d"%{"#":fr_num}+"/scaler_order_parameter_%(#)07d"%{"#":p_comp}]
    phir = np.zeros(phi.shape[0])
    phii = np.zeros(phi.shape[0])

    for j in range(0,phi.shape[0]):
        phir[j] = phi[j][0]
        phii[j] = phi[j][1]
    print np.mean(abs(phir))
    print np.mean(abs(phii))
    print np.mean((phir))
    print np.mean((phii))

    
    # check interactive plotting and turn it off
    istatus = plt.isinteractive();
    print istatus
    if istatus:plt.ioff()

    leg_hands = []
    leg_strs = []

    fig = plt.figure()
    ax = fig.add_axes([.1,.1,.8,.8])
    ax.hold(True)
    ax.set_aspect('equal')

    ax.quiver(x,y,phir,phii,scale = 100,headlength = 2,headwidth= 1.5,headaxislength=2)
    ax.set_title(sname + " temp: " + str(temp))
#    plt.plot(x,y,'ro')


    
    if not fnameg == None:
        f_path = '/home/tcaswell/python/figures/' + sname + '/'
        if not os.path.isdir(f_path):
            os.makedirs(f_path)
        
        fnameg = f_path + str(key) + '_phi6.png'
        fig.savefig(fnameg,dpi = 500)
        
        
    if istatus:
        print "displaying figure"
        plt.ion()
        plt.show()
    else:
        print "closing figure"
        plt.close(fig)


def _plot_file_frame_nsize(key,conn,fr_num,fnameg=None):
    '''

    '''
    (fname,p_comp) = conn.execute("select fout,comp_key from comps where function = 'phi6' and dset_key = ?;",(key,)).fetchone()
    (comp_number,) = conn.execute("select comp_key from comps where function = 'Iden' and fout = ?",(fname,)).fetchone()

    (sname,stype,temp) = conn.execute("select sname,dtype,temp from dsets where key = ?",(key,)).fetchone()

    
    f = h5py.File(fname,'r')
    
    x = f["/frame%(#)06d"%{"#":fr_num}+"/x_%(#)07d"%{"#":comp_number}]
    y = f["/frame%(#)06d"%{"#":fr_num}+"/y_%(#)07d"%{"#":comp_number}]
    ns = f["/frame%(#)06d"%{"#":fr_num}+"/neighborhood_size_%(#)07d"%{"#":p_comp}]

    print np.mean(ns)
    
    # check interactive plotting and turn it off
    istatus = plt.isinteractive();
    print istatus
    if istatus:plt.ioff()

    leg_hands = []
    leg_strs = []

    fig = plt.figure()
    ax = fig.add_axes([.1,.1,.8,.8])
    ax.hold(True)
    #ax.set_aspect('equal')

    #    ax.scatter(x,y,c=ns)
    ax.hist(ns,bins=range(0,10))
    ax.set_title(sname + " temp: " + str(temp))
#    plt.plot(x,y,'ro')


    
    if not fnameg == None:
        f_path = '/home/tcaswell/python/figures/' + sname + '/'
        if not os.path.isdir(f_path):
            os.makedirs(f_path)
        
        fnameg = f_path + str(key) + '_phi6.eps'
        fig.savefig(fnameg)
        
        
    if istatus:
        print "displaying figure"
        plt.ion()
        plt.show()
    else:
        print "closing figure"
        plt.close(fig)




def plot_file_frame_pos(key,conn, fr_num):
    (fname,p_comp) = conn.execute("select fout,comp_key from comps where function = 'phi6' and dset_key = ? and date = '2010-04-19';",(key,)).fetchone()
    (comp_number,) = conn.execute("select comp_key from comps where function = 'Iden' and fout = ?",(fname,)).fetchone()

    (sname,stype,temp) = conn.execute("select sname,dtype,temp from dsets where key = ?",(key,)).fetchone()

    
    f = h5py.File(fname,'r')

    x = f["/frame%(#)06d"%{"#":fr_num}+"/x_%(#)07d"%{"#":comp_number}]
    y = f["/frame%(#)06d"%{"#":fr_num}+"/y_%(#)07d"%{"#":comp_number}]
    ns = f["/frame%(#)06d"%{"#":fr_num}+"/neighborhood_size_%(#)07d"%{"#":p_comp}][:]
    istatus = plt.isinteractive();
    print istatus
    if istatus:plt.ioff()

    
    fig = plt.figure()
    ax = fig.add_axes([.1,.1,.8,.8])
    ax.set_aspect('equal')
    sc = ax.scatter(x,y,c=ns)
    plt.colorbar(sc)
        
    if istatus:
        print "displaying figure"
        plt.ion()
        plt.show()
    else:
        print "closing figure"
        plt.close(fig)




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
    istatus = plt.isinteractive();
    if istatus:plt.ioff()
    
    dset_names = ['bin_count', 'bin_edges']
    fig = plt.figure()
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
        plt.ion()
        plt.show()
    else:
        plt.close(fig)
        
        


def make_gofr_tmp_series(sname,conn,fnameg=None,fnamegn=None,gtype='gofr',date = None):
    '''Takes in a sample name and plots all of the g(r) for it '''
    r_scale = 6.45/60
    dset_names = ['bin_count', 'bin_edges']
    
    res = gen.get_list_gofr(sname,conn,date = date,gtype = gtype)
    print "there are " + str(len(res)) + " entries found"
    # check interactive plotting and turn it off
    istatus = plt.isinteractive();
    print istatus
    if istatus:plt.ioff()

    leg_hands = []
    leg_strs = []

    fig = plt.figure()
    ax = fig.add_axes([.1,.1,.8,.8])
    ax.hold(True)
    ax.grid(True)
    ax.set_color_cycle([cm.jet(x/len(res)) for x in range(len(res))])
    gn_g = []
    gn_t = []
    gn_p = []

    
    gn_fig = plt.figure()
    gn_ax = gn_fig.add_axes([.1,.1,.8,.8])


    
    for r in res:


        print r[3]
        g = gen.get_gofr_group(r[1],gtype,r[0])
        leg_hands.append(ax.plot(g[dset_names[1]][:]*r_scale,g[dset_names[0]]))
        leg_strs.append(str(r[2]))
        g1 = np.max(g[dset_names[0]])
        if not np.isnan(g1):
            try:
                gn_p.append((float(r[2]),g1))
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
    fig.legend(leg_hands,leg_strs)
    ax.set_title(sname)
    ax.set_xlabel(r'r [$\mu m$]')
    ax.set_ylabel(r'G(r)')
    
    gn_ax.plot(gn_t,np.array(gn_g)-1,'x-')
    gn_ax.grid(True)
    gn_ax.set_title(sname + r' $g_1(T)$')
    gn_ax.set_xlabel('T')
    gn_ax.set_ylabel('$g_1$')


    if not fnameg == None:
        fig.savefig(fnameg)

    if not fnamegn == None:
        gn_fig.savefig(fnamegn)

        
    if istatus:
        print "displaying figure"
        plt.ion()
        plt.show()
    else:
        print "closing figure"
        plt.close(fig)
        plt.close(gn_fig)

    return zip(gn_t,gn_g)

def make_sofq_3D_plot(key,conn,Q):
    '''From the key plots s(q) as computed from the 3D g(r)'''
    res = conn.execute("select comp_key from comps where dset_key=? and function='gofr3D'",(key,)).fetchall()
    if not len(res)==1:
        raise util.dbase_error("can't find 3D gofr")

    plt.figure()
    g = gen.get_gofr3D(res[0][0],conn)

    S = sofQ(g,Q)

    res2 = conn.execute("select comp_key from comps where dset_key=? and function='gofr'",(key,)).fetchall()
    if not len(res2)==1:
        raise util.dbase_error("can't find gofr")
    
    g2 = gen.get_gofr2D(res2[0][0],conn)
    S2 = gen.sofQ(g2,Q)


    istatus = plt.isinteractive();
    if istatus:plt.ioff()

    # plot s(q)
    leg_hands = []
    leg_strs = []
    
    fig = plt.figure()
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
        plt.show()
        plt.ion()
    
    else:
        print "closing figure"
        plt.close(fig)
        plt.close(gn_fig)




def make_2d_gofr_plot(key,conn,fname = None):
    # add error handling on all of these calls
    
    # get comp_number of gofr
    res = conn.execute("select comp_key,fout from comps where dset_key == ? and function == 'gofr'",(key,)).fetchall()
    (g_ck,g_fname) = res[0]
        
    # get dset name
    (sname, temp) = conn.execute("select sname,temp from dsets where key == ? ",(key,)).fetchone()

    print sname + " " + str(temp)
    group = gen.get_gofr_group(g_fname,'gofr',g_ck)



    # make plot
    istatus = plt.isinteractive();
    if istatus:plt.ioff()
    
    dset_names = ['bin_count', 'bin_edges']
    fig = plt.figure()
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
        plt.ion()
        plt.show()
    else:
        plt.close(fig)
        
        

def plot_file_nsize_hist(key,conn,fnameg=None):
    '''

    '''
    (fname,p_comp) = conn.execute("select fout,comp_key from comps where function = 'phi6' and dset_key = ? and date = '2010-04-19';",(key,)).fetchone()
    (comp_number,) = conn.execute("select comp_key from comps where function = 'Iden' and fout = ?",(fname,)).fetchone()

    (sname,stype,temp) = conn.execute("select sname,dtype,temp from dsets where key = ?",(key,)).fetchone()

    print p_comp
    
    f = h5py.File(fname,'r')

    bins = np.array(range(0,14))-.5
    hist_cum = np.zeros((1,len(bins)-1))
    nf = f.attrs["number-of-planes"]
    for fr in range(0,nf):
        ns = f["/frame%(#)06d"%{"#":fr}+"/neighborhood_size_%(#)07d"%{"#":p_comp}][:]
        tmp = np.histogram(ns,bins,new=True)
        hist_cum += tmp[0]

    print hist_cum
    

    # check interactive plotting and turn it off
    istatus = plt.isinteractive();
    print istatus
    if istatus:plt.ioff()

    leg_hands = []
    leg_strs = []

    fig = plt.figure()
    ax = fig.add_axes([.1,.1,.8,.8])
    ax.hold(True)
    ax.grid(True)
    #ax.set_aspect('equal')

    #    ax.scatter(x,y,c=ns)
    ax.step(bins[:-1],hist_cum.T)
    ax.set_title(sname + " temp: " + str(temp))
#    plt.plot(x,y,'ro')


    
    if not fnameg == None:
        f_path = '/home/tcaswell/python/figures/' + sname + '/'
        if not os.path.isdir(f_path):
            os.makedirs(f_path)
        
        fnameg = f_path + str(key) + '_nn_hist.eps'
        fig.savefig(fnameg)
        
        
    if istatus:
        print "displaying figure"
        plt.ion()
        plt.show()
    else:
        print "closing figure"
        plt.close(fig)



def try_fits(dset_key,conn):
    """Try's a couple of fits and plots the results """

    # get out the computation number
    res = conn.execute("select comp_key from comps where function = 'gofr3D' and dset_key = ?",(dset_key,)).fetchall()
    if not len(res) == 1:
        raise "die"

    # get gofr
    gofr = gen.get_gofr3D(res[0][0],conn)
    gofr = fitting._trim_gofr(gofr,.2)
    
    # fits
   
    (p_out1_2,cov1_2,err1_2) = fitting.fit_gofr(gofr,2,fitting.fun_decay_exp_inv,(2,7.35,1.5,0,0,0))
    (p_out2_2,cov2_2,err2_2) = fitting.fit_gofr(gofr,2,fitting.fun_decay_exp,(1.5,7.35,1.5,0,0,0))

    
    # plots
    

    
    # check interactive plotting and turn it off
    istatus = plt.isinteractive();
    print istatus
    if istatus:plt.ioff()

    leg_hands = []
    leg_str = []

    fig = plt.figure()
    ax = fig.add_axes([.1,.1,.8,.8])
    ax.hold(True)
    ax.grid(True)
    #ax.set_aspect('equal')
    leg_hands.append(ax.step(gofr.x,gofr.y-1))
    leg_str.append("g(r)")


    leg_hands.append(ax.step(gofr.x,fitting.fun_decay_exp_inv(p_out1_2,gofr.x)))
    leg_str.append("exp inv 2")


    leg_hands.append(ax.step(gofr.x,fitting.fun_decay_exp(p_out2_2,gofr.x)))
    leg_str.append("exp 2")



    print p_out1_2
    print "exp inv 2 err: " + str(err1_2)
    print p_out2_2
    print "exp 2 err: " + str(err2_2)



    ax.legend(leg_hands,leg_str)
    ax.set_title('g(r) fitting')
    ax.set_xlabel(r' r [$\mu$m]')
    ax.set_ylabel('g(r)')
    
            
    if istatus:
        print "displaying figure"
        plt.ion()
        plt.show()
    else:
        print "closing figure"
        plt.close(fig)


def tmp_series_gn(sname,conn,dtype = 't',gtype = 'gofr'):
    """Makes plots for g_1 to g_n """

    res = conn.execute("select comps.comp_key,dsets.temp\
    from comps,dsets where comps.dset_key = dsets.key and comps.function=? and\
    dsets.sname = ? and dtype = ?",(gtype,sname,dtype,)
                       ).fetchall()


    temps = [r[1] for r in res]
    gofrs = [gen.get_gofr2D(r[0],conn) for r in res]
    fits = [fitting.fit_gofr2(g,2.1,fitting.fun_decay_exp_inv) for g in gofrs]
    peaks = [gen.find_peaks_fit(g,fitting.fun_decay_exp_inv_dr,p.beta)
             for g,p in zip(gofrs,fits)]


    
    istatus = non_i_plot_start()
    
    # make g_n plots for peaks
    fig,ax = set_up_plot()
    leg_strs = []
    leg_hands = []
    for j in range(min([len(p[0]) for p in peaks])):
        leg_hands.append(ax.plot(temps,np.array([p[0][j][1] for p in peaks])-1,'x-'))
        leg_strs.append('$g_'+str(j)+'$' )

    ax.legend(leg_hands,leg_strs)
    ax.set_title('$g_n$ maximums')
    ax.set_xlabel(r'T [C]')
    ax.set_ylabel('g(peak) -1')


    
    # make g_n plots for troughs
    fig,ax = set_up_plot()
    leg_strs = []
    leg_hands = []
    for j in range(min([len(p[1]) for p in peaks])):
        leg_hands.append(ax.plot(temps,np.array([p[1][j][1] for p in peaks])-1,'x-'))
        leg_strs.append('$g_'+str(j)+'$' )

    ax.legend(leg_hands,leg_strs)
    ax.set_title('$g_n$ minimums')
    ax.set_xlabel(r'T [C]')
    ax.set_ylabel('g(peak) -1')


    # make plot for trough locations
    fig,ax = set_up_plot()
    leg_strs = []
    leg_hands = []
    for j in range(min([len(p[1]) for p in peaks])):
        leg_hands.append(ax.plot(temps,np.array([p[1][j][0] for p in peaks]),'x-'))
        leg_strs.append('$g_'+str(j)+'$' )

    ax.legend(leg_hands,leg_strs)
    ax.set_title('minimum locations')
    ax.set_xlabel(r'T [C]')
    ax.set_ylabel('r [$\mu m$]')


    # make plot for peak locations
    fig,ax = set_up_plot()
    leg_strs = []
    leg_hands = []
    for j in range(min([len(p[0]) for p in peaks])):
        leg_hands.append(ax.plot(temps,np.array([p[0][j][0] for p in peaks]),'x-'))
        leg_strs.append('$g_'+str(j)+'$' )

    ax.legend(leg_hands,leg_strs)
    ax.set_title('peak locations')
    ax.set_xlabel(r'T [C]')
    ax.set_ylabel('r [$\mu m$]')

    
    non_i_plot_stop(istatus)
    
    return peaks

def non_i_plot_start():
    istatus = plt.isinteractive();
    print istatus
    if istatus:plt.ioff()
    return istatus

def non_i_plot_stop(istatus):
    if istatus:
        print "displaying figure"
        plt.ion()
        plt.show()
    else:
        print "closing all figures"
        plt.close('all')

def set_up_plot():
    fig = plt.figure()
    ax = fig.add_axes([.1,.1,.8,.8])
    ax.hold(True)
    ax.grid(True)
    
    return fig,ax



def make_gofr_by_plane_plots(comp,conn):
    istatus = non_i_plot_start()
    (fig,ax) = set_up_plot()
    
    gc = gen.get_gofr_by_plane_cps(comp,conn)
    
    dset = conn.execute("select dset_key from comps where comp_key = ?",(comp,)).fetchone()[0]
    (temp,dtype) = conn.execute("select temp,dtype from dsets where key = ?",(dset,)).fetchone()

    ax.plot([np.max(g.y) for g in gc])
    ax.set_title("dset: " + str(dset) + " temp: " + str(temp) + "C  dtype:" + dtype)
    
    
    non_i_plot_stop(istatus)
    


def make_gofr_by_plane(d_lst,conn):
    istatus = non_i_plot_start()
    (fig,ax) = set_up_plot()
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
    non_i_plot_stop(istatus)


class color_mapper:
    def __init__(self,mn,mx,name = 'jet'):
        self._mn = mn
        self._mx = mx
        self.cmp = cm.get_cmap(name)
    def get_color(self,t):
        return self.cmp((t-self._mn)/(self._mx-self._mn))
    


def gn_type_plots(sname,conn):

    istatus = non_i_plot_start()
    
    # make g_n plots for peaks
    fig,ax = set_up_plot()
    leg_strs = []
    leg_hands = []


    for t in ['t','z']:
        res = conn.execute("select comps.comp_key,dsets.temp\
        from comps,dsets where comps.dset_key = dsets.key and comps.function=? and\
        dsets.sname = ? and dtype = ?",('gofr',sname,t,)
                       ).fetchall()


        temps = [r[1] for r in res]
        gofrs = [gen.get_gofr2D(r[0],conn) for r in res]
        fits = [fitting.fit_gofr2(g,2.1,fitting.fun_decay_exp_inv) for g in gofrs]
        peaks = [gen.find_peaks_fit(g,fitting.fun_decay_exp_inv_dr,p.beta)
                 for g,p in zip(gofrs,fits)]


    

        for j in range(min([len(p[0]) for p in peaks])):
            leg_hands.append(ax.plot(temps,np.array([p[0][j][1] for p in peaks])-1,'x-'))
            leg_strs.append('$g_'+str(j)+'$ ' + t )



    ax.legend(leg_hands,leg_strs)
    ax.set_title('$g_n$ maximums')
    ax.set_xlabel(r'T [C]')
    ax.set_ylabel('g(peak) -1')


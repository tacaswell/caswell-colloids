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
import trackpy.cpp_wrapper as cpp_wrapper
import matplotlib
import matplotlib.pyplot as plt
import random
import itertools
import h5py
import numpy
import math

# change to take 
def _plot_file_frame_phi6(f,comp_number,fr_num):
    '''
    Takes an open hdf_file handle, f, 
    '''
    x = f["/frame%(#)06d"%{"#":fr_num}+"/x_%(#)07d"%{"#":comp_number}]
    y = f["/frame%(#)06d"%{"#":fr_num}+"/y_%(#)07d"%{"#":comp_number}]
    phi = f["/frame%(#)06d"%{"#":fr_num}+"/scaler_order_parameter_%(#)07d"%{"#":comp_number}]
    phir = zeros(phi.shape[0])
    phii = zeros(phi.shape[0])
    for j in range(0,phi.shape[0]):
        phir[j] = phi[j][0]
        phii[j] = phi[j][1]
        

    plt.quiver(x,y,phir,phii)
    plt.plot(x,y,'ro')
    plt.show()


def _plot_file_frame_pos(f,comp_number, fr_num):
    x = f["/frame%(#)06d"%{"#":fr_num}+"/x_%(#)07d"%{"#":comp_number}]
    y = f["/frame%(#)06d"%{"#":fr_num}+"/y_%(#)07d"%{"#":comp_number}]
    plt.plot(x,y,'ro')
    plt.show()


def _draw_gofr_hex_lines(ax,r0):
    '''This function will draw on a graph vertical linse where peaks in g(r) are expected
    for a hex packing'''

    lin_scale = .85#((4+2*math.sqrt(3))/2 -2)/2
    irr_pos = [2*math.sqrt(3),  math.sqrt(28)]
    irr_pos_txt = [r'$2\, \sqrt{3}$',r'$\sqrt{28}$']
    for s in range(2,9):
        ax.plot(s*r0*numpy.ones(2),[0 ,3],'r')
        ax.annotate(str(s),xy=(s*r0,2.5),xycoords='data',
                    xytext=(-1,0),textcoords='offset points')
    for s,t in zip(irr_pos,irr_pos_txt):
        ax.annotate(t,xy=(s*r0,2.75),xycoords='data',
                    xytext=(-1,0),textcoords='offset points')
        ax.plot(s*r0*numpy.ones(2),[0 ,3],'k')
    for s in range(1,6):
        ax.plot((1+ s*lin_scale)*2*r0*numpy.ones(2),[0 ,3],'m')
        ax.annotate(str(s),xy=(2*r0*(1+s*lin_scale),2.25),xycoords='data',
                    xytext=(-1,0),textcoords='offset points')
    

def _get_gofr_group(fname,prefix,comp_num):
    '''Returns the h5py group that specified '''
    return  h5py.File(fname,'r')[prefix + "_%(#)07d"%{"#":comp_num}]
    
def make_2dv3d_plot(key,conn,fname = None):
    # add error handling on all of these calls
    
    # get comp_number of gofr
    res = conn.execute("select comp_key,fout from comps where dset_key == ? and function == 'gofr'",(key,)).fetchall()
    (g_ck,g_fname) = res[0]
        
        
    # get comp_number of 3D gofr
    (g3D_ck,g3D_fname) = conn.execute("select comp_key,fout from comps where dset_key == ? and function == 'gofr3D'",(key,)).fetchone()


    # get dset name
    (sample_name, temp) = conn.execute("select sname,temp from dsets where key == ? ",(key,)).fetchone()

    print sample_name + " " + str(temp)
    group = _get_gofr_group(g_fname,'gofr',g_ck)
    group3D = _get_gofr_group(g3D_fname,'gofr3D',g3D_ck)


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
    d0 = group3D[dset_names[1]][numpy.argmax(group3D[dset_names[0]])]
    print numpy.argmax(group3D[dset_names[0]])
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
        plt.draw()
        plt.ion()
    else:
        close(fig)
        
        


def make_gofr_tmp_series(sname,conn,fnameg=None,fnamegn=None):
    '''Takes in a sample name and plots all of the g(r) for it '''
    dset_names = ['bin_count', 'bin_edges']
    
    res = conn.execute("select comps.comp_key,comps.fout,dsets.temp from comps,dsets where comps.dset_key = dsets.key and comps.function='gofr3D' and dsets.sname = ?",(sname,)).fetchall()

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

    gn_g = []
    gn_t = []
    gn_p = []
    gn_fig = plt.figure()
    gn_ax = gn_fig.add_axes([.1,.1,.8,.8])
    
    for r in res:
        g = _get_gofr_group(r[1],'gofr3D',r[0])
        leg_hands.append(ax.plot(g[dset_names[1]],g[dset_names[0]]))
        leg_strs.append(str(r[2]))
        try:
            gn_p.append((float(r[2]),numpy.max(g[dset_names[0]])))
            # gn_t.append((float(r[2]))
            # gn_g.append(numpy.max(g[dset_names[0]]))
        except (ValueError,TypeError ) :
            gn_p.append((25,numpy.max(g[dset_names[0]])))
            # gn_t.append(25)
            # gn_g.append(numpy.max(g[dset_names[0]]))
            pass

    gn_p.sort(lambda x,y: int(numpy.sign(x[0]-y[0])))
    for p in gn_p:
        gn_t.append(p[0])
        gn_g.append(p[1])
    print gn_t
    print gn_g
    fig.legend(leg_hands,leg_strs)
    ax.set_title(sname)
    ax.set_xlabel(r'r [$\mu m$]')
    ax.set_ylabel(r'G(r)')
    
    gn_ax.plot(gn_t,gn_g,'x-')
    gn_ax.set_title(r'$g_1(T)$')
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

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
import fitting


def get_gofr3D(comp_num,conn):
    '''Takes in computation number and database connection and extracts the given g(r) and returns it as a
    numpy array'''

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




def make_gofr_2Dv3D(conn):
    keys = conn.execute("select key from dsets where dtype = 'z'").fetchall()
    for k in keys:
        plots.make_2dv3d_plot(k[0],conn,'figures/2Dv3D/%(#)02d.png'%{"#":k[0]})

    

def get_gofr2D(comp_num,conn):
    '''Takes in computation number and database connection and extracts the given g(r) and returns it as a
    numpy array'''

    dset_names = ['bin_count', 'bin_edges']
    
    res = conn.execute("select fout from comps where comp_key = ?",(comp_num,)).fetchall()
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

    
def find_peaks(gofr_pairs,expt_spc = None):
    '''Takes in a pairs object and returns a list of tuples with (x[indx],y[indx],indx)

    expt_spc is the expected spacing of the peaks

    Assumes that the largest peak comes first
    '''

    # find the largest (assumed to be first) peak
    
    
        
    if expt_spc is  None:
        # if no guess given, make a guess based on the location of the first peak
        pass


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
    pfit = fit_peak(gp.x[indx-wind:indx+wind],gp.y[indx-wind:indx+wind])

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
        indx = val_to_indx(gp.x,crit_p)
        # pick out window around that box
        pfit = fit_peak(gp.x[indx-wind:indx+wind],gp.y[indx-wind:indx+wind])
        

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



def fit_peak(x,y):
    """Fits a quadratic to the data points handed in """
    def quad(B,x):
        return B[0] *(x -B[1]) ** 2 + B[2]

    beta = (0,np.mean(x),y[val_to_indx(x,np.mean(x))])

    data = sodr.Data(x,y)
    model = sodr.Model(quad)
    worker = sodr.ODR(data,model,beta)
    out = worker.run()


    
    ## plts.figure()
    ## plts.plot(x,y)
    ## plts.plot(x,quad(out.beta,x))
    ## plts.title(out.beta[1])
    return out

def get_list_gofr (sname,conn,date = None,gtype = 'gofr'):
    """Takes in a  sample name and an optional string specifying
    the date of the computations to use"""

    print date
    print gtype
    if date is None:
        return conn.execute("select comps.comp_key,comps.fout,dsets.temp,comps.fin \
        from comps,dsets where comps.dset_key = dsets.key and comps.function=? and\
        dsets.sname = ?",(gtype,sname,)).fetchall()
    else:
        return conn.execute("select comps.comp_key,comps.fout,dsets.temp,comps.fin \
        from comps,dsets where comps.dset_key = dsets.key and comps.function=? and\
        dsets.sname = ? and dsets.dtype = 't' and comps.date = ?",(gtype,sname,date)).fetchall()

    
def get_gofr_by_plane_cps(comp_num,conn):
    fname = conn.execute("select fout from comps where comp_key = ?",(comp_num,)).fetchall()
    fname = fname[-1][0]
    F =  h5py.File(fname)
    g = F['gofr_by_plane_%(#)07d'%{'#':comp_num}]

    g_l = [cord_pairs(g[c]['bin_edges'][:],g[c]['bin_count'][:]) for c in g]

    
    del g
    F.close()
    return g_l

def get_gofr_by_plane_tmp(comp_num,conn):
    (fname,) = conn.execute("select fout from comps where comp_key = ?",(comp_num,)).fetchone()
    
    F =  h5py.File(fname)
    g = F['gofr_by_plane_%(#)07d'%{'#':comp_num}]

    temps = [g[c].attrs['temperature'] for c in g]

    
    del g
    F.close()
    return temps


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
        
        

def make_gofr_tmp_series(comp_list,conn,fnameg=None,fnamegn=None,gtype='gofr',date = None):
    '''Takes in a sample name and plots all of the g(r) for it '''
    r_scale = 6.45/60
    dset_names = ['bin_count', 'bin_edges']

    
    res = [conn.execute("select comps.comp_key,comps.fout,dsets.temp,comps.fin \
    from comps,dsets where comps.comp_key = ? and dsets.key = comps.dset_key",c).fetchone()
           for c in comp_list]

    res.sort(key=lambda x:x[2])

    
    (sname,) = conn.execute("select sname from dsets where key in (select dset_key from comps where\
    comp_key = ?)",(res[0][0],)).fetchone()

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
    
    gn_g = []
    gn_t = []
    gn_p = []

    
    gn_fig = plt.figure()
    gn_ax = gn_fig.add_axes([.1,.1,.8,.8])

    def get_gofr_tmp(fname,comp_num):
        F = h5py.File(fname,'r')
        t = F[gtype + "_%(#)07d"%{"#":comp_num}].attrs['temperature']
        F.close()
        return t
    
    temps = [get_gofr_tmp(r[1],r[0]) for r in res]
    tmax = max(temps)
    tmin = min(temps)
    ax.set_color_cycle([cm.jet((t-tmin)/(tmax-tmin)) for t in temps])
    for r,temperature in zip(res,temps):


        
        F = h5py.File(r[1],'r')
        g = F[gtype + "_%(#)07d"%{"#":r[0]}]
        
        
        leg_hands.append(ax.plot(g[dset_names[1]][:]*r_scale,g[dset_names[0]]))
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
        F.close()
    gn_p.sort(lambda x,y: int(np.sign(x[0]-y[0])))
    for p in gn_p:
        gn_t.append(p[0])
        gn_g.append(p[1])
    print gn_t
    print gn_g
    fig.legend(leg_hands,leg_strs,loc=0)
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

def make_2d_gofr_plot(comp_key,conn,fname = None):
    # add error handling on all of these calls
    
    # get comp_number of gofr
    (key,g_fname) = conn.execute("select dset_key,fout from comps where comp_key == ? and function == 'gofr'",(comp_key,)).fetchone()
    
        
    # get dset name
    (sname, temp) = conn.execute("select sname,temp from dsets where key == ? ",(key,)).fetchone()

    print sname + " " + str(temp)
    group = gen.get_gofr_group(g_fname,'gofr',comp_key)



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


def gn2D(comp_lst,conn,fig,txt):
    

    res = [conn.execute("select comps.comp_key,dsets.temp\
    from comps,dsets where comps.comp_key = ? and dsets.key = comps.dset_key",c).fetchone()
           for c in comp_lst]

    
    
    temps = [r[1] for r in res]
    gofrs = [gen.get_gofr2D(r[0],conn) for r in res]
    fits = [fitting.fit_gofr2(g,2.1,fitting.fun_decay_exp_inv) for g in gofrs]
    peaks = [gen.find_peaks_fit(g,fitting.fun_decay_exp_inv_dr,p.beta)
             for g,p in zip(gofrs,fits)]


    fig.plot(temps,np.array([p[0][0][1] for p in peaks])-1,txt)

def tmp_series_gn2D(comp_list,conn):
    """Makes plots for g_1 to g_n
    
    comp_list : list of 1 element tuple that contain the computation
    number to be plotted
    """

    res = [conn.execute("select comp_key,fout\
    from comps where comps.comp_key = ?",c).fetchone()
           for c in comp_list]

    res.sort(key=lambda x:x[1])
    
    def get_gofr_tmp(fname,comp_num):
        F = h5py.File(fname,'r')
        t = F["gofr_%(#)07d"%{"#":comp_num}].attrs['temperature']
        F.close()
        return t
    
    temps = [get_gofr_tmp(r[1],r[0]) for r in res]


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

    #ax.legend(leg_hands,leg_strs,loc=0)
    add_labels(ax,'$g_n$ maximums',r'T [C]','g(peak) -1')


    
    # make g_n plots for troughs
    fig,ax = set_up_plot()
    leg_strs = []
    leg_hands = []
    for j in range(min([len(p[1]) for p in peaks])):
        leg_hands.append(ax.plot(temps,np.array([p[1][j][1] for p in peaks])-1,'x-'))
        leg_strs.append('$g_'+str(j)+'$' )

    #ax.legend(leg_hands,leg_strs,loc=0)
    add_labels(ax,'$g_n$ minimums',r'T [C]','g(peak) -1')
    

    # make plot for trough locations
    fig,ax = set_up_plot()
    leg_strs = []
    leg_hands = []
    for j in range(min([len(p[1]) for p in peaks])):
        leg_hands.append(ax.plot(temps,np.array([p[1][j][0] for p in peaks]),'x-'))
        leg_strs.append('$g_'+str(j)+'$' )

    #ax.legend(leg_hands,leg_strs,loc=0)
    add_labels(ax,'minimum locations',r'T [C]','r [$\mu m$]')


    # make plot for peak locations
    fig,ax = set_up_plot()
    leg_strs = []
    leg_hands = []
    for j in range(min([len(p[0]) for p in peaks])):
        leg_hands.append(ax.plot(temps,np.array([p[0][j][0] for p in peaks]),'x-'))
        leg_strs.append('$g_'+str(j)+'$' )

    #ax.legend(leg_hands,leg_strs,loc=0)
    add_labels(ax,'peak locations',r'T [C]','r [$\mu m$]')

    
    non_i_plot_stop(istatus)
    
    return peaks

def make_gofr_by_plane_plots(comp,conn):
    istatus = non_i_plot_start()
    (fig,ax) = set_up_plot()
    
    gc = gen.get_gofr_by_plane_cps(comp,conn)
    temps = gen.get_gofr_by_plane_tmp(comp,conn)
    dset = conn.execute("select dset_key from comps where comp_key = ?",(comp,)).fetchone()[0]
    (temp,dtype,fname) = conn.execute("select temp,dtype,fname from dsets where key = ?",(dset,)).fetchone()

    wind = 30
    max_indx = [np.argmax(g.y[5:])+5 for g in gc]
    betas = [gen.fit_peak(g.x[m-wind:m+wind],g.y[m-wind:m+wind]) for m,g in zip(max_indx,gc)]
    ax.step(temps,[b.beta[2] for b in betas])
    ax.step(temps,[np.max(g.y) for g in gc])
    #    ax.set_title("dset: " + str(dset) + " temp: " + str(temp) + "C  dtype:" + dtype)
    ax.set_title("dset: " + str(dset) + " " + fname.split('/')[-1])
    
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
    add_labels(ax,'$g_n$ maximums',r'T [C]','g(peak) -1')


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

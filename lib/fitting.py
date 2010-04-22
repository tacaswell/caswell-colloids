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

from util import cord_pairs
import numpy as np
import scipy.optimize
import general

def _trim_gofr(gofr,r0):
    """     Takes in gofr as a cord_pairs and r0 in
    units of 1/d, d as determined by the location
    of the first peak.

    Returns a cords_pair object with g(r>r0) 
    """
    # find the first peak
    d0 = gofr.x[np.argmax(gofr.y)]

    print "the Diameter is " + str(d0)

    # convert r0 to real units
    r0 *=d0

    
    
    # cut the data at r0
    r_arg = np.argmax([v for v in gofr.x if v<= r0])
    gofr_trim = cord_pairs(gofr.x[r_arg:],gofr.y[r_arg:])



    return gofr_trim

def fun_decay_exp_inv(r,p):
    """Returns C/r exp(- r/a) cos(K(r)+phi_0)
    evaluated at r.  p = (a,K,C,phi_0)"""
    return (p[2] / r) * np.exp(-r/p[0]  ) * np.cos(p[1]*r + p[3])+ r*p[4] + p[5]


def fun_decay_exp(r,p):
    """Returns C exp(- r/a) cos(K(r)+phi_0)
    evaluated at r.  p = (a,K,C,phi_0)"""
    return (p[2] ) * np.exp(-r/p[0]  ) * np.cos(p[1]*r + p[3]) + (r)*p[4] + p[5]


def fun_decay_inv(r,p):
    """Returns C/r cos(K(r)+phi_0)
    evaluated at r.  p = (K,C,phi_0)"""
    return (p[1] / r) * np.cos(p[0]*r + p[2])


def fit_gofr(gofr,r0,func,p0 =(1,12,.5,0)):
    """
    Fits to 
    Takes in gofr as a cord_pairs and r0 in
    units of 1/d, d as determined by the location
    of the first peak.

    Returns the tuple that is the argument for the function
    used for the fitting, the covariance matrix, and the
    sum of the residual squared
    """
    # trim the g(r) data
    gofr = _trim_gofr(gofr,r0)
    

    def local_fun(p):
        return func(gofr.x,p) - gofr.y+1

    (p_out,p_cov,id,m,flg) = scipy.optimize.leastsq(local_fun,p0,full_output=1)
    print flg
    return p_out,p_cov,(sum(local_fun(p_out))**2)/len(gofr.x)
    # shove into the numpy fitting code


def fit_tmp_series_gofr2D(sname,conn):
    '''Takes in a sample name and fits all of the 2D gofrs'''
    
    
    res = conn.execute("select comps.comp_key,comps.fout,dsets.temp from comps,dsets where comps.dset_key = dsets.key and comps.function='gofr3D' and dsets.sname = ?",(sname,)).fetchall()


    temps = []
    fits = []
    for r in res:
        gofr = general.get_gofr3D(r[0],conn)
        fits.append(fit_gofr(gofr,2,fun_decay_exp_inv,(2,7.35,1.5,0,0,0)))
        try:
            temps.append(float(r[2]))
        except (ValueError,TypeError ) :
            temps.append(25)

            
    return zip(fits,temps)

def make_sofq_3D_plot(key,conn,Q):
    '''From the key plots s(q) as computed from the 3D g(r)'''
    res = conn.execute("select comp_key from comps where dset_key=? and function='gofr3D'",(key,)).fetchall()
    if not len(res)==1:
        raise util.dbase_error("can't find 3D gofr")

    plt.figure()
    g = general.get_gofr3D(res[0][0],conn)

    S = sofQ(g,Q)

    res2 = conn.execute("select comp_key from comps where dset_key=? and function='gofr'",(key,)).fetchall()
    if not len(res2)==1:
        raise util.dbase_error("can't find gofr")
    
    g2 = general.get_gofr2D(res2[0][0],conn)
    S2 = general.sofQ(g2,Q)


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
    group = general.get_gofr_group(g_fname,'gofr',g_ck)



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



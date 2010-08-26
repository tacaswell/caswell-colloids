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

def open_conn():
    '''Opens the data base at the standard location and returns the connection'''
    conn =sqlite3.connect('/home/tcaswell/colloids/processed/processed_data.db')
    conn.execute('PRAGMA foreign_keys = ON;')
    return conn

def delete_comp(comp_num,conn):
    """Deletes a computation from the database, does not touch the actual data """
    r = conn.execute("select * from comps where comp_key = ?",(comp_num,)).fetchone();
    print "Preparing to delete the row " + str(comp_num)
    print r
    if util.query_yes_no("are you sure you want to delete this entry? ",default='no')=='yes':
        conn.execute("delete from comps where comp_key = ?",(comp_num,))
        conn.commit()
    pass

def get_fout_comp(key,conn,func):
    '''Returns the fout name and the computation number for the iden on the dataset given by key'''
    res = conn.execute("select fout,comp_key from comps where dset_key == ? and function == ?",
                       (key,func)).fetchall()
    if not len(res) == 1:
        raise lib.util.dbase_error("either no or too many iden operations found")

    return res[0]
    # 

def get_xyzI_f(h5file,frame,comp_num):
    '''Extracts the x,y,z,I from the given frame and comp_num from the open file handed in '''

    cords = lib.util.Cords3D()
    cords.x = h5file["/frame%(#)06d"%{"#":frame}+"/x_%(#)07d"%{"#":comp_num}][:]
    cords.y = h5file["/frame%(#)06d"%{"#":frame}+"/y_%(#)07d"%{"#":comp_num}][:]
    cords.z = h5file["/frame%(#)06d"%{"#":frame}+"/z_%(#)07d"%{"#":comp_num}][:]
    cords.I = h5file["/frame%(#)06d"%{"#":frame}+"/intensity_%(#)07d"%{"#":comp_num}][:]

    return cords


def get_xyzI(key,conn,frame):
    '''Returns four vectors for x,y,z,I for the data set given by key'''

    (fname, comp_number) = get_fout_comp(key,conn,'link3D')

    f = h5py.File(fname,'r')

    cord = get_xyzI_f(f,frame,comp_number)
    f.close()
    
    return cord

    pass
        

def make_z_slices_series(conn,key,step,base_dir,base_name):
    """Takes in a data base connection, a dset key, a directory to dump to, and a basename """

    # check to make sure base_dir exists
    cords = get_xyzI(key,conn,0)

    lib.pov.make_z_slices(cords,step,base_dir+'/'+base_name)
    pass


def make_phi6(conn,key,save=None):
    # get fname from

    lib.plots._plot_file_frame_phi6(key,conn,15,save)

def make_nsize(conn,key,save=None):
    # get fname from

    lib.plots._plot_file_frame_nsize(key,conn,5,save)

def mean_n_size_frame(key,conn,fr_num):

    (fname,p_comp) = conn.execute("select fout,comp_key from comps where function = 'phi6' and dset_key = ?;",(key,)).fetchone()


    (sname,stype,temp) = conn.execute("select sname,dtype,temp from dsets where key = ?",(key,)).fetchone()

    
    f = h5py.File(fname,'r')
   

    ns = f["/frame%(#)06d"%{"#":fr_num}+"/neighborhood_size_%(#)07d"%{"#":p_comp}]

    nmean = np.mean(ns)
    print nmean
    return nmean
    

    

def sofQ(c_pair,Q ):
    '''computes the Fourier transform of the c_pair passed in at all q in Q.'''


    tmp_y = np.hstack((np.flipud(c_pair.y),c_pair.y))-1
    tmp_x = np.hstack((np.flipud(c_pair.x),c_pair.x))
    tmp_y = tmp_y*tmp_x
    tmp_x = np.hstack((np.flipud(-c_pair.x),c_pair.x))

    plt.plot(tmp_x,tmp_y)
    
    dx = np.mean(np.diff(c_pair.x))
    sq = lambda q: abs(dx*(1j/q)*sum(tmp_y*np.exp(1j*tmp_x*q*2*np.pi)))**2
    S = map(sq,Q)
       
    
    return S



def val_to_indx(x,x0):
    """Returns the index of the first value in x greater than x0, assumes the list is sorted"""
    return bisect.bisect_left(x,x0)


def crazy_s(g):
    """Computes the excess structural entropy ala a number of papers """
    return sum([ (y*np.log(y) + y-1)*(x**2)
                 for (x,y) in itertools.izip(g.x,g.y)
                 if y != 0])*np.mean(np.diff(g.x))*np.pi
    
    



def estimate_phi3D(key,conn):
    """Estimates the volume fraction from the number of particles, the
    position of the first peak and the dimensions of the data"""

    # check to make sure that both the link and 3D g(r) exist and get
    # logistic information
    link_res = conn.execute("select fout,comp_key from comps where dset_key = ? and function = 'link3D'",(key,)).fetchall()
    if len(link_res) == 0:
        raise "link computation not found"
    link_res = link_res[-1]

    g_res = conn.execute("select comp_key from comps where dset_key = ? and function = 'gofr3D'",(key,)).fetchall()
    if len(g_res) == 0:
        raise dbase_error("g(r) 3D computation not found")
    g_res = g_res[-1]



    # get g(r) data and fit first peak
    g = get_gofr3D(g_res[0],conn)
    pram = fitting.fit_gofr2(g,2.1,fitting.fun_decay_exp_inv) 
    peaks = find_peaks_fit(g,fitting.fun_flipper(fitting.fun_decay_exp_inv_dr),pram.beta)
    r = peaks[0][0][0]/2
    # get dimension from data set
    F = h5py.File(link_res[0],'r')
    dim = F.attrs.get('dims',0);

    
    # get number of entries in data set
    p_count = len(F["/frame000000/x_%(#)07d"%{"#":link_res[1]}])
    # compute
    phi = (np.pi * r**3 * 4/3) * p_count / (np.prod(dim))
    # clean up everything
    F.close()
    return (phi,p_count/(.63*np.prod(dim)/(np.pi*(.98/2)**3*4/3)))
    pass

## def open_gofr_by_plane_group(comp_num,conn):
##     fname = conn.execute("select fout from comps where comp_key = ?",(comp_num,)).fetchall()
##     fname = fname[-1][0]
##     return  h5py.File(fname)['gofr_by_plane_%(#)07d'%{'#':comp_num}]
    

def get_exp_time(comp_num,conn):
    """
    Given a set number returns the effective exposure of the first frame in an Iden
    """
    (func,fname) = conn.execute("select function,fout from comps where comp_key = ? and " +
                                " function like 'Iden%'",(comp_num,)).fetchone()

    
    F = h5py.File(fname,'r')
    if 'Exposure' in F.attrs.keys():
        exp_time = float(F.attrs['Exposure'].strip().split(' ')[0])
    elif 'Exposure' in  F['/frame000000'].attrs.keys():
        exp_time = F['/frame000000'].attrs['Exposure']

    F.close()
    
    return exp_time
    

def ff(n):
    """Formats frame names """
    return "frame%(#)06d"%{"#":n}

def fd(str_,n):
    """ formats dset names"""
    return str_ + "_%(#)07d"%{"#":n}

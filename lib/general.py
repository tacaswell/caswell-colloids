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
import lib.pov 
import lib.plots 
import lib.util 
import numpy as np
from util import cord_pairs

def open_conn():
    '''Opens the data base at the standard location and returns the connection'''
    conn =sqlite3.connect('/home/tcaswell/colloids/processed/processed_data.db')
    conn.execute('PRAGMA foreign_keys = ON;')
    return conn

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


def make_gofr_2Dv3D(conn):
    keys = conn.execute("select key from dsets where dtype = 'z'").fetchall()
    for k in keys:
        plots.make_2dv3d_plot(k[0],conn,'figures/2Dv3D/%(#)02d.png'%{"#":k[0]})
        

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
    

def get_gofr_group(fname,prefix,comp_num):
    '''Returns the h5py group that is  specified '''
    return  h5py.File(fname,'r')[prefix + "_%(#)07d"%{"#":comp_num}]



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
    

def get_gofr2D(comp_num,conn):
    '''Takes in computation number and database connection and extracts the given g(r) and returns it as a
    numpy array'''

    dset_names = ['bin_count', 'bin_edges']
    
    res = conn.execute("select fout from comps where comp_key = ?",(comp_num,)).fetchall()
    if not len(res) == 1:
        print len(res)
        raise util.dbase_error("error looking up computation")

    
    g = get_gofr_group(res[0][0],'gofr',comp_num)
    gofr = np.array(g[dset_names[0]])
    bins = np.array(g[dset_names[1]])*6.45/60
    return cord_pairs(bins,gofr)
    

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

    
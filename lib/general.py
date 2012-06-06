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
import bisect
import datetime
import itertools
import sqlite3

import h5py
import numpy as np

import fitting
import util 


def open_conn(fname=None):
    '''Opens the data base at the standard location and returns the connection'''
    if fname is None:
        fname = "/home/tcaswell/colloids/proc_db.db"
    conn =sqlite3.connect(fname)
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
    
def get_acq_time(dset_key,conn):
    """
    Given a dset key, gets the local acquisition time of the first plane.  Needs to have
    had iden run on the data set to work
    """

    (fname,iden_key) = conn.execute("select fout,comp_key from iden where dset_key = ?"
                                    ,(dset_key,)).fetchone()
    F = h5py.File(fname,'r')

    acq_time = F['/frame000000'].attrs['acquisition-time-local']
    

    F.close()
    del F

    return acq_time

def ff(n):
    """Formats frame names """
    return "frame%(#)06d"%{"#":n}

def fd(str_,n):
    """ formats dset names"""
    return str_ + "_%(#)07d"%{"#":n}

def extract_spatial_calibration(group):
    """group : a open valid hdf group pointing in to a frame

    returns the spatial calibration
    """

    if not group.attrs['spatial-calibration-state']:
        raise Error("Spatial calibration not enabled")

    

    return group.attrs['spatial-calibration-x']

def get_z_pos(comp_num,conn):
    """
    given an iden like computation number returns the (z,y,x) cords of the first frame
    """
    (func,fname) = conn.execute("select function,fout from comps where comp_key = ? and " +
                                " function like 'Iden%'",(comp_num,)).fetchone()

    pos = None
    F = h5py.File(fname,'r')
    if 'z-position' in  F['/frame000000'].attrs.keys():
        pos = (F['/frame000000'].attrs['z-position'],
               F['/frame000000'].attrs['stage-position-y'],
               F['/frame000000'].attrs['stage-position-x'])

    F.close()
    del F
    return pos


def get_iden_attrs_fr(comp_num,conn,attr_lst):
    """
    given an iden like computation number returns the (z,y,x) cords of the first frame
    """
    (fname,) = conn.execute("select fout from iden where comp_key = ?",(comp_num,)).fetchone()

    out = None
    F = h5py.File(fname,'r')
    if all([at in F['/frame000000'].attrs for at in attr_lst]):
        out = [ [F[g].attrs[at] for g in F if not g=='parameters' ]for at in attr_lst]

    F.close()
    del F
    return out

def get_iden_attrs_tl(comp_num,conn,attr_lst):
    """
    given an iden like computation number returns the (z,y,x) cords of the first frame
    """
    (fname,) = conn.execute("select fout from iden where comp_key = ?",(comp_num,)).fetchone()

    out = None
    F = h5py.File(fname,'r')
    if all([at in F.attrs for at in attr_lst]):
        out = [ F.attrs[at] for at in attr_lst]

    F.close()
    
    return out

def set_plane_temp_constant(iden_key,conn,temp):
    """Takes an Iden_key, a database connection and a temperature and
    sets the value of temperature for all of the planes in the iden
    object to temp.  This is for making older iden files work with
    code that expects the temperature to be in the meta-data"""

    (fname,) = conn.execute("select fout from comps where comp_key =? and function like 'Iden%'",(iden_key,)).fetchone()
    F = h5py.File(fname,'r+')
    
    for g in F:
        if g == 'parameters':
            break
        grp = F[g]
        if 'temperature' in grp.attrs.keys():
            raise Exception("This file already has temperatures for the planes")
        grp.attrs['temperature'] = temp
    F.close()

def list_iden_attrs(comp_num,conn):
    """Returns the list of top level attributes and a list of the per from attributes"""
    (func,fname) = conn.execute("select function,fout from comps where comp_key = ? and " +
                                " function like 'Iden%'",(comp_num,)).fetchone()

    out = None
    F = h5py.File(fname,'r')
    out = [F.attrs.keys(), F['/frame000000'].attrs.keys()]
           
    F.close()
    del F
    return out



def set_plane_exposure(iden_key,conn):
    """Takes an Iden_key, a database connection.  This takes the
    Exposure in the top level data, parses it an puts it in the frame
    level meta-data where the file spec says it should be.  This is
    for making older iden files work with code that expects the
    exposure to be in the meta-data"""

    (fname,) = conn.execute("select fout from comps where comp_key =? and function like 'Iden%'",(iden_key,)).fetchone()
    F = h5py.File(fname,'r+')
    exp = F.attrs['Exposure']
    (exp_time,exp_units) = exp.strip().split(' ')
    exp_time = int(exp_time)
    for g in F:
        if g == 'parameters':
            break
        grp = F[g]
        if 'Exposure' in grp.attrs.keys():
            raise Exception("This file already has exposure for the planes")
        grp.attrs['Exposure'] = exp_time
        grp.attrs['Exposure units'] = exp_units
    F.close()
    del F



def fix_dtime(iden_key,conn):
    """This fixes the dtime values which are apparently screwed up in
    old Iden files"""
    atr = 'dtime'
    (fname,) = conn.execute("select fout from comps where comp_key =? and function like 'Iden%'",(iden_key,)).fetchone()
    F = h5py.File(fname,'r+')
    grp = F[ff(0)]
    del grp.attrs[atr]
    grp.attrs.create(atr,0,None,'int32')
    tmp = grp.attrs["acquisition-time-local"]
    init_t = datetime.datetime.strptime(tmp[:17],"%Y%m%d %H:%M:%S")
    init_t.replace(microsecond=int(tmp[18:])*1000)
    del grp
    for j in range(1,F.attrs['number-of-planes']):
        
        grp = F[ff(j)]
        del grp.attrs[atr]
        tmp = grp.attrs["acquisition-time-local"]
        cur_t = datetime.datetime.strptime(tmp[:17],"%Y%m%d %H:%M:%S")
        cur_t.replace(microsecond=int(tmp[18:])*1000)
        dt = cur_t - init_t
        dtime = dt.seconds*(10**3) + dt.microseconds/(10**3)
        grp.attrs.create(atr,dtime,None,'int32')
        del grp
        init_t = cur_t
    F.close()
    del F

def avg_dtime(iden_key,conn):
    '''Returns the average dtime of an iden set '''
    atr = 'dtime'
    (fname,) = conn.execute("select fout from comps where comp_key =?",(iden_key,)).fetchone()
    F = h5py.File(fname,'r')
    dtime_sum = 0
    frame_count = F.attrs['number-of-planes']
    for j in range(1,frame_count):
        grp = F[ff(j)]
        dtime_sum += grp.attrs[atr]
        del grp
        
    F.close()
    del F

    return dtime_sum/frame_count



def np_to_org(a):
    for b in a:
        print "|" + '|'.join(['%.2f'%c for c in b]) +'|'


def get_plane_temps(dset_key,conn):
    """Returns an array with the temperatures of the planes in the iden file """

    (fname,iden_key) = conn.execute("select fout,comp_key from iden where dset_key = ?"
                                    ,(dset_key,)).fetchone()
    try:
        F = h5py.File(fname,'r')
        temps = [F[grp].attrs['temperature'] for grp in F if not grp =='parameters']
    finally:
        F.close()
        del F

    return temps

def get_bin_centers(left_edges):
    """returns an array of the same size which represents the center
    of the bin if the input is the left edges"""

    bin_diffs = np.diff(left_edges)/2
    bin_diffs = np.append(bin_diffs,np.mean(bin_diffs))

    return left_edges + bin_diffs

def second_deriv(edges,values):
    ''' smooth the values and takes the second derivative

    Smoothing code taken from http://www.scipy.org/Cookbook/SignalSmooth

    '''
    window_len = 3
    # this mirrors the last few points out so that the convolution is happy
    s=np.r_[values[window_len-1:0:-1],values,values[-1:-window_len:-1]]
    w = np.ones(window_len,'d')
    values = np.convolve(w/w.sum(),s,mode='valid')[(window_len//2):-(window_len//2)]
    
    h = np.mean(np.diff(edges))
    dr2 = (values[2:] - 2*values[1:-1] + values[:-2])/(h*h)

    window_len = 3
    s=np.r_[dr2[window_len-1:0:-1],dr2,dr2[-1:-window_len:-1]]
    w = np.ones(window_len,'d')
    
    dr2 = np.convolve(w/w.sum(),s,mode='valid')[(window_len//2):-(window_len//2)]
    
    return edges[1:-1],dr2

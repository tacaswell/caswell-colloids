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

import plots as plt
import h5py
import numpy
import pdb
import matplotlib

def _ff(n):
    """Formats frame names """
    return "frame%(#)06d"%{"#":n}

def _fd(str_,n):
    """ formats dset names"""
    return str_ + "_%(#)07d"%{"#":n}


def write_delta_T(comp,conn):
    """Goes back and adds the delta T md to msd results"""
    (fin,fout) = conn.execute("select fin,fout from comps where comp_key = ?",(comp,)).fetchone()
    (iden_comp,) = conn.execute("select comp_key from comps where fout = ? and function = 'Iden'",(fin,)).fetchone()
    Fin = h5py.File(fin,'r')
    # get delta T
    dt = Fin['frame000006'].attrs['dtime']
    Fin.close()
    del Fin
    # set it
    Fout = h5py.File(fout,'r+')
    Fout[_fd('msd_squared',comp)].attrs['dtime'] = dt
    Fout[_fd('mean_squared_disp',comp)].attrs['dtime'] = dt
    Fout[_fd('mean_disp',comp)].attrs['dtime'] = dt
    Fout.close()
    del Fout
    
def plot_msd(comp_key,conn,fig=None):
    if fig is None:
        fig = plt.Figure('t[s]','msd','msd')
        pass
    
    (fin,dset) = conn.execute("select fout,dset_key from comps where comp_key = ?",(comp_key,)).fetchone()
    (temp,) = conn.execute("select temp from dsets where key = ?",(dset,)).fetchone()
    
    
    
    Fin = h5py.File(fin,'r')
    g_name = _fd('mean_squared_disp',comp_key)
    msd = Fin[g_name]
    msd = Fin[g_name]['data'][:]
    dt = Fin[g_name].attrs['dtime']
    print 
    print 'the delta is ',  dt, 'for comp ' ,comp_key
    t = (numpy.arange(len(msd))+1)*dt

    cm = plt.color_mapper(27,33)
    
    fig.plot(t,msd,label=str(temp),color=cm.get_color(temp))
    Fin.close()
    del Fin
    return fig

def plot_msd_series(comp_key_lst,conn,sname=None):
    """ Takes in a list of msd computations and makes a nice plot"""

    tltstr = 'mean squared displacement'
    if not sname is None:
        tltstr = tltstr + " " + sname
        


    cmplst = [c + conn.execute("select temp from dsets where key in (select dset_key from msd_prams where comp_key = ?) ",c).fetchone() for c in comp_key_lst]

    print cmplst
    cmplst.sort(key=lambda x:x[1])
    
    fig = plt.Figure('t[s]',r'$\langle \Delta \rangle ^2$',tltstr,count=len(comp_key_lst),func=matplotlib.axes.Axes.loglog)
    for c in cmplst:
        plot_msd(c[0],conn,fig)

def delete_msd(comp_key,conn):
    """Deletes all record of a computation from both the data files and the data base """
    # get data file name
    (fname,) = conn.execute("select fout from comps where comp_key = ?",(comp_key,)).fetchone()
    # open file
    F = h5py.File(fname)
    # delete data sets
    del F[_fd('mean_squared_disp',comp_key)]
    del F[_fd('mean_disp',comp_key)]
    del F[_fd('msd_squared',comp_key)]
    # make sure changes take and clean up file
    F.flush()
    F.close()
    del F
    # delete from database
    conn.execute("delete from comps where comp_key = ?",(comp_key,))
    conn.execute("delete from msd_prams where comp_key = ?",(comp_key,))
    conn.commit()
                 

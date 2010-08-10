#Copyright 2009 Thomas A Caswell
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


import h5py
import os.path
import math
import re
import numpy as np

def _get_real_dims(h5f,comp_num):
    real_dims = np.zeros(3)
    real_dims[0] = math.ceil(h5f["/frame000000/x_%(#)07d"%{"#":comp_num}][:].max())
    real_dims[1] = math.ceil(h5f["/frame000000/y_%(#)07d"%{"#":comp_num}][:].max())
    real_dims[2] = math.ceil(h5f["/frame000000/z_%(#)07d"%{"#":comp_num}][:].max())
    return real_dims


def _get_real_2d_dims(h5f,comp_num):
    real_dims = np.zeros(2)
    real_dims[0] = math.ceil(h5f["/frame000000/x_%(#)07d"%{"#":comp_num}][:].max())
    real_dims[1] = math.ceil(h5f["/frame000000/y_%(#)07d"%{"#":comp_num}][:].max())
    return real_dims


def _get_attr_dims(h5f):
    return h5f.attrs.get('dims',0);


def fix_all_dims(comp_num,conn):
    
    fname = conn.execute("select fout from comps where comp_key = ?",(comp_num,)).fetchone()[0]
    F = h5py.File(fname,'r+')
 

    real_dims = _get_real_dims(F,comp_num)

    
    del F.attrs['dims']
    F.attrs.create('dims',real_dims,None,'float32')
    F.close()

def fix_z_dim(comp_num,conn):
    
    fname = conn.execute("select fout from comps where comp_key = ?",(comp_num,)).fetchone()[0]
    F = h5py.File(fname,'r+')
    
    attr_dims = _get_attr_dims(F)
    real_dims = _get_real_dims(F,comp_num)
    attr_dims[2] = real_dims[2]
    
    F.attrs.modify('dims',attr_dims)
    F.close()


def compare_dims(comp_num,conn):
    fname = conn.execute("select fout from comps where comp_key = ?",(comp_num,)).fetchone()[0]
    F = h5py.File(fname,'r')
    
    attr_dims = _get_attr_dims(F)
    print attr_dims
    real_dims = _get_real_dims(F,comp_num)
    print real_dims
    F.close()
    return attr_dims,real_dims        


def fix_xy_dims(comp_num,conn):
    
    fname = conn.execute("select fout from comps where comp_key = ?",(comp_num,)).fetchone()[0]
    F = h5py.File(fname,'r+')
 

    real_dims = _get_real_2d_dims(F,comp_num)

    
    #del F.attrs['dims']
    F.attrs.create('dims',real_dims,None,'float32')
    F.close()



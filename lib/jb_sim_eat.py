#Copyright 2011 Thomas A Caswell
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

# this module is for eating and parsing data files from Justin Burton

# the expected format of the text file is:
#x1
#y1
#z1
#x2
#y2
#z2

import h5py
import numpy as np
import datetime

def eat_file(fname):
    """Expects a absolute path name.  Returns three lists for x,y,z"""
    F = open(fname)
    data = F.readlines()
    F.close()
    
    if len(data)%3 != 0:
        raise Exception("the data has the wrong number of particles")

    x,y,z = zip(* np.reshape(np.array([float(d) for d in data]),(-1,3)))
    
    return x,y,z


def generate_file(x,y,z,c_key,phi,fname):
    F = h5py.File(fname,'w')
    g = F.create_group(ff(0))

    x_ds = g.create_dataset(fd('x',c_key),x.shape,x.dtype,compression='szip',chunks=(1,2000))
    x_ds[:] = x
    y_ds = g.create_dataset(fd('y',c_key),x.shape,x.dtype,compression='szip',chunks=(1,2000))
    y_ds[:] = y
    z_ds = g.create_dataset(fd('z',c_key),x.shape,x.dtype,compression='szip',chunks=(1,2000))
    z_ds[:] = z

    g.attrs['name'] = 'sim data'
    g.attrs['Exposure'] = 0
    g.attrs['Exposure units'] = 'na'
    g.attrs['dtime'] = 0
    g.attrs['spatial-calibration-state'] = 1
    g.attrs['spatial-calibration-x'] = 1
    g.attrs['spatial-calibration-y'] = 1 
    g.attrs['spatial-calibration-units'] = 'arb'
    tmp_date = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')
    g.attrs['acquisition-time-local'] = tmp_date
    g.attrs['modification-time-local'] = tmp_date
    g.attrs['stage-position-x'] = 0
    g.attrs['stage-position-y'] = 0
    g.attrs['z-position'] = 0
    g.attrs['pixel-size-x'] = 1
    g.attrs['pixel-size-y'] = 1
    g.attrs['temperature'] = 0
    

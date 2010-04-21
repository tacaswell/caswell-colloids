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
import lib.general as general
import trackpy.cpp_wrapper

def main_func():
    conn = general.open_conn()

    cmps = conn.execute("select comp_key from comps where function = 'phi6'").fetchall();
    for c in cmps:
        (fname,) = conn.execute("select fout from comps where comp_key = ? ;",c).fetchone()
        print fname
        print c

        f = h5py.File(fname,'r+')
        nf = f.attrs["number-of-planes"]
        del f["/parameters/scaler_order_parameter_%(#)07d"%{"#":c[0]}]
        del f["/parameters/neighborhood_size_%(#)07d"%{"#":c[0]}]
        for fr in range(0,nf):
            del f["/frame%(#)06d"%{"#":fr}+"/neighborhood_size_%(#)07d"%{"#":c[0]}]
            del f["/frame%(#)06d"%{"#":fr}+"/scaler_order_parameter_%(#)07d"%{"#":c[0]}]
            
            pass

    

  



#main_func()


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
import general
import appwrappy.cpp_wrapper
def main_func():
    conn = general.open_conn()
    res = conn.execute('select key from dsets;').fetchall()
    for r in res:
        print r[0]
        
        if len(conn.execute("select comp_key from comps where function = 'phi6' and dset_key = ?;",r).fetchall()) ==0:
            appwrappy.cpp_wrapper.do_phi6(r[0],conn)
    conn.close();


main_func()


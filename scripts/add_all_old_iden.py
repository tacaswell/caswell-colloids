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


_params = ['threshold',
 'p_rad',
 'hwhm',
 'd_rad',
 'mask_rad',
 'shift_cut',
 'rg_cut',
 'e_cut']

def extract_iden_prams(comp,conn):
    """Takes a computation ID, extracts the paramters from the hdf file,
    Returns tuple (dataset,comp_key,[params])"""


    (fin,dsk) = conn.execute("select fout,dset_key from comps where comp_key = ?",(comp,)).fetchone()
    F = h5py.File(fin,'r')

    g = F['parameters']['x_%(#)07d'%{'#':comp}]
    pr = [g.attrs[p] for p in _params]
    
    return (dsk,comp) + tuple(pr)

def add_iden_pram_to_db(params,conn):
    """Adds a set of Iden parameters to the database.  Take a tuple of the form
    (dset_key,comp_key,threshold,p_rad,hwhm,d_rad,mask_rad,shift_cut,rg_cut,e_cut) """

    conn.execute("insert into Iden_prams" +
                 " (dset_key,comp_key,threshold,p_rad,hwhm,d_rad,mask_rad,shift_cut,rg_cut,e_cut) " +
                 "values (?,?,?,?,?,?,?,?,?,?);",params)
    conn.commit()


def main_loop(conn):
    comps = conn.execute("select comp_key from comps where function = 'Iden'").fetchall()
    for c in comps:
        p = extract_iden_prams(c[0],conn)
        add_iden_pram_to_db(p,conn)
    pass
                 

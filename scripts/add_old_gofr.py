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
import lib.general
import numpy as np
_params = [
    'nbins',
    'max_range',
    'comp_count'
    ]

def add_table(conn):
    conn.execute('CREATE TABLE gofr_prams (comp_key INTEGER, dset_key INTEGER, iden_key INTEGER,  nbins INTEGER, max_range FLOAT,FOREIGN KEY(dset_key) REFERENCES dsets(key),FOREIGN KEY(comp_key) REFERENCES comps(comp_key),FOREIGN KEY(iden_key) REFERENCES Iden_prams(comp_key))')
    conn.commit()

def add_gofr_pram_to_db(c,conn):
    """Adds a set of gofr_by_plane parameters to the database.  Take a tuple of the form
    (dset_key,comp_key,'nbins','max_range','comp_count') 
    """
    (dset_key ,comp_key ,fname) = c
    print c
    if comp_key == 513 or comp_key == 512:
        return
    F = h5py.File(fname,'r')

    nbins = len(F['gofr' + "_%(#)07d"%{"#":comp_key}]['bin_count'])
    max_rng = np.round(np.mean(np.diff(F['gofr' + "_%(#)07d"%{"#":comp_key}]['bin_edges']))*nbins)
    
    if 'data_file' in F['gofr' + "_%(#)07d"%{"#":comp_key}].attrs.keys():
        fin = F['gofr' + "_%(#)07d"%{"#":comp_key}].attrs['data_file']
        (iden_key,) = conn.execute("select comp_key from comps where fout = ? and function = 'Iden'",(fin,)).fetchone()

    params = (dset_key,comp_key,nbins,max_rng,iden_key)
    conn.execute("insert into gofr_prams" +
                  " (dset_key,comp_key,nbins,max_range,iden_key) " +
                  "values (?,?,?,?,?);",params)

    conn.commit()
    print params
    F.close()

def main_loop(conn):
    comps = conn.execute("select dset_key,comp_key,fout from comps where function = 'gofr'").fetchall()
    for c in comps:
        add_gofr_pram_to_db(c ,conn)
    pass
                 
if __name__ == "__main__":
    conn = lib.general.open_conn()
    main_loop(conn)
    conn.close()

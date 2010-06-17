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


_params = [
    'nbins',
    'max_range',
    'comp_count'
    ]

def add_table(conn):
    conn.execute('CREATE TABLE gofr_by_plane_prams (comp_key INTEGER, dset_key INTEGER,  nbins INTEGER, max_range FLOAT,comp_count INTEGER,FOREIGN KEY(dset_key) REFERENCES dsets(key),FOREIGN KEY(comp_key) REFERENCES comps(comp_key))')
    conn.commit()

def add_gbp_pram_to_db(params,conn):
    """Adds a set of gofr_by_plane parameters to the database.  Take a tuple of the form
    (dset_key,comp_key,'nbins','max_range','comp_count') 
    """

    conn.execute("insert into gofr_by_plane_prams" +
                 " (dset_key,comp_key,'nbins','max_range','comp_count') " +
                 "values (?,?,?,?,?);",params)

    conn.commit()


def main_loop(conn):
    comps = conn.execute("select dset_key,comp_key from comps where function = 'gofr_by_plane'").fetchall()
    for c in comps:
        add_gbp_pram_to_db(c + (1000,100,25),conn)
    pass
                 

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


def add_col(conn):
    conn.execute("alter table dsets add frame_count INTEGER ")
    conn.commit()

def update_frame_dset(conn):
    dlst = conn.execute("select fout,dset_key from comps where function = 'Iden'").fetchall()
    for (fin,key) in dlst:
        F = h5py.File(fin,'r')
        frame_count = len(F.keys())-1

        conn.execute("update dsets set frame_count = ? where key = ?",(frame_count,key,))
        conn.commit()
        F.close()

    

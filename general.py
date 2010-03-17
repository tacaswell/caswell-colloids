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
import trackpy.vis.pov as pov
import trackpy.vis.plots as plots

class cords:
    def __init__(self):
        self.x = None
        self.y = None
        self.z = None
        self.I = None
        self._ind = 0

def open_conn():
    '''Opens the data base at the standard location and returns the connection'''
    return sqlite3.connect('/home/tcaswell/colloids/processed/processed_data.db')

def get_iden_fout(key,conn):
    '''Returns the fout name and the computation number for the iden on the dataset given by key'''

    # 

def get_xyzI_f(h5file,frame,comp_num):
    '''Extracts the x,y,z,I from the given frame and comp_num from the open file handed in '''
    
    pass

def get_xyzI(key,conn,frame):
    '''Returns four vectors for x,y,z,I for the data set given by key'''

    (fname, comp_number) = get_iden_fout(key,conn)

    f = h5py.File(fname,'r')

    cord = get_xyzI_f(f,frame,comp_number)
    f.close()
    
    return cord

    pass

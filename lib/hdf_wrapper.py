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

import h5py
import general as gen

class read_parts:
    def __init__(self,fname,comp_num):
        self.F = h5py.File(fname,'r')
        self.cnum = comp_num

    def __del__(self):
        self.F.close()
        del self.F

    
    def get_cords(self,frame_num):
        return (self.F[gen.ff(frame_num)][gen.fd('x',self.cnum)],
                self.F[gen.ff(frame_num)][gen.fd('y',self.cnum)])
    

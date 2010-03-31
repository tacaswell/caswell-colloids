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


class Cords3D:
    """A class for carrying around sets of 3D data"""
    def __init__(self):
        self.x = None
        self.y = None
        self.z = None
        self.I = None
        
    def __iter__(self):
        return _quad_objects(self.x,self.y,self.z,self.I)

class Cords2D:
    """A class for carrying around sets of 2D data"""
    def __init__(self):
        self.x = None
        self.y = None
        self.I = None
        
    def __iter__(self):
        return _triple_objects(self.x,self.y,self.I)

    
def _triple_objects(a,b,c):
    """ generator for triples"""
    for i in range(len(a)):
        yield (a[i],b[i],c[i])
def _quad_objects(a,b,c,d):
    """ generator for quads"""
    for i in range(len(a)):
        yield (a[i],b[i],c[i],d[i])

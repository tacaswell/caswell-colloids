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
from __future__ import division

import numpy as np


class gofr_comp:
    def __init__(self,nbins,max_r):
        self.nbins = nbins
        self.vals = np.zeros(self.nbins)
        self.max_r = max_r
        self.min_r = 0
        self.bin_edges = np.linspace(self.min_r,self.max_r,self.nbins+1)
        self.add_count = 0
        self.scaleing = self.nbins /(self.max_r - self.min_r)
    def add_particle(self,cord,list):
        """Adds a particle to the g(r) computation by taking distances
        against list"""

        cord = np.array(cord)
        
        for p in list:
            p = np.array(p)
            r = np.sqrt(np.sum((p-cord)**2))
            
            if r >= self.max_r or r < self.min_r or r ==0:
                continue
            r -= self.min_r
            
            self.vals[np.floor(r*self.scaleing)]+=1

        self.add_count+=1


def indx_to_cord_2D(indx,hash_size):

    if indx <0 or (not indx < np.prod(hash_size) ):
        raise "Index out of range"
    
    return (indx%hash_size[0],indx//hash_size[0])
    
def cord_to_indx_2D(cord,hash_size):
    """Converts coordinate position to an index.  """
    assert len(cord) == len(hash_size)
    cord = np.array(cord)
    if any(cord  >= hash_size) or any(cord < 0):
        print cord
        print hash_size
        raise "cord out of range"
    return int(cord[0] + cord[1] * hash_size[0])

def hash_fun(cord,box_size,hash_size):
    return cord_to_indx_2D(np.floor(np.array(cord)/box_size),hash_size)



class hash_case:
    def __init__(self,dims,box_size):
        self.dims = dims
        self.hash_dims = np.ceil(np.array(dims)/box_size)
        self.box_size = box_size
        self.hash_table = [[] for j in range(int(np.prod(self.hash_dims)))]
    def get_region(self,center, rrange):
        center = np.array(center)
        region = []
        shifts = [np.array([j,k]) for j in range(-rrange,rrange + 1) for k in range(-rrange,rrange + 1)]
        boxes =  [region.extend(self.hash_table[cord_to_indx_2D(center + s,self.hash_dims)]) for s in shifts]
        return region
    def add_particle(self,cord):
        self.hash_table[hash_fun(cord,self.box_size,self.hash_dims)].append(cord)
        
    def compute_corr(self,buf,fun):
        centers = [(j,k)
                   for j in range(buf,int(self.hash_dims[0]-buf))
                   for k in range(buf,int(self.hash_dims[1]-buf))]
        for c in centers:
            print c
            cur_box = self.hash_table[cord_to_indx_2D(c,self.hash_dims)]
            cur_region = self.get_region(c,buf)
            for p in cur_box:
                fun(p,cur_region)

def link_dumb(h1,h2,m_range):
    """This function takes in a pair of hash case (which in this code
    are single planes) and does linking This does not check to make
    sure links are 1-1"""

    pairs = []
    buf = 1
    # loop over boxes in h1
    centers = [(j,k)
               for j in range(buf,int(h1.hash_dims[0]-buf))
               for k in range(buf,int(h1.hash_dims[1]-buf))]
    for c in centers:
        # get list of particles in current box
        cur_box_lst = h1.hash_table[cord_to_indx_2D(c,h1.hash_dims)]
        # get list of particles in current region in h1
        h2_region = h2.get_region(c,1)
        for sp in cur_box_lst:
            canidate_part = h2_region[np.argmin([dist(sp,p) for p in h2_region])]
            
            if dist(sp,canidate_part) < m_range:
                pairs.append((sp,canidate_part))

    return pairs

def dist(a,b):
    return np.sum(np.sqrt((a-b)**2))


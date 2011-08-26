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

import h5py
import lib.img as li
import lib.vanHove as lvh
import numpy as np


class gofr_comp:
    def __init__(self,nbins,max_r):
        self.nbins = nbins               # number of data points to take
        self.vals = np.zeros(self.nbins) # data array
        self.max_r = max_r              # maximum r to be added
        self.min_r = 0                  # minimum r to be added
        self.bin_edges = np.linspace(self.min_r,self.max_r,self.nbins+1) # locations of bin edges
        self.add_count = 0                                    # number of particles added
        self.scaleing = self.nbins /(self.max_r - self.min_r) # scaling for binning
        self.spat_dim = 2               # number of spatial dimensions
    def add_particle(self,cord,lst):
        """Adds a particle to the g(r) computation by taking distances
        against list"""

        cord = np.array(cord)
        
        for p in lst:
            p = np.array(p)
            r = np.sqrt(np.sum((p-cord)**2))
            
            if r >= self.max_r or r < self.min_r or r ==0:
                continue
            r -= self.min_r
            
            self.vals[np.floor(r*self.scaleing)]+=1

        self.add_count+=1
    def normalize(self):
        """ returns the left bin edges and the normalized values for
        g(r)"""

        
        if self.spat_dim == 2:
            area_norm = np.pi*np.diff(self.bin_edges**2)
            rho = (np.sum(self.vals) + self.add_count)/(np.pi*(self.max_r**2))
        elif self.spat_dim == 3:
            area_norm = 4*np.pi*np.diff(self.bin_edges**3)/3
            rho = (np.sum(self.vals) + self.add_count)/(4*np.pi*(self.max_r**3)/3)
        else:
            raise Exception("not a valid dimension")
        print rho
        print self.add_count
        return self.bin_edges[:-1], self.vals/(area_norm*rho)

def indx_to_cord(indx,hash_size):

    if indx <0 or (not indx < np.prod(hash_size) ):
        raise "Index out of range"

    if len(hash_size)==2:
        cord = (indx%hash_size[0],indx//hash_size[0])
    elif len(hash_size)==3:
        cord = (indx%hash_size[0],(indx//hash_size[0])%hash_size[1],indx//(hash_size[0]*hash_size[1]))
    return cord

    
def cord_to_indx(cord,hash_size):
    """Converts coordinate position to an index.  """
    assert len(cord) == len(hash_size)
    cord = np.array(cord)
    if any(cord  >= hash_size) or any(cord < 0):
        print cord
        print hash_size
        raise Exception("cord out of range")
    if len(cord)==2:
        indx = int(cord[0] + cord[1] * hash_size[0])
    elif len(cord)==3:
        indx = int(cord[0] + cord[1] * hash_size[0] + cord[2] * hash_size[0]*hash_size[1] )
    return indx

def hash_fun(cord,box_size,hash_size):
    return cord_to_indx(np.floor(np.array(cord)/box_size),hash_size)



class hash_case:
    def __init__(self,dims,box_size):
        self.dims = dims
        self.hash_dims = np.ceil(np.array(dims)/box_size)
        self.box_size = box_size
        self.hash_table = [[] for j in range(int(np.prod(self.hash_dims)))]
        self.spat_dims = len(dims)
    def get_region(self,center, rrange):
        center = np.array(center)
        region = []
        if self.spat_dims == 2:
            shifts = [np.array([j,k]) for j in range(-rrange,rrange + 1) for k in range(-rrange,rrange + 1)]
            boxes =  [region.extend(self.hash_table[cord_to_indx(center + s,self.hash_dims)]) for s in shifts]
        elif self.spat_dims ==3:
            shifts = [np.array([j,k,m])
                      for j in range(-rrange,rrange + 1)
                      for k in range(-rrange,rrange + 1)
                      for m in range(-rrange,rrange + 1)]
            boxes =  [region.extend(self.hash_table[cord_to_indx(center + s,self.hash_dims)]) for s in shifts]
        return region
    def add_particle(self,cord):
        self.hash_table[hash_fun(cord,self.box_size,self.hash_dims)].append(cord)
        
    def compute_corr(self,buf,fun):
        if self.spat_dims ==2:
            centers = [(j,k)
                       for j in range(buf,int(self.hash_dims[0]-buf))
                       for k in range(buf,int(self.hash_dims[1]-buf))]
        elif self.spat_dims ==3:
            centers = [(j,k,m)
                       for j in range(buf,int(self.hash_dims[0]-buf))
                       for k in range(buf,int(self.hash_dims[1]-buf))
                       for m in range(buf,int(self.hash_dims[2]-buf))]
        for c in centers:
            cur_box = self.hash_table[cord_to_indx(c,self.hash_dims)]
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
        cur_box_lst = h1.hash_table[cord_to_indx(c,h1.hash_dims)]
        # get list of particles in current region in h1
        h2_region = h2.get_region(c,1)
        if len(h2_region)>0:
            for sp in cur_box_lst:
                canidate_part = h2_region[np.argmin([dist(sp,p) for p in h2_region])]
                
                if dist(sp,canidate_part) < m_range:
                    pairs.append((sp,canidate_part))

    return pairs

def dist(a,b):
    return np.sum(np.sqrt((a-b)**2))





def link_frame_avg(dset_key,conn,frame=25):
    """ Takes in a dset_key a connection, and a frame in the single frame data set.

    It then links all possible averaging combinations against each
    other using link_dumb and returns arrays of the statistics on the
    distribution.
    """
    ce = conn.execute
    iden_lst = ce("select comp_key,frames_avged,fout from Iden where dset_key = ?",dset_key).fetchall()
    F_lst = [h5py.File(i[2],'r') for i in iden_lst]
    cut_pram = {'e_cut':.5,'rg_cut':15,'shift_cut':1.5}
    xy_lst = [li.extract_centers(Fi,frame//c[1],c[0],cut_pram) for (Fi,c) in zip(F_lst,iden_lst)]

    H_lst = []
    # fill list of hash objects
    for j in range(0,len(F_lst)):
        H_lst.append(hash_case([1390,1040],10))

    # fill hash object with particles
    for h,xy in zip(H_lst,xy_lst):
        for x,y in zip(*xy):
           h.add_particle(np.array([x,y]))

    res_std = np.zeros([len(H_lst),len(H_lst)])
    res_mean = np.zeros([len(H_lst),len(H_lst)])
    res_hwhm = np.zeros([len(H_lst),len(H_lst)])
    res_msd = np.zeros([len(H_lst),len(H_lst)])

    # does the commutation
    for j in range(0,len(H_lst)):
        for k in range(j+1,len(H_lst)):
            pairs = link_dumb(H_lst[j],H_lst[k],1)
            res_mean[j][k] = np.mean([s[0]-e[0] for (s,e) in pairs] + [s[1]-e[1] for (s,e) in pairs])
            res_std[j][k] = np.std([s[0]-e[0] for (s,e) in pairs] + [s[1]-e[1] for (s,e) in pairs])
            res_hwhm[j][k] = lvh._vh_hwhm(*(
                np.histogram(
                    [s[0]-e[0] for (s,e) in pairs] +
                    [s[1]-e[1] for (s,e) in pairs],
                    1000)[::-1]))
            res_msd[j][k] = lvh._vh_msd(*(
                np.histogram(
                    [s[0]-e[0] for (s,e) in pairs] +
                    [s[1]-e[1] for (s,e) in pairs],
                    1000)[::-1]))
            

    for f in F_lst:
        f.close()
        del f

    del F_lst

    return res_std,res_mean,res_hwhm,res_msd

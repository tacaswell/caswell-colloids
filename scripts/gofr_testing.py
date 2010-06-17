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


import lib.static_calc as sc
import itertools
from lib.util import cord_pairs
import h5py
import numpy as np
import numpy.random.mtrand as nrm


def make_fuzz_fun(fz_per_step):
    def add_fuzz(p,step):
        return np.array(p) + fz_per_step*step* nrm.standard_normal(2) 
    return add_fuzz

def rnd_pos(p,step):
    return np.around(np.array(p),decimals = step)

def shift_pos(p,step):
    p = np.array(p)
    return (step* np.around(p) + p)/(step + 1)

def test_change(dset_key,conn,fun):
    """Takes in a data set key and computes g(r) with different amounts of noise added
    to the positions returns a list of the g(r)"""
    max_rng = 100
    nbins = 1000    
    buff = 1
    cull_rate = .5
    
    (fname,comp_num) = conn.execute("select fout,comp_key from comps where dset_key =?" +
                                    " and function = 'Iden'",(dset_key,)).fetchone()
    F = h5py.File(fname,'r')
    frame = 10
    
    x = F["/frame%(#)06d"%{"#":frame}+"/x_%(#)07d"%{"#":comp_num}]
    y = F["/frame%(#)06d"%{"#":frame}+"/y_%(#)07d"%{"#":comp_num}]

    parts = [a for a in itertools.izip(x,y) if nrm.random_sample() < cull_rate]

    
    steps = 5
    

    gs = []

    for j in range(0,steps):
        hc = sc.hash_case((np.ceil(np.max(x)),np.ceil(np.max(y))),max_rng)
        gofr = sc.gofr_comp(nbins,max_rng)
        for p in parts:
            hc.add_particle(fun(p,j))
            pass
        hc.compute_corr(buff,gofr.add_particle)

        gs.append(cord_pairs(gofr.bin_edges, gofr.vals))
        
    

    return gs
        
        

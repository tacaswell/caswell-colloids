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
import numpy.random as nr


def run_steps(start,max_accepted_steps,max_steps,generate_step,log_p_fun,T):
    '''Runs MC until there are max_accepted_steps accepted steps or until
    there are max_steps total tried.

    The simulation is started at the point start

    the function log_p_fun returns the log of the probability density.

    Proposed steps are generate using generate step'''

    accpted_steps = 0
    steps = 0
    cur_pos = start
    cur_log_p = log_p_fun(start,T)
    step_chain = [(cur_pos,cur_log_p)]
    
    # main loop
    while (max_accepted_steps > accpted_steps) and (max_steps > steps):
        # proposed step
        prop_pos = generate_step(cur_pos)

        # get the ratio of probability for the two positions 
        prop_log_p = log_p_fun(prop_pos,T)
        
        ## if accpted_steps > max_accepted_steps/2:
        ##     T = T*.5
        ##     print accpted_steps,steps
            
        if prop_log_p > cur_log_p:
            p = 1
        else:
            p = np.exp(prop_log_p - cur_log_p)
        
        # get a uniform random number
        u = nr.uniform(0,1)
        
        # if to accept or not
        if u<p:
            step_chain.append((prop_pos,prop_log_p))
            cur_pos = prop_pos
            cur_log_p = prop_log_p
            accpted_steps +=1
        # increment steps
        steps += 1

    print accpted_steps,steps
    return step_chain

def p_fun_generator(p_list,T,ip_fun,g_scale):
    '''Generates probability densities for a single particle being
    added to the cluster, operates with k = 1 and temperature T.

    uses ip_fun for inter particle potential and a linear potential to
    attract everything to the origin.

    the value of r^2 is passed in to ip_fun to same some calculation time.

    '''

    def log_p_fun(pos):
        r_s = [np.sum((p-pos)**2) for p in p_list]
        return -(np.sum([ip_fun(r) for r in r_s ])) /T

    return log_p_fun



def lj_generator(ep,rmin,r_max):
    '''
    Generates Lenard-Jones penitential functions for use in p_fun_generator

    ep is the over all energy scale, rmin is the location of the
    minimum, rmax is the maximum range of the interaction
    
    '''
    def lj(r):
        if r>r_max:
            return 0
        r_tmp = np.power(rmin**2/r,3)
        r_max_tmp = np.power(rmin**2/r_max,3)
        return ep*( (np.power(r_tmp,2) - 2 * r_tmp) - (np.power(r_max_tmp,2) - 2 * r_max_tmp))
    

    return lj

def gen_step(pos):
    '''Generates the next step '''

    return pos + nr.normal(0,.05,np.shape(pos))


def list_loop(p_list,n,lj):
    chain_stack = []
    for j in range(0,n):
        #lpf = p_fun_generator(p_list,.1,lj,1);
        lpf = p_fun_hash(p_list,2.5,lj)
        prop_point = nr.uniform(-1,1,3)
        prop_point = prop_point/np.sqrt(np.sum(prop_point**2))*(np.power(len(p_list)/(4*.74),1.0/3))
        chain_stack.append(run_steps(prop_point,5000,500000,gen_step,lpf,.1))
        p_list.append(np.array([np.mean([c[0][0] for c in chain_stack[-1][1000:]]),
                                np.mean([c[0][1] for c in chain_stack[-1][1000:]]),
                                np.mean([c[0][2] for c in chain_stack[-1][1000:]])]))
                                                                        
    return p_list,chain_stack
                                                                        
 

class p_fun_hash:
    def __init__(self,p_list,ip_range, ip_fun):
        self.p_list = p_list
        self.ip_range = ip_range
        self.ip_fun = ip_fun
        
        
        self.cur_center_hash = None
        self.cur_center = None
        self.hash_list = None
    def __call__(self,pos,T):
        # check the hash of the position
        tmp_hash = self.compute_hash(pos)
        # if not the same as the current hash value, reload
        if self.cur_center_hash is None or len(self.hash_list) == 0 or (
            (tmp_hash != self.cur_center_hash).any() and np.sum((self.cur_center - pos)**2)>.25 ):
            print self.cur_center_hash
            print tmp_hash
            print pos
            print np.sqrt(np.sum(pos**2))
            self.cur_center_hash  = tmp_hash
            self.cur_center = pos
            self.fill_hash()
            
        # if there are no particles near by, force to center aggressively
        if len(self.hash_list) ==0:
            print 'not in contact with anything'
            return -np.sum(pos**2)*1500/T

        # compute the distances needed
        r_s = [np.sum((p-pos)**2) for p in self.hash_list]
               
        return -(np.sum([self.ip_fun(r) for r in r_s ])) /T


    def fill_hash(self):
        
        self.hash_list = [p for p in self.p_list
                          if ((self.compute_hash(p)>= (self.cur_center_hash-np.ones(self.cur_center_hash.shape)))*
                              (self.compute_hash(p)<= (self.cur_center_hash+np.ones(self.cur_center_hash.shape))))
                          .all()]
        print 'refilled hash_list',len(self.hash_list)
        
        
    def compute_hash(self,pos):
        return np.floor(pos/self.ip_range)


#Copyright 2011 Thomas A Caswell
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
import scipy.odr as sodr
import scipy.optimize as sopt
import numpy.lib.polynomial   as nlp
from scipy import interpolate

T_cuts = (30.459495617262299, 31.902701811433076)
lr_1 = (-0.0107074083705739, 0.83490661832579649)
lr_2 = (-0.077343942915963221, 2.9286950745111477)


def T_to_phi_factory(T_c,T_to_r_func):
    '''factory for generating T->phi/phi* conversion functions'''
    p_c = T_to_r_func(T_c)**3
    def tmp(T):
        '''Converts from T to \phi/\phi* '''
        return (T_to_r_func(T)**3)/p_c
    return tmp


def T_to_phi_err_factory(T_c,T_to_r_func,tem,tep=None):
    '''factory for generating T->phi/phi* conversion functions'''
    if tep is None:
        tep = tem
    p_c = T_to_r_func(T_c)**3
    def tmp(T):
        '''Converts from T to \phi/\phi* '''
        return np.abs((T_to_r_func(T)**3)/p_c - (T_to_r_func(T+tep)**3)/p_c),np.abs((T_to_r_func(T-tem)**3)/p_c - (T_to_r_func(T)**3)/p_c)
    return tmp

def smooth_connect(T1,T2,m1,b1,m2,b2):
    """Generates the coefficients for a cubic fit to smoothly connect two lines

    returns a tuple to be feed in to poly1d"""
    
    return ((b2*T1**3 - 3*b2*T1**2*T2 + 3*b1*T1*T2**2 + 2*m1*T1**2*T2**2 - 
             2*m2*T1**2*T2**2 - b1*T2**3)/(T1 - T2)**3,
            (m2*T1 - m1*T2)/(T1 - T2) + 
            (3*T1*T2*(-2*b1 + 2*b2 - m1*T1 + m2*T1 - m1*T2 + m2*T2))/
            (T1 - T2)**3,(3*b1*T1 - 3*b2*T1 + 2*m1*T1**2 - 2*m2*T1**2 + 
                          3*b1*T2 - 3*b2*T2 + 2*m1*T1*T2 - 2*m2*T1*T2 +
                          2*m1*T2**2 - 2*m2*T2**2)/(T1 - T2)**3,
            (-2*b1 + 2*b2 - m1*T1 + m2*T1 - m1*T2 + m2*T2)/(T1 - T2)**3)[::-1]
       

def bi_linear_T_to_r_factory(T_cuts,lr_1,lr_2):
    '''Factory that returns a function to maps T->r assuming a linear fit'''
    p1 = nlp.poly1d(lr_1)
    p2 = nlp.poly1d(lr_2)
    pc = nlp.poly1d(smooth_connect(*(T_cuts + lr_1 + lr_2)))
    def tmp(T):
        '''Map T -> r using a bi linear fit with a cubic interpolation
        between them'''
        def local(T):
            if T < T_cuts[0]:
                return p1(T)
            elif(T<T_cuts[1]):
                return pc(T)
            else:
                return p2(T)
        
        T = np.array(T)
        if len(T.shape) == 0:
            return local(T)
        else:
            return np.array([local(t) for t in T])
    return tmp


def T_to_phi_factory(T_c,T_to_r_func):
    '''factory for generating T->phi/phi* conversion functions'''
    p_c = T_to_r_func(T_c)**3
    def tmp(T):
        '''Converts from T to \phi/\phi* '''
        return (T_to_r_func(T)**3)/p_c
    return tmp

def linear_T_to_r_factory(m,b):
    '''Factory that returns a function to maps T->r assuming a linear fit'''
    def tmp(T):
        '''Map T -> r using a linear fit'''
        return m*T + b
    return tmp


    

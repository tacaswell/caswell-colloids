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

_alpha = 1.4                            # mm^2/s

def heat_kernel(x,t):
    '''Returns the heat kernel '''
    return 1/(np.sqrt(4*np.pi*_alpha*t)) * np.exp(-x**2/(4*_alpha*t ))

def iterate_heat(u0,x,t):
    '''Given a distribution u0 at points x returns the distribution
    evolved forward by t seconds'''
    dy = np.mean(np.diff(x))
    ut = np.array([np.sum(u0 * heat_kernel(np.linspace(x[0]-y,x[-1]-y,len(x)),t)*dy) for
                   y in x])
    return ut
    
              

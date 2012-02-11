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

# taken from 10.1103/PhysRevLett.10.321 PRL, vol 10,num 8, pg 321 (1963)
#from __future__ import division

def Al(l):
    '''from after equation 10'''
    A = 0
    for j in range(0,2):
        A += H(m,j(m,l)) 

    return A/3

def H(m):
    ''' '''
    pass

def f(e):
    '''Definition of f (paragraph under (10) '''
    return (3 + 3*e - e**2)/(4*e**2)

def chi_p(f):
    '''Definition of f (paragraph under (10) '''
    return     (f + np.sqrt(f**2 + 1/8))**(1/3)
    
def chi_m(f):
    '''Definition of f paragraph under (10) '''
    chi = (f - np.sqrt(f**2 + 1/8))
    sg = np.sign(chi)
    return sg*(np.abs(chi)**(1/3))
    
def tl(l,e):
    ''' equation 10'''
    J = np.exp(np.pi*2j/3)
    return 2*e/(1-e) * (-1 + chi_p(f(e))*(J**l) + chi_m(f(e))*(J**(-l)))
    

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


    

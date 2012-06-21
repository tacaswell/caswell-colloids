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



import appwrappy.cpp_wrapper as cw

def main(min_step,step_incr,step_count,comps,conn,nbins,max_range):
    '''Computes  a sweep of step max_step and min_track_length of van Hove calculations

    expects comps to be a list of tuples where the first element is a tracking computation key
    '''

    

    # set up loop over steps
    steps = [min_step + j * step_incr for j in range(0,step_count)]
    for s in steps:
        # loop over computations
        for c in comps:
            # build dictionary
            pram_i = {'nbins':nbins,'max_step':s,'min_track_length':s}
            pram_f = {'max_range':max_range}
            # do computation
            cw.do_vanHove(c[0],conn,pram_i,pram_f)

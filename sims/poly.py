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
import lib.plots as lp
import matplotlib

def sim_4_pt(rads):
    '''Takes in 4 radiuses and returns a list of the possible
    permutations of the distance between the tips of the 4 particle
    layout'''

    perms = [[rads[0],rads[1],rads[2],rads[3]],
             [rads[0],rads[2],rads[1],rads[3]],
             [rads[0],rads[3],rads[1],rads[2]],
             [rads[1],rads[2],rads[0],rads[3]],
             [rads[1],rads[3],rads[0],rads[2]],
             [rads[2],rads[3],rads[0],rads[1]]]

    
    
    return [compute_dist_4(p) for p in perms]

def compute_dist_4(rads):
    '''Computes the actual distance '''

    a = compute_height([rads[0],rads[2],rads[3]])
    b = compute_height([rads[1],rads[2],rads[3]])
    return np.sqrt((a[0] + b[0])**2 + (a[1] - b[1])**2)

def compute_height(r):
    ''' computes the height of a triangle with the given sides '''
    a = r[0] + r[1]
    b = r[0] + r[2]
    c = r[1] + r[2]

    d = (b**2 + c**2 - a**2)/(2*c)
    
    return (np.sqrt(b**2 - d**2),d)

def sim_3_pt(rads):
    '''Takes in 3 radius, returns the all the permutations of the 3 arrangement'''

    perms = [rads,
             [rads[1],rads[2],rads[0]],
             [rads[2],rads[0],rads[1]]
             ]

    return [compute_dist_3(p) for p in perms]

def compute_dist_3(rads):
    '''computes actual distance'''
    return rads[0] + 2*rads[1] + rads[2]




def run_sim(scale,count):
    '''generates distances and sigma for the given scale using count
    sets of radii '''
    a = []
    b = []
    if scale == 0:
        for j in range(0,count):
            a.extend( sim_3_pt([0.5]*3))
            a.extend( sim_3_pt([0.5]*3))                        
            b.extend( sim_4_pt([0.5]*4))
    else:
        for j in range(0,count):
            a.extend( sim_3_pt(nr.normal(.5,scale,3)))
            a.extend( sim_3_pt(nr.normal(.5,scale,3))) 
            b.extend( sim_4_pt(nr.normal(.5,scale,4)))

        
    sig = scale/.5

    return (sig,a,b)

def gen_sims(start,stop,steps,count):
    return [run_sim(s,count) for s in np.linspace(start,stop,steps)]


def compute_hists(data):
    return [(sig,np.histogram(ga+gb,np.linspace(.5,3.5,500),new=True),
             np.histogram(ga,np.linspace(.5,3.5,500),new=True),
             np.histogram(gb,np.linspace(.5,3.5,500),new=True))
            for (sig,ga,gb) in data]


def plot_sweep_hist(hists):
    fig = lp.Figure('r','N','second peak sim',
                    func=matplotlib.axes.Axes.step,count=len(hists))
    for (sig,h,ha,hb) in hists:
        fig.plot(h[1][:-1],h[0],label='%.2f'%(sig*100))

def plot_partials(hist):
    (sig,(h,e),(ha,ea),(hb,eb)) = hist
    tmp_fig = lp.Figure('r','N','second peak sim %.2f'%(sig*100),
                func=matplotlib.axes.Axes.step,count=3)
    tmp_fig.plot(e[:-1],h,label='total')
    tmp_fig.plot(ea[:-1],ha,label='3pt')
    tmp_fig.plot(eb[:-1],hb,label='4pt')

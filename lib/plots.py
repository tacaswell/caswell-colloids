#Copyright 2010,2011 Thomas A Caswell
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

import datetime
import os.path

import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt


class cord_pairs:
    def __init__(self,x,y):
        self.x = x
        self.y = y
    def __iter__(self):
        return itertool.izip(x,y)


def set_up_axis(xlabel,ylabel,title,*args,**kwargs):
    fig,ax = set_up_plot()
    add_labels(ax,title,xlabel,ylabel)
    return ax

def set_up_plot():
    fig = plt.figure()
    ax = fig.add_axes([.1,.1,.8,.8])
    ax.hold(True)
    ax.grid(True)
    return fig,ax

def add_labels(ax,title_txt,xlabel_txt,ylabel_txt):
    """Sets the title, xlabel, and ylabel texts on the axis ax """
    ax.set_title(title_txt)
    ax.set_xlabel(xlabel_txt)
    ax.set_ylabel(ylabel_txt)


def save_figure(fname,fig):
    """Saves a figure to the proper date folder """
    spath = '/home/tcaswell/colloids/figures/' + str(datetime.date.today()) + '/'
    if not  os.path.isdir(spath):
        os.makedirs(spath,0755)
    if os.path.isfile(spath+fname+'.png') or os.path.isfile(spath+fname+'.eps'):
        raise Exception("file name exists: " + spath + fname)

    
    fig.savefig(spath + fname+'.eps',format='eps',dpi=400)
    print '[[' + spath + fname+'.eps]]'
    fig.savefig(spath + fname+'.png',format='png')
    print '[[' + spath + fname+'.png]]'
    return spath + fname+'.eps',spath + fname+'.png'


class tac_figure:
    def __init__(self,xlabel,ylabel,title,*args,**kwargs):
        
        self.fig,self.ax = set_up_plot()
        self.leg_hands = []
        self.leg_strs = []
        self.cmap = None
        if 'cmap' in kwargs:
            self.cmap = kwargs['cmap']
        if 'func' in kwargs:
            self.func = kwargs['func']
        else:
            self.func = matplotlib.axes.Axes.plot
        if 'count' in kwargs:
            count = kwargs['count']
            if self.cmap is None:
                self.cmap = cm.get_cmap('winter')
            self.ax.set_color_cycle([self.cmap(j/(count-1)) for j in range(count)] )
        
        
        add_labels(self.ax,title,xlabel,ylabel)
        
    def draw_line(self,x,y,*args,**kwargs):
        if 'label' in kwargs:
            txt = kwargs['label']
            del kwargs['label']
        else:
            txt = str(len(self.leg_hands))
        if not ('lw' in kwargs or 'linewidth' in kwargs):
            kwargs['lw'] = 2

        self.leg_hands.append(self.func(self.ax,x,y,*args,**kwargs))
        self.leg_strs.append(txt)
        self.ax.legend(self.leg_hands,self.leg_strs,loc=0)
        plt.draw()
        
    def err_plot(self,x,y,y_err,x_err,*args,**kwargs):
        if 'label' in kwargs:
            txt = kwargs['label']
            del kwargs['label']
        else:
            txt = str(len(self.leg_hands))
        if not ('lw' in kwargs or 'linewidth' in kwargs):
            kwargs['lw'] = 3
        self.leg_hands.append(self.ax.errorbar(x,y,y_err,x_err,*args,**kwargs))
        self.leg_strs.append(txt)
        self.ax.legend(self.leg_hands,self.leg_strs,loc=0)
        plt.draw()
    def axis(self,x_lim,y_lim):
        self.ax.set_xlim(x_lim)
        self.ax.set_ylim(y_lim)
        plt.draw()


    def save(self,fname):
        save_figure(fname,self.fig)
        
class color_mapper:
    def __init__(self,mn,mx,name = 'jet'):
        self._mn = mn
        self._mx = mx
        self.cmp = cm.get_cmap(name)
    def get_color(self,t):
        return self.cmp((t-self._mn)/(self._mx-self._mn))

def non_i_plot_start():
    istatus = plt.isinteractive();
    print istatus
    if istatus:plt.ioff()
    return istatus

def non_i_plot_stop(istatus):
    if istatus:
        print "displaying figure"
        plt.ion()
        plt.draw()
        plt.show()
        
    else:
        print "closing all figures"
        plt.close('all')


    


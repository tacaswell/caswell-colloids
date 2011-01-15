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

# needed for the wrapper classes
import PIL
import numpy as np

# needed for plotting
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm



class Stack_wrapper:
    def __init__(self,fname):
        '''fname is the full path '''
        self.im  = PIL.Image.open(fname)
        
        self.im.seek(0)
        # get image dimensions from the meta data the order is flipped
        # due to row major v col major ordering in tiffs and numpy
        self.im_sz = [self.im.tag[0x101][0],
                      self.im.tag[0x100][0]]
    
    def get_frame(self,j):
        '''Extracts the jth frame from the image sequence.
        if the frame does not exist return None'''
        try:
            self.im.seek(j)
        except EOFError:
            return None
        return np.reshape(self.im.getdata(),self.im_sz)
        

class Series_wrapper:
    def __init__(self,base_name,padding,ext):
        '''base name includes the full path up to the numbering
        assumes that the numbering is left padded with 0 to padding places
        the extension does not include the .'''
        self.base_name = base_name + '%0' + str(padding) + 'd.' + ext
        

    def get_frame(self,j):
        '''Extracts the jth frame from the image sequence.
        if the frame does not exist return None'''
        try:
            print self.base_name%j
            im = PIL.Image.open(self.base_name%j)
        except IOError:
            return None
        img_sz = im.size[::-1]
        return np.reshape(im.getdata(),img_sz)



        
def plot_centers(img,x,y):
    '''img_wrapper is assumed to be any object that implements
    get_frame(j).

    centers is a sequence of 2 element sequences
    '''
    # local helper functions
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

            
    istatus = non_i_plot_start()
    # this turns off interactive plotting (makes it go a bit faster,
    # not really important here, but habit from other code
    
    fig = plt.figure()
    ax = fig.add_axes([.1,.1,.8,.8])
    ax.hold(True)
    # this flipud is here as an artifact of a mis-matched broken
    # symmetry and because for reason's I do not understand tiff's
    # store images last row first, and my center finding code reads
    # the lines in in order.  I should change my code, but I have too
    # much other stuff written around it now for it to be worth 'fixing'
    ax.imshow(np.flipud(img),interpolation='nearest',cmap=cm.gray)
    ax.plot(x,y,'rx')

    non_i_plot_stop(istatus)
    # turns interactive plotting back on.
    
    pass


def extract_centers(F,frame,comp_num):
    """F is an open hdf object that is under my file scheme (don't ask
    unless you really want to know)"""

    # local helper functions
    def ff(n):
        """Formats frame names """
        return "frame%(#)06d"%{"#":n}

    def fd(str_,n):
        """ formats dset names"""
        return str_ + "_%(#)07d"%{"#":n}

    x = F[ff(frame)][fd('x',comp_num)]
    y = F[ff(frame)][fd('y',comp_num)]

    return x,y

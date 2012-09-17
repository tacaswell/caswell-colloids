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

import os.path

import PIL.Image
import h5py
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np


# needed for the wrapper classes
def parse_mm_xml_string(xml_str):

    def _write(md_dict,name,val):
        if name == "acquisition-time-local" or name == "modification-time-local":
            tmp = int(val[18:])
            val = val[:18] + "%(#)03d"%{"#":tmp}
            val = datetime.datetime.strptime(val,'%Y%m%d %H:%M:%S.%f')
        md_dict[name] = val



    def _parse_attr(file_obj,dom_obj):
        if dom_obj.getAttribute("id") =="Description":
            _parse_des(file_obj,dom_obj)
        elif dom_obj.getAttribute("type") =="int":
            _write(file_obj,dom_obj.getAttribute("id"),int(dom_obj.getAttribute("value")))
        elif  dom_obj.getAttribute("type") =="float":
            _write(file_obj,dom_obj.getAttribute("id"),float(dom_obj.getAttribute("value")))
        else: 
            _write(file_obj,dom_obj.getAttribute("id"), dom_obj.getAttribute("value").encode('ascii'))

    def _parse_des(file_obj,des_obj):
        des_string = des_obj.getAttribute("value")
        des_split = des_string.split("&#13;&#10;")

        for x in des_split:
            tmp_split = x.split(":")
            if len(tmp_split) ==2:
                _write(file_obj,tmp_split[0],tmp_split[1].encode('ascii'))

    dom = xml.dom.minidom.parseString(xml_str)

    props = dom.getElementsByTagName("prop")
    f = dict()
    for p in props:
        _parse_attr(f,p)

    return f


# needed for plotting

class Stack_wrapper:
    def __init__(self,fname):
        '''fname is the full path '''
        self.im  = PIL.Image.open(fname)

        self.im.seek(0)
        # get image dimensions from the meta data the order is flipped
        # due to row major v col major ordering in tiffs and numpy
        self.im_sz = [self.im.tag[0x101][0],
                      self.im.tag[0x100][0]]
        self.cur = self.im.tell()

    def get_frame(self,j):
        '''Extracts the jth frame from the image sequence.
        if the frame does not exist return None'''
        try:
            self.im.seek(j)
        except EOFError:
            return None
        
        self.cur = self.im.tell()
        return np.reshape(self.im.getdata(),self.im_sz).astype('uint16')
    def __iter__(self):
        self.im.seek(0)
        self.old = self.cur
        self.cur = self.im.tell()
        return self

    def next(self):
        try:
            self.im.seek(self.cur)
            self.cur = self.im.tell()+1
        except EOFError:
            self.im.seek(self.old)
            self.cur = self.im.tell()
            raise StopIteration
        return np.reshape(self.im.getdata(),self.im_sz)

    def get_meta(self,j):
        cur = self.im.tell()
        if cur != j:
            self.im.seek(j)
            xml_str = im_wrap.im.tag[270]
            self.im.seek(cur)
        else: 
            xml_str = im_wrap.im.tag[270]
        return parse_xml_string(xml_str)
class Series_wrapper:
    def __init__(self,base_name,ext,padding = None):
        '''base name includes the full path up to the numbering
        assumes that the numbering is left padded with 0 to padding places
        the extension does not include the .'''
        if padding is None:
            im_guess = len(os.listdir(os.path.dirname(base_name)))
            padding = int(np.log10(im_guess))+1
            
        self.base_name = base_name + '%0' + str(padding) + 'd.' + ext
        

    def get_frame(self,j):
        '''Extracts the jth frame from the image sequence.
        if the frame does not exist return None'''
        try:
            print self.base_name%j
            # added extra +1 to cope with mm starting naming at 1, not the god-given 0
            im = PIL.Image.open(self.base_name%(j+1))
        except IOError:
            print "didn't find the file"
            return None
        img_sz = im.size[::-1]
        return np.reshape(im.getdata(),img_sz).astype('uint16')

def return_img_wrap(dset_key,conn):
    (fname,ftype) = conn.execute("select fname,ftype from dsets where dset_key = ?",dset_key).fetchone()
    if ftype == 1:
        return Stack_wrapper(fname)
    elif ftype ==2:
        return Series_wrapper(fname,'tif')
    raise Exception("not a valid format type")
def plot_centers_simple(iden_key,conn,frame,ax=None,alpha=None):
    ''' Function that does all of the look up for you '''
    
    (iden_fname,dset_key,
     sname,img_fname,ftype) = conn.execute("select fout,iden.dset_key,sname,fname,ftype "+
                                           " from iden inner join dsets on "+
                                           "dsets.dset_key = iden.dset_key where comp_key = ?",
                                           (iden_key,)).fetchone()
    
    cpn =  ['e_cut','rg_cut','shift_cut'] 
    cp = [.5,15,1.5]
    cut_pram = dict(zip(cpn,cp))
    cut_pram = None
    if ftype == 1:
        im_wrap = Stack_wrapper(img_fname)
        pass
    elif ftype ==2:
        im_wrap = Series_wrapper(img_fname,'tif')
        pass
    F = h5py.File(iden_fname,'r')

    plot_centers(F,im_wrap,frame,iden_key,ax,alpha,cut_pram)

    F.close()
    del F


def plot_centers(F,s_wrapper,frame,comp_key,ax=None,alpha=None,cut_prams=None):
    '''Function that automates some of the extraction steps '''

    x,y = extract_centers(F,frame-1,comp_key,cut_prams)
    img = s_wrapper.get_frame(frame)

    ax = _plot_centers(img,x,y)
    
    ## x_trim,y_trim = [],[]
    ## thresh = 2
    
    ## for x1,y1 in zip(x,y):
    ##     fail_flag = False
    ##     n_count = 0
    ##     for x2,y2 in zip(x,y):
            
    ##         d = np.sqrt((x1-x2)**2 + (y1-y2)**2)
    ##         if d< 25:
    ##             n_count +=1
    ##         if n_count > (thresh+1):
    ##             fail_flag = True
    ##             continue
    ##     if not fail_flag:
    ##         x_trim.append(x1)
    ##         y_trim.append(y1)
    ## print len(x_trim),len(y_trim)
    ## ax.plot(x_trim,y_trim,'r.')
    
def _plot_centers(img,x,y,ax=None,alpha=None):
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
    if ax is None:
        fig = plt.figure()
        ax = fig.add_axes([.1,.1,.8,.8])
        
    ax.hold(True)
    # this flipud is here as an artifact of a mis-matched broken
    # symmetry and because for reason's I do not understand tiff's
    # store images last row first, and my center finding code reads
    # the lines in in order.  I should change my code, but I have too
    # much other stuff written around it now for it to be worth 'fixing'
    ax.imshow(np.flipud(img),interpolation='nearest',cmap=cm.gray,alpha=alpha)
    ax.plot(x,y,'r.')

    non_i_plot_stop(istatus)
    # turns interactive plotting back on.
    
    return ax


def extract_centers(F,frame,comp_num,cut_pram = None):
    """F is an open hdf object that is under my file scheme (don't ask
    unless you really want to know)"""

    # local helper functions
    def ff(n):
        """Formats frame names """
        return "frame%(#)06d"%{"#":n}

    def fd(str_,n):
        """ formats dset names"""
        return str_ + "_%(#)07d"%{"#":n}
    
    x = F[ff(frame)][fd('x',comp_num)][:]
    y = F[ff(frame)][fd('y',comp_num)][:]

    if cut_pram is not None:
        R =  F[ff(frame)][fd('R2',comp_num)][:]
        e =  F[ff(frame)][fd('eccentricity',comp_num)][:]
        dx =  F[ff(frame)][fd('x_shift',comp_num)][:]
        dy =  F[ff(frame)][fd('y_shift',comp_num)][:]
        x,y = filter_centers(x,y,R,e,dx,dy,cut_pram)
    return x,y

def filter_centers(x,y,R,e,dx,dy,cut_pram):
    """Filters centers using the standard cut limits """

    indx = np.array([True] * len(x))

    if 'shift_cut' in cut_pram:
        tmp_indx = np.sqrt((dx)**2 + (dy)**2) < cut_pram['shift_cut']
        indx = np.logical_and(indx,tmp_indx)
    if 'rg_cut' in cut_pram:
        tmp_indx = R < cut_pram['rg_cut']
        indx = np.logical_and(indx,tmp_indx)
    if 'e_cut' in cut_pram:
        tmp_indx = e < cut_pram['e_cut']
        indx = np.logical_and(indx,tmp_indx)


    return x[indx],y[indx]

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

import h5py
import numpy as np


class Gofr_wrapper (object):
    """ This is a class to wrap around the gofr data files to make the
    keeping track of the hdf files simpler and to facilitate
    generalization of the mess of functions that are currently in
    gofr.py

    Everything interesting is stored in local variables.

    """
    def __init__ (comp_key,conn):
        """Does all the extraction

        Argos
        comp_key -- 1-tuple with the key of the gofr computation to be extracted
        conn -- db object
        """

        # initial set up
        self.gofr = None
        self.bins = None
        self.T = None
        self.scale_factor = None
        self.rho = None

        # constants
        self.scale_factor = 6.45/60   # TODO this should be extracted 
        dset_names = ['bin_count', 'bin_edges']
        # do sanity checks on comp_key
        
        res = conn.execute("select fout from gofr where comp_key = ?",
                           comp_key).fetchall()
        if not len(res) == 1:
            print len(res)
            raise util.dbase_error("error looking up computation")
        
        try:
            F = h5py.File(res[0][0],'r')
            g = F["gofr_%(#)07d"%{"#":comp_num}]
            # extract the bin values
            self.gofr = np.array(g[dset_names[0]])

            # extract the left bin edges and self.scale_factor
            self.bins = np.array(g[dset_names[1]])*self.scale_factor
            # this is to get the value of the middle of the bins, not
            # the left edges
            bin_steps = np.diff(self.bins)
            self.bins += np.concatenate((bin_steps,[np.mean(bin_steps),]))

            # extract rho
            self.rho = g.attrs['rho']/(self.scale_factor**2)

            # temperature
            if 'temperature' in  g.attrs:
                t = g.attrs['temperature']
            else:
                # hack to get the file format from 2010-05-24 to work
                (data_fname,) = conn.execute(
                    "select fname from dsets where dset_key in " +
                    "(select dset_key from gofr where comp_key = ?)",
                    comp_key).fetchone()

                pos_temps = re.findall('[23][0-9]-[0-9]',data_fname)
                if len(pos_temps) > 0:
                    t = float(pos_temps[0].replace('-','.'))
                elif data_fname.find('room')>0:
                    t = 27.0

                else:
                    pos_temps = re.findall('[23][0-9]\.[0-9]',data_fname)
                    if len(pos_temps) > 0:
                        t = float(pos_temps[0].replace('-','.'))
                    else:
                        print data_fname
                        raise Exception("can't figure out temperature")

        finally:
            # clean up the h5 file
            F.close()
            del F
        
        # do extraction

        pass

class Gofr_bp_wrapper (object):
    """ Clone of above to deal with by plane files use a setter
    function to select which plane to extract, init takes argument for
    initial loaded plane

    make the internals look identical to Gofr_wrapper
    """


##################
#helper functions#
##################
def find_peaks_fit(bin_centers,gofr,dfun,p):
    """Looks for peaks based on hints from the derivative of the expected function and
    the fit parameters
    returns 4 lists, one entry per peak found
    lmax [(location, height),]
    lmin [(location, depth),]
    diffs [expected-actual, ] (min and max interleaved)
    sz [coeff in quad,]
    """
    
    dfun = fitting.fun_flipper(dfun)
    wind = 20;

    lmax = []
    lmin = []
    diffs = []
    sz = []
    # find the first max
    indx = np.argmax(gp.y[15:])+15
    pfit = fit_quad_to_peak(bin_centers[indx-wind:indx+wind],gp.y[indx-wind:indx+wind])

    lmax.append((pfit.beta[1],pfit.beta[2]))
    diffs.append(bin_centers[indx]-pfit.beta[1])
    # get expected spacing from
    e_spc = (2*np.pi)/p[1]
    
    print e_spc
    
    cur_state = np.sign(pfit.beta[0])
    sz.append(pfit.beta[0])
    #start at first peak + 1/4 e_spc
    cur_pos = pfit.beta[1] + e_spc/4
    # step over x range in e_spc/2 steps
    max_pos = np.max(bin_centers)
    print max_pos
    print cur_pos
    while cur_pos < max_pos:
        # get zero crossing from dfun
        try:
            crit_p = sopt.brentq(dfun,cur_pos,cur_pos + e_spc/2,args=p)
            #print 'diff from expected ' + str(crit_p - (cur_pos + e_spc/4))
        except ValueError:
            print "no zero found"
            break
        # convert value into indx
        indx = gen.val_to_indx(bin_centers,crit_p)
        # pick out window around that box
        pfit = fit_quad_to_peak(bin_centers[indx-wind:indx+wind],gofr[indx-wind:indx+wind])
        

        # determine if max or min
        # add center/value to max or min output
        if np.abs(pfit.beta[0])<.05:
            print "peak too small"
            break

        print 'diff from crossing and min/max ' + str( crit_p - pfit.beta[1])
        diffs.append(crit_p - pfit.beta[1])
        sz.append(pfit.beta[0])
        
        if pfit.beta[0] >0:
            if cur_state != -1:
                print "found two peaks in a row"
                break
            lmin.append((pfit.beta[1],pfit.beta[2]))
            cur_state = np.sign(pfit.beta[0])
        elif pfit.beta[0]<0:
            if cur_state != 1:
                print "found two troughs in a row"
                break
            lmax.append((pfit.beta[1],pfit.beta[2]))
            cur_state = np.sign(pfit.beta[0])

            
        # increment cur_pos
        cur_pos = crit_p +  e_spc/4
        wind +=1
        pass
    
    return lmax,lmin,diffs,sz

def fit_quad_to_peak(x,y):
    """Fits a quadratic to the data points handed in """
    def quad(B,x):
        return B[0] *(x -B[1]) ** 2 + B[2]

    beta = (0,np.mean(x),y[gen.val_to_indx(x,np.mean(x))])

    data = sodr.Data(x,y)
    model = sodr.Model(quad)
    worker = sodr.ODR(data,model,beta)
    out = worker.run()
    return out


# Copyright 2011, Vinothan N. Manoharan, Thomas G. Dimiduk, Rebecca W. Perry,
# Jerome Fung, and Ryan McGorty, Guangnan Meng
#
# This file is part of Holopy.
#
# Holopy is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Holopy is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Holopy.  If not, see <http://www.gnu.org/licenses/>.
"""
Functions for tracking particles through reconstructions.



.. moduleauthor:: Ryan McGorty <mcgorty@fas.harvard.edu>
.. moduleauthor:: Guangnan Meng
"""

import numpy as np
import scipy as sp
from scipy import ndimage
from scipy.ndimage import sobel
from scipy.ndimage import maximum_position, minimum_position
import os
import holopy.io
from holopy.io import load

class Particles:
    def __init__(self, lambda_noise=1.0, particle_diameter=3.0, neighborcutoff=0.25):
        self.lambda_noise = lambda_noise
        self.particle_diameter = particle_diameter
        self.neighborcutoff = neighborcutoff

        self.pixelsize = np.array([1.,1.,1.])
        
        self.boxcar = np.array(self.particle_diameter/self.pixelsize, dtype=np.int32)
        odd_offset = np.mod(self.boxcar, 2) + 2
        self.boxcar=self.boxcar+odd_offset + 1

        print self.boxcar
        
        [x, y, z] = np.mgrid[-self.boxcar[0]/2+1:self.boxcar[0]/2+1, 
                                 -self.boxcar[1]/2+1:self.boxcar[1]/2+1,
                                 -self.boxcar[2]/2+1:self.boxcar[2]/2+1]
        
        #convolution kernal:
        contrast_kernal=np.ones(self.boxcar)/np.prod(self.boxcar)
        gaussian_kernal=np.exp(-(x**2+y**2+z**2)/(4*self.lambda_noise**2))
        b=gaussian_kernal.sum()
        #Still need to get the kernal right.
        '''
        Cx=numpy.exp(-(numpy.arange(-self.boxcar[0]/2+1,self.boxcar[0]/2+1)**2)/(2*self.lambda_noise**2))
        Cy=numpy.exp(-(numpy.arange(-self.boxcar[1]/2+1,self.boxcar[1]/2+1)**2)/(2*self.lambda_noise**2))
        Cz=numpy.exp(-(numpy.arange(-self.boxcar[2]/2+1,self.boxcar[2]/2+1)**2)/(2*self.lambda_noise**2))
        C=Cx.sum()*Cy.sum()*Cz.sum()
        K0=1./(B*C)-B/numpy.prod(self.boxcar)
        '''
        
        self.bandandoupassttpassfilter=gaussian_kernal/b
        #centroid mask kernel, assume spherical particle shape:
        self.distance=np.array((x*self.pixelsize[0])**2
                                  +(y*self.pixelsize[1])**2
                                  +(z*self.pixelsize[2])**2)
        self.centroidmask= np.array(self.distance<(self.particle_diameter/2.)**2, 
                                       dtype=np.int32)

def localmx3d(stack, diameter=None, separation=None, masscut = 0.05, fix=False):
    """
    Finds local maximums.

    :param stack: 3D array
    :param diameter: particle diameter
    :param separation: distance between particles
    :param masscut: defines threshold, throws away pixels below this
        fraction of the maximum value
        

    """

    pos = np.zeros((1, 3))
    if stack.max()==0:
        print "Warning: Zero intensity."
        return pos
    
    if diameter is None:
        diameter = np.array((5, 5, 5))
  
    diameter = np.ceil(diameter)
    # in case diameter is not set as odd numbers.
    offset = np.mod(diameter+1., 2.)
    diameter = np.array(diameter + offset, dtype=int)
    w = diameter/2
    
    # matrix array to setup masks.
    [x, y, z] = np.mgrid[-w[0]:w[0]+1.,-w[1]:w[1]+1.,-w[2]:w[2]+1.]
    # make up distance matrix and centroid mask
    # act_distance to determine the centroid mask
    act_distance=np.array((x/diameter[0])**2+(y/diameter[1])**2+(z/diameter[2])**2) 
    # centroid_mask will be used in filtering data and finding features.
    centroid_mask = np.array(act_distance<0.25, dtype=np.int32)
    
    if separation is None:
        # the minimum allowable separation between features,
        # default as diameter + 1.
        separation = diameter + 1
    separation = np.ceil(separation)

    if fix:
        stack = np.fix(stack) # throw out decimal points to save computational time??
    
    mc = masscut*stack.max()
    b = np.where(stack<mc)
    gd = ndimage.morphology.grey_dilation(stack, size = diameter)
    gd[b] = -1.
    can_loc = np.where(gd == stack)
    can_loc = np.array(can_loc).transpose()
    candidate_loc = can_loc + w #shift boundary
    
    # Pad zeros around boundaries, to prevent mask out-of-bounds.
    a = np.zeros((stack.shape+ 2*w))
    a[w[0]:-w[0], w[1]:-w[1], w[2]:-w[2]] = stack
    stack = a
    gd = 0; a = 0; b=0; x = 0; y = 0; z = 0; act_distance = 0 # to clear memory
    '''
    # leave at least 3 pixels around center
    mph_msk = np.maximum((diameter-4)/2, 1)
    # print mph_msk
    mph_msk = np.ones(mph_msk)

    # Apply morphology analysis
    # binary erosion approach.
    bi_eros = ndimage.morphology.binary_erosion(stack, structure = mph_msk, iterations=1)
    mph_stk = stack*bi_eros
    candidate_loc = np.where(bi_eros == True)
    
    # sort candidate_loc by pixel intensity.
    sorted_index = np.argsort(stack[candidate_loc])
    candidate_loc = np.array(candidate_loc).transpose()
    '''
    if candidate_loc.shape[0] == 0:
        print 'No pixel meets the criteria!'
        return pos

    # Go through the candidate pixels and locate local maxima
    while candidate_loc.shape[0]>1: 
        loc = candidate_loc[0, :]
        # load around pixels into sampling array
        sample = centroid_mask*stack[loc[0]-w[0]:loc[0]+w[0]+1,
                                     loc[1]-w[1]:loc[1]+w[1]+1,
                                     loc[2]-w[2]:loc[2]+w[2]+1]
        # if the pixel is not the maximum within the centroid mask,
        # kick the loc out and move to next loc
        if stack[loc[0], loc[1], loc[2]] < sample.max():
            candidate_loc = np.delete(candidate_loc, 0, 0)
        # if the pixel is the maximum within the centroid mask, 
        # select the pixel location and kick out neigbors.
        elif stack[loc[0], loc[1], loc[2]] == sample.max():
            current_pos = (loc-w) # coordinates in pixels, minus boundary due to the padding layer.
            pos=np.vstack((pos, current_pos))
            
            # to kick  surround pixels out in candidat_loc
            dis = (candidate_loc - loc + 0.0)/separation
            dis = np.sqrt(np.sum(dis**2, axis=1))
            nuke_ind = np.array(np.where(dis<0.5))
            nuke_ind = np.reshape(nuke_ind, nuke_ind.size)
            candidate_loc = np.delete(candidate_loc, nuke_ind, 0)
            # print candidate_loc.shape[0]
    
    return pos[1:, :]


def findparticles(stack, diameter=None, separation = None, masscut = 0.05, threshold = 0.05, fix=False):
    
    particles = np.zeros((1, 5))
    if stack.max()==0:
        print "Warning: Zero intensity."
        particles=np.vstack((particles, particles))
        return particles
    
    if diameter is None:
        diameter = np.array((5, 5, 5))
    if sp.isscalar(diameter):
        diameter = diameter*sp.ones(3)
        
    # in case the diameter is not set as odd numbers.
    offset = np.mod(diameter+1., 2.)
    diameter = np.array(diameter + offset, dtype=int)
    w = diameter/2
    
    if separation is None:
        # the minimum allowable separation between features,
        # default as diameter + 1.
        separation = diameter + 1
    if sp.isscalar(separation):
        separation = separation*sp.ones(3)

    if fix:
        stack = np.fix(stack) # throw out decimal points to save computational time??        

    # set up masks.  
    # matrix array to setup masks.
    [x, y, z] = np.mgrid[-w[0]:w[0]+1.,-w[1]:w[1]+1.,-w[2]:w[2]+1.]
    # distance will be used in calculating second moment.
    distance=np.array(x**2+y**2+z**2)
    # act_distance to determine the centroid mask
    act_distance=np.array((x/diameter[0])**2+(y/diameter[1])**2+(z/diameter[2])**2) 
    # centroid_mask will be used in filtering data and finding features.
    centroid_mask = np.array(act_distance<0.25, dtype=np.int32)
    
    # find local maxima
    max_loc = localmx3d(stack=stack, diameter = diameter, separation = separation, masscut=masscut)
    max_loc = max_loc + w #shift boundary
    # Pad zeros around boundaries, to prevent mask out-of-bounds.
    a = np.zeros((stack.shape+ 2*w))
    a[w[0]:-w[0], w[1]:-w[1], w[2]:-w[2]] = stack
    stack=a
    a=0; x=0; y=0; z=0;  act_distance=0 # to clear memory
    
    if max_loc.all() == 0:
        print 'No local maximum found!'
        return particles
    
    # Go through local maxima and centroid.
    for loc in max_loc[:, ]:
        #print loc
        # centroid_mask
        sample = stack[loc[0]-w[0]:loc[0]+w[0]+1,
                                     loc[1]-w[1]:loc[1]+w[1]+1,
                                     loc[2]-w[2]:loc[2]+w[2]+1]
        m0 = sample.sum() # integrated brightness.
        # centroiding shift, see J.Col.Int.Sci. 179, 298, (1996), Eq.5.
        shift=(np.dot(sample.sum(axis=2).sum(axis=1), np.arange(-w[0], w[0]+1.)),
               np.dot(sample.sum(axis=0).sum(axis=1), np.arange(-w[1], w[1]+1.)), 
               np.dot(sample.sum(axis=0).sum(axis=0), np.arange(-w[2], w[2]+1.)))/m0
        
        # TO-DO: if the controiding shift is too large
        '''
        counter = 0
        while abs(shift.max())>1.0:
            loc= np.array(loc + shift, dtype=int)
            # centroid_mask
            sample = stack[loc[0]-w[0]:loc[0]+w[0]+1,
                           loc[1]-w[1]:loc[1]+w[1]+1,
                           loc[2]-w[2]:loc[2]+w[2]+1]
            m0 = sample.sum() # integrated brightness.
            # centroiding shift, see J.Col.Int.Sci. 179, 298, (1996), Eq.5.
            shift=(np.dot(sample.sum(axis=2).sum(axis=1), np.arange(-w[0], w[0]+1.)),
                   np.dot(sample.sum(axis=0).sum(axis=1), np.arange(-w[1], w[1]+1.)), 
                   np.dot(sample.sum(axis=0).sum(axis=0), np.arange(-w[2], w[2]+1.)))/m0
            print counter, shift
            counter = counter +1
        '''
        # move along the centroiding shift and optimize around new center.
        m2 = sample*distance
        m2 = m2.sum()/m0 # second momentum
        pos=(loc+shift-w) # coordinates in real units, minus boundary due to the padding layer.
        particles=np.vstack((particles, (pos[0], pos[1], pos[2], m0, m2)))

    # throw out the features with integrated intensity less than threshold.
    m0_cutoff=threshold*particles[:,3].max()
    loc=np.where(particles[:, 3]>m0_cutoff)
    particles = np.reshape(particles[loc, :], (-1, 5))
    '''
    # throw out the features with m2 out of range.
    m2_cutoff = np.average(particles[:, 4]) - 1.5*np.std(particles[:,4])
    loc=np.where(particles[:, 4]>m2_cutoff)
    particles = np.reshape(particles[loc, :], (-1, 5))
    '''
    print str(particles.shape[0])+" particles found."
    return particles

def _centroid_region(volume, threshold, absolute_thresh=None):
    """
    Takes a reconstructed volume and centroids
    the values above a certain threshold.

    :param volume: reconstructed volume
    :type volume: ndarray

    :param threshold: keep only values above this fraction of max
    :type threshold: float

    :param absolute_thresh: keep values about this (not relative)
    :type absolute_thresh: float

    :return: centroid
    :rtype: ndarray
    """
    if absolute_thresh is None:
        max_values_to_keep = volume.max() * threshold
    else:
        max_values_to_keep = absolute_thresh
    volume = volume*(volume>max_values_to_keep)
    return ndimage.center_of_mass(volume)

def _create_mask(sz, is2d=True):
    """
    Create a "mask"

    :param sz: size of mask
    """
    if sp.mod(sz,2)==0:
        print "Size cannot be even."
        return 
    mask_size = int(sz/2)
    if is2d:
        n,m = sp.ogrid[-1*mask_size:mask_size+1, -1*mask_size:mask_size+1]
        mask = n.repeat(sz,1)**2 + m.repeat(sz,0)**2
    else:
        n,m,o = sp.ogrid[-1*mask_size:mask_size+1, -1*mask_size:mask_size+1, -1*mask_size:mask_size+1]
        mask = n.repeat(sz,1)**2 + m.repeat(sz,0)**2 + o**2
    circ_mask = mask.copy()
    circ_mask[sp.where(circ_mask > mask_size**2)] = 0
    return circ_mask
    
def _get_new_volume_bounds(centrd, xy_shape, z_shape, zstep):
    """
    Takes x,y,z coordinate result of centroid
    """
    x_min = centrd[0] - (xy_shape/2)
    x_max = centrd[0] + (xy_shape/2)
    y_min = centrd[1] - (xy_shape/2)
    y_max = centrd[1] + (xy_shape/2)
    z_min = centrd[2] - ((z_shape/2)*zstep)
    z_max = centrd[2] + ((z_shape/2)*zstep)
    return x_min,x_max,y_min,y_max,z_min,z_max

def radius_of_gyration(recon_roi, sz, isTwoDim=True, sumOverZ=True):
    """
    Computes radius of gyration of region of the reconstructed
    area (by default) or volume.
    """
    mask = _create_mask(sz, is2d=isTwoDim)
    low = int(sp.floor(sz/2.))
    high = int(sp.ceil(sz/2.))
    xc,yc,zc = recon_roi[:,:,:].shape
    xc=xc/2; yc=yc/2; zc=zc/2

    #For testing:
    print mask.shape
    print recon_roi.shape
    
    if isTwoDim:
        if sumOverZ:
            #Sum over z-slices
            roi = recon_roi[xc-low:xc+high, yc-low:yc+high,:].sum(axis=2)
        else:
            roi = recon_roi[xc-low:xc+high,yc-low:yc+high,:]
    else:
        if sz>min(roi.shape):
            print "Mask must be smaller than region."
            return False
        roi = recon_roi[xc-low:xc+high, yc-low:yc+high, zc-low:zc+high]
    if isTwoDim:
        if sumOverZ:
            mass = sp.sum(sp.sum(roi,axis=1),axis=0)
            rad_gyr = sp.sum(sp.sum(roi*mask,axis=1),axis=0)
        else:
            roi = roi.swapaxes(0,2)
            mass = sp.sum(sp.sum(roi,axis=1),axis=1)
            rad_gyr = sp.sum(sp.sum(roi*mask,axis=1),axis=1)
    else:
        mass = sp.sum(sp.sum(sp.sum(roi,axis=1),axis=1),axis=0)
        rad_gyr = sp.sum(sp.sum(sp.sum(roi*mask,axis=1),axis=1),axis=0)
    return rad_gyr/mass

def get_stats_of_recons(holo, z_min, z_max, zsteps,
                           x_min, x_max,
                           y_min, y_max,
                           threshold, sz=None, all_at_once=True):
    """
    Reconstructs a volume. Finds centroid and radius of gyration.
    """
    zs,zstep = sp.linspace(z_min, z_max, zsteps, endpoint=True, retstep=True)
    if all_at_once:
        recs = reconstruct(holo, zs)
        recs = recs[:,:,:,0]
    else:
        recs = sp.zeros((holo.shape[0], holo.shape[1], len(zs)))
        for i in range(0,len(zs)):
            rec = reconstruct(holo, zs[i])
            recs[:,:,i] = rec[:,:,:,0]
    recs = abs(recs**2)
    center_relative = _centroid_region(recs[x_min:x_max,y_min:y_max,:], threshold)
    center = sp.array([x_min, y_min, z_min]) + (center_relative*sp.array([1.,1.,zstep]))
    if sz is None:
        sz = (x_max-x_min)/3. #Completely arbitrary
    print center_relative
    radgyr = radius_of_gyration(recs[center[0].round()-sz:center[0].round()+sz+1,
                                     center[1].round()-sz:center[1].round()+sz+1,
                                     :],2*sz+1)
    return sp.hstack((center, radgyr))

def guess_location(holo, approx_x, approx_y, zdists,
                   z_min = 3e-6, z_max = 100e-6, region_size = 50, thresh = 0.6):
    """
    Takes hologram and approximate coordinates of something of interest. Finds
    a better guess of the x,y and z location.
    """

    #If zdists is just a number, just take a guess given reasonable range.
    if sp.isscalar(zdists):
        zdists = sp.linspace(z_min,z_max,num=zdists)

    #If zdists is scalar then, zsteps will be zero
    if sp.isscalar(zdists):
        zstep=0.
    else:
        zstep = zdists[1]-zdists[0]

    rec = reconstruct(holo, zdists)
    x_min = max(0,approx_x-region_size)
    y_min = max(0,approx_y-region_size)
    x_max = min(holo.shape[0], approx_x+region_size)
    y_max = min(holo.shape[1], approx_y+region_size)
    center_relative = _centroid_region(abs(rec[x_min:x_max,y_min:y_max,:,0]**2), threshold=thresh)
    center = sp.array([x_min, y_min, z_min]) + (center_relative*sp.array([1.,1.,zstep]))
    return center

def fast_focus_detection(holo, approx_x, approx_y, zdists, z_min=3e-6, z_max=100e-6, xy_strip_size=5,
                         gradient_filter=True, find_max_int=True):
    """
    Find the z-distance of an object by reconstructing
    slices.
    Guesses for x and y should be better than in function
    guess_location. 
    
    """
    if sp.isscalar(zdists):
        zdists = sp.linspace(z_min, z_max, num=zdists)

    rec_xz = reconstruct(holo[:,approx_y-xy_strip_size:approx_y+xy_strip_size],
                              zdists, recon_along_xz=True)
    rec_yz = reconstruct(holo[approx_x-xy_strip_size:approx_x+xy_strip_size,:],
                              zdists, recon_along_xz=True)

    #Apply gradient filter
    if gradient_filter:
        rec_xz = sobel(abs(rec_xz),axis=1)
        rec_yz = sobel(abs(rec_yz),axis=1)

    if find_max_int:
        location_xz = maximum_position(abs(rec_xz))[1]
        location_yz = maximum_position(abs(rec_yz))[1]
    else:
        location_xz = minimum_position(abs(rec_xz))
        location_yz = minimum_position(abs(rec_yz))
        
    average_z_of_max = sp.mean([zdists[location_xz.round()], zdists[location_yz.round()]])
    return average_z_of_max


def track_listfiles(filelist, background, opticsclass, approx_x, approx_y, approx_z,
                    zstep, number_recs, xy_box_size, 
                    z_minimum=3e-6, z_maximum=100e-6, threshold=0.5, rad_sz=5,
                    use_absolute_threshold=False, invert=True):
    """
    Tracks something from a list of hologram files.
    """
    tracking_data = sp.zeros((len(filelist),5))
    
    zdists = sp.arange(approx_z-(zstep*number_recs/2.),
                       approx_z+(zstep*number_recs/2.),
                       zstep)

    tracking_data[0,0:3] = sp.array([approx_x, approx_y, approx_z])
    
    for i in range(1, len(filelist)):
        x_min,x_max,y_min,y_max,z_min,z_max = _get_new_volume_bounds(tracking_data[i-1,0:3], xy_box_size, number_recs, zstep)
        holo = load(filelist[i], bg=background, optics=opticsclass)
        if invert:
            holo = holo.max()-holo
        stats = get_stats_of_recons(holo, z_min, z_max, number_recs,
                           x_min, x_max,
                           y_min, y_max,
                           threshold, sz=rad_sz, all_at_once=True)
        tracking_data[i,0:4] = stats
        
        
    return tracking_data                         
    
    
        
    
    
    
    
    

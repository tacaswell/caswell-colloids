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

import itertools

import h5py
import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import EllipseCollection

import general as gen
import img
import plots
from general import ff, fd
import static_calc as lsc

class track:
    def __init__(self,start_frame,positions = None):
        if positions is None:
            self.positions = []
        else:
            self.positions = positions 
        self.start_frame = start_frame
        pass
    def __len__(self):
        return len(self.positions)
    def append(self,pos):
        self.positions.append(np.array(pos))
    def msd(self,max_steps):
        max_steps = int(max_steps)
        msd_vec = np.zeros(max_steps)
        for j in np.arange(max_steps)+1:
            msd_vec[j-1] = np.mean(np.sum(np.diff(np.vstack(self.positions[::j])[:,:2],1,0)**2,1))
        return msd_vec

    def trim_track(self,thresh):
        '''Trims any points with I<thresh off either end of the track
        returns via changing the object, returns none if the track is emptied
        returns a copy for good measure'''

        poss = self.positions

        while len(poss)>0 and poss[0][2] < thresh:
            poss.pop(0)
            t.start_frame +=1

        while len(poss)>0 and poss[-1][2] < thresh:
            poss.pop()

        if len(poss) ==0:
            return None

        return self

    def split_track(self,thresh):
        '''Splits a track on positions that have I<thresh

        returns a list of the new tracks
        '''

        t = trim_track(self,thresh)
        if t is None:
            return []
        
        t_list = [self]
        poss = self.positions
        for j in range(len(poss)):
            if poss[j][2] < thresh:
                t.positions = poss[:j]

                t2 = lt.track(t.start_frame + j+1,poss[j+1:])
                t_list.extend(split_track(t2,thresh))
                break

        return t_list
            


class plane_wrapper:
    def __init__(self,g,iden_key,track_key):
        self.x = g[fd('x',iden_key)][:]
        self.y = g[fd('y',iden_key)][:]
        self.I = g[fd('intensity',iden_key)][:]
        self.next_part = g[fd('next_part',track_key)][:]

    def get_particle(self,p_id):
        return np.array((self.x[p_id],self.y[p_id],self.I[p_id])),self.next_part[p_id]
    
class data_wrapper:
    def __init__(self,F,iden_key,track_key):
        self.iden_key = iden_key
        self.track_key = track_key
        self.planes = [plane_wrapper(F[ff(f)],iden_key,track_key) for f in range(0,F.attrs['number-of-planes'])]
        pass

    def get_particle(self,frame,p_id):
        return self.planes[frame].get_particle(p_id)


def D_hist(trk_lst,max_t):
    msd_vecs = [t.msd(max_t) for t in trk_lst]
    X = matrix(arange(1,max_t+1)).T

    fts = np.array([inv(X.T*X)*(X.T*matrix(msd_vecs).T) for v in msd_vecs])[0]
    
    print mean(fts)
    print std(fts)

    return mean(fts),std(fts)

    
def extract_track(DW,frame_num,part_num,trk):
    """starting with particle part_num in frame_num extracts the path
    going forwards"""
    pos,nxt = DW.get_particle(frame_num,part_num)
    trk.append(pos)
        
    
    if  nxt != -1:
        extract_track(DW,frame_num+1,nxt,trk)
    return trk

def extract_all_tracks(track_comp_key,conn,min_track_len):
    '''Extracts all of the tracks in a file'''
    # make sure track id is valid
    # extract comp_id
    (fin,dset_key,iden_key) = conn.execute("select fout,dset_key,iden_key from tracking where " +
                                  "comp_key=?",
                                  (track_comp_key,)).fetchone()
    
    # open file
    F = h5py.File(fin,'r')
    try:
        DW = data_wrapper(F,iden_key,track_comp_key)
        # get list of particles that are of the minimum length or longer
        tracks = [extract_track(DW,p,i,track(p)) for (i,p,l) in
                  itertools.izip( F[fd('tracking',track_comp_key)]['start_particle'],
                                  F[fd('tracking',track_comp_key)]['start_plane'],
                                  F[fd('tracking',track_comp_key)]['length'])
                                  if l > min_track_len  ]
    finally:
        # close hdf file and clean up
        F.close()
        del F

    return tracks

def print_info(F,frame_num,part_num,comp_num):
    """Prints returns (prev_part,next_part,track_id) for the given
    particle and computation"""
    return (F[ff(frame_num)][fd('prev_part',comp_num)][part_num],
            F[ff(frame_num)][fd('next_part',comp_num)][part_num],
            F[ff(frame_num)][fd('track_id',comp_num)][part_num])


def plot_tracks(track_comp_key,region,init_frame,conn):
    """Takes in a track_comp key and a region of the image and plots
    the tracks going forward of all the particles that are in the
    initial frame

    region = [x_offset y_offset x_len y_len]
    """

    # make sure track id is valid
    # extract comp_id
    (fin,dset_key,iden_key) = conn.execute("select fout,dset_key,iden_key from tracking where" +
                                  " comp_key=?",
                                  (track_comp_key,)).fetchone()
    # open file
    F = h5py.File(fin,'r')

    try:
        
        # extract list of particles in ROI
        x = F[ff(init_frame)][fd('x',iden_key)][:]
        y = F[ff(init_frame)][fd('y',iden_key)][:]

        dtime = F[ff(init_frame)].attrs['dtime']

        sp_scale = gen.extract_spatial_calibration(F[ff(init_frame)])

        ind = matplotlib.mlab.find((x>region[0]) * (y>region[1]) \
              * (x<(region[0] + region[2] ) )*( y<(region[1] + region[3])))

        
        
        
        # loop over particles in region and extract tracks
        tracks = [extract_track(F,init_frame,j,track_comp_key,iden_key,[]) for j in ind]
        # set up figure
        (fig,ax) = plots.set_up_plot()

        ax.set_title(' frame: ' + str(init_frame)
                     + ' dtime: ' + str(dtime) + 'ms')

        def trk_len_hash(trk):
            return len(trk)
        def trk_disp_hash(trk):
            return np.sqrt((np.sum(np.array(trk[-1]) - np.array(trk[0]))**2))

        trk_hash = trk_len_hash
        
        t_len = [trk_hash(trk) for trk in tracks]
        cm = plots.color_mapper(min(t_len),max(t_len))
        print (min(t_len),max(t_len))
        
        # loop over tracks and plot
        [ax.plot(np.array([t[0] for t in trk])*sp_scale,
                 np.array([t[1] for t in trk])*sp_scale,
                 '-',
                 color=cm.get_color(trk_hash(trk)))
         for trk in tracks]

        # plot the starting points
        ax.plot(np.array(x[ind])*sp_scale,
                np.array(y[ind])*sp_scale,'xk')
        
        ax.set_aspect('equal')
        plt.draw()
        
    finally:
        # close hdf file and clean up
        F.close()
        del F



def _extract_tracks(F,start_plane,start_indx,iden_key,track_key):

    frame_num = start_plane

    g_x = F[ff(frame_num)][fd('x',iden_key)]
    g_y = F[ff(frame_num)][fd('y',iden_key)]
    g_nxt = F[ff(frame_num)][fd('next_part',track_key)]

    trks = [[(g_x[s],g_y[s],g_nxt[s])] for s in start_indx]

    frame_num += 1

    while ff(frame_num) in F.keys():
        g_x = F[ff(frame_num)][fd('x',iden_key)]
        g_y = F[ff(frame_num)][fd('y',iden_key)]
        g_nxt = F[ff(frame_num)][fd('next_part',track_key)]

        for t in trks:
            s = t[-1][-1]
            if s != -1:
                t.append((g_x[s],g_y[s],g_nxt[s]))
        frame_num += 1
        
    return trks


    
def population_sort_tracks(track_key,conn,frame_start,frame_step, cut_off):
    """This function takes in a track_key, a connection, and a displacement cut off.


    Three lists are return: two lists of tuples that pair initial
    position, final position, and displacement of each track that is
    present in frame frame_start and lasts until at least frame_start
    + frame_step corresponding to above and below the cut off.

    A list of the particles that are present in frame_start and do not
    last long enough is ruterns"""
    
    # defs
    short_list = []
    long_list = []
    die_list = []
    
    # sql stuff
    (iden_key,fout) = conn.execute("select iden_key,fout from tracking inner join comps on tracking.comp_key = comps.comp_key where " +
                                  "tracking.comp_key = ?",    
                                   track_key).fetchone()
    
    # open file
    F = h5py.File(fout,'r')
    try:
        # extract the initial position data
        frame_tmp = F[ff(frame_start)]
        init_pos = zip(frame_tmp[fd('x',iden_key)][:],
                      frame_tmp[fd('y',iden_key)][:],
                      )
        track_id = frame_tmp[fd('track_id',track_key[0])][:]
        init_next_part = frame_tmp[fd('next_part',track_key[0])][:]
        del frame_tmp
        frame_tmp = F[ff(frame_start + frame_step)]
        # extract the final position data
        final_pos = zip(frame_tmp[fd('x',iden_key)][:],
                        frame_tmp[fd('y',iden_key)][:])
        del frame_tmp
        
        # extract the next_particle data from all the frames between the two
        next_part = []
        for j in range(0,frame_step-1):
            next_part.append(F[ff(frame_start + j + 1)][fd('next_part',track_key[0])][:])
            
            
        pass
    finally:
        # make sure we clean up the hdf stuff
        F.close()
        print 'cleaned up'
        del F
    # walk along tracks to sort in to lists
    for pos,nxt_p,tid in zip(init_pos,init_next_part,track_id):
        die_flag = False
        if tid == -1:
            continue
        for j in range(0,frame_step-1):
            if nxt_p == -1:
                die_flag = True
                break
            nxt_p = next_part[j][nxt_p]
        if die_flag or nxt_p == -1:
            die_list.append(pos)
            continue

        fn_pos = final_pos[nxt_p]
        disp = np.array(pos) - np.array(fn_pos)
        
        if np.abs(disp[0]) > cut_off or np.abs(disp[1]) > cut_off:
            long_list.append((pos,fn_pos))
        else:
            short_list.append((pos,fn_pos))
            
    return short_list,long_list,die_list

def plot_population(short_list,long_list,dead_list,ax=None):
    """function do deal with plotting the results of population_sort_tracks """
    if ax is None:
        fig = plt.figure()
        ax = fig.gca()
        
    ax.set_aspect('equal')

    # plot the lost particles 
    ec = EllipseCollection(
        9,
        9,
        0,
        units='x',
        offsets=np.vstack(dead_list),
        transOffset=ax.transData,
        facecolors='',
        edgecolors='m'
        )
      
    ax.add_collection(ec)
    
    # plot less mobile particles
    i,j = zip(*short_list)
    ec = EllipseCollection(
        9,
        9,
        0,
        units='x',
        offsets=np.vstack(i),
        transOffset=ax.transData,
        facecolors='',
        edgecolors='b'
        )
      
    ax.add_collection(ec)
    ax.quiver(*zip(*[t[0]+t[1] for t in  zip(i,[(jj[0]-ii[0],jj[1]-ii[1]) for jj,ii in zip(i,j)])]),
              color='r',scale_units='xy',scale=1,angles='xy')
    
    
    # plot more mobile particles 
    i,j = zip(*long_list)
    ec = EllipseCollection(
        9,
        9,
        0,
        units='x',
        offsets=np.vstack(i),
        transOffset=ax.transData,
        facecolors='',
        edgecolors='g')
    
    ax.add_collection(ec)
    
    ax.quiver(*zip(*[t[0]+t[1] for t in  zip(i,[(jj[0]-ii[0],jj[1]-ii[1]) for jj,ii in zip(i,j)])]),
              color='k',scale_units='xy',scale=1,angles='xy')
    return ax

def plot_population_img(conn,frame,track_key,short_list,long_list,dead_list,ax=None):
    """function do deal with plotting the results of population_sort_tracks """
    ax = plot_population(short_list,long_list,dead_list)
    
    
    
    (iden_fname,dset_key,
     sname,img_fname,ftype) = conn.execute("select fout,dsets.dset_key,sname,fname,ftype "+
                                           " from tracking inner join dsets on "+
                                           "dsets.dset_key = tracking.dset_key where comp_key = ?",
                                           (track_key,)).fetchone()
    
    
    
    
    if ftype == 1:
        im_wrap = img.Stack_wrapper(img_fname)
        pass
    elif ftype ==2:
        im_wrap = img.Series_wrapper(img_fname,'tif')
        pass
    img_f = im_wrap.get_frame(frame)
    ax.imshow(np.flipud(img_f),interpolation='nearest',cmap=cm.gray,alpha=1)

def coarse_grain(track_key,conn,frame_start,forward_step,grid_size):
    """Coarse grains tracking data to detect anisotropic drift """
    
    # sql stuff
    (iden_key,fout) = conn.execute("select iden_key,fout from tracking inner join comps on tracking.comp_key = comps.comp_key where " +
                                  "tracking.comp_key = ?",    
                                   track_key).fetchone()
    
    
    
    # open file
    F = h5py.File(fout,'r')
    
    # extract the data dimensions
    dims = F.attrs['dims'][:]
    hash_dims = dims//grid_size +1
    # define a local has function
    def hash_fun(x,y):
        '''Local hash function '''
        return int(x//grid_size + (y//grid_size )*hash_dims[0])


    # set up data structures
        
    cg_vectors = [np.array([0,0]) for i in range(0,hash_dims.prod())]
    cg_count = np.zeros(hash_dims.prod())
    
    try:
        # extract all of the relevant data
        start_frame = F[ff(frame_start)]

        init_pos = zip(start_frame[fd('x',iden_key)][:],
                       start_frame[fd('y',iden_key)][:],
                       )

        init_next_part = start_frame[fd('next_part',track_key[0])][:]

        end_frame = F[ff(frame_start+forward_step)]
        final_pos = zip(end_frame[fd('x',iden_key)][:],
                       end_frame[fd('y',iden_key)][:],
                       )
        
        next_part = []
        for j in range(0,forward_step-1):
            frame_tmp = F[ff(frame_start + j + 1)]
            next_part.append(frame_tmp[fd('next_part',track_key[0])][:])
            del frame_tmp
        drift_val = (start_frame.attrs[fd('cumulative_displacement',track_key[0])][:]
                     - end_frame.attrs[fd('cumulative_displacement',track_key[0])][:])
        del start_frame
        del end_frame
    finally:
        # clean up hdf
        F.close()
        del F

    # loop over particles in start frame
    for pos,nxt_p in zip(init_pos,init_next_part):
        die_flag = False
        if nxt_p == -1:
            continue
        for j in range(0,forward_step-1):
            if nxt_p == -1:
                die_flag = True
                break
            nxt_p = next_part[j][nxt_p]
        if die_flag or nxt_p == -1:
            continue

        fn_pos = final_pos[nxt_p]
        disp =  np.array(fn_pos) - np.array(pos)
        n = hash_fun(*pos)
        try:
            cg_vectors[n] = cg_vectors[n]  + disp
            cg_count[n] += 1
        except Exception:
            print n
            
            
            
    # average the vectors
    cg_avg = [v for (v,c) in zip( cg_vectors,cg_count)]
    print drift_val
    print cg_vectors[:10]
    cg_avg = [(v )/c - drift_val for (v,c) in zip( cg_vectors,cg_count)]
    print cg_avg[:10]
    # split up the results
    u,v = zip(*cg_avg)
    x,y = np.meshgrid(grid_size*(np.arange(0,hash_dims[0]) + .5),
                      grid_size*(np.arange(0,hash_dims[1]) + .5))
    return x,y,u,v,cg_count


def coarse_grain_dedrift(track_key,conn,frame_start,forward_step,grid_size):
    """Coarse grains tracking data to detect anisotropic drift.  This
    version is slower because it computes the cumulative displacement
    as it goes instead of using the MD in the tracking file.  This
    needs to be done because the was a bug that did not compute the
    displacement."""
    
    # sql stuff
    (iden_key,fout) = conn.execute("select iden_key,fout from tracking inner join comps on tracking.comp_key = comps.comp_key where " +
                                  "tracking.comp_key = ?",    
                                   track_key).fetchone()
    
    
    try:
        # open file
        F = h5py.File(fout,'r')

        # extract the data dimensions
        dims = F.attrs['dims'][:]
        hash_dims = dims//grid_size +1
        # define a local has function
        def hash_fun(x,y):
            '''Local hash function '''
            return int(x//grid_size + (y//grid_size )*hash_dims[0])


        # set up data structures

        cg_vectors = [np.array([0,0]) for i in range(0,hash_dims.prod())]
        cg_count = np.zeros(hash_dims.prod())

    
        # extract all of the relevant data
        try:
            start_frame = F[ff(frame_start)]
            end_frame = F[ff(frame_start+forward_step)]
            
            init_pos = zip(start_frame[fd('x',iden_key)][:],
                           start_frame[fd('y',iden_key)][:],
                           )

            init_next_part = start_frame[fd('next_part',track_key[0])][:]

        
            final_pos = zip(end_frame[fd('x',iden_key)][:],
                           end_frame[fd('y',iden_key)][:],
                           )
            drift_val = (start_frame.attrs[fd('cumulative_displacement',track_key[0])][:]
                         - end_frame.attrs[fd('cumulative_displacement',track_key[0])][:])
        except Exception:
            raise
        finally:
            del start_frame
            del end_frame
        print F
        
        next_part = []
        poses = []
        for j in range(0,forward_step):
            frame_tmp = None
            try:
                frame_tmp = F[ff(frame_start + j + 1)]
                next_part.append(frame_tmp[fd('next_part',track_key[0])][:])
                poses.append(zip(frame_tmp[fd('x',iden_key)][:],
                                 frame_tmp[fd('y',iden_key)][:],
                                 ))
                del frame_tmp
            except Exception:
                print ff(frame_start + j + 1)
                print F
                print F.keys()
                raise
            
                
        
    finally:
        # clean up hdf
        F.close()
        del F

    # loop over particles in start frame
    disp_sum = [(0,0) for j in range(0,forward_step-1)]
    plane_count = np.zeros(forward_step-1)

    for pos,nxt_p in zip(init_pos,init_next_part):
        die_flag = False
        prev_pos = pos
        if nxt_p == -1:
            continue
        for j in range(0,forward_step-1):
            if nxt_p == -1:
                die_flag = True
                break
            try:
                tp = poses[j]
                tmp_pos = tp[nxt_p]
                disp_sum[j] = disp_sum[j] + ( np.array(tmp_pos) -np.array(prev_pos))
                plane_count[j] +=1
                prev_pos = tmp_pos
            except Exception:
                print j,nxt_p
            nxt_p = next_part[j][nxt_p]
        if die_flag or nxt_p == -1:
            continue

        fn_pos = final_pos[nxt_p]
        disp =  np.array(fn_pos) - np.array(pos)
        n = hash_fun(*pos)
        try:
            cg_vectors[n] = cg_vectors[n]  + disp
            cg_count[n] += 1
        except Exception:
            print n

    drift_val = np.sum([v/c for (v,c) in zip(disp_sum,plane_count)],0)
    print drift_val
    print np.mean([(v)/c for (v,c) in zip( cg_vectors,cg_count)],0)
    # average the vectors
    cg_avg = [(v)/c-drift_val for (v,c) in zip( cg_vectors,cg_count)]

    # split up the results
    u,v = zip(*cg_avg)
    x,y = np.meshgrid(grid_size*(np.arange(0,hash_dims[0]) + .5),grid_size*(np.arange(0,hash_dims[1]) + .5))
    return x,y,np.array(u),np.array(v),cg_count


def remove_track_comp(comp_key,conn):
    '''Completely removes a tracking computation'''
    
    (fin,) = conn.execute("select fout from tracking where comp_key=?",
                                  (comp_key,)).fetchone()

    # the order is important to keep the foreign constraints happy
    conn.execute("delete from tracking where comp_key = ?",(comp_key,))
    conn.execute("delete from comps where comp_key = ?",(comp_key,))
    # commit to db, commit before deleting the data as unmarked data is less irritating
    # than non-existing data
    conn.commit()

    
    # open file
    F = h5py.File(fin,'r+')
    
    for group in F.keys():
        if group[:5] == 'frame':
            del F[group][fd('prev_part',comp_key)]
            del F[group][fd('next_part',comp_key)]
            del F[group][fd('track_id',comp_key)]
    for pram in F['parameters'].keys():
        if int(pram.split('_')[-1]) == comp_key:
            del F['parameters'][pram]
    # delete the top level 
    del F[fd('tracking',comp_key)]
    
    F.close()
    del F

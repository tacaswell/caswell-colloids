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

import h5py
import plots
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import itertools
import general as gen
from general import ff
from general import fd



def extract_track(F,frame_num,part_num,track_group,iden_group,trk):
    """starting with particle part_num in frame_num extracts the path
    going forwards"""
    trk.append((F[ff(frame_num)][fd('x',iden_group)][part_num],
                F[ff(frame_num)][fd('y',iden_group)][part_num]))
    
    nxt = F[ff(frame_num)][fd('next_part',track_group)][part_num]
    
    if  nxt != -1:
        extract_track(F,frame_num+1,nxt,track_group,iden_group,trk)
    return trk

def extract_all_tracks(track_comp_key,conn,min_track_len):
    '''Extracts all of the tracks in a file'''
    # make sure track id is valid
    # extract comp_id
    (fin,dset_key) = conn.execute("select fout,dset_key from comps where\
    function='tracking' and comp_key=?",eu
                                  (track_comp_key,)).fetchone()
    (d_temp,) = conn.execute("select temp from dsets where key = ?",
                             (dset_key,)).fetchone()
    (iden_key, ) = conn.execute("select iden_key from tracking_prams where comp_key=?",
                                (track_comp_key,)).fetchone()

    # open file
    F = h5py.File(fin,'r')
    try:
        # get list of particles that are of the minimum length or longer
        t_lst = [(i,p,l) for (i,p,l) in
                 itertools.izip( F[fd('tracking',track_comp_key)]['start_particle'],
                                 F[fd('tracking',track_comp_key)]['start_plane'],
                                 F[fd('tracking',track_comp_key)]['length'])
                 if l > min_track_len]
        tracks = [extract_track(F,p,i,track_comp_key,iden_key,[]) for (i,p,l) in t_lst]
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
    (fin,dset_key) = conn.execute("select fout,dset_key from comps where\
    function='tracking' and comp_key=?",
                                  (track_comp_key,)).fetchone()
    (d_temp,) = conn.execute("select temp from dsets where key = ?",
                             (dset_key,)).fetchone()
    (iden_key, ) = conn.execute("select iden_key from tracking_prams where comp_key=?",
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

        ax.set_title('T: ' + str(d_temp) + ' frame: ' + str(init_frame)
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
    for pos,nxt_p in zip(init_pos,init_next_part):
        die_flag = False
        for j in range(0,frame_step-1):
            if nxt_p == -1:
                die_flag = True
                break
            nxt_p = next_part[j][nxt_p]
        if die_flag or nxt_p ==-1:
            die_list.append(pos)
            continue

        fn_pos = final_pos[nxt_p]
        sq_disp = np.sum((np.array(pos) - np.array(fn_pos))**2)
        if sq_disp > cut_off:
            long_list.append((pos,fn_pos,sq_disp))
        else:
            short_list.append((pos,fn_pos,sq_disp))
        if sq_disp > 100000:
            print pos,fn_pos,nxt_p
    
    return short_list,long_list,die_list

def plot_population(short_list,long_list,dead_list):
    """function do deal with plotting the results of population_sort_tracks """
    fig = plt.figure()
    i,j,d = zip(*short_list)
    ax = fig.gca()
    ax.scatter(*zip(*i),s=10,c='r')
    
    i,j,d = zip(*long_list)
    ax.scatter(*zip(*i))

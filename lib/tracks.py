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

def _ff(n):
    """Formats frame names """
    return "frame%(#)06d"%{"#":n}

def _fd(str_,n):
    """ formats dset names"""
    return str_ + "_%(#)07d"%{"#":n}
def extract_track(F,frame_num,part_num,track_group,iden_group,trk):
    """starting with particle part_num in frame_num extracts the path going forwards """
    trk.append((F[_ff(frame_num)][_fd('x',iden_group)][part_num],F[_ff(frame_num)][_fd('y',iden_group)][part_num]))
    nxt = F[_ff(frame_num)][_fd('next_part',track_group)][part_num]
    if  nxt != -1:
        extract_track(F,frame_num+1,nxt,track_group,iden_group,trk)


def print_info(F,frame_num,part_num,comp_num):
    """Prints returns (prev_part,next_part,track_id) for the given particle and computation """
    return (F[_ff(frame_num)][_fd('prev_part',comp_num)][part_num],
            F[_ff(frame_num)][_fd('next_part',comp_num)][part_num],
            F[_ff(frame_num)][_fd('track_id',comp_num)][part_num])

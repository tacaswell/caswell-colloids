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

# this module adds a column to the diet table and adds a table for format names

import sqlite3


def fill_format_names(conn_new):
    fill_data = [ (1,'mm_stack'),(2,'mm_series')]
    for c in fill_data:
        print c
        conn_new.execute("insert into format_types (ftype_key,ftype_name) values (?,?)",c)
    conn_new.commit()

def fill_dsets(conn_old,conn_new):
    # get data from old database
    dsets_old = conn_old.execute("select dset_key,fname,dtype,ddate,sname from dsets")

    # shove in to new
    for d in dsets_old:
        conn_new.execute("insert into dsets (dset_key,fname,dtype,ddate,sname,ftype) " +
                         "values (?,?,?,?,?,?) ",
                         d + (1,))
    # commit
    conn_new.commit()
    pass


    


def _fill_fun(conn_old,conn_new,table,fields):
    

    val_old = conn_old.execute("select " + ','.join(fields) + " from " + table).fetchall()

    
    for v in val_old:
        conn_new.execute("insert into "+ table +" (" + ','.join(fields) + ") values (" +
                         ','.join(['?']*len(fields)) + ")",v)
    
    conn_new.commit()



def fill_func_names(conn_old,conn_new):
    tbl =    'func_names'
    flds = ('func_key',
            'func_name')
    _fill_fun(conn_old,conn_new,tbl,flds)



def fill_comps(conn_old,conn_new):
    tbl = 'comps'
    flds = ('comp_key',
            'func_key',
            'dset_key')
    
    _fill_fun(conn_old,conn_new,tbl,flds)
    
    
    
def fill_iden(conn_old,conn_new):
    tbl = 'iden'
    flds =('comp_key',
           'dset_key',
           'threshold',
           'hwhm',
           'p_rad',
           'd_rad',
           'mask_rad',
           'top_cut',
           'frames_avged',
           'fout',
           'date')
    _fill_fun(conn_old,conn_new,tbl,flds)

def fill_gofr(conn_old,conn_new):
    tbl = 'gofr'
    flds = ('comp_key',
            'dset_key',
            'iden_key',
            'nbins',
            'max_range',
            'shift_cut',
            'rg_cut',
            'e_cut',
            'fin',
            'fout',
            'date')
    _fill_fun(conn_old,conn_new,tbl,flds)
    
def fill_gofr_by_plane(conn_old,conn_new):
    tbl = 'gofr_by_plane'
    flds = ('comp_key',
            'dset_key',
            'iden_key',
            'nbins',
            'max_range',
            'comp_count',
            'shift_cut',
            'rg_cut',
            'e_cut',
            'fin',
            'fout',
            'date')
    _fill_fun(conn_old,conn_new,tbl,flds)
    
def fill_tracking(conn_old,conn_new):
    tbl = 'tracking'
    flds = ('comp_key',
            'iden_key',
            'dset_key',
            'search_range',
            'shift_cut',
            'rg_cut',
            'e_cut',
            'fin',
            'fout',
            'date')
    _fill_fun(conn_old,conn_new,tbl,flds)
    
def fill_msd_old(conn_old,conn_new):
    tbl = 'msd_old'
    flds = ('comp_key',
            'iden_key',
            'dset_key',
            'search_range',
            'msd_steps',
            'min_track_length',
            'shift_cut',
            'rg_cut',
            'e_cut',
            'fin',
            'fout',
            'date')
    _fill_fun(conn_old,conn_new,tbl,flds)

def fill_msd(conn_old,conn_new):
    tbl = 'msd'
    flds = ('comp_key',
            'iden_key',
            'track_key',
            'dset_key',
            'msd_steps',
            'min_track_length',
            'fin',
            'fout',
            'date')
    _fill_fun(conn_old,conn_new,tbl,flds)
    
def fill_trk_stat(conn_old,conn_new):
    tbl = 'trk_stat'
    flds = ('comp_key',
            'iden_key',
            'dset_key',
            'search_range',
            'steps',
            'hist_bins',
            'hist_range',
            'shift_cut',
            'rg_cut',
            'e_cut',
            'fin',
            'fout',
            'date')
    _fill_fun(conn_old,conn_new,tbl,flds)


def fill_trk_stat_trk(conn_old,conn_new):
    tbl = 'trk_stat_trk'
    flds = ('comp_key',
            'trk_key',
            'dset_key',
            'search_range',
            'steps',
            'hist_binsINTEGER',
            'hist_range',
            'fin',
            'fout',
            'date')
    _fill_fun(conn_old,conn_new,tbl,flds)

def fill_vanHove(conn_old,conn_new):
    tbl = 'vanHove'
    flds = ('comp_key',
            'dset_key',
            'track_key',
            'min_track_length',
            'max_step',
            'max_range',
            'nbins',
            'fin',
            'fout',
            'date')
    _fill_fun(conn_old,conn_new,tbl,flds)
    
def main(conn,conn_new):
    fill_format_names(conn_new)
    fill_func_names(conn,conn_new)
    fill_dsets(conn,conn_new)
    fill_comps(conn,conn_new)
    fill_iden(conn,conn_new)
    fill_gofr(conn,conn_new)
    fill_gofr_by_plane(conn,conn_new)
    fill_tracking(conn,conn_new)
    fill_msd(conn,conn_new)
    fill_vanHove(conn,conn_new)
    fill_trk_stat(conn,conn_new)
    

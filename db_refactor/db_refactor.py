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

import sqlite3


def fill_dsets(conn_old,conn_new):
    # get data from old database
    dsets_old = conn_old.execute("select key,fname,dtype,date,sname from dsets")
    # shove in to new
    for d in dsets_old:
        conn_new.execute("insert into dsets (dset_key,fname,dtype,ddate,sname) values (?,?,?,?,?) ",d)
    # commit
    conn_new.commit()
    pass

def fill_comps(conn_old,conn_new):
    # get data from old database
    comps_old = conn_old.execute("select comp_key,function,dset_key from comps")
    # shove in to new
    for c in comps_old:
        conn_new.execute("insert into comps (comp_key,func_name,dset_key) values (?,?,?)",c)
    # commit
    conn_new.commit()
    

def fill_iden(conn_old,conn_new):

    fields= ['comp_key',
                'dset_key',
                'threshold',
                'hwhm',
                'p_rad',
                'd_rad',
                'mask_rad',
                'top_cut']
    
    fields_in = fields + ['avg_count']
    # get data from the strait iden table
    iden_old = conn_old.execute("select " + ','.join(fields) + " from Iden_prams ")

    fields.append('frames_avged')
    for io in iden_old:
        if( io[0] <580 and io[0] >568):
            conn_new.execute("insert into iden (" +
                             ','.join(fields) + ") values ("
                             + ','.join(['?']*len(fields)) + ")",
                             io + (0,))
        else:
            conn_new.execute("insert into iden (" +
                             ','.join(fields) + ") values ("
                             + ','.join(['?']*len(fields)) + ")",
                             io + (0,))
    # get data from iden_avg
    iden_old = conn_old.execute("select " + ','.join(fields_in) + " from Iden_avg_prams ")
    for io in iden_old:
        conn_new.execute("insert into iden (" +
                         ','.join(fields) + ") values ("
                         + ','.join(['?']*len(fields)) + ")",
                         io)
    
    conn_new.commit()

    pass

def _fill_fun(conn_old,conn_new,t_old,t_new,f_flds,i_flds):

    p_old = conn_old.execute("select " + ','.join(f_flds) + " from " + t_old).fetchall()
    for p in p_old:
        # get the iden parameters, this is only valid because of how
        # the code has been run, always using the values selected when
        # Iden was run, in general this may not work
        i_pram = ()
        iden_key = p[2:3]
        
        if( iden_key[0] <580 and iden_key[0] >568):
            i_type = 'Iden'
        else:
            (i_type,) = conn_old.execute("select function from comps where comp_key = ?",iden_key).fetchone()
        if i_type == 'Iden':

            i_pram += conn_old.execute("select " + ','.join(i_flds) + " from Iden_prams where comp_key = ?",iden_key).fetchone()
        elif i_type =='Iden_avg':
            i_pram += conn_old.execute("select " + ','.join(i_flds) + " from Iden_avg_prams where comp_key = ?",iden_key).fetchone()
        else:
            print p


        conn_new.execute("insert into "+ t_new +" (" + ','.join(f_flds+i_flds) + ") values (" +
                         ','.join(['?']*len(f_flds + i_flds)) + ")",p + i_pram)
    
    conn_new.commit()

def fill_gofr(conn_old,conn_new):
    
    g_flds = ['comp_key',
              'dset_key',
              'iden_key',
              'nbins',
              'max_range',]
    i_flds = ['shift_cut',
              'rg_cut',
              'e_cut']              

    _fill_fun(conn_old,conn_new,'gofr_prams','gofr',g_flds,i_flds)
        
    pass

def fill_gofr_by_plane(conn_old,conn_new):
    
    g_flds = ['comp_key',
              'dset_key',
              'iden_key',
              'nbins',
              'max_range',
              'comp_count']
    i_flds = ['shift_cut',
              'rg_cut',
              'e_cut']              
    _fill_fun(conn_old,conn_new,'gofr_by_plane_prams','gofr_by_plane',g_flds,i_flds)
    pass

def fill_tracking(conn_old,conn_new):
    
    g_flds = ['comp_key',
              'dset_key',
              'iden_key',
              'search_range',]
    i_flds = ['shift_cut',
              'rg_cut',
              'e_cut']        
    _fill_fun(conn_old,conn_new,'tracking_prams','tracking',g_flds,i_flds)
    pass

def fill_msd(conn_old,conn_new):
    
    g_flds = ['comp_key',
              'dset_key',
              'iden_key',
              'search_range',
              'msd_steps',
              'min_track_length']
    i_flds = ['shift_cut',
              'rg_cut',
              'e_cut']
    _fill_fun(conn_old,conn_new,'msd_prams','msd',g_flds,i_flds)
    
    pass

def fill_trk_stat(conn_old,conn_new):
    
    g_flds = ['comp_key',
              'dset_key',
              'iden_key',
              'search_range',
              'steps',
              'hist_bins',
              'hist_range']
    i_flds = ['shift_cut',
              'rg_cut',
              'e_cut']              
    _fill_fun(conn_old,conn_new,'trk_stat_prams','trk_stat',g_flds,i_flds)
    pass

def fill_vanHove(conn_old,conn_new):
    vh_f = ['comp_key',
            'dset_key',
            'track_key',
            'min_track_length',
            'max_step',
            'max_range',
            'nbins']
    
    
    # get data from the strait iden table
    vh_old = conn_old.execute("select " + ','.join(vh_f) + " from vanHove_prams ")
    for vh in vh_old:
        conn_new.execute("insert into vanHove (" + ','.join(vh_f) +
                         ") values ("+ ','.join(['?']*len(vh_f)) + ")",vh)
        
def main(conn,conn_new):
    dr.fill_dsets(conn,conn_new)
    dr.fill_comps(conn,conn_new)
    dr.fill_iden(conn,conn_new)
    dr.fill_gofr(conn,conn_new)
    dr.fill_gofr_by_plane(conn,conn_new)
    dr.fill_tracking(conn,conn_new)
    dr.fill_msd(conn,conn_new)
    dr.fill_vanHove(conn,conn_new)
    dr.fill_trk_stat(conn,conn_new)
    

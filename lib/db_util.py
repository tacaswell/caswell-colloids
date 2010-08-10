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

from util import dbase_error
import sqlite3

def get_fname(key,conn):
    return conn.execute("select fname from dsets where key = ?",(key,)).fetchone()[0]


_table_lookup = {'gofr':'gofr_prams',
                 'Iden':"select * from Iden_prams where comp_key = ?",
                 'gofr_by_plane':'gofr_by_plane_prams'
                 }

def display_prams(comp_key,pram_list,conn):
    """Looks up and prints the parameters for the computation and returns a dictionary
    with the parameters

    comp_key : the computation key

    pram_list : list of the parameters to be read out

    conn : sqlite database connection"""

    tmp_rf = conn.row_factory
    conn.row_factory = sqlite3.Row
    
    fun = conn.execute("select function from comps where comp_key = ?",(comp_key,)).fetchone()
    if fun is None:
        raise dbase_error("comp_key invalid")

    fun_table = _table_lookup[fun[0]]

    row = conn.execute(fun_table,(comp_key,)).fetchone()
    
    pram_dic = {}
    for p in pram_list:
        pram_dic[p] = row[p]

    for k in pram_dic:
        print k + ": " + str(pram_dic[k])

    conn.row_factory = tmp_rf
    return pram_dic
     
    

    

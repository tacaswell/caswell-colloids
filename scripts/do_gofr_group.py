#!/usr/bin/env python

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



import lib.general as gen
import trackpy.cpp_wrapper as cw

def gofr_group(conn,sname):

    res_t = conn.execute("select dsets.key from dsets,comps where dsets.sname = ? and dsets.key = comps.dset_key and comps.function = 'Iden' and dsets.dtype = 't'",(sname,)).fetchall()

    res_z = conn.execute("select dsets.key from dsets,comps where dsets.sname = ? and dsets.key = comps.dset_key and comps.function = 'Iden' and dsets.dtype = 'z'",(sname,)).fetchall()
    for r in res_t:
        print r[0]
        cw.do_gofr(r[0],conn)
    for r in res_z:
        print r[0]
        cw.do_gofr(r[0],conn)
        cw.do_link3D(r[0],conn)
        cw.do_gofr3D(r[0],conn)

def main_loop():
    conn = gen.open_conn();
    gofr_group(conn,'2010-04-26-2')
    gofr_group(conn,'2010-04-26-6-thin')
    gofr_group(conn,'2010-04-26-6')


if __name__ == "__main__":
    main_loop()


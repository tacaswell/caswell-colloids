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
import appwrappy.cpp_wrapper as cw

def gofr_date(conn,date):

    ## res_t = conn.execute("select key from dsets where date = ? and dtype = 't'",(date,)).fetchall()


    ## for r in res_t:
    ##     print r[0]
    ##     cw.do_gofr(r[0],conn)

    res_z = conn.execute("select key from dsets where date = ? and dtype = 'z'",(date,)).fetchall()
    for r in res_z:
        print r[0]
        #cw.do_gofr(r[0],conn)
        #cw.do_link3D(r[0],conn)
        cw.do_gofr3D(r[0],conn)

def main_loop():
    conn = gen.open_conn();
    gofr_date(conn,'2010-05-24')

    conn.close()

if __name__ == "__main__":
    main_loop()


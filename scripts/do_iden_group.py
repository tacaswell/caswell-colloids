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
import sqlite3 as sq
import trackpy.cpp_wrapper as cw


def iden_group(conn,sname):
    res = conn.execute("select key from dsets where sname = ?",(sname,)).fetchall()
    for r in res:
        print r[0]
        cw.do_Iden(r[0],conn)

def main_loop():
    conn = gen.open_conn();
    iden_group(conn,'2010-04-26-2')
    #iden_group(conn,'2010-04-26-6-thin')
    #iden_group(conn,'6b')

if __name__ == "__main__":
    main_loop()



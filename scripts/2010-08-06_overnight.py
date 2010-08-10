import lib.general as gen
import trackpy.cpp_wrapper as cw

conn = gen.open_conn()
ce = conn.execute

## set_lst = ce("select key from dsets where date = '2009-04-07'").fetchall()
## for s in set_lst:
##     cw.do_Iden(s[0],conn)
##     cw.do_Iden_avg(s[0],conn,3)

pram_f = {'max_range':100,'search_range':5,'box_side_len':5,'hist_range':5}
pram_i = {'nbins':2000,'comp_count':200,'min_track_length':20,'msd_steps':300,'comp_count':300,'steps':300,'hist_bins':250}
cpp
iden_lst = ce("select comp_key from comps where dset_key in (select key from dsets where date = '2009-04-07')").fetchall()
for c in iden_lst:
    cw.do_gofr_by_plane(c[0],conn,pram_i,pram_f)
    cw.do_msd(c[0],conn,pram_i,pram_f)
    

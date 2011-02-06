import lib.general as gen
import lib.plots as plts
import lib.msd as lm

conn = gen.open_conn("/home/tcaswell/colloids/proc_db.db")
ce = conn.execute

steps = [5,10,15,20,25]
styles = ['x','^','o','<','>']

msds = [ce("select comp_key from msd where date = '2010-12-08' and msd_steps = ?",(c,)).fetchall()
        for c in steps]

fits = [ [lm.fit_msd(c[0],conn) for c in m]
         for m in msds]



fig = plts.Figure('T [C]',r'r [$\mu$m]','NIPAM radius vs Temperature')

for fts,sty,lab in zip(fits,styles,steps):
    fig.plot([f[1] for f in fts],[f[2] for f in fts],sty,label=str(lab) + ' steps')


fig.save('test')

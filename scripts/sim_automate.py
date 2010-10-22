from __future__ import division

import sims.mc_helper as sm
import h5py
import numpy as np
import numpy.random as nr
import sys

def save_to_hdf(fname,p,chains,run,count):
    Fout = h5py.File(fname,'r+')
    run_name = "run_%(#)04d"%{"#":run}
    if run_name not in Fout.keys():
        raise Exception("file not inited properly")

    grp = Fout[run_name]

    plist_name = "plist_%(#)04d"%{"#":count}
    pos_name = "pos_%(#)04d"%{"#":count}
    log_name = "log_p_%(#)04d"%{"#":count}

    grp.create_dataset(plist_name,data=p)

    logs = [c[1] for c in chains]
    pos = [c[0] for c in chains]
    pos = np.vstack(pos)
    
    grp.create_dataset(pos_name,data=pos)
    grp.create_dataset(log_name,data=logs)

    Fout.close()
    del Fout


def create_hdf(fname,run,T,ij_e,rmin,r_max):
    Fout = h5py.File(fname,'a')
    run_name = "run_%(#)04d"%{"#":run}
    grp = Fout.create_group(run_name)
    grp.attrs['T'] = T
    grp.attrs['rmin'] = rmin
    grp.attrs['r_max'] = r_max
    grp.attrs['ij_e'] = ij_e
    Fout.close()
    del Fout

def main(fname,rmin,r_max,ij_e,T,part_to_add):
    lj = sm.lj_generator(ij_e,rmin,r_max)
    
    
    for run in range(0,50):
        create_hdf(fname, run,T,ij_e,rmin,r_max)
        p_list = [np.array([0,0,0]),np.array([0,1,0]),np.array([np.sqrt(3)/2,.5,0])]

        count = 0
        #save_to_hdf(fname,p_list,[],run,count)

        for j in range(0,part_to_add):
            print len(p_list)
            count += 1
            lpf = sm.p_fun_hash(p_list,r_max,lj)
            prop_point = nr.uniform(-1,1,3)
            prop_point = prop_point/np.sqrt(np.sum(prop_point**2))*(np.power(len(p_list)/(4*.74),1.0/3))

            chain = sm.run_steps(prop_point,5000,1000000,sm.gen_step,lpf,T)

            save_to_hdf(fname,p_list,chain,run,count)
            
            p_list.append(np.array([np.mean([c[0][0] for c in chain[1000:]]),
                                np.mean([c[0][1] for c in chain[1000:]]),
                                np.mean([c[0][2] for c in chain[1000:]])]))

        print 'finished run'

if __name__ == '__main__':
    fname = sys.argv[1]
    T = float(sys.argv[2])
    main(fname,1,2.5,10,T,150)

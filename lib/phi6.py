import h5py
import numpy as np
import img as li

import os

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from general import ff,fd

# change to take 
def _plot_file_frame_phi6(key,conn,fr_num,fnameg=None):
    '''

    '''
    (fname,p_comp) = conn.execute("select fout,comp_key from comps where function = 'phi6' and dset_key = ?;",(key,)).fetchone()
    (comp_number,) = conn.execute("select comp_key from comps where function = 'Iden' and fout = ?",(fname,)).fetchone()

    (sname,stype,temp) = conn.execute("select sname,dtype,temp from dsets where key = ?",(key,)).fetchone()

    
    f = h5py.File(fname,'r')

    fc = f.attrs['number-of-planes']
    if fr_num <fc:
        fr_num = fc-1
        

    x = f[ff(fr_num)][fd('x',iden_key)][:]
    x = f[ff(fr_num)][fd('y',iden_key)][:]
    phi = f[ff(fr_num)][fd('scaler_order_parameter',iden_key)]

    
    phir = np.zeros(phi.shape[0])
    phii = np.zeros(phi.shape[0])

    for j in range(0,phi.shape[0]):
        phir[j] = phi[j][0]
        phii[j] = phi[j][1]
    print np.mean(abs(phir))
    print np.mean(abs(phii))
    print np.mean((phir))
    print np.mean((phii))

    
    # check interactive plotting and turn it off
    istatus = plt.isinteractive();
    print istatus
    if istatus:plt.ioff()

    leg_hands = []
    leg_strs = []

    fig = plt.figure()
    ax = fig.add_axes([.1,.1,.8,.8])
    ax.hold(True)
    ax.set_aspect('equal')

    ax.quiver(x,y,phir,phii,scale = 100,headlength = 2,headwidth= 1.5,headaxislength=2)
    ax.set_title(sname + " temp: " + str(temp))
#    plt.plot(x,y,'ro')


    
    if not fnameg == None:
        f_path = '/home/tcaswell/python/figures/' + sname + '/'
        if not os.path.isdir(f_path):
            os.makedirs(f_path)
        
        fnameg = f_path + str(key) + '_phi6.png'
        fig.savefig(fnameg,dpi = 500)
        
        
    if istatus:
        print "displaying figure"
        plt.ion()
        plt.show()
    else:
        print "closing figure"
        plt.close(fig)


def _plot_file_frame_nsize(key,conn,fr_num,fnameg=None):
    '''

    '''
    (fname,p_comp) = conn.execute("select fout,comp_key from comps where function = 'phi6' and dset_key = ?;",(key,)).fetchone()
    (comp_number,) = conn.execute("select comp_key from comps where function = 'Iden' and fout = ?",(fname,)).fetchone()

    (sname,stype,temp) = conn.execute("select sname,dtype,temp from dsets where key = ?",(key,)).fetchone()

    
    f = h5py.File(fname,'r')
    
    x = f[ff(fr_num)][fd('x',iden_key)]
    x = f[ff(fr_num)][fd('y',iden_key)]
    ns = f[ff(fr_num)][fd('neighborhood_size',iden_key)]

    print np.mean(ns)
    
    # check interactive plotting and turn it off
    istatus = plt.isinteractive();
    print istatus
    if istatus:plt.ioff()

    leg_hands = []
    leg_strs = []

    fig = plt.figure()
    ax = fig.add_axes([.1,.1,.8,.8])
    ax.hold(True)
    #ax.set_aspect('equal')

    #    ax.scatter(x,y,c=ns)
    ax.hist(ns,bins=range(0,10))
    ax.set_title(sname + " temp: " + str(temp))
#    plt.plot(x,y,'ro')


    
    if not fnameg == None:
        f_path = '/home/tcaswell/python/figures/' + sname + '/'
        if not os.path.isdir(f_path):
            os.makedirs(f_path)
        
        fnameg = f_path + str(key) + '_phi6.eps'
        fig.savefig(fnameg)
        
        
    if istatus:
        print "displaying figure"
        plt.ion()
        plt.show()
    else:
        print "closing figure"
        plt.close(fig)


def _extract_phi6_values(phi6_key,conn,fr_num):
    '''Takes in a phi6 computation number, assumed to be a length one tuple
    as comes out of the sql function the sql connection, and a frame number'''
    
    (fname,iden_key,dset_key) = conn.execute("select fout,iden_key,dset_key from phi6" +
                                             " where  comp_key = ?",
                                             phi6_key).fetchone()
    

    f = h5py.File(fname,'r')

    ns = f[ff(fr_num)][fd('neighborhood_size',phi6_key[0])][:]
    
    x = f[ff(fr_num)][fd('x',iden_key)][:]
    # hack to deal with bug in size of output vectors, now fixed in the c++
    if(x.shape != ns.shape):
        if(x.shape > ns.shape and len(x.shape) ==1):
            ns_full = np.append(ns,-1*np.ones(np.array(x.shape) - np.array(ns.shape)))
    else:
        ns_full = ns
    x = x[ns>0]
    y = f[ff(fr_num)][fd('y',iden_key)][:]

    y = y[ns_full>0]
    
    phi = f[ff(fr_num)][fd('scaler_order_parameter',phi6_key[0])]

    phi = phi[ns>0]
    f.close()
    del f


    return x,y,phi,ns[ns>0]



def testing(phi6_key,conn,fr_num):
    '''Takes in a phi6 computation number, assumed to be a length one tuple
    as comes out of the sql function the sql connection, and a frame number'''
    
    (fname,iden_key,dset_key) = conn.execute("select fout,iden_key,dset_key from phi6" +
                                             " where  comp_key = ?",
                                             phi6_key).fetchone()
    

    f = h5py.File(fname,'r')
    ns = f[ff(fr_num)][fd('neighborhood_size',phi6_key[0])][:]
    x = f[ff(fr_num)][fd('x',iden_key)][:]
    
    f.close()
    del f

    return (x.shape,ns.shape)
def mean_phi(phi6_key,conn,fr_num):
    ''' '''
    print phi6_key
    x,y,phi,ns = _extract_phi6_values(phi6_key,conn,fr_num)
    (pr,pi) = zip(*phi)

    bar_p = np.array([np.sqrt(p[0]**2 + p[1]**2) for p in phi])

    
    
    return np.mean(bar_p), (np.mean(pr),np.mean(pi))
    
def plot_file_frame_pos(phi6_key,conn, fr_num,img = None):
    
    (sname,stype) = conn.execute("select sname,dtype from dsets where dset_key in (select dset_key from phi6 where comp_key = ?)",
                                      phi6_key).fetchone()
    

    x,y,phi,ns = _extract_phi6_values(phi6_key,conn,fr_num)
    
    (pr,pi) = zip(*phi)
    
    istatus = plt.isinteractive();
    print istatus
    if istatus:plt.ioff()

    
    fig = plt.figure()
        
    ax = fig.add_axes([.1,.1,.8,.8])
    ax.set_aspect('equal')
    bar_p = np.array([np.sqrt(p[0]**2 + p[1]**2) for p in phi])
    if img is not None:
        ax.imshow(np.flipud(img),interpolation='nearest',cmap=cm.gray)
    print np.max(bar_p)
    print np.min(bar_p)
    print np.mean(bar_p)
    #sc = ax.scatter(x,y,c=bar_p)
    ax.quiver(x,y,pr,pi,ns)
    #plt.colorbar(sc)
    
    if istatus:
        print "displaying figure"
        plt.ion()
        plt.show()
    else:
        print "closing figure"
        plt.close(fig)

    
def make_phi6_movie(phi6_key,conn,path,x_range,y_range,frame_count ):
    
    (sname,img_fname,ftype) = conn.execute("select sname,fname,ftype from dsets where dset_key in " +
                                           "(select dset_key from phi6 where comp_key = ?)",
                                           phi6_key).fetchone()
    (iden_fname ) = conn.execute("select fout from iden where comp_key in "+
                                  "(select iden_key from phi6 where comp_key = ?)",
                                  phi6_key).fetchone()

    range_str = str(x_range[0]) + '-' + str(x_range[1]) \
                + '_' + str(y_range[0]) + '-' + str(y_range[1])


    if ftype == 1:
        im_wrap = li.Stack_wrapper(img_fname)
        pass
    elif ftype ==2:
        im_wrap = li.Series_wrapper(img_fname)
        pass
    
    if not os.path.exists(path):
        os.makedirs(path,0751)



    istatus = plt.isinteractive();
    print istatus
    if istatus:plt.ioff()
    
    

    for fr in range(0,frame_count):
        
        x,y,phi,ns = _extract_phi6_values(phi6_key,conn,fr)
        (pr,pi) = zip(*phi)
        img = im_wrap.get_frame(fr)

        fig = plt.figure()
        
        ax = fig.add_axes([.1,.1,.8,.8])
        ax.set_aspect('equal')
        ax.axis('off')
        ax.set_xlim(x_range)
        ax.set_ylim(y_range)
        bar_p = np.array([np.sqrt(p[0]**2 + p[1]**2) for p in phi])
        
        ax.imshow(np.flipud(img),interpolation='nearest',cmap=cm.gray)
    
    
        ax.quiver(x,y,pr,pi,ns,clim=[0,9])
        
        plt.savefig(path + '/' + sname + '_' + range_str + '_%05d.png'%fr ,format='png')
        
        plt.close(fig)

    if istatus:
    
        plt.ion()
    

        
        
def make_2d_gofr_plot(comp_key,conn,fname = None):
    # add error handling on all of these calls
    
    # get comp_number of gofr
    (key,g_fname) = conn.execute("select dset_key,fout from comps where comp_key == ? and function == 'gofr'",(comp_key,)).fetchone()
    
        
    # get dset name
    (sname, temp) = conn.execute("select sname,temp from dsets where key == ? ",(key,)).fetchone()

    print sname + " " + str(temp)
    group = gen.get_gofr_group(g_fname,'gofr',comp_key)



    # make plot
    istatus = plt.isinteractive();
    if istatus:plt.ioff()
    
    dset_names = ['bin_count', 'bin_edges']
    fig = plt.figure()
    ax = fig.add_axes([.1,.1,.8,.8])
    ax.hold(True)
    ax.plot(group[dset_names[1]][:],group[dset_names[0]])


    # finds the location of the maximum, assume to be the first peak
    d0 = group[dset_names[1]][np.argmax(group[dset_names[0]])]
    print np.argmax(group[dset_names[0]])
    print d0
    #_draw_gofr_hex_lines(ax,d0/2)
    ax.set_title(sname + " temp: " + str(temp) + ' d0' + str(d0))
    ax.set_ylim(0,3)

    ax.grid(True)

    # save figure
    if not fname == None:
        f_path = '/home/tcaswell/python/figures/' + sname + '/'
        if not os.path.isdir(f_path):
            os.makedirs(f_path)
        
        fname = f_path + str(key) + '_gofr2D.png'
        fig.savefig(fname)


     
        
    if istatus:
        plt.ion()
        plt.show()
    else:
        plt.close(fig)

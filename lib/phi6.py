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
    
    x = f["/frame%(#)06d"%{"#":fr_num}+"/x_%(#)07d"%{"#":comp_number}]
    y = f["/frame%(#)06d"%{"#":fr_num}+"/y_%(#)07d"%{"#":comp_number}]
    phi = f["/frame%(#)06d"%{"#":fr_num}+"/scaler_order_parameter_%(#)07d"%{"#":p_comp}]
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
    
    x = f["/frame%(#)06d"%{"#":fr_num}+"/x_%(#)07d"%{"#":comp_number}]
    y = f["/frame%(#)06d"%{"#":fr_num}+"/y_%(#)07d"%{"#":comp_number}]
    ns = f["/frame%(#)06d"%{"#":fr_num}+"/neighborhood_size_%(#)07d"%{"#":p_comp}]

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


def plot_file_frame_pos(key,conn, fr_num):
    (fname,p_comp) = conn.execute("select fout,comp_key from comps where function = 'phi6' and dset_key = ? and date = '2010-04-19';",(key,)).fetchone()
    (comp_number,) = conn.execute("select comp_key from comps where function = 'Iden' and fout = ?",(fname,)).fetchone()

    (sname,stype,temp) = conn.execute("select sname,dtype,temp from dsets where key = ?",(key,)).fetchone()

    
    f = h5py.File(fname,'r')

    x = f["/frame%(#)06d"%{"#":fr_num}+"/x_%(#)07d"%{"#":comp_number}]
    y = f["/frame%(#)06d"%{"#":fr_num}+"/y_%(#)07d"%{"#":comp_number}]
    ns = f["/frame%(#)06d"%{"#":fr_num}+"/neighborhood_size_%(#)07d"%{"#":p_comp}][:]
    istatus = plt.isinteractive();
    print istatus
    if istatus:plt.ioff()

    
    fig = plt.figure()
    ax = fig.add_axes([.1,.1,.8,.8])
    ax.set_aspect('equal')
    sc = ax.scatter(x,y,c=ns)
    plt.colorbar(sc)
        
    if istatus:
        print "displaying figure"
        plt.ion()
        plt.show()
    else:
        print "closing figure"
        plt.close(fig)

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

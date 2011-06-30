import img
import numpy as np

import openpiv.tools
import openpiv.process



def tac_piv(f_base,start,end):
    # load images
    sw = img.Series_wrapper(f_base,'tif')
    frame_a = sw.get_frame(start).astype(np.int32)
    frame_b = sw.get_frame(end).astype(np.int32)

    # copy tutorial
    u, v, sig2noise = openpiv.process.extended_search_area_piv( frame_a, frame_b, window_size=24, overlap=12, dt=0.02, search_area_size=64, sig2noise_method='peak2peak' )

    x, y = openpiv.process.get_coordinates( image_size=frame_a.shape, window_size=24, overlap=12 )

    u, v, mask = openpiv.validation.sig2noise_val( u, v, sig2noise, threshold = 1.3 )

    u, v = openpiv.filters.replace_outliers( u, v, method='localmean', max_iter=10, kernel_size=2)

    return x,y,u,v,sig2noise

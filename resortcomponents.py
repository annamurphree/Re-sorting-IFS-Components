# -*- coding: utf-8 -*-
"""
Created on Fri Mar  6 16:52:52 2020

@author: Anna Murphree

This is an algorithm that will re-sort two components of IFU data based on velocity arrays.
Specifically, this code sorts the sides of the arrays: it places spaxels with greater velocities 
in the left side of component A and the right side of component B, and spaxels with lower velocities 
in the right side of component A and the left side of component B. 

Inputs: c1, c2, & velocity arrays for c1 and c2 (the original components)
Outputs: cA & cB (new sorted components A and B)

"""

def re_sort(c1, c2, vel1, vel2):
    import numpy as np
    from astropy.io import fits
    import matplotlib.pyplot as plt
    import matplotlib
    cmap = matplotlib.cm.seismic
    #cmap = matplotlib.cm.afmhot
    cmap.set_bad('black', 1.)
    
    hdu1 = fits.open(c1)
    data1 = hdu1[0].data
    hdu2 = fits.open(c2)
    data2 = hdu2[0].data
    hdu3 = fits.open(vel1)
    v1 = hdu3[0].data
    hdu4 = fits.open(vel2)
    v2 = hdu4[0].data
    
    # this will plot the original c1 and c2 arrays
    plt.imshow(data1, cmap=cmap, origin='lower')
    plt.colorbar(shrink=0.7)
    plt.show()
    plt.imshow(data2, cmap=cmap, origin='lower')
    plt.colorbar(shrink=0.7)
    
    #data1[data1 > 1e98] = np.nan
    #data2[data2 > 1e98] = np.nan

    # to make the velocity masks:
    # create holders the size of the whole array
    big1 = np.empty((22, 33))
    big1[:] = np.nan
    big2 = np.empty((22, 33))
    big2[:] = np.nan
    
    # this sorts the velocity arrays
    for i in range(big1.shape[0]):
        for j in range(big1.shape[1]):
            # if one side is nan and the other is not, make them both nan
            if data1[i, j] == np.nan and data2[i, j] != np.nan:
                data2[i, j] == np.nan
            elif data2[i, j] == np.nan and data1[i, j] != np.nan:
                data1[i, j] == np.nan
            
            # if neither of them are nan, we can sort them:
            elif v1[i,j] != np.nan and v2[i,j] != np.nan:
                #if v1 is bigger, put it into big1
                if v1[i,j] > v2[i,j]:
                    big1[i,j] = v1[i,j]
                #if v2 is bigger, put it into big2    
                elif v2[i,j] > v1[i,j]:
                    big2[i,j] = v2[i,j]
                #else: 
                #    big1[i,j] = data1[i,j]
                #    big2[i,j] = data1[i,j]
    
    # make them into masks
    big1[np.isnan(big1)] = 0
    big2[np.isnan(big2)] = 0
    big1[big1 != 0] = 1
    big2[big2 != 0] = 1
    # now, big1 is a mask of all the spaxels where v1 is higher
    # and big2 is a mask of all the spaxels where v2 is higher

    # apply the velocity masks:
    # define the left, middle, and right sides of the flux arrays
    f1left = data1[:, 0:10]
    f1mid = data1[:, 10:24] #the narrow component middle
    f1right = data1[:, 24:33]
    f2left = data2[:, 0:10]
    f2mid = data2[:, 10:24] #the broad component middle
    f2right = data2[:, 24:33]
    
    # define the left, middle, and right sides of the velocity masks
    b1left = big1[:, 0:10]
    #b1mid = big1[:, 10:24]
    b1right = big1[:, 24:33]
    b2left = big2[:, 0:10]
    #b2mid = big2[:, 10:24]
    b2right = big2[:, 24:33]
    
    ''
    # define the left, middle, and right sides of the velocity arrays
    v1left = v1[:, 0:10]
    v1mid = v1[:, 10:24]
    v1right = v1[:, 24:33]
    v2left = v2[:, 0:10]
    v2mid = v2[:, 10:24]
    v2right = v2[:, 24:33]
    ''
    # multiply each part of the flux arrays by the velocity masks
    # component A:
    A1left = np.multiply(f1left, b1left)
    A2left = np.multiply(f2left, b2left)
    A1right = np.multiply(f1right, b2right)
    A2right = np.multiply(f2right, b1right)
    # add the parts together into the final component arrays
    Aleft = np.array(np.add(A1left, A2left))
    Aright = np.array(np.add(A1right, A2right))
    A = np.concatenate((Aleft, f1mid, Aright), axis=1)
    # component B:
    B1left = np.multiply(f1left, b2left)
    B2left = np.multiply(f2left, b1left)
    B1right = np.multiply(f1right, b1right)
    B2right = np.multiply(f2right, b2right)
    # add the parts together into the final component arrays
    Bleft = np.add(B1left, B2left)
    Bright = np.add(B1right, B2right)
    B = np.concatenate((Bleft, f2mid, Bright),axis=1)
    
    '' #this will sort the velocities into A and B arrays too, if you want that
    #velocities:
    vA1left = np.multiply(v1left, b1left)
    vA2left = np.multiply(v2left, b2left)
    vA1right = np.multiply(v1right, b2right)
    vA2right = np.multiply(v2right, b1right)
    # add the parts together into the final component arrays
    vAleft = np.array(np.add(vA1left, vA2left))
    vAright = np.array(np.add(vA1right, vA2right))
    vA = np.concatenate((vAleft, v1mid, vAright), axis=1)
    # B:
    vB1left = np.multiply(v1left, b2left)
    vB2left = np.multiply(v2left, b1left)
    vB1right = np.multiply(v1right, b1right)
    vB2right = np.multiply(v2right, b2right)
    # add the parts together into the final component arrays
    vBleft = np.array(np.add(vB1left, vB2left))
    vBright = np.array(np.add(vB1right, vB2right))
    vB = np.concatenate((vBleft, v2mid, vBright), axis=1)
    ''
    
    A[A == 0] = np.nan
    B[B == 0] = np.nan
    vA[A == 0] = np.nan
    vB[B == 0] = np.nan    
    
    # make sure the arrays are the same shape you started with
    np.reshape(A, (22, 33))
    np.reshape(B, (22, 33))
    np.reshape(vA, (22, 33))
    np.reshape(vB, (22, 33))

    ''
    # this shows the sorted A and B arrays
    plt.imshow(A, cmap=cmap, origin='lower')
    plt.colorbar(shrink=0.7)
    plt.show()
    plt.imshow(B, cmap=cmap, origin='lower')
    plt.colorbar(shrink=0.7)
    '' # write the new arrays to fits files
    fits.writeto('comp_A.fits', A, overwrite = True)
    fits.writeto('comp_B.fits', B, overwrite = True)
    fits.writeto('v_A.fits', vA, overwrite = True)
    fits.writeto('v_B.fits', vB, overwrite = True)

re_sort('Halphafc1_err.fits', 'Halphafc2_err.fits', 'v50c1err.fits', 'v50c2err.fits')
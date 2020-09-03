"""
Created on Thu Sep  3 16:28:55 2020

@author: Anna Murphree, annammurphree@gmail.com

This code re-sorts two flux components of IFU data based on velocity arrays.
Specifically, it sorts the sides of the arrays: it places spaxels with greater velocities 
in the left side of component A and the right side of component B, and spaxels with lower velocities 
in the right side of component A and the left side of component B. 

Inputs: flux components c1 & c2, velocity arrays v1 & c2, emission line (for file names)
Outputs: re-sorted components A and B
"""

# inputs: original flux components c1 & c2, original velocity components v1 & v2, emission line (for file name)
# outputs: sorted flux componented A & B
def re_sort(c1, c2, vel1, vel2, eline, save='no'):
    import numpy as np
    from astropy.io import fits
    import matplotlib.pyplot as plt
    import matplotlib
    
    # open the data files
    hdu1 = fits.open(c1)
    data1 = hdu1[0].data
    hdu2 = fits.open(c2)
    data2 = hdu2[0].data
    hdu3 = fits.open(vel1)
    v1 = hdu3[0].data
    hdu4 = fits.open(vel2)
    v2 = hdu4[0].data
    
    # call all the bad data nan (bad = 1e98, but anything over 50 is bad too)
    data1[data1 > 50] = np.nan
    data2[data2 > 50] = np.nan

    # to make the velocity masks:
    # create holders the size of the whole array
    big1 = np.empty((22, 33))
    big1[:] = np.nan
    big2 = np.empty((22, 33))
    big2[:] = np.nan
    # loop through each spaxel in the velocity arrays:
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
    
    #make them into masks:
    big1[np.isnan(big1)] = 0
    big2[np.isnan(big2)] = 0
    big1[big1 != 0] = 1
    big2[big2 != 0] = 1
    # now, big1 is a mask of all the spaxels where v1 is higher
    # and big2 is a mask of all the spaxels where v2 is higher
    
    # apply the velocity masks:
    # define the left, middle, and right sides of the flux arrays
    f1left = data1[:, 0:10]
    f1mid = data1[:, 10:24]   #the narrow component middle
    f1right = data1[:, 24:33]
    f2left = data2[:, 0:10]
    f2mid = data2[:, 10:24]   #the broad component middle
    f2right = data2[:, 24:33]
    
    # define the left and right sides of the velocity masks (we'll leave the middle unsorted)
    b1left = big1[:, 0:10]
    #b1mid = big1[:, 10:24]
    b1right = big1[:, 24:33]
    b2left = big2[:, 0:10]
    #b2mid = big2[:, 10:24]
    b2right = big2[:, 24:33]
    
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
    
    # make sure the final arrays are the right shape:
    np.reshape(A, (22, 33))
    np.reshape(B, (22, 33))

    # plot the two original and new components:
    fig, axs = plt.subplots(2,2, sharey=True, sharex=True)
    plt.suptitle(f'{eline}', fontsize=14)
    cmap = matplotlib.cm.rainbow
    cmap.set_bad('black', 1.)
    img = axs[0,0].imshow(data1, cmap=cmap, origin='lower')
    axs[0,0].set_title('C1')
    axs[0,1].imshow(data2, cmap=cmap, origin='lower')
    axs[0,1].set_title('C2')
    axs[1,0].imshow(A, cmap=cmap, origin='lower')
    axs[1,0].set_title('CA')
    axs[1,1].imshow(B, cmap=cmap, origin='lower')
    axs[1,1].set_title('CB')
    cbaxes = fig.add_axes([0.95, 0.1, 0.03, 0.8])
    plt.colorbar(img, cax = cbaxes)
    
    if save == 'yes':
    # write the new components to new fits files:
        fits.writeto(f'{eline}_A.fits', A, overwrite = True)
        fits.writeto(f'{eline}_B.fits', B, overwrite = True)

re_sort('OI6300fc1.fits', 'OI6300fc2.fits', 'v50c1.fits', 'v50c2.fits', 'OI6300')
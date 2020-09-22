"""
Created on Mon Sep 14 14:58:39 2020

@author: annamurphree

This code opens a datacube of IFS data, extracts emission line arrays for a 
given header, and saves them into .fits files. You can specify which header,
and then if you want a specific emission line's files. 
"""

def get_that_data(cube):
    from astropy.io import fits
    import idlsave
    import numpy as np
    
    data = idlsave.read(f'{cube}')
    print('Data headers: ', data.keys())
    header = input("Which header? ")
    head = data[f'{header}']        # a recarray
    tabledata = head.table_data     # a ndarray
    print('Components: ', tabledata[0].pkey)
    
    # the different headers are organized differently:
    if header == 'emlwav' or header == 'emlwaverr' or header == 'emlsig' or header == 'emlsigerr':
        comp1 = tabledata[0][2][1][0][8]
        comp2 = tabledata[0][1][1][0][8]
        print('Emission lines: ', comp1.pkey)
        eline = input("Which line? ")
        comps = {'c2':comp2, 'c1':comp1}
        
    elif header == 'emlflx' or header == 'emlflxerr':
        ftot = tabledata[0][0][1][0][8]
        comp2 = tabledata[0][1][1][0][8]
        comp2pk = tabledata[0][4][1][0][8]
        comp1 = tabledata[0][5][1][0][8]
        comp1pk = tabledata[0][6][1][0][8]
        print('Emission lines: ', ftot.pkey)
        eline = input("Which line? ")
        comps = {'ftot':ftot, 'c2':comp2, 'c2pk':comp2pk, 'c1':comp1, 'c1pk':comp1pk}
        
    elif header == 'emlcvdf':
        fluxerr = tabledata[0][0][1][0][8]
        vel = tabledata[0][1]
        cumfluxnorm = tabledata[0][2][1][0][8]
        cumfluxnormerr = tabledata[0][3][1][0][8]
        flux = tabledata[0][5][1][0][8]
        print('Emission lines: ', flux.pkey)
        eline = input("Which line? ")
        comps = {'fluxerr':fluxerr, 'flux':flux, 'vel':vel, 'cfn':cumfluxnorm, 'cfnerr':cumfluxnormerr}
    else: 
        print('Header does not exist')
    
    # loop through each component in the header:
    for c in comps: 
        for key in range(len(comps[c].pkey)):
            line = comps[c].pkey[key]
            if line != None:
                line = line.decode('utf-8')
                if line == eline:
                    print('Line: ', line)
                    data = comps[c].pvalue[key]
                    data[data > 1e90] = np.nan   # call all bad data nan
                    fits.writeto(f'{header}_{c}_{line}.fits', data)
                elif eline == 'all':
                    print('Line: ', line)
                    data = comps[c].pvalue[key]
                    data[data > 1e90] = np.nan   # call all bad data nan
                    fits.writeto(f'{header}_{c}_{line}.fits', data)            
    
get_that_data('3zw35.lin.xdr')
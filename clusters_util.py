import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import ascii
from astropy.table import Table
import os
DTD_path = os.getenv('DTD') + '/'
cluster_path = DTD_path + 'LMC_Star_Clusters/'


def write_cluster_rrl_table(cluster_name, ra, dec, log_age, log_M, rad, n_RRL, n_back, filename = ''):
    """
    Write Cluster information and RRL counts per cluster
    """

    rrlTable = Table()
    rrlTable['Name'] = cluster_name
    rrlTable['RA'] = np.around(ra, decimals=2)
    rrlTable['Dec'] = np.around(dec, decimals=2)
    rrlTable['LogAge'] = log_age
    rrlTable['logM'] = log_M
    rrlTable['Rad'] = np.around(rad, decimals=2)
    rrlTable['N_RRL'] = n_RRL
    rrlTable['N_back'] = np.around(np.mean(n_back, axis=1), decimals=2)
    rrlTable['s_back'] = np.around(np.sqrt(rrlTable['N_back']), decimals=2)
    rrlTable['Det_Sig'] = np.around(np.absolute((rrlTable['N_RRL'] - rrlTable['N_back']))/rrlTable['s_back'],\
                                    decimals=2)
    print 'Clusters with a significant RRL population\n\n'
    print 'Names - ', rrlTable['Name'][rrlTable['Det_Sig']>2.0]
    print 'Number of RRLs - ', rrlTable['N_RRL'][rrlTable['Det_Sig']>2.0]
    print 'Significance (Nsigma) - ', rrlTable['Det_Sig'][rrlTable['Det_Sig']>2.0]

    rrlTable.write(filename, overwrite=True, format='ascii.tab')

def size_from_Bica_catalog(objName, unit='deg'):
    """
    Returns location and information for the requested cluster
    
    Parameters:
        objName: Name of the object (str)

        unit: Unit of size ['arcsec'|'arcmin'|'deg']
    Returns:
        size of the object (in arcmin)
    """
    
    conv_factor = {'arcsec':60., 'arcmin':1., 'deg':1./60.}

    bima_Cl = ascii.read(cluster_path + 'LMC_Bica_Catalog.tsv', data_start=3)
    cl_names = [cl_name.replace(' ', '').split(',') for cl_name in bima_Cl['Names']]
    wh_inbima = [i for i,name in enumerate(cl_names) if objName in name]
    if len(wh_inbima)==0:
        return 0
    else:
        return bima_Cl['amaj'][wh_inbima[0]]*conv_factor[unit]

def cluster_ra_dec(ra_str, dec_str):
    """
    Converts RA, DEC of clusters from string format (vv:vv:vv) to degrees

    Parameters:
        ra_str - RA (in hh:mm:ss)

        dec_str - DEC (in dd:mm:ss)

    Returns:
        RA and DEC in degrees
    """

    coords = [ra + ' ' + dec for (ra, dec) in zip(ra_str, dec_str)]
    coords_deg = SkyCoord(coords, unit=('hourangle', 'deg'), frame='icrs')
    return coords_deg.ra.deg, coords_deg.dec.deg

def pc_to_pixel_radius(rad, unit='arcsec'):
    """
    Convert the radius 'rad' in pc to radius in pixels

    Parameters:
       rad - radius (in parsecs)
       
       unit - units of the returned radius 
              ['arcsec'|'arcmin'|'deg']

    Returns:
       Radius in requested unit.

    """
    rad_dict = {'arcsec':rad/0.24, 'arcmin':rad/(0.24*60), 'deg':rad/(0.24*60*60)}
    return rad_dict[unit]

def num_rrl_within_box(ra_cen, dec_cen, size, rrl_ra, rrl_dec):
    """
    Returns the number of RRLs per box defined by the coordinates 
    and size.

    Parameters:
        ra_cen - RA of box center (in degrees)
        
        dec_cen - DEC of box center (in degrees)
        
        size - size of box (in degrees)

        rrl_ra - RAs of RRLyrae (in degrees)

        rrl_dec - DECs of RRLyrae (in degrees)

    Returns:
        Number of RRLs in box
    """
    rad_to_deg = 3.14/180.
#    ras_in_box = (rrl_ra <= ra_cen + (size/np.sin(-dec_cen*rad_to_deg))) & (rrl_ra >= ra_cen - (size/np.sin(-dec_cen*rad_to_deg)))
    ras_in_box = (rrl_ra <= ra_cen + size) & (rrl_ra >= ra_cen - size)
    rrl_ras_temp = rrl_ra[ras_in_box]
    rrl_dec_temp = rrl_dec[ras_in_box]
    
    decs_in_box = (rrl_dec_temp <= dec_cen + size) & (rrl_dec_temp >= dec_cen - size)
    return rrl_dec_temp[decs_in_box].size

def num_rrl_around_box(ra_cen, dec_cen, size, rrl_ra, rrl_dec):
    """
    Returns the number of RRLs per box in the surrounding area of 
    specified box

    Parameters:
        ra_cen - RA of box center (in degrees)
        
        dec_cen - DEC of box center (in degrees)
        
        size - size of box (in degrees)

        rrl_ra - RAs of RRLyrae (in degrees)

        rrl_dec - DECs of RRLyrae (in degrees)

    Returns:
        Number of RRLs in surrounding boxes
    """    
    return np.array([num_rrl_within_box(ra_cen + (2.0*size), dec_cen, size, rrl_ra, rrl_dec),\
            num_rrl_within_box(ra_cen - (2.0*size), dec_cen, size, rrl_ra, rrl_dec),\
            num_rrl_within_box(ra_cen, dec_cen + (2.0*size), size, rrl_ra, rrl_dec),\
            num_rrl_within_box(ra_cen, dec_cen - (2.0*size), size, rrl_ra, rrl_dec),\
            num_rrl_within_box(ra_cen + (2.0*size), dec_cen + (2.0*size), size, rrl_ra, rrl_dec),\
            num_rrl_within_box(ra_cen + (2.0*size), dec_cen - (2.0*size), size, rrl_ra, rrl_dec),\
            num_rrl_within_box(ra_cen - (2.0*size), dec_cen + (2.0*size), size, rrl_ra, rrl_dec),\
            num_rrl_within_box(ra_cen - (2.0*size), dec_cen - (2.0*size), size, rrl_ra, rrl_dec)])

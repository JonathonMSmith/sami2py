""" Utilities for archiving data in a netCDF file format
    
    loads the .dat files from the trchive and reorganizes them into a netCDF
    including meta data
"""
import os
import xarray as xr
import numpy as np
from sami2py import __version__
from .utils import get_unformatted_data


def _archive_model(path, info, clean, test):
    """Moves the model output files to a common archive

    Parameters
    ----------
    path : (string)
        full path of file destination
    clean : (boolean)
        If True, then delete dat files locally
    fejer : (boolean)
        Specifies whether Fejer-Scherliess model is used
        If False, then 'exb.inp' is also archived
    """
    import shutil
    import subprocess

    if info['fmtout']:
        filelist = ['glonf.dat', 'glatf.dat', 'zaltf.dat',
                    'denif.dat', 'vsif.dat', 'tif.dat', 'tef.dat',
                    'time.dat', 'sami2py-1.00.namelist']
        if info['outn']:
            filelist.append('dennf.dat')
            filelist.append('u4f.dat')
    else:
        filelist = ['glonu.dat', 'glatu.dat', 'zaltu.dat',
                    'deniu.dat', 'vsiu.dat', 'tiu.dat', 'teu.dat',
                    'time.dat', 'sami2py-1.00.namelist']
        if info['outn']:
            filelist.append('dennu.dat')
            filelist.append('u4u.dat')

    if os.path.isfile(filelist[0]):
        try:
            os.stat(path)
        except FileNotFoundError:
            os.makedirs(path)

        hash = subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD'])
        with open(os.path.join(path, 'version.txt'), 'w+') as f:
            f.write('sami2py v' + __version__ + '\n')
            f.write('short hash ' + hash.decode("utf-8"))

        # load dat files into xarray
        if test:
            model = xr.Dataset()
        else:
            model = _load_model(info)

        if not info['fejer']:
            model.attrs['ExB model'] = 'Fouerier Series'
            model.attrs['Fouerier Coeffs'] = np.loadtxt('exb.inp')

        if clean:
            """Delete dat files if clean is True."""
            for list_file in filelist[:-1]:
                os.remove(list_file)
        # export model with metadata to the archive directory as a netCDF
        if not test:
            filename = os.path.join(path, 'sami.nc')
            model.to_netcdf(filename)
    else:
        print('No files to move!')



def _calculate_slt(ut, lon0, day):
    """Calculates Solar Local Time for reference point of model

    Returns
    -------
    self.slt : (float)
        Solar Local Time in hours

    """

    local_time = np.mod((ut * 60 + lon0 * 4), 1440) / 60.0
    mean_anomaly = 2 * np.pi * day / 365.242
    delta_t = (-7.657 * np.sin(mean_anomaly) +
               9.862 * np.sin(2 * mean_anomaly + 3.599))
    slt = local_time - delta_t / 60.0
    return slt

def _load_model(info):
    """Loads model results
    """

    nf = 98
    nz = 101
    ni = 7
    lon0 = info['lon']
    day = info['day']

    # Get times
    time = np.loadtxt('time.dat')
    ut = time[:, 1] + time[:, 2] / 60 + time[:, 3] / 3600

    slt = _calculate_slt(ut, lon0, day)
    nt = len(ut)

    if info['fmtout']:
        # Get Location
        glat = np.loadtxt('glatf.dat')
        glon = np.loadtxt('glonf.dat')
        zalt = np.loadtxt('zaltf.dat')

        # Get plasma values
        deni = np.loadtxt('denif.dat')
        vsi = np.loadtxt('vsif.dat')
        ti = np.loadtxt('tif.dat')
        te = np.loadtxt('tef.dat')

        # get neutral values
        if info['outn']:
            denn = np.loadtxt('dennf.dat')
            u4 = np.loadtxt('u4f.dat')
    else:
        # Get Location
        glat = get_unformatted_data('', 'glat')
        glon = get_unformatted_data('', 'glon')
        zalt = get_unformatted_data('', 'zalt')

        # Get plasma values
        dim0 = nz*nf*ni + 2
        dim1 = nt
        deni = get_unformatted_data('', 'deni',
                                    dim0=dim0, dim1=dim1, reshape=True)
        vsi = get_unformatted_data('', 'vsi',
                                   dim0=dim0, dim1=dim1, reshape=True)
        ti = get_unformatted_data('', 'ti',
                                  dim0=dim0, dim1=dim1, reshape=True)
        if info['outn']:
            denn = get_unformatted_data('', 'denn',
                                        dim0=dim0, dim1=dim1, reshape=True)

        # Electron Temperatures and neutral wind have only one species
        dim0 = nz*nf + 2
        te = get_unformatted_data('', 'te',
                                  dim0=dim0, dim1=dim1, reshape=True)
        if info['outn']:
            u4 = get_unformatted_data('', 'u4',
                                      dim0=dim0, dim1=dim1, reshape=True)

    glat = np.reshape(glat, (nz, nf), order="F")
    glon = np.reshape(glon, (nz, nf), order="F")
    zalt = np.reshape(zalt, (nz, nf), order="F")
    deni = np.reshape(deni, (nz, nf, ni, nt), order="F")
    vsi = np.reshape(vsi, (nz, nf, ni, nt), order="F")
    ti = np.reshape(ti, (nz, nf, ni, nt), order="F")
    te = np.reshape(te, (nz, nf, nt), order="F")
    data = xr.Dataset({'deni': (['z', 'f', 'ion', 'ut'], deni),
                       'vsi': (['z', 'f', 'ion', 'ut'], vsi),
                       'ti': (['z', 'f', 'ion', 'ut'], ti),
                       'te': (['z', 'f', 'ut'], te),
                       'slt': (['ut'], slt)},
                      coords={'glat': (['z', 'f'], glat),
                              'glon': (['z', 'f'], glon),
                              'zalt': (['z', 'f'], zalt),
                              'ut': ut})
    data.attrs = info
    if info['outn']:
        print(info['outn'])
        denn = np.reshape(denn, (nz, nf, ni, nt), order="F")
        data['denn'] = (('z', 'f', 'ion', 'ut'), denn)
        u4 = np.reshape(u4, (nz, nf, nt), order="F")
        data['u4'] = (('z', 'f', 'ut'), u4)
    return data

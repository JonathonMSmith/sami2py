import os
import numpy as np
import xarray
from netCFD4 import Dataset


def _read_fortran_data(path, clean, fejer, fmtout, info):
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
    if fmtout:
        filelist = ['glonf.dat', 'glatf.dat', 'zaltf.dat',
                    'denif.dat', 'vsif.dat', 'tif.dat', 'tef.dat',
                    'time.dat', 'sami2py-1.00.namelist']
    else:
        filelist = ['glonu.dat', 'glatu.dat', 'zaltu.dat',
                    'deniu.dat', 'vsiu.dat', 'tiu.dat', 'teu.dat',
                    'time.dat', 'sami2py-1.00.namelist']
    # set up the netcdf dataset dimensions
    time = np.loadtxt('time.dat')
    ut = time[:, 1] + time[:, 2] / 60 + time[:, 3] / 3600
    dataset = Dataset('model_data.nc', 'w', format='NETCDF4_CLASSIC')
    nf = dataset.createDimension('nf', 98)
    nz = dataset.createDimension('nz', 101)
    nt = dataset.createDimension('nt', len(ut))
    ni = dataset.createDimension('ni', 7)
    # set global attributes
    dataset.setncatts(info)
    # grab all of the data output from sami2 fortran and store in the netcdf
    # dataset with the proper dimensions
    tmp_var = dataset.createVariable('time', np.float, nt)
    tmp_var[:] = ut
    if os.path.isfile(filelist[0]):
        try:
            os.stat(path)
        except FileNotFoundError:
            os.makedirs(path)
        # loop through model output
        for list_file in filelist:
            if list_file == 'time.dat' or list_file == 'sami2py-1.00.namelist':
                continue
            # if it's formatted data then simply load it with numpy otherwise
            # open the binary dat file and then use numpy to format it
            if fmtout:
                f_dat = np.loadtxt(list_file)
            else: 
                f = open(list_file, 'rb')
                f_dat = np.fromfile(f, dtype='float')[1:-1]
            # depending on which model output it is determine what the 
            # appropriate dimensions are for storing
            if np.max(f_dat.size()) == len(nz)*len(nf):
                dim = ('nz', 'nf')
            elif np.max(f_dat.size()) == len(nz)*len(nf)*len(nt):
                dim = ('nz', 'nf', 'nt')
            else:
                dim = ('nz', 'nf', 'ni', 'nt')
            # store the model variable in the dataset
            tmp_var = dataset.createVariable(list_file[:-5], np.float, dim)
            tmp_var[:] = f_dat

        if clean:
            for list_file in filelist[:-1]:
                os.remove(list_file)
        if not fejer:
            exb_inp = np.loadtxt('exb.inp')
            func = dataset.createDimension('func', 2)
            term = dataset.createDimension('term', 10)
            exb = dataset.createVariable('exb', np.float, ('term', 'func'))
            exb[:] = exb_inp
    else:
        raise FileNotFoundError('No model output found to archive')

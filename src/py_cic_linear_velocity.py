

import os, sys
import numpy as np
from numpy.lib.index_tricks import nd_grid
import re

sys.path.append("/home/xiaodongli/software/pythonlib")
import stdA_py3
from struct import unpack
from scipy import fftpack

a = 1

f = 0.519609

h = 0.7
H = stdA_py3.Hz(0.3071, -1, h, 1/a-1) / h

printstr ='''Usage: py_cic_linear_velocity inputpath
    -boxsize boxsize

Example:

Bash example:
    '''

if len(sys.argv) <=1:
    print(printstr)
    sys.exit()

inputpath = sys.argv[1]
print('\t set inputpath as: ', inputpath)

boxsize = ngrid  = None

if len(sys.argv) > 2:
    for iarg in range(2,len(sys.argv),2):
        str1, str2 = sys.argv[iarg:iarg+2]
        if str1 in ['-boxsize', '-BOXSIZE']:
            boxsize = int(str2)
            print('\t set boxsize as: ', boxsize)
        else:
            print('Unknown option: ', str1)
            print(printstr)
            sys.exit()

def rho_to_v(rho,boxsize,a_f_H):
    a,H,f = a_f_H
    rho_mean = np.mean(rho)
    delta = rho / rho_mean
    delta_k = fftpack.fftn(delta)
    n = np.shape(delta_k)[0]
    k_vector = np.empty((n,n,n,3),dtype=float)
    k_1d = fftpack.fftfreq(n,d=boxsize/(2*np.pi)/n)
    k_vector_source = np.meshgrid(k_1d, k_1d, k_1d)
    for i in range(3):
        k_vector[:,:,:,i] =  k_vector_source[i]

    k_norm = np.clip(np.linalg.norm(k_vector,axis=3),1e-10,None,None)
    v_x = np.empty(k_vector.shape,dtype=float)
    for i in range(3):
        v_x[:,:,:,i] = np.real(fftpack.ifftn(1j*k_vector[:,:,:,i]/k_norm**2 * a*H*f*delta_k))
    return v_x

ngrid = int(re.search('(?<=nc)(.*)(?=_)',inputpath).group())
nsplit = int(re.search('(?<=nsplit)(.*)',inputpath).group())
if nsplit>1:
    # ngrid_part = int(ngrid/nsplit)
    # rhos_grid = np.empty((nsplit,nsplit,nsplit),dtype=object)
    # extra_numbers = np.empty((int(nsplit**3),2),dtype=float)
    # for i in range(nsplit):
    #     for j in range(nsplit):
    #         for k in range(nsplit):
    #             rhofile = "subgrids.ifile{0:d}.rhogrid".format(i*nsplit**2+j*nsplit+k+1)
    #             rhofile = inputpath+"/"+rhofile
    #             print("Now read {0:s}".format(rhofile))
    #             rho_source = np.fromfile(rhofile,dtype='float32')
    #             first_number, last_number = rho_source[[0,-1]]
    #             rhos_grid[i,j,k]=rho_source[1:-1].reshape(ngrid_part,ngrid_part,ngrid_part)
    #             extra_numbers[i,j,k]=(first_number,last_number)
    #     rhos_grid = rho_source[1:-1].reshape(ngrid,ngrid,ngrid)
    pass
elif nsplit==1:
    rhofile = "subgrids.rhogrid"
    rhofile = inputpath+"/"+rhofile
    print("Now read {0:s}".format(rhofile))
    rho_source = np.fromfile(rhofile,dtype='float32')
    first_number, last_number = rho_source[[0,-1]]
    rho = rho_source[1:-1].reshape(ngrid,ngrid,ngrid)
    v_x = rho_to_v(rho,boxsize,(a,f,H))

    file_id = ['vx','vy','vz']
    for i in range(len(file_id)):
        output_file = rhofile.replace(".rhogrid",".{0:s}grid".format(file_id[i]))
        v = v_x[:,:,:,i].ravel()
        v = np.insert(v,0,first_number)
        v = np.insert(v,-1,last_number)
        v.tofile(output_file)
        print("{0:s} save done".format(output_file))
    print()
else:
    raise ValueError("nsplit must be more than 1")



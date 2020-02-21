import numpy as np
import pynbody
import os , sys



#nc = 512; xyzmin = 0; xyzmax = 512;

print(' CIC rho/px/py/pz grid calculation, multi files ')
print(' Usage: py_cic filename nc xyzmin xyzmax outputname')

if len(sys.argv)<5:
    sys.exit()

filename = sys.argv[1]
nc, xyzmin, xyzmax = [float(xx) for xx in sys.argv[2:5]]; nc = int(nc)
dx = (xyzmax-xyzmin)/float(nc); dx_inv = 1./dx

print('#######################################\n', ' * Start cic calculation of rho/px/py/pz fields...')

if True:
    files = os.popen("ls "+filename).read().split()
    nowf = open(files[0]+'.'+str(len(files))+'files.info','w')
    try:
        outputname = sys.argv[5]
    except: 
        outputname = files[0]+'.'+str(len(files))+'files'
        print('outputname not available; use inputfile\'s name: \n\t',outputname)
    nowf.write('nc, xyzmin, xyzmax, outputname = '+' '.join([str(xx) for xx in [nc,xyzmin,xyzmax,outputname,'\n']]))
    print('nc, xyzmin, xyzmax, outputname = '+' '.join([str(xx) for xx in [nc,xyzmin,xyzmax,'\n']]))
    print('Preparing for processing ', len(files), 'files:')
    nowf.write('Will measuring rho/px/py/pz filed for '+str(len(files))+' files:\n')
    for nowfile in files: 
        print('\t',nowfile)
        nowf.write('\t'+nowfile+'\n')
    nowf.close()
    rhogrid, vxgrid, vygrid, vzgrid = np.zeros((nc,nc,nc)).astype('float32'), np.zeros((nc,nc,nc)).astype('float32'), np.zeros((nc,nc,nc)).astype('float32'), np.zeros((nc,nc,nc)).astype('float32')
    for nowfile in files:
      print(' open ', nowfile, '...')
      simdata = pynbody.load(nowfile)
      pos, vel = simdata['pos'], simdata['vel']
      ndat = len(simdata)
      print(' \tIn total ', ndat, 'particles.')
      for i in range(ndat):
        X, V = pos[i,:], vel[i,:]
        ix = (X * dx_inv).astype('int32')
        w = 1. - (X*dx_inv - ix)
        ix0 = (ix+nc)%nc; ix1 = (ix+1+nc)%nc 
        if i%50000 == 0: print('\t\t',i,'-th particle...') 
        w1, w2, w3, w4, w5, w6, w7, w8 = w[0]*w[1]*w[2], w[0]*(1-w[1])*w[2], w[0]*w[1]*(1-w[2]), w[0]*(1-w[1])*(1-w[2]), \
            (1-w[0])*w[1]*w[2], (1-w[0])*(1-w[1])*w[2], (1-w[0])*w[1]*(1-w[2]), (1-w[0])*(1-w[1])*(1-w[2])
        for grid, wei in [
                [rhogrid, 1],
                [vxgrid, V[0]],
                [vygrid, V[1]],
                [vzgrid, V[2]],
                ]:
          grid[ix0[0],  ix0[1], ix0[2]] += w1*wei
          grid[ix0[0],  ix1[1], ix0[2]] += w2*wei
          grid[ix0[0],  ix0[1], ix1[2]] += w3*wei
          grid[ix0[0],  ix1[1], ix1[2]] += w4*wei
          grid[ix1[0],  ix0[1], ix0[2]] += w5*wei
          grid[ix1[0],  ix1[1], ix0[2]] += w6*wei
          grid[ix1[0],  ix0[1], ix1[2]] += w7*wei
          grid[ix1[0],  ix1[1], ix1[2]] += w8*wei

      print(' Finishing processing ', i, 'particles.')
print(' write data...')
for outfile, grid in [
        [outputname+'.rhogrid',rhogrid],
        [outputname+'.pxgrid',vxgrid],
        [outputname+'.pygrid',vygrid],
        [outputname+'.pzgrid',vzgrid],
        ]:
    print('\t write to file: ', outfile)
    grid.reshape(-1).tofile(outfile, format='')


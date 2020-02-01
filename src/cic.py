import numpy as np
import pynbody
import os , sys


simid = int(sys.argv[1])
snpstr = sys.argv[2]

nc = 512; xyzmin = 0; xyzmax = 512;
dx = (xyzmax-xyzmin)/float(nc); dx_inv = 1./dx


print('#######################################\n', ' * simid, snpstr = ', simid, snpstr)

if True:
    #files = os.popen("ls /home/xiaodongli/data/colas/cola_multiverse/BigMD_3000_mocks/"+str(simid)+"/snp0"+str(simid)+snpstr+'*').read().split()
    files = os.popen("ls /home/xiaodongli/data/colas/cola_multiverse/BigMD_3000_mocks/"+str(simid)+"/snp0"+str(simid)+snpstr).read().split()
    nowf = open(str(simid)+'.'+snpstr+'.info','w')
    nowf.write('nc, xyzmin, xyzmax = '+' '.join([str(xx) for xx in [nc,xyzmin,xyzmax,'\n']]))
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
        for grid, wei in [
                [rhogrid, 1],
                [vxgrid, V[0]],
                [vygrid, V[1]],
                [vzgrid, V[2]],
                ]:
          grid[ix0[0],  ix0[1], ix0[2]] += w[0]*w[1]*w[2]*wei; 
          grid[ix0[0],  ix1[1], ix0[2]] += w[0]*(1-w[1])*w[2]*wei;  
          grid[ix0[0],  ix0[1], ix1[2]] += w[0]*w[1]*(1-w[2]) *wei; 
          grid[ix0[0],  ix1[1], ix1[2]] += w[0]*(1-w[1])*(1-w[2]) *wei; 
          grid[ix1[0],  ix0[1], ix0[2]] += (1-w[0])*w[1]*w[2] *wei; 
          grid[ix1[0],  ix1[1], ix0[2]] += (1-w[0])*(1-w[1])*w[2] *wei; 
          grid[ix1[0],  ix0[1], ix1[2]] += (1-w[0])*w[1]*(1-w[2]) *wei; 
          grid[ix1[0],  ix1[1], ix1[2]] += (1-w[0])*(1-w[1])*(1-w[2])*wei; 

      print(' Finishing processing ', i, 'particles.')
print(' write data...')
for outfile, grid in [
        [str(simid)+str(snpstr)+'.rhogrid',rhogrid],
        [str(simid)+str(snpstr)+'.pxgrid',vxgrid],
        [str(simid)+str(snpstr)+'.pygrid',vygrid],
        [str(simid)+str(snpstr)+'.pzgrid',vzgrid],
        ]:
    grid.reshape(-1).tofile(outfile, format='')


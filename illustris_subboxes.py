
import numpy as np
import os
import numpy as np
import struct
make_plot = True
if make_plot:
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
import pandas as pd

def conv_to_log(A):
    return np.sign(A) * np.log(np.abs(A) + 1.)

#siparttype, nsplit, overlaprat = 0, 2, 0.
class illustris_subbox:

    def __init__(self, sim='Illustris-3', snapdir='Snapshot', snap=135, parttype=0,
                 nsplit=2, overlaprat=0.0, path=None):
        self.output_columns = {'PartType0': ['Coordinates', 'Velocities', 'Masses', 'StarFormationRate',
                   'GFM_Metallicity', 'NeutralHydrogenAbundance',
                   'StarFormationRate', 'InternalEnergy'],
         'PartType1': ['Coordinates', 'Velocities'],
         'PartType4': ['Coordinates', 'Velocities', 'Masses'],
         'PartType5': ['Coordinates', 'Velocities', 'BH_Mass', 'BH_Progs',
                       'BH_Mdot', 'Masses', 'HostHaloMass' ]
        }
        self.sim = sim; self.snapdir = snapdir; self.nsplit = nsplit; self.snap=str(snap)
        self.parttype = parttype; self.nsplit = nsplit; self.overlaprat=overlaprat
        self.feature_names = self.output_columns['PartType'+str(parttype)]
        self.cic_feature_names = (' '.join(self.feature_names).replace('Coordinates', '').\
            replace('Velocities', ' vx vy vz ')+' number').split()
        self.n_cic_feature = len(self.cic_feature_names)
        if path == None:
            if self.sim == 'Illustris-3':
                self.path = '/mnt/xiaodong/1/Illustris_xiaodong/Illustris_cic/'+self.sim+'/'+self.snapdir+\
                    '/snap_'+self.snap+'.PartType'+str(self.parttype)
            else:
                self.path = '/mnt/xiaodong/1/Illustris_xiaodong/Illustris_cic/'+self.sim+'/'+self.snapdir+\
                    '/snap_'+self.snap+'.PartType'+str(self.parttype)+\
                    '/snap_'+self.snap+'.PartType'+str(self.parttype)
        else:
            self.path = path
        self.infofile = self.path+'/nsplit'+str(self.nsplit)+'_overlaprat%.3f'%(self.overlaprat)+'.info'
        self.filename_npar_dict = {}
        for nowstr in open(self.infofile, 'r').readlines():
            if nowstr[0] == '#': continue
            filename, npar, nfeature, boxsize = nowstr.split()[:4]; self.boxsize = float(boxsize)
            self.filename_npar_dict[filename] = [int(npar), int(nfeature), float(boxsize)]
        self.overlap_distance = self.overlaprat * self.boxsize / float(self.nsplit)

    def subbox_data_file(self, i1, i2, i3):
        return self.path+'/nsplit'+str(self.nsplit)+'_overlaprat%.3f'%self.overlaprat+\
            '_subbox'+str(i1+1)+'_'+str(i2+1)+'_'+str(i3+1)

    def subbox_data_par(self, i1, i2, i3):
        for nowkey in self.filename_npar_dict:
            if nowkey.split('/')[-1] == self.subbox_data_file(i1,i2,i3).split('/')[-1]:
                return self.filename_npar_dict[nowkey]

    def load_subbox_data(self, i1, i2, i3, datatype='DataFrame'):
        '''load in subbox data as "array" or "DataFrame"'''
        npar, nfeature, boxsize = self.subbox_data_par(i1, i2, i3)
        nowf = open(self.subbox_data_file(i1,i2,i3), 'rb')
        data = np.array(struct.unpack(str(npar*nfeature)+'f', nowf.read(npar*nfeature*4))).reshape(npar, nfeature)
        cols = ' '.join(self.output_columns['PartType'+str(self.parttype)]).replace('Coordinates', 'x y z ').replace('Velocities', 'vx vy vz ').split()
        nowf.close()
        if datatype == 'array':
            return data
        elif datatype == 'DataFrame':
            return pd.DataFrame(data, columns=cols)
        else:
            print('Error (load_subbox_data): datatype should be "array" or "DataFrame"! We get: ', datatype)

    def check_files(self, i1=None, i2=None, i3=None):
        if i1 == None or i2 == None or i3 == None:
            for i1 in range( self.nsplit):
                for i2 in range(self.nsplit):
                    for i3 in range( self.nsplit):
                        subbox_file = self.subbox_data_file(i1,i2,i3)
                        print(' * subbox ',i1,i2,i3,':\n\t',os.popen('ls -alh '+subbox_file).read()[:-1], )
                        print('\tnpar, nfeature, boxsize = ', self.subbox_data_par(i1, i2, i3))
        else:
            subbox_file = self.subbox_data_file(i1,i2,i3)
            print(' * subbox ',i1,i2,i3,':\n\t',os.popen('ls -alh '+subbox_file).read()[:-1], )
            print('\tnpar, nfeature, boxsize = ', self.subbox_data_par(i1, i2, i3))
    def scatter_2d(self, i1_range = [0], i2_range = [0], i3_range = [0], zmin = 0., zmax = 300,
                   fig=None, ax=None, pointsize=0.01):
        if not make_plot:
            print(' scatter_2d exist since make_plot = ', make_plot); return
        if ax == None:
            fig, ax = plt.subplots(1, 1, figsize=(6,6));
        for i1 in i1_range:
            for i2 in i2_range:
                for i3 in i3_range:
                    df = self.load_subbox_data(i1, i2, i3); 
                    rows = np.where((df['z'] < df['z'].min()+zmax)&(df['z'].min()+zmin<df['z']))[0]
                    ax.scatter(df['x'][rows], df['y'][rows],
                            s=pointsize, label='_'.join(['subbox']+[str(xx) for xx in [i1,i2,i3]]))
            #axs[iparttype].legend(loc='upper left');
        ax.set_title(self.sim+', parttype = '+str(self.parttype)+', z<'+str(zmax)); ax.grid()
        return fig, ax

    def scatter_3d(self, i1=0, i2=0, i3=0,
                   fig=None, ax=None, pointsize=0.01):
        if not make_plot:
            print(' scatter_3d exist since make_plot = ', make_plot); return
        if ax == None:
            fig = plt.figure(figsize=(6,6)); ax = fig.add_subplot(111, projection='3d')
        df = self.load_subbox_data(i1, i2, i3)
        #print(df.describe())
        #rows = np.where((df['x']<15000) & (df['y']<15000) & (df['z']<15000))[0]
        ax.scatter(df['x'], df['y'], df['z'], alpha=0.3, s=pointsize, c=df['z'],
                   label='_'.join(['subbox']+[str(xx) for xx in [i1,i2,i3]]))
        #axs[iparttype].legend(loc='upper left');
        ax.set_title(self.sim+', parttype = '+str(self.parttype)); ax.grid()
        return fig, ax

    def do_cic(self, nc, i1s=None,i2s=None,i3s=None, create_bash=True):
        '''Generate bash script for cic filed calculation for i1, i2, i3 in range of i1s, i2s, i3s.
        if i1s,i2s,i3s not given, then compute all fields
        nc: number of cells in each direction
        '''
        if i1s == None: i1s = range(self.nsplit);
        if i2s == None: i2s = range(self.nsplit);
        if i3s == None: i3s = range(self.nsplit);
        dxyz = self.boxsize / float(self.nsplit)
        if create_bash:
            bashfile = self.path+'/nsplit'+str(self.nsplit)+'_overlaprat%.3f'%(self.overlaprat)+'_'+str(nc)+'cic.sh'
            bashf = open(bashfile, 'w')
        for i1 in i1s:
            for i2 in i2s:
                for i3 in i3s:
                    xmin = dxyz * i1 - self.overlap_distance*0.5
                    ymin = dxyz * i2 - self.overlap_distance*0.5
                    zmin = dxyz * i3 - self.overlap_distance*0.5
                    boxsize = dxyz + self.overlap_distance
                    npar, nfeature, totalbd = self.subbox_data_par(i1,i2,i3)
                    cmd = ' '.join([str(xx) for xx in ['LSS_illustris_cic  -input  ',self.subbox_data_file(i1,i2,i3),
                     '-nc ', nc, '-xmin',xmin,'-ymin',ymin,'-zmin',zmin,'-boxsize','%.2f'%boxsize,
                                '-npar',npar,'-nfeature',nfeature,'-output',
                                    self.subbox_data_file(i1,i2,i3)+'_cic', '-nsplit 1']])
                    if create_bash:
                        bashf.write(cmd+'\n')
                    else:
                        print(cmd)
        if create_bash:
            print('commands written to : ', bashfile)
            return bashfile
    def subbox_cic_files(self, i1,i2,i3,nc,do_check=False):
        cic_dict = {}
        for ifeature, feature in enumerate(self.cic_feature_names):
            nowfile = self.subbox_data_file(i1,i2,i3)+'_cic_'+str(nc)+'grid_feature'+str(ifeature+1)
            cic_dict[feature] = nowfile
            if do_check:
                print(feature,':\t', os.popen('ls -alh '+nowfile).read())
        return cic_dict
    def cic_file(self, i1,i2,i3,nc, ifeature):
        return self.subbox_data_file(i1,i2,i3)+'_cic_'+str(nc)+'grid_feature'+str(ifeature+1)
    def load_cic_file(self, i1,i2,i3,nc, ifeature):
        gridfile = self.cic_file(i1,i2,i3,nc,ifeature)
        nowf = open(gridfile, 'rb')
        grid = struct.unpack('1i'+str(nc**3)+'f1i', nowf.read())
        return np.array(grid[1:-1]).reshape(nc,nc,nc)
    def imshow_cic(self, i1,i2,i3,nc, ifeature, logfield=True, iz=0, fig = None, ax = None):
        if not make_plot:
            print(' imshow_cic exit since make_plot = ', make_plot); return
        grid = self.load_cic_file(i1,i2,i3,nc, ifeature)
        if fig == None or ax == None:
            fig, ax = plt.subplots()
        xy = grid[:,::-1,iz].T
        if logfield: xy = conv_to_log(xy)
        ax.imshow(xy); ax.set_title(', '.join([self.sim, 'parttype='+str(self.parttype), '\n', 
            str(self.cic_feature_names[ifeature]), str(nc)+'-grid']))
        return fig, ax


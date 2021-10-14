import numpy as np
import os
import xitools

##### Settings

##### 
#--------------------
# 1. path
mainpath='/home/xiaodongli/projects/cosmology_grid_related/snap_analysis/'
print('Set mainpath as: \n\t', mainpath)

#--------------------
# 2. massbin edges (of bigmd)
class massbin:
    def __init__(self):
        #self.mock = 'bigmd'
        self.mock_redshift = 0.5053
        self.info = '''
        diff_massbin mocks are made from on bigmd snapshot 33 
            (a=0.6643, z=0.5053)
        '''
        print(self.info)
        self.massbin_edges = np.loadtxt('/home/xiaodongli/projects/cosmology_grid_related/snap_analysis/MCF_diff_massbin_bigmd/masscut_scheme2/'+
                              '33_nbar1e-3.shiftz.allcolumns.nbar0.6e-4.mass_ranges.txt')
        self.massbin_edges = dict(zip([x for x in self.massbin_edges[:,0]], self.massbin_edges[:,[2,1]]))
        print(' (massbin init) Read-in massbin_edges as: ', )
        for key in self.massbin_edges:
            print('\t', key, self.massbin_edges[key])
    def printinfo(self):
        print(self.info)
        for key in massbin_edges:
            print('(massbin)\t', key, massbin_edges[key])
print('Set up massbin information of bigmd snap')
massbin = massbin()   
massbin_edges = massbin.massbin_edges

#--------------------
# 3. bigmd snapid-redshift relation

class bigmd:
    def __init__(self):
        data = np.loadtxt('/home/xiaodongli/data/BigMDPL/snaps/Redshifts.csv', delimiter=',')
        self.snapz = dict(zip(data[:,1].astype('int'), data[:,3]))
        self.zsnap = dict(zip(data[:,3], data[:,1].astype('int'), ))
print('Set up bigmd snapz, zsnap information')
bigmd = bigmd()

#--------------------
# 4. SDSS data information
class sdss_info:
    def __init__(self):
        #self.*** = ***
        pass
    ## TBA
    
#---------------------
def number_to_str(x):
    if x == int(x):
        return str(int(x))
    else:
        return str(x)


###cat = mockfile, ranfile, inifile, 2pcffile, rrfile
class file_system_class:
    '''
    Example:
        file_system = file_system_class()

        file_system.sdss('data', ibin= 3, weipow=0.2)
        file_system.diff_cosmo(om=0.3071, w=-1, weipow=0.2) ##  red can be 0.4, 0.5, 0.6; nbar can be 3e-4, 3.5e-4, 4e-4
        file_system.diff_massbin(massbin=0.5)
        
        file_system.print_available()
        file_system.check_file(printinfo=False, return_if_found_missing=True);
    '''
    def __init__(self):
        self.diff_massbin_nbars = ['0.6e-4']
        self.diff_massbin_snapids = [33]
        self.diff_massbin_numNB_margins = [[2,30], [6,30],[20,100] ]              
        
        self.diff_cosmo_omws = [
             [0.2671, -1.0],
             [0.2871, -1.2], [0.2871, -1.0], [0.2871, -0.8], 
             [0.3071, -1.2], [0.3071, -1.0], [0.3071, -0.8], [0.3071, -0.6],
             [0.3271, -1.2], [0.3271, -1.0], [0.3271, -0.8],
             [0.3471, -1.0]
                    ]
        
        #[[om,w] for om in [0.2671, 0.2871, 0.3071, 0.3271, 0.3471] for w in [-0.8, -1, -1.2]]
        self.diff_cosmo_reds = [0.4, 0.5, 0.6]
        self.diff_cosmo_numNB_margins = [[10,30], [30,30],[100,100] ]              
        self.diff_cosmo_nbars = ['3e-4', '3.5e-4', '4e-4']
        self.diff_cosmo_shiftstrs= ['shiftz', 'shiftr']
        
        self.sdss_ibins = [1,2,3]
        self.sdss_sample_nsample = {'data':1, 'J08':4, 'Patchy':10 }
        
        self.weipows = self.diff_cosmo_weipows = self.diff_massbin_weipows = \
            [-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1]
    def print_available(self):
        names = '''weipows, diff_massbin_nbars, diff_massbin_snapids, diff_massbin_numNB_margins, diff_cosmo_omws, diff_cosmo_reds, diff_cosmo_numNB_margins, diff_cosmo_nbars,  sdss_sample_nsample, sdss_ibins'''
        lists = [self.weipows, self.diff_massbin_nbars, self.diff_massbin_snapids, 
            self.diff_massbin_numNB_margins,
            self.diff_cosmo_omws, self.diff_cosmo_reds, self.diff_cosmo_numNB_margins, self.diff_cosmo_nbars,
            self.sdss_sample_nsample, self.sdss_ibins,  ]
        print(' (file_system.print_available) Availables:\n')
        for name, nowlist in zip(names.split(','), lists):
            print('\t', name.strip(' '), ':\n\t\t', nowlist)
            if name.strip(' ') in 'weipows diff_massbin_numNB_margins diff_cosmo_nbars'.split():
                print()
    ### -------------------------------------
    ### file_system.sdss
    ### -------------------------------------
    def sdss(self, sample='data', imock=0, ibin=1, weipow=0.2, numNB=30,  sbin=150,smax=150,mubin=120, 
             sdsscat='DR12v4-CMASS'):
        '''
        self.sdss_ibins = [1,2,3]
        self.sdss_sample_nsample = {'data':1, 'J08':4, 'Patchy':10 }
        available weipows: [-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1,   -3, -2, 2, 3]
        available numNB: 30 
        '''
        path= '/home/xiaodongli/data/boss2pcf/data/'+sdsscat+'/xyzw.binsplitted/datafile/'
        if sample == 'data':
            return path+'data/data.xyzw.'+str(ibin)+'of3.numNB'+str(numNB)+'.wei'+number_to_str(weipow)+\
                '.weipow1.%ds0to%d.%dmu.2pcf'%(sbin,smax,mubin)
        elif sample == 'J08':
            return path+'J08/J08.RSD.%03i'%(imock)+'.xyzw.'+str(ibin)+'of3.numNB'+str(numNB)+'.wei'+number_to_str(weipow)+\
                '.weipow1.%ds0to%d.%dmu.2pcf'%(sbin,smax,mubin)
        elif sample == 'Patchy': # warning: Pathcy is using a different path!!! # BOSSLi
            return '/home/xiaodongli/data/mnt_1/boss2pcf/data/'+sdsscat+'/'+\
                'xyzw.binsplitted/No0PatchyV6C.RSD.xyzw/PatchyV6C.RSD.%04i'%imock+'.xyzw.'+str(ibin)+'of3.numNB'+str(numNB)+'.wei'+number_to_str(weipow)+'.weipow1.%ds0to%d.%dmu.2pcf'%(sbin,smax,mubin)
    ### -------------------------------------
    ### file_system.diff_massbin
    ### -------------------------------------
    def diff_massbin(self, massbin=0.5, nbarstr='0.6e-4', numNB=6, margin=30, snapid=33,
                     weipow=None, shiftstr='shiftz',cat='2pcffile',xran=10,
                     sbin=150,smax=150,mubin=120):
        '''
            available massbin: 0, 0.5, 1.0, 1.5, .. 14.5, 15.0, 15.5 
            available nbar:  0.6e-4
            available numNB margin:  [2,30], [6,30],[20,100]   # default: 6, 30
            available weipow: [-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1]
            
            path: /home/xiaodongli/projects/cosmology_grid_related/snap_analysis/MCF_diff_massbin_bigmd/masscut_scheme2
        '''
        path = mainpath+'/MCF_diff_massbin_bigmd/masscut_scheme2/'
        mockfile = path + str(snapid)+'_nbar1e-3.'+shiftstr+'.allcolumns.nbar'+str(nbarstr)+\
                    '.massbin'+number_to_str(massbin)+\
                    '.txt.numNB'+str(numNB)+'.margin'+str(margin)+'.00.boxsize2500.00'
        ranfile = path + 'nbar'+nbarstr+'.box2500.x'+str(xran)+'ran.numNB'+str(numNB*xran)+'.margin'+str(margin)+'.00.boxsize2500.00'
        if cat == 'mockfile':
            return mockfile
        elif cat == 'ranfile':
            return ranfile
        elif cat =='inifile':
            inifile = mockfile + '.zplus1000000000.weipow'+str(weipow)+'.%d'%sbin+'s0to%d'%smax+'.%d'%mubin+'mu.ini'
            return inifile
        elif cat =='2pcffile':
            pcffile = mockfile + '.zplus1000000000.weipow'+str(weipow)+'.%d'%sbin+'s0to%d'%smax+'.%d'%mubin+'mu.2pcf'
            return pcffile
        elif cat == 'rrfile':
            rrfile = ranfile + '.weipow'+str(weipow)+'.%d'%sbin+'s0to%d'%smax+'.%d'%mubin+'mu.rr'
            return rrfile
        else:
            print('ERROR!')
    ### -------------------------------------
    ### file_system.diff_cosmo
    ### -------------------------------------
    def diff_cosmo(self, om=0.3071,w=-1.0000,redshift=0.4,nbarstr='3e-4',weipow=None,
                   numNB=30,margin=30,cat='2pcffile',
                   shiftstr='shiftz',xran=10,sbin=150,smax=150,mubin=120):
        '''
        self.diff_cosmo_reds = [0.4, 0.5, 0.6]
        self.diff_cosmo_numNB_margins = [[10,30], [30,30],[100,100] ]   # default: 30,30           
        self.diff_cosmo_nbars = ['3e-4', '3.5e-4', '4e-4']
        self.diff_cosmo_shiftstrs= ['shiftz', 'shiftr']
        available weipows: [-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1,  ]
        '''
        red_snap_dict = {0.4:'t', 0.5:'s', 0.6:'r' }
        snap_red_dict = {'t':0.4, 's':0.5, 'r':0.6 }

        path = mainpath
        mockfile = path + 'diff_cosmo_1/om%.4f'%om+'_omb0.0482_oml%.4f'%(1-om)+'_w%.4f'%w+'_Hubble0.68_sig80.8228_ns0.9600_simsnap01810'+\
                    red_snap_dict[redshift]+'_rockstar_halo_xyzvxvyvz_mvir_vmax.'+str(shiftstr)+'.nbar'+nbarstr+'.numNB'+str(numNB)+\
                    '.margin'+str(margin)+'.00.boxsize680.00'
        ranfile = path + 'random_file/nbar'+nbarstr+'.x'+str(xran)+'ran.numNB'+str(numNB*xran)+'.margin'+str(margin)+'.00.boxsize680.00'

        if cat == 'mockfile':
            return mockfile
        elif cat == 'ranfile':
            return ranfile
        elif cat == 'rrfile':
            rrfile = ranfile + '.weipow%.1f'%weipow+'.%d'%sbin+'s0to%d'%smax+'.%d'%mubin+'mu.rr'
            return rrfile 
        else:
            if shiftstr =='shiftz':
                if cat =='inifile':
                    inifile = mockfile + '.zplus1000000000.weipow%.1f'%weipow+'.%d'%sbin+'s0to%d'%smax+'.%d'%mubin+'mu.ini'
                    return inifile
                elif cat =='2pcffile':
                    pcffile = mockfile + '.zplus1000000000.weipow%.1f'%weipow+'.%d'%sbin+'s0to%d'%smax+'.%d'%mubin+'mu.2pcf'
                    return pcffile
                else:
                    print('ERROR!')
            else:
                if cat =='inifile':
                    inifile = mockfile + 'weipow%.1f'%weipow+'.%d'%sbin+'s0to%d'%smax+'.%d'%mubin+'mu.ini'
                    return inifile
                elif cat =='2pcffile':
                    pcffile = mockfile + 'weipow%.1f'%weipow+'.%d'%sbin+'s0to%d'%smax+'.%d'%mubin+'mu.2pcf'
                    return pcffile
                else:
                    print('ERROR!')
    ### -------------------------------------
    ### file_system.check_filess
    ### -------------------------------------                           
    def check_file(self, printinfo=False, return_if_found_missing=False):
        missing_files =[]
        print(' (file_system.check_file) check all available files...')
        
        ####--------------------------------------
        ## 1. check diff_mass files
        ####--------------------------------------
        #print('             (Step 1) check diff_massbin files')
        #print('                available massbins= ',massbin_edges.keys() )
        #print('                available nbars   = ',self.diff_massbin_nbars )
        #print('                available weipows = ',self.diff_massbin_weipows )
        #print('                available snapids = ',self.diff_massbin_snapids )
        #print('                available numNB/margins = ',self.diff_massbin_numNB_margins )
        
        flag = True; ifile = 0; imissing = 0
        for snapid in self.diff_massbin_snapids:
        #    print('                * snpid = ', snapid, ', checking different massbin...', end='\n\t\t')
            for massbin in massbin_edges.keys():
        #        print(massbin, end=' ==> ')
                for nbar in self.diff_massbin_nbars:
                    for weipow in self.diff_massbin_weipows:
                        for numNB, margin in self.diff_massbin_numNB_margins:
                            nowfile = self.diff_massbin(massbin=massbin, snapid=snapid, nbarstr=nbar, 
                                    numNB=numNB, margin=margin, weipow=weipow,); 
                            ifile += 1
                            if printinfo: print(os.popen('ls -alh '+nowfile).read())
                            if not os.path.isfile(nowfile):
                                print('WARNING! (file_system.checkfile): ', nowfile, 'not found!')
                                flag = False; imissing+=1; missing_files.append(nowfile)
                                if return_if_found_missing: return
        if flag: 
            print('                 Finish checking ', ifile, ' diff_massbin files (all pass!!!!)')
        else:
            print('                 Finish checking ', ifile, ' diff_massbin files (',imissing,' missing!!!!)')
        
        ####--------------------------------------
        ## check 2. diff_cosmo files
        ####--------------------------------------
        #print('             (Step 2). check diff_cosmo files')
        #print('                available omws          = ',self.diff_cosmo_omws )
        #print('                available reds          = ',self.diff_cosmo_reds )
        #print('                available numNB_margins = ',self.diff_cosmo_numNB_margins )
        flag = True; ifile = 0; imissing = 0
        for om, w in self.diff_cosmo_omws:
            for red in self.diff_cosmo_reds:
                for nbar in self.diff_cosmo_nbars:
                    for numNB, margin in self.diff_cosmo_numNB_margins:
                        for weipow in self.diff_cosmo_weipows:
                         for shiftstr in ['shiftz']:
                         #for shiftstr in self.diff_cosmo_shiftstrs:
                            nowfile = self.diff_cosmo(om=om, w=w, redshift=red,
                                            numNB=numNB, margin=margin,  weipow=weipow, shiftstr =shiftstr)
                            ifile += 1
                            if printinfo: print(os.popen('ls -alh '+nowfile).read())
                            if not os.path.isfile(nowfile):
                                print('WARNING! (file_system.checkfile): ', nowfile, 'not found!')
                                flag = False; imissing+=1; missing_files.append(nowfile)
                                if return_if_found_missing: return
        if flag: 
            print('                 Finish checking ', ifile, ' diff_cosmo files (all pass!!!!)')
        else:
            print('                 Finish checking ', ifile, ' diff_cosmo files (',imissing,' missing!!!!)')  
        
        ####--------------------------------------
        ### 3. check sdss files
        ####--------------------------------------
        flag = True; ifile = 0; imissing = 0
        for sample in self.sdss_sample_nsample.keys():
            nsample = self.sdss_sample_nsample[sample]
            print(' * sdss sample ', sample, ' (%d) samples'%nsample)
            for isample in range(nsample):
                for ibin in self.sdss_ibins:
                    print(' ** sdss ibin ', ibin)
                    for numNB, margin in self.diff_cosmo_numNB_margins:
                        print(' *** sdss numNB, margin ', numNB, margin, 'weipow = ', self.weipows)
                        for weipow in self.weipows:
                            nowfile = self.sdss(ibin=ibin, numNB=numNB, weipow=weipow, 
                                               sample=sample, imock=isample)
                            ifile += 1
                            if printinfo: print(os.popen('ls -alh '+nowfile).read())
                            if not os.path.isfile(nowfile):
                                print('WARNING! (file_system.checkfile): ', nowfile, 'not found!')
                                flag = False; imissing+=1; missing_files.append(nowfile)
                                if return_if_found_missing: return
        if flag: 
            print('                 Finish checking ', ifile, ' sdss data files (all pass!!!!)')
        else:
            print('                 Finish checking ', ifile, ' sdss data files (',imissing,' missing!!!!)') 
        print(' (file_system.check_file) Done.')
        return missing_files

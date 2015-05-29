#!/home/xiaodongli/software/anaconda/bin/python
# Filename: stdA.py 

import commands
import matplotlib.pyplot as plt
import numpy as np
import mpl_toolkits.mplot3d.axes3d as Axes3D
import commands, os, random, scipy.integrate, copy, sys, time, datetime
import scipy.ndimage as ndimage
import glob as glob
import textwrap

pythonlibPATH = commands.getoutput('echo $pythonlibPATH')
pyfile=pythonlibPATH+'/'+'stdA.py'

#################################
# Useful Definiations
#################################

FIG_X_SIZE = 18; FIG_Y_SIZE = 6; FIG_FONT_SIZE = 18

Figlabelsize = 16; plt.rc('xtick', labelsize=Figlabelsize); plt.rc('ytick', labelsize=Figlabelsize)

CONST_C = 299792.458 #unit: km/s


PLOT_COLOR_ARRAY = ['b', 'r', 'g', 'k', 'c', 'y', 'gray']
PLOT_STYLE_ARRAY = ['-', '--', '-.', '..']
NUM_PLOT_COLOR = len(PLOT_COLOR_ARRAY)
NUM_PLOT_STYLE = len(PLOT_STYLE_ARRAY)
NOW_PLOT_COLOR = 0; 
NOW_PLOT_STYLE = 0

import itertools


itcolors = itertools.cycle(['b', 'g', 'r', 'c', 'm', 'y', 'k', 'gray'])

def get_color_style(sepreturn=False,init=False):
    global NOW_PLOT_COLOR, NOW_PLOT_STYLE
    i = mod(NOW_PLOT_COLOR, NUM_PLOT_COLOR); j = mod(NOW_PLOT_COLOR, NUM_PLOT_STYLE)
    color_style = PLOT_COLOR_ARRAY[i]+PLOT_STYLE_ARRAY[j]
    NOW_PLOT_COLOR = NOW_PLOT_COLOR + 1; NOW_PLOT_STYLE = NOW_PLOT_STYLE + 1
    if init:
	NOW_PLOT_COLOR = 0; NOW_PLOT_STYLE = 0;
    if sepreturn:
	return PLOT_COLOR_ARRAY[i], PLOT_STYLE_ARRAY[j]
    else:
	return color_style


############################################################
### General useful functions
############################################################

def callsys(cmd): # This file only used in .py python code; not for ipython notebook (you can not get any screen print!)
	os.system(cmd)

def isfile(filename, exit_if_not_exit=False):
	if not exit_if_not_exit:
		return os.path.isfile(filename)
	else:
		if not os.path.isfile(filename):
			print '(isfile) (Exit now...) File not found:\n\t', filename
			sys.exit()
		else:
			return True

def loadtxt_rand(filename, rat=0.1, printinfo=False):
	nowf = open(filename, 'r') 
	rlt = []
	nlines = 0
	nlines_read = 0
	while True:
		nowstr = nowf.readline()
		if nowstr == '':
			break
		nlines += 1
		x = random.uniform(0,1)
		if x <= rat:
			wordlist = nowstr.split()
			if wordlist[0] == '#':
				continue
			else:
				rlt.append([float(x) for x in wordlist])
				nlines_read += 1
	nowf.close()
	if printinfo:
		print nlines_read, ' lines read from ', nlines, ' lines; file = ', filename
	return rlt
	
def get_funame(function):
	totstr = str(function)
	nowstr = ''
	for i in range(10, len(totstr)):
		if totstr[i] == ' ':
			break
		else:
			nowstr += totstr[i]
	return nowstr

def isfiles(files, printinfo=True, exit_if_not_exit=False, noprintatall=False):
    nmissfile = 0
    isfilelist = []
    missfilelist = []
    for nowfile in files:
	if isfile(nowfile, exit_if_not_exit=exit_if_not_exit):
		isfilelist.append(nowfile)
	else:
		missfilelist.append(nowfile)
		nmissfile += 1
    if printinfo:
	print ' (isfiles) Finishing checking ', len(files), 'files; #-missing = ', nmissfile
	print '\tList of existing files ('+str(len(isfilelist))+'):'
	for nowfile in isfilelist[0:min(10,len(isfilelist))]:
		print '\t\t', nowfile
	if len(isfilelist) > 10:
		print '\t\t......'
	if nmissfile != 0:
		print '\tList of missing files('+str(len(missfilelist))+'):'
		for nowfile in missfilelist[0:min(10,len(missfilelist))]:
			print '\t\t', nowfile
		if len(missfilelist) > 10:
			print '\t\t......'
    if nmissfile == 0:
	if printinfo: 
		print ' (isfiles) Pass!!!'
	return True
    else:
	if not noprintatall:
		print ' (isfiles) Not pass!!! nmissfile = ', nmissfile
	return False


def roundint(x):
	return int(x+0.5)

def FindRoot( f, f0, xl, xr, tol = 1e-6, maxstep=100 ): # Find Root: f(x) == f0
    fl = f(xl); fr = f(xr);
    if abs(fl-f0) < tol:
        return xl, fl-f0, 0, 'S'
    if abs(fr-f0) < tol:
        return xr, fr-f0, 0, 'S'
    if fl > f0 and fr > f0:
        print "Error! f(xl)/f(xr) = ", fl, fr, "both larger than ", f0, "!"
        return 0, 0, 0, 'E' # E means error
    if fl < f0 and fr < f0:
        print "Error! f(xl)/f(xr) = ", fl, fr, "both smaller than ", f0, "!"
        return 0, 0, 0, 'E'
    
    xm = xl+(xr-xl)*(f0-fl)/(fr-fl); fm = f(xm)
    step = 0
    while abs(fm-f0) > tol and step <= maxstep:
        if (fm > f0 and fl > f0) or (fm < f0 and fl < f0):
            xl = xm; fl = fm
        else:
            xr = xm; fr = fm
        xm = xl+(xr-xl)*(f0-fl)/(fr-fl); fm = f(xm)
        step +=1
    return xm, fm-f0, step, 'S' # S means success 

def shfile_baekdu(shfilename, jobname, ncpu, cmdstr, printinfo=True):
	nowf = open(shfilename, 'w')
	nowf.write('#!/bin/bash\n#$ -V\n#$ -cwd\n#$ -S /bin/bash\n#$ -N '+jobname+'\n\n### Use Parallel Environment mpi,1320 CPUs:\n#$ -pe openmpi '+str(ncpu)+'\n#$ -q normal\n#$ -R yes\n'+cmdstr)
	nowf.close()
	if printinfo:
		print ' (shfile_baekdu) create baekdu fmt shfile: \n\t', shfilename

def dimavg(x1, x2, avgdim=1):
    """ 
	Average of two numbers in some dimension
		(x1**avgdim + x2**avgdim) / 2.0)**(1.0/avgdim)	
    """
    return ((x1**avgdim + x2**avgdim) / 2.0)**(1.0/avgdim)	

############################################################
#### Geometry
############################################################

def Sphere_Vol(r, sky_coverage=1.0):
	return r**3.0*4.0*np.pi/3.0*sky_coverage

def getangle(x, y): # get theta, phi according to x,y. 2D
    theta = np.arcsin(y/np.sqrt(x**2+y**2))
    if x<0:
        theta = np.pi-theta
    if x>0 and y<0:
        theta = 2*np.pi+theta
    return theta

def getSC(x, y, z): # get r, theta, phi given x,y,z. 3D
    r=np.sqrt(x**2+y**2+z**2)
    theta=getangle(z,np.sqrt(x**2+y**2))
    phi=getangle(x,y)
    return r, theta, phi

def frommutotheta(oneminusmu, returnstr=False): # convert 1-mu to theta; mu \equiv cos(theta)
        mu = 1.0 - oneminusmu
        theta = np.arccos(mu)
        theta *= (180/np.pi)
        if returnstr:
            return '%.1f'%theta
        else:
            return theta

def XYofCircle(r=1.0, num=100):
	theta = np.linspace(0.0, 2.0*np.pi, num)
	X = r*np.cos(theta)
	Y = r*np.sin(theta)
	return X, Y

def radecr_to_xyz(ra, dec, r, degreefac = np.pi/180.0):
    ra = ra * degreefac
    dec = dec * degreefac
    x = r*np.cos(dec) * np.cos(ra)
    y = r*np.cos(dec) * np.sin(ra)
    z = r*np.sin(dec)
    return x,y,z
def xyz_to_radecr(x,y,z,degreefac = np.pi/180.0):
    r = np.sqrt(x*x+y*y+z*z)
    dec = (np.pi/2.0-getangle(z,np.sqrt(x*x+y*y))) / degreefac # dec = pi/2 - theta
    ra = getangle(x,y) / degreefac # phi
    return ra,dec, r
def list_xyz_to_radecr(X,Y,Z,degreefac = np.pi/180.0):
    ndat = len(X)
    if len(X) != len(Y) or len(X) != len(Z):
	print 'Error! length of X,Y,Z not match!: ', len(X), len(Y), len(Z)
	ndat = min(len(X), len(Y), len(Z))
    RA, DEC, R = range(ndat), range(ndat), range(ndat)
    for row in range(ndat):
	RA[row], DEC[row], R[row] = xyz_to_radecr(X[row], Y[row], Z[row],degreefac = degreefac)
    return RA, DEC, R
def list_xyz_to_r(X,Y,Z):
    ndat = len(X)
    if len(X) != len(Y) or len(X) != len(Z):
	print 'Error! length of X,Y,Z not match!: ', len(X), len(Y), len(Z)
	ndat = min(len(X), len(Y), len(Z))
    R = range(ndat)
    for row in range(ndat):
	R[row] = np.sqrt(X[row]**2.0 + Y[row]**2.0 + Z[row]**2.0)
    return R
def ReviseRA(RA):
	rlt = [ra for ra in RA]
	for row in range(len(rlt)):
            if rlt[row] > 180:
                rlt[row] = rlt[row] - 360
	return rlt

Area_AllSky = 4*np.pi*(180.0/np.pi)**2.0 # Total degree squares of the whole sky
def area_from_radec(ramin, ramax, decmin, decmax, degreefac = np.pi/180.0):
    return (np.sin(decmax*degreefac)-np.sin(decmin*degreefac)) * (ramax-ramin) / degreefac

############################################################
### Cosmology
############################################################

def Hz(omegam, w, h, z):
    return 100*h*np.sqrt(omegam*(1.0+z)**3.0 + (1.0-omegam)*(1.0+z)**(3.0*(1+w)))

def comov_r(omegam, w, h, z):
    global CONST_C
    x, y = scipy.integrate.quad(lambda x: 1.0/Hz(omegam, w, h, x), 0, z)
    return CONST_C * x * h 

def comov_r_dft(z, omegam=0.26, w=-1.0, h=0.7):
    return comov_r(omegam, w, h, z)
    
def get_z(omegam, w, h, given_comov_r, zl=0.0, zr=5.0):
    z, deltar, num_step, flag = FindRoot(lambda z: comov_r(omegam, w, h, z), given_comov_r, zl, zr)
    return z

def get_z_dft(given_comov_r, omegam=0.26, w=-1.0, h=0.7):
	return get_z(omegam, w, h, given_comov_r)

def CosmoConvert(x, y, z, om1, w1, om2, w2):
	r = np.sqrt(x*x+y*y+z*z)
	red = get_z(om1,w1,0.7,r)
	r2 = comov_r(om2, w2, 0.7, red)
	return x*r2/r, y*r2/r, z*r2/r

def CosmoConvert2D(x, y, om1, w1, om2, w2):
	r = np.sqrt(x*x+y*y)
	red = get_z(om1,w1,0.7,r)
	r2 = comov_r(om2, w2, 0.7, red)
	return x*r2/r, y*r2/r


def km_to_Mpctoh(v, h):
    return v / (3.08567758e19) * h

def omwstr(om, w, fmt='%.4f', sep1='', sep2='_'):
    return 'om'+sep1+str(fmt)%om+sep2+'w'+sep1+str(fmt)%w
def cosmoconvertedfilename(filename, om, w): ### filename
	return filename + '.cosmo-converted.' + omwstr(om,w)


### very specialized cosmological functions

execfile(pythonlibPATH+'/stdA_cosmo.py')

############################################################
### Arrays
############################################################

execfile(pythonlibPATH+'/stdA_array.py')

############################################################
### Files, strings
############################################################

def getfilelist(path):
	"""	List of files from ls"""
	return commands.getoutput('ls '+path).split()

def shfile__cmdlist(shfilename, cmdlist, binlocation='#!/bin/bash'):
	f0 = open(shfilename, 'w')
	f0.write(binlocation+'\n')
	for cmd in cmdlist:
		f0.write(cmd + '\n')
	f0.close()

def shfile__qsubjoblist(shfilename, jobshfilelist, binlocation='#!/bin/bash'):
	f0 = open(shfilename, 'w')
	f0.write(binlocation+'\n')
	for jobshfile in jobshfilelist:
		f0.write('qsub '+jobshfile + '\n')
	f0.close()

def loadfile(datafile, noreturn=True):
    rlt = []
    nowfile = open(datafile, 'r')
    while True:
        nowstr = nowfile.readline()
        if nowstr == '':
            break
	if not noreturn:
	        rlt.append(nowstr)
	else:
		rlt.append(nowstr[0:len(nowstr)-1])
    return rlt

def load_strarray(datafile, exitcode='#'):
    data = []
    nowfile = open(datafile, 'r')
    while True:
        nowstr = nowfile.readline()
        if nowstr == '':
            break
        data.append(str_to_numbers(nowstr, exitcode=exitcode, do_float_conver=False))
    return data

def cp_except(file1, file2, except_lines = [], except_strs = []):
	if len(except_lines) != len(except_strs):
		print 'ERROR (cp_except)!!! length mismatch: ', len(except_lines), len(except_strs)
		return 0
	for row in range(1,len(except_lines)):
		if except_lines[row] < except_lines[row-1]:
			print 'ERROR (cp_except)!!! except_lines must arranged from small to big but we find misordering: ', \
				row-1, except_lines[row-1], ' and ', row, except_lines[row]
			return 0
	if len(except_lines) == 0:
		commands.getoutput('cp '+file1+' '+file2)
		return 1
	nowiexcept = 0; nowlexcept = except_lines[nowiexcept]; nowstrexcept = except_strs[nowiexcept]
	f1 = open(file1, 'r')
	f2 = open(file2, 'w')
	nowl = 1
	while True:
		nowstr = f1.readline()
		if nowstr == '':
			break
		if nowl == nowlexcept:
			nowstr = nowstrexcept
			if nowiexcept != len(except_lines) -1:
				nowiexcept += 1
				nowlexcept = except_lines[nowiexcept]
				nowstrexcept = except_strs[nowiexcept]
		f2.write(nowstr)
		nowl += 1
	return 1

def indice_of_character(string,char):
	indice = []
	for i in range(len(string)):
		if string[i] == char:
			indice.append(i)
	return indice

def str_merge(strlist, div=' '):
	rlt = ''
	for i in range(len(strlist)-1):
		rlt += (strlist[i] + div)
	rlt += strlist[i+1]
	return rlt

def fmtstrlist(X, fmtstr='', div='', endstr=''):
	nowstr = ''
	for row in range(len(X)-1):
		x = X[row]
		nowstr += (fmtstr%x+div)
	x = X[len(X)-1]
	nowstr += (fmtstr%x+endstr)
	return nowstr

## separate the component of path and filename from a combination of path/file. 
## e.g., for /home/xiaodongli/haha.txt, will return ('/home/xiaodongli/', 'haha.txt')
def separate_path_file(path):
	""" Separate the component of path and filename from a combination of path/file. 
		e.g., for /home/xiaodongli/haha.txt, will return ('/home/xiaodongli/', 'haha.txt')
	"""
	if path == '':
		return '', ''
	else:
		indice = indice_of_character(path, '/')
		if indice == []:
			return '', path
		else:
			nowi = indice[len(indice)-1]
			return path[0:nowi+1], path[nowi+1:len(path)]

############################################################
### Statistics
############################################################

### CLs of 1-7 sigma
sig1=0.683;sig2=0.954;sig3=0.9973;sig4=0.999937;sig5=0.99999943;sig6=0.999999998;sig7=0.9999999999974

# Get the mean & variance, with errors
def get_stat_from_list(x,getvarer = False):
    n = len(x)
    x_av    = sum(x)/len(x)
    x_var   = (sum(x[row]**2.0 for row in range(n)) - n * x_av**2.0) /(n-1.0) #correct the bias using n-1
    x_av_er = np.sqrt(x_var) / np.sqrt(n-1.0) # correct the bias using n-1 
    if getvarer:
        x_var_er = np.sqrt(2.0/(n-1.0)) * x_var
        return x_av, x_var, x_av_er, x_var_er
    else:
        return x_av, x_var, x_av_er

### Advanced statistics

execfile(pythonlibPATH+'/stdA_stat.py')

############################################################
### Plottings
############################################################

def figax(figxsize=8, figysize=6, xlim=[], ylim=[], xlabel='', ylabel='', labelfontsize=22,
		tl=0.5,tr=0.5,tx='', tha='center', tva='center', tfs=16, title='', titlefontsize=18, noax=False,
		inputax=None):
	nowfig = plt.figure(figsize=(figxsize, figysize))
	if noax:
		nowax = None
	if not noax or inputax != None:
		if inputax != None:
			nowax = inputax
		else:
			nowax = nowfig.add_subplot(111)
		if xlim != []:
			nowax.set_xlim(xlim[0], xlim[1])
		if ylim != []:
			nowax.set_ylim(ylim[0], ylim[1])
		if xlabel != '':
			nowax.set_xlabel(xlabel, fontsize=labelfontsize)
		if ylabel != '':
			nowax.set_ylabel(ylabel, fontsize=labelfontsize)
		if tx != '':
			nowax.text(tl, tr, tx, horizontalalignment=tha, verticalalignment=tva,transform=ax.transAxes, fontsize=tfs)
		if title != '':
			nowax.set_title(title, fontsize=titlefontsize)
	return nowfig, nowax

def axtext(nowax, tl=0.5,tr=0.5,tx='', tha='center', tva='center', axtransform=True, 
		fd={'style':'normal', 'weight':'bold', 'size':13, 'color':'c'}):
	if tx != '':
		if axtransform:
			nowax.text(tl, tr, tx, horizontalalignment=tha, verticalalignment=tva,transform=ax.transAxes, fontdict=fd)
		else:
			nowax.text(tl, tr, tx, horizontalalignment=tha, verticalalignment=tva, fontdict=fd)

def ppt_firstpage(filename='sample', title='Please input title', datestr=None, name = 'Xiao-Dong Li', 
                  tiltfs=26, datafs=22, namefs=20, figxsize=8, figysize=6, fmt='pdf'):
    """ Create the first page of ppt"""
    fig, ax = figax(figxsize=figxsize, figysize=figysize)
    ax.text(0.5, 0.6, title, horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=tiltfs)
    if datestr==None:
        datestr = datetime.date.today()
    ax.text(0.5, 0.3, datestr, horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=datafs)
    
    ax.text(0.5, 0.2, name, horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=namefs)
    ax.set_xticks([])
    ax.set_yticks([])

    figname = filename+'.'+fmt
    print 'Saving to file: ', figname
    fig.savefig(figname, format=fmt)
    plt.show()

def ppt_textpage(filename='sample', inputtext='Please input text',  textfs=18, wrapsize=35,
                  figxsize=8, figysize=6, fmt='pdf'):
    """ Create a text page of ppt"""
    fig, ax = figax(figxsize=figxsize, figysize=figysize)
    if wrapsize != None:
        inputtextlist = textwrap.wrap(inputtext, width=wrapsize)
        inputtext=''
        for word in inputtextlist:
            inputtext += (word+'\n')
    ax.text(0.5, 0.5, inputtext, horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=textfs)
    ax.set_xticks([])
    ax.set_yticks([])
    
    figname = filename+'.'+fmt
    print 'Saving to file: ', figname
    fig.savefig(figname, format=fmt)
    plt.show()

def invert_ax(ax, axes='x'):
	if axes=='x':
		x1,x2=ax.get_xlim();
		ax.set_xlim(x2,x1)
	elif axes=='y':
		y1,y2=ax.get_ylim();
		ax.set_ylim(y2,y1)
	elif axes=='z':
		z1,z2=ax.get_zlim();
		ax.set_zlim(z2,z1)
	else:
		print 'ERROR (invert_ax): wrong axes (must be x,y,z)! ', axes


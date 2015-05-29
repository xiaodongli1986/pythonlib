
import commands
import numpy as np
import sys
import stdA


plotfunlist = ['plot', 'scatter', 'hist', 'errorbar', 'scatter3d','plot3d','scatterradec', 'histr', 'histrcube', 'scatterrw']
plotfun3dlist = ['scatter3d', 'plot3d']

printstr = 'Usage: EXE plotfun filename [-xcol xcol] [-ycol ycol] [-zcol zcol] [options...]'+\
	'\n\tInputs: '+\
	'\n\t\tplotfun:   '+stdA.str_merge(plotfunlist)+''+\
	'\n\t\tfilename:  name of data file '+\
	'\n\t\txcol, ycol, zcol:  columns as x, y, z '+\
	'\n\tOptions:'+\
	'\n\t\t[-savefig savefig]  T or F; save the plot as an eps file;  '+\
	'\n\t\t[-randrat randrat]  Randomly selecly a portion of file to plot; must be positive; can >1;'+\
	'\n\t\t[-figfmt figfmt]  format of figure; by default png'+\
	'\n\t\t[-showfig showfig]  whether show the figure; by default True; set as False if creating many files'+\
	'\n\t\t[-titlefs/-titlefontsize, -labelfs/-labelfontsize, -fs/-fontsize ]  size of font'+\
	'\n\t\t[-figtitle figtitle]  title of figure; by default name of the data file'+\
	'\n\t\t[-xlabel, -ylabel, -zlabel, -bins, -linecolor/-lc, -color/-c, -linewidth/-lw, -pointsize/-ps  -figxsize, -figysize ] detail plot settings '+\
	'\n\t\t[-singlecolfile/T,F]    set as true for single column file'+\
	'\n\t\t[-singleplot/T,F]       if T, all files plotted in one plot'+\
	'\n\t\t[-automaticcolor,autoc,autocolor,automaticc/T,F]   if T, automatically assign color'+\
	'\n\t\t[-showleg/T,F]          if T, show legend'+\
	'\n\t\t[-setxmin/-xmin, -setxmax/-xmax, -setymin/-ymin, -setymax/-ymax]   set the range of x,y in figure'+\
	'\n\t\t[-histrange]     set the range of histogram; in form of e.g. 0-100 '
	

cmdargs = sys.argv

if len(cmdargs) < 3:
	print printstr; sys.exit()

plotfun = cmdargs[1]


if not (plotfun in plotfunlist):
	print 'ERROR (PlotIJ)!: unknown plotfun: ', plotfun
	print printstr; sys.exit()

icmd = 2

filenames = commands.getoutput('ls '+cmdargs[2]).split()

print plotfun, len(filenames), 'files.'

xcol = 1
ycol = 2
zcol = 3
wcol = 4



colors = ['k', 'b', 'g', 'y', 'r', 'c', 'gray', 'm' ]

randrat = 1.1
figfmt = 'png'
titlefs = 18
labelfs = 18
legfs = 18
xlabel = 'x'
ylabel = 'y'
zlabel = 'z'
linecolor = 'k'
color = 'k'
bins = 50
origfigtitle = ''
singlecolfile = False
showfig = True
figxsize =8 
figysize =6
singleplot=False
automaticcolor=False
showleg=False
histrange=None

xmin=None
xmax=None
ymin=None
ymax=None

		
if plotfun in ['scatter', 'scatter3d', 'scatterradec', 'scatterrw']:
	linewidth=0
	pointsize = 2
else:
	linewidth=1
	pointsize = 0

if plotfun == 'scatterradec':
	xlabel = 'ra'; ylabel='dec'
if plotfun in ['hist', 'histr', 'histrcube']:
	ylabel = 'N'
if plotfun == 'histr':
	xlabel = 'r'	
if plotfun == 'histrcube':
	xlabel = '$r^3$'
if plotfun == 'scatterrw':
	xlabel = 'r'

if len(filenames) == 1:
	savefig = False
	showfig = True
else:
	savefig = True
	showfig = False


if len(cmdargs) >=4:

	optioncmds = cmdargs[3:len(cmdargs)]

#	optionstart = False

	for icmd in range(len(optioncmds)):

		if np.mod(icmd, 2) == 1:
			continue

		else:

			opt1 = optioncmds[icmd]
			if icmd+1 > len(optioncmds) -1:
				print 'ERROR (PlotIJ)!: failed to find value for arg ', opt1
				print printstr; sys.exit()

			opt2 = optioncmds[icmd+1]

			if opt1 == '-xcol':
				xcol = int(opt2)
			elif opt1 == '-ycol':
				ycol = int(opt2)
			elif opt1 == '-zcol':
				zcol = int(opt2)
			elif opt1 == '-savefig':
				if opt2[0] == 'T':
					savefig = True
				elif opt2[0] == 'F':
					savefig = False
				else:
					print 'ERROR (PlotIJ)!: wrong savefig; must start with T or F: ', opt2
					sys.exit()
			elif opt1 == '-randrat':
				randratstr = opt2
				randrat = float(opt2)
				if randrat < 0.0:
					print 'ERROR (PlotIJ)!: wrong option randrat; shall be positive!'
					sys.exit()
			elif opt1 == '-figfmt':
				figfmt=opt2
			elif opt1 == '-showfig':
				if opt2[0] == 'T':
					showfig = True
				elif opt2[0] == 'F':
					showfig = False
				else:
					print 'ERROR (PlotIJ)!: wrong showfig; must start with T or F: ', opt2
					sys.exit()
			elif opt1 in ['-titlefs', '-titlefontsize']:
				titlefs=int(opt2)
			elif opt1 in ['-labelfs', '-labelfontsize']:
				labelfs=int(opt2)
			elif opt1 in ['-legfs', '-legfontsize', '-legendfs', '-legendfontsize']:
				legfs=int(opt2)
			elif opt1 in ['-fs', '-fontsize']:
				titlefs=int(opt2)
				labelfs = titlefs
				legfs = titlefs
			elif opt1 == '-xlabel':
				xlabel=opt2
			elif opt1 == '-ylabel':
				ylabel=opt2
			elif opt1 == '-zlabel':
				zlabel=opt2
			elif opt1 in ['-color', '-c']:
				color=opt2
			elif opt1 in ['-linecolor', '-lc']:
				linecolor=opt2
			elif opt1 == '-bins':
				bins=int(opt2)
			elif opt1 == '-figtitle':
				origfigtitle=opt2
			elif opt1 == '-figxsize':
				figxsize=float(opt2)
			elif opt1 == '-figysize':
				figysize=float(opt2)
			elif opt1 in ['-linewidth', '-lw']:
				linewidth=float(opt2)
			elif opt1 in ['-pointsize', '-ps']:
				pointsize=float(opt2)	
			elif opt1 == '-singlecolfile':
				if opt2[0] == 'T':
					singlecolfile = True
				elif opt2[0] == 'F':
					singlecolfile = False
				else:
					print 'ERROR (PlotIJ)!: wrong singlecolfile; must start with T or F: ', singlecolfile
					sys.exit()
			elif opt1 == '-singleplot':
				if opt2[0] == 'T':
					singleplot = True
				elif opt2[0] == 'F':
					singleplot = False
				else:
					print 'ERROR (PlotIJ)!: wrong singleplot; must start with T or F: ', opt2
					sys.exit()
			elif opt1 in ['-automaticcolor', '-autoc', '-autocolor', '-automaticc']:
				if opt2[0] == 'T':
					automaticcolor = True
				elif opt2[0] == 'F':
					automaticcolor = False
				else:
					print 'ERROR (PlotIJ)!: wrong automaticcolor; must start with T or F: ', opt2
					sys.exit()
			elif opt1 in ['-setxmin', '-xmin']:
				xmin = float(opt2)
			elif opt1 in ['-setxmax', '-xmax']:
				xmax = float(opt2)
			elif opt1 in ['-setymin', '-ymin']:
				ymin = float(opt2)
			elif opt1 in ['-setymax', '-ymax']:
				ymax = float(opt2)
			elif opt1 == '-showleg':
				if opt2[0] == 'T':
					showleg = True
				elif opt2[0] == 'F':
					showleg = False
				else:
					print 'ERROR (PlotIJ)!: wrong showleg; must start with T or F: ', opt2
					sys.exit()
			elif opt1 == '-histrange':
				nowi = stdA.indice_of_character(opt2, '-')[0]
				histrange = (float(opt2[0:nowi]), float(opt2[nowi+1:len(opt2)]))
			else:
				print 'ERROR (PlotIJ)!: wrong option: ', opt1
				sys.exit()

			
### Make the plot

if singleplot:
	fig, ax = stdA.figax(figxsize=figxsize, figysize=figysize)

ifile = 0
for filename in filenames:
	ifile += 1
	if automaticcolor:
		linecolor = colors[np.mod(ifile,len(colors))]
		color = colors[np.mod(ifile,len(colors))]
	if not singleplot:
		fig, ax = stdA.figax(figxsize=figxsize, figysize=figysize)

	if randrat >= 1.0:
		data = np.loadtxt(filename)
	else:
		data = stdA.loadtxt_rand(filename, rat=randrat, printinfo=True)

	if singlecolfile:
		data = [[x] for x in data]

	if plotfun == 'plot':
		colstr = str(xcol)+'-'+str(ycol)
		X, Y = stdA.XYfromdata(data, xcol-1, ycol-1)
		ax.plot(X,Y,lw=linewidth,markersize=pointsize,c=linecolor,label=stdA.separate_path_file(filename)[1])
	elif  plotfun == 'scatter':
		colstr = str(xcol)+'-'+str(ycol)	
		X, Y = stdA.XYfromdata(data, xcol-1, ycol-1)
		ax.scatter(X,Y,lw=linewidth,s=pointsize,c=color,label=stdA.separate_path_file(filename)[1])
	elif plotfun == 'scatterradec':
		colstr = str(xcol)+'-'+str(ycol)+'-'+str(zcol)
		X, Y, Z = stdA.XYZfromdata(data, xcol-1, ycol-1, zcol-1)
		RA, DEC, R = stdA.list_xyz_to_radecr(X, Y, Z)	
		ax.scatter(RA,DEC,lw=linewidth,s=pointsize,c=color,label=stdA.separate_path_file(filename)[1])

	elif plotfun == 'scatterrw':
		colstr = str(xcol)+'-'+str(ycol)+'-'+str(zcol)
		X, Y, Z = stdA.XYZfromdata(data, xcol-1, ycol-1, zcol-1)
		if xcol == ycol == zcol:
			R = X
		else:
			R = [np.sqrt(X[row]**2.0 + Y[row]**2 + Z[row]**2) for row in range(len(X))]
		W = stdA.Xfromdata(data, wcol-1)
		ax.scatter(R,W,lw=linewidth,s=pointsize,c=color,label=stdA.separate_path_file(filename)[1])

	elif  plotfun == 'hist':
		colstr = str(xcol)
		X = stdA.Xfromdata(data, xcol-1)
		ax.hist(X,bins=bins,label=stdA.separate_path_file(filename)[1], range=histrange)
	elif plotfun == 'histr':
		colstr = str(xcol)+'-'+str(ycol)+'-'+str(zcol)
		X, Y, Z = stdA.XYZfromdata(data, xcol-1, ycol-1, zcol-1)
		if xcol == ycol == zcol:
			R = X
		else:
			R = [np.sqrt(X[row]**2.0 + Y[row]**2 + Z[row]**2) for row in range(len(X))]
		ax.hist(R,bins=bins,label=stdA.separate_path_file(filename)[1], range=histrange)
	elif plotfun == 'histrcube':
		colstr = str(xcol)+'-'+str(ycol)+'-'+str(zcol)
		X, Y, Z = stdA.XYZfromdata(data, xcol-1, ycol-1, zcol-1)
		if xcol == ycol == zcol:
			R = X
		else:
			R = [(X[row]**2.0 + Y[row]**2 + Z[row]**2)**1.5 for row in range(len(X))]
		ax.hist(R,bins=bins,label=stdA.separate_path_file(filename)[1], range=histrange)
	elif  plotfun == 'errorbar':
		colstr = str(xcol)+'-'+str(ycol)+'-'+str(zcol)
		X, Y, Z = stdA.XYZfromdata(data, xcol-1, ycol-1, zcol-1)
		ax.errorbar(X,Y,Z,label=stdA.separate_path_file(filename)[1], range=histrange)
	elif  plotfun == 'scatter3d':
		colstr = str(xcol)+'-'+str(ycol)+'-'+str(zcol)
		X, Y, Z = stdA.XYZfromdata(data, xcol-1, ycol-1, zcol-1)
		ax = fig.add_subplot(111, projection='3d')
		ax.scatter(X,Y,Z,lw=linewidth,s=pointsize,c=color,label=stdA.separate_path_file(filename)[1])

	elif  plotfun == 'plot3d':
		colstr = str(xcol)+'-'+str(ycol)+'-'+str(zcol)
		X, Y, Z = stdA.XYZfromdata(data, xcol-1, ycol-1, zcol-1)
		ax = fig.add_subplot(111, projection='3d')
		ax.plot(X,Y,Z,lw=linewidth,markersize=pointsize,c=linecolor,label=stdA.separate_path_file(filename)[1])

	ax.set_xlabel(xlabel, fontsize=labelfs)
	ax.set_ylabel(ylabel, fontsize=labelfs)		
	if plotfun in plotfun3dlist:
		ax.set_zlabel(zlabel, fontsize=labelfs)

	print '*** ',plotfun,' column ', colstr, ' of file ', filename
	ax.grid()
	if origfigtitle == '':
		figtitle = stdA.separate_path_file(filename)[1]
	else:
		figtitle = origfigtitle
	ax.set_title(figtitle, fontsize=titlefs)
	
	if xmin != None:
		ax.set_xlim(left=xmin)
	if xmax != None:
		ax.set_xlim(right=xmax)
	if ymin != None:
		ax.set_ylim(bottom=ymin)
	if ymax != None:
		ax.set_ylim(top=ymax)


	if showleg:
		ax.legend(loc='best', frameon=False,fontsize=legfs)

	if not singleplot:
		if savefig:
			figname = filename
			if randrat <1:
				figname += ('.randrat'+randratstr)
			figname += ('.col-'+colstr+'.'+plotfun+'.'+figfmt)
			print '\tSaving figure to: ', figname,  '...'
			fig.savefig(figname, format=figfmt)
		if showfig: stdA.plt.show()
if singleplot:
	if savefig:
		figtitle += ( ' et al ('+str(len(filenames))+' files)')
		ax.set_title(figtitle, fontsize=titlefs)
		figname = filename+'.'+str(len(filenames))+'files'
		if randrat <1:
			figname += ('.randrat'+randratstr)
		figname += ('.col-'+colstr+'.'+plotfun+'.'+figfmt)
		print '\tSaving figure to: ', figname,  '...'
		fig.savefig(figname, format=figfmt)
	if showfig: stdA.plt.show()


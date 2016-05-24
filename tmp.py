def smu_xis_loadinmock(imumin, 
		print_xis=False):
		''' Load in the xis computed from all mocks; returning the dictonary.
			format: catname, catname2, RSDstr, ibin, imock'''
	 	### 1. Load in xis from mock 
		xis_mock = {}
		xisfile = smu_xis_xisdir + 'MockResult.imumin%03i'%imumin+'.xis'
		for catname in catnamelist:
		  xis_mock[catname] = {}
		  for catname2 in catname2list:
			xis_mock[catname][catname2] = {}
			for RSDstr in ['noRSD', 'RSD']:
				xis_mock[catname][catname2][RSDstr] = {}
		   		for ibin in range(totbin):
				    xis_mock[catname][catname2][RSDstr][ibin]= []
		nowf = open(xisfile, 'r')
		print '  Open and read: ', xisfile
		while True:
			nowstr = nowf.readline()
			if nowstr == '': break
			nowstr = nowstr.split()
			if nowstr[0][0] == '#': continue
			catname, catname2, RSDstr,  imock, ibin = \
				nowstr[0] , nowstr[1], nowstr[2], int(nowstr[3]), int(nowstr[4])
			xis_mock[catname][catname2][RSDstr][ibin].append([float(nowstr[row]) for row in range(5,len(nowstr))])
			if False and random.random() > 0.99: ### Checking the correctness of xi loaded from file...
				X1 = [float(nowstr[row]) for row in range(5,len(nowstr))]
                                galfile = binsplittedfilename(mockfile(catname, catname2, imock, RSDstr, ), ibin+1,totbin)
                             	Tpcfrltfile = Tpcfrltfilename(galfile, smax_mock, smax_mock, nummubins)
				X2 = smu_xis(Tpcfrltfile, smax=smax_mock, sfact=0, make_plot=False, imumin=imumin);
				print len(X1), len(X2)
				sys.exit()
		nowf.close()
		print '  Load in Done. # of xis loaded: '
		if print_nummock:
		 for catname in catnamelist:
		  print '#################'
                  for catname2 in catname2list:
			print '######'
                        for RSDstr in ['noRSD', 'RSD']:
				print 
                                for ibin in range(totbin):
					print '    ', catname, catname2, RSDstr, ibin,':', \
						 len(xis_mock[catname][catname2][RSDstr][ibin])
		return xis_mock


execfile('/home/xiaodongli/software/pythonlib/Tpcftools.py')

filename = "data"

filename2 = "data.radecr"

data = np.loadtxt(filename)

X,Y,Z = XYZfromdata(data)
RA,DEC,R = list_xyz_to_radecr(X,Y,Z)

f2 = open(filename2,'w')

for i in range(len(R)):
	f2.write(str(RA[i])+' '+str(DEC[i])+' '+str(R[i])+'\n')
f2.close()



import numpy as np
#import matplotlib.pyplot as plt
import sys
from sklearn import neighbors

printstr = 'Usage: py_neighbors -inputfile inputfile -boxsize boxsize -num_NB num_NB -buffer buffer \n\tbuffer: add buffer region before neighbor search (for peoriodic boundary condition sample)'

cmdargs = sys.argv
if len(cmdargs) <7:
    print(printstr)
    sys.exit()
buffer = 0

for iarg in range(1, len(cmdargs), 2):
    str1, str2 = cmdargs[iarg], cmdargs[iarg+1]
    print(str1, str2)
    if str1 in ['-inputfile']:
        inputfile = str2
    elif str1 in ['-boxsize']:
        boxsize = float(str2)
    elif str1 in ['-num_NB']:
        num_NB = int(str2)
    elif str1 in ['-buffer']:
        buffer = float(str2)
    else:
        print('unknown option:', str1)
        print(printstr)
        sys.exit()

print(' (py_neighbors) options:')
print('\tinputfile = ', inputfile)
print('\tboxsize= ', boxsize)
print('\tnum_NB= ', num_NB)
print('\tbuffer= ', buffer)

if True:
    
    
    outputfile = inputfile+'.row_of_'+str(num_NB) + 'neighbors.peo_bound_condition_buffer'+str(buffer)
    outputfile2 = inputfile+'.row_of_'+str(num_NB) + 'neighbors.peo_bound_condition_buffer'+str(buffer)+'.distances'
    
    data = np.loadtxt(inputfile)
    cols = np.arange(0, len(data))
    data = np.column_stack([data, cols.reshape(-1, 1)])
    
    ndat = len(data)
    
    print(' (py_neighbors)  load in data from ', inputfile)
    print(' (py_neighbors)  len(data) = ', len(data))
    print(' (py_neighbors)  adding buffer:  boxsize, buffer= ', boxsize, buffer)
    
    def add_buffer(data, boxsize=boxsize, buffer=buffer):
        newdata = []
        for row in range(len(data)):
            x, y, z = data[row,:3]
            for xshift in [-1,0,1]:
                for yshift in [-1,0,1]:
                    for zshift in [-1,0,1]:
                        if xshift==0 and yshift==0 and zshift==0: continue
                        nowx = x+xshift*boxsize
                        nowy = y+yshift*boxsize
                        nowz = z+zshift*boxsize
                        
                        if -buffer < nowx < boxsize+buffer and \
                            -buffer < nowy < boxsize+buffer and \
                            -buffer < nowz < boxsize+buffer:
                                newdata.append(
                                    np.concatenate( [ np.array([nowx,nowy,nowz]), data[row,3:]  ])
                                )
        return np.row_stack( [data, np.array(newdata)] )
    
    data = add_buffer(data, boxsize, buffer)
    print(' (py_neighbors)  after adding buffer len(data) = ', len(data))
    
    print(' (py_neighbors)  start searching for neighbors...')
    
    model = neighbors.NearestNeighbors(n_neighbors = num_NB+1)
    model.fit(X=data[:,:3])
    
    distance_of_neighbors, py_neighbors = [], []
    for row in range(ndat):
        if row %10000 == 0:
            print('        searching for ', row, '-th object... (%.2f'%(row/ndat*100)+'%)')
        X = data[row, :3]
        distance, index = model.kneighbors(X.reshape(-1,3),   )
        distance = distance[0]; index = index[0]
        rows = data[:,-1][index].astype('int')
        distance_of_neighbors.append(distance[1:])
        py_neighbors.append(rows[1:])

    print(' (py_neighbors) saving results to: \n\t\t', outputfile, '\n\t\t', outputfile2 )
    
    np.savetxt(outputfile, np.array(py_neighbors), fmt='%d')
    np.savetxt(outputfile2, np.array(distance_of_neighbors), fmt='%.7f')
    

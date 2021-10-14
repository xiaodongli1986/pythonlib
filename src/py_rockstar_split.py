import os, sys
import numpy as np

#nsnap = 37
#ncore = 12   # ncore <= 24

printstr = '''# Example: py_rockstar_split  -nsnap 37  -ncore 12 
                  ### for 1 Node of TianHe, ncore<=24. ### '''

args = sys.argv

if len(args)<=1 or len(args)%2!=1:
    print(printstr)
    sys.exit()
print('len(args)= ', len(args))
for iarg in range(1, len(args), 2):
    arg1, arg2 = args[iarg], args[iarg+1]
    if arg1[0] != '-':
        print('Unkown arg1 :', arg1)
        print(printstr)
        sys.exit()
        
    arg1 = arg1[1:]
    
    if arg1 == 'nsnap':
        nsnap = int(arg2)
        print('\t nsnap= ', nsnap)
    elif arg1 == 'ncore':
        ncore = int(arg2)
        if ncore > 24:
            print('For 1 Node of TianHe, ncore<=24.')
            sys.exit()
        
        print('\t ncore= ', ncore)
        


remainder = nsnap%ncore
idx_list  = np.arange(0, nsnap-remainder, int(nsnap/ncore))

with open('1_py_rockstar_multisnap_all.sh', mode='w') as ff:
    ff.write('#!/bin/sh ' +'\n')
    
    for i in idx_list[:-(remainder)]:
        
        with open('1_py_rockstar_multisnap_'+str(i)+'_'+str(i+2)+'.sh', mode='w') as f:
            f.write('#!/bin/sh '     +'\n')
            f.write('ranseed=04000 ' +'\n')
            f.write('basename=sim '  +'\n')
            for j in range(int(nsnap/ncore)):
                f.write('py_rockstar ${basename}\*snap${ranseed}.'+str(i+j).zfill(3)+'.\*'+'\n')
            f.close()
            
            ff.write('yhrun -N 1 -n 1 '+'1_py_rockstar_multisnap_'+str(i)+'_'+str(i+2)+'.sh'+' & \n')

    for i in idx_list[-(remainder):]:
        
        with open('1_py_rockstar_multisnap_'+str(i)+'_'+str(i+3)+'.sh', mode='w') as f:
            f.write('#!/bin/sh '     +'\n')
            f.write('ranseed=04000 ' +'\n')
            f.write('basename=sim '  +'\n')
            for j in range(int(nsnap/ncore+1)):
                f.write('py_rockstar ${basename}\*snap${ranseed}.'+str(i+j).zfill(3)+'.\*'+'\n')
            f.close()

            ff.write('yhrun -N 1 -n 1 '+'1_py_rockstar_multisnap_'+str(i)+'_'+str(i+3)+'.sh'+' & \n')
            
    ff.write('wait')
ff.close()







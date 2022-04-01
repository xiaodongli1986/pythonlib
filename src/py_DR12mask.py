# Developed by jiacheng

import os,sys

printstr = 'Example: \n'+"py_DR12mask -radecfile absolute_path/radecfile.txt -skyregion NGC"+"\n skyregion: NGC, SGC, or TOT "

args = sys.argv
    
if len(args)<=1 or len(args)%2!=1:
    print(printstr)
    sys.exit()

for iarg in range(1, len(args), 2):

    arg1, arg2 = args[iarg], args[iarg+1]

    if arg1[0] != '-':
        print('Unkown arg1 :', arg1)
        print(printstr)
        sys.exit()

    arg1=arg1[1:]

    if arg1 == 'radecfile':
        radecfile = arg2
        print('\t radecfile= ', radecfile)

    elif arg1 == 'skyregion':
        skyregion = arg2

        if skyregion in ['NGC', 'SGC', 'TOT']:
            print('\t skyregion= ', skyregion)
        else:
            print('skyregion should be NGC, SGC, or TOT !')
            sys.exit()                    
    else:
        print('printstr')
        sys.exit()

print('jsub job ...')

cmd = 'jsub -o AddMask.out -n 9 mpirun -n 9 python /home/xiaodongli/software/pythonlib/DR12_mask/MPI_create_DR12_CMASSmask.py -radecfile '+radecfile+' -skyregion '+skyregion

print('cmd=\n'+cmd)

os.popen(cmd)











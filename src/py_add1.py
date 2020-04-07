
import sys
print('Usage: py_add1 inputfile')

inputfile = sys.argv[1]
outputfile = inputfile+'.add1'

print("open ", inputfile, 'for read...')
print("open ", outputfile, 'for write ...')

nowf0 = open(inputfile,'r')
nowf = open(outputfile, 'w')
iline = 0
while True:
    nowstr = nowf0.readline()
    if nowstr == '': break
    nowf.write(nowstr[:len(nowstr)-1]+' 1\n')
    iline += 1
print (iline, 'lines written to ', outputfile)



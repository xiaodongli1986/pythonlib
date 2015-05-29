#!/usr/bin/python
import sys
import commands
import os

def str_to_numbers(str1,exitcode='#', do_float_conver=True):
    floatlist = []
    str2 = ''
    numberunderconstruct = False
    for i in range(len(str1)):
        duru = str1[i]
        if duru == exitcode:
            break
        elif duru == '\n' or duru == ' ' or duru == '\t' or duru == ',':
            if numberunderconstruct == True:
                floatlist.append(str2); str2 = ''
                numberunderconstruct = False
        else:
            if numberunderconstruct == False:
                numberunderconstruct = True
                str2 = duru
            else:
                str2 += duru
    if str2 != '':
        floatlist.append(str2);
    if do_float_conver:
        floatlist = [float(floatlist[row]) for row in range(len(floatlist))]
    return floatlist


### Declare the function of this code
print 
print 'Add a new main programe to the py codes...'
print 'Usage: EXE Name-of-the-main-program '
print

### Get the command: which programme you want to add?
cmdargs = sys.argv
if len(cmdargs) != 2:
	print 'Error: len(cmdargs) = ', len(cmdargs)
	sys.exit()
nowEXEname = cmdargs[1]


print
print 'We will add a new EXE named as  ', nowEXEname

### Make a save of the Makefile... Very headache code
output = commands.getoutput('ls')
output = str_to_numbers(output, do_float_conver = False)
MFs = []
MFIs = [-1]
for MF in output:
	if MF[0:13] =='Makefile.SAVE':
		MFs.append(MF)
		MFIs.append(int(MF[13:len(MF)]))
saveindex = max(MFIs) + 1
tmpstr = 'Makefile.SAVE'+str(saveindex)
print 'Copy of old Makefile saved to :', tmpstr
os.system('cp Makefile '+tmpstr)
file1 = tmpstr
os.system('rm Makefile')
file2 = 'Makefile'

### Generate new Makefile...
nowf1 = open(file1, 'r')
nowf2 = open(file2, 'w')

while True:
	nowstr = nowf1.readline()
	if nowstr == '':
		break
	nowstr0 = nowstr[0:len(nowstr)-1]
	now_str_array = str_to_numbers(nowstr, do_float_conver = False)
	if len(now_str_array) == 0:
		nowf2.write(nowstr)
		continue
	elif now_str_array[0] == 'EXEs':
		nowf2.write(nowstr0 + ' ../bin/py_'+nowEXEname+'\n')
		continue
	else:
		nowf2.write(nowstr)
		continue
nowf1.close()
nowf2.close()

pyname = 'py_'+nowEXEname+'.py'
nowf3 = open(pyname, 'w')
nowf3.write('import commands\n')
nowf3.write('pythonlibPATH=commands.getoutput(\'echo $pythonlibPATH\')\n')
nowf3.write('execfile(pythonlibPATH+\'/stdA.py\')\n\n')
nowf3.write('print \'This is an empty python programme: add something you want!\'\n')
nowf3.close()

print '\nNew file created: ', pyname, '\n\n'


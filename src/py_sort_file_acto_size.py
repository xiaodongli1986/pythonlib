
import sys, random, numpy, os

print('Usage: \n\tEXE inputfilelist (this file containing the names of files you want to sort)')
print('Example:\n\tpy_ran file_names.txt')
args = sys.argv

if len(args) <2:
	print('ERROR (ran.py)! Must have at least 2 arguments.')
	sys.exit()

inputfilelist = args[1]
files = open(inputfilelist).readlines()

print('Checking ', len(files), 'files...')

print('ls -alh for all files...')
file_infos = []

for nowfile in files:
    file_infos.append(os.popen('ls -l '+nowfile).read().split())

for file_info in file_infos:
    file_info[4] = int(file_info[4])

print('Sorting according to sizes of files...')
def get_val(file_info):
    return file_info[4]
file_infos.sort(key = get_val)

outputfilelist = inputfilelist+'.sorted'
print('Write sorted result to ', outputfilelist)
nowf = open(outputfilelist, 'w')
for file_info in file_infos:
    nowf.write(file_info[len(file_info)-1]+'\n')
    print(os.popen('ls -l '+file_info[len(file_info)-1]).read())
nowf.close()

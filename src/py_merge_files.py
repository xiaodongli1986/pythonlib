
import sys, random, numpy, os

print 'Usage: \n\tEXE filename delimiter'
print 'Example:\n\tpy_merge_files tmpfile.ibox . (This will merge tmpfile.ibox0.XX into tmpfile.ibox0, merge tmpfile.ibox1.XX into tmpfile.ibox1, ...)'
args = sys.argv

if len(args) <2:
	print 'ERROR (py_merge_files.py)! Must have at least 2 arguments.'
	sys.exit()

filename = args[1]
sep = args[2]
files = os.popen("ls "+filename+'*').read().split()
print 'filename, sep = ', filename, sep

suffixes = []
for nowfile in files:
    nowstr = nowfile[len(filename):]
    #print 'nowfile, nowstr = ', nowfile, nowstr
    for ichar, nowchar in enumerate(nowstr):
        if nowchar == sep:
            break
    suffixes.append(nowstr[:ichar])

suffixes_set = list(set(suffixes))

inputfile_dict = {}
outputfile_dict = {}
print 'In total ', len(suffixes_set), 'files will be generated:'
for suffix in suffixes_set:
    inputfile_dict [suffix] = []
    #print '\t', filename+suffix
    outputfile_dict [suffix] = filename+suffix

ifile = 0
#print files
for nowfile in files:
    #print ifile, nowfile, suffixes[ifile]
    inputfile_dict[suffixes[ifile]].append(nowfile)
    ifile += 1

for key in inputfile_dict.keys():
#    continue
    print '\t', len(inputfile_dict[key]), ' inputfiles will be merged...'
    #print '\t\t', inputfile_dict[key]
    for nowinputfile in inputfile_dict[key]:
        print '\t\t', os.popen("ls -alh "+nowinputfile).read()
    print '\tNow merge them into', outputfile_dict[key], '...'
    for ifile, nowfile in enumerate(inputfile_dict[key]):
        print '\t\tmerging', ifile,'- th file', nowfile, '...' 
        if ifile == 0: 
            os.popen('cp '+nowfile+' '+outputfile_dict[key])
        else:
            os.popen('cat '+nowfile+' >> '+outputfile_dict[key])
    print '\tMerged file generated:'
    print '\t\t', os.popen("ls -alh "+nowinputfile).read(), '\n'


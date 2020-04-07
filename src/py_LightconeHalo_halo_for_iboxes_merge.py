
import os, numpy, sys

print 'Usage: EXE num-of-iboxes overlap xyzmin xyzmax dir-name file-name output_1-eighth-shell(T or F) in-one-dir (T or F)'
print '## The files with names "dir-name$i/file-name", where i goes from 1 to num-of-iboxes, will be merged, with their overlapping region properly treated (so that no redundant particlse are included) '
print '## The last argument tells whether just outputing 1/8 shell (xyz>0). By default False (output all sphere)'

#filename = 'test_lightcone.nbox4_overlap10.0_xyz-1800.0to1800.0.ibox'


args = sys.argv

if len(args) < 7:
    sys.exit()

nbox, overlap, xyzmin, xyzmax, dirname, filename  = args[1:7]

shell_output = False
in_one_dir = False

if len(args) >= 8:
    shell_output = args[7] 
    if shell_output[0] in ['T', 't']:
        shell_output = True
    elif shell_output[0] in ['F', 'f']:
        shell_output = False
    print 'set according to input:  args[7], shell_output = ', args[7], shell_output
if len(args) >= 9:
    nowstr = args[8] 
    if nowstr[0] in ['T', 't']:
        in_one_dir = True
    elif nowstr[0] in ['F', 'f']:
        in_one_dir = False
    print 'set according to input:  args[8], in_one_dir   = ', args[8], in_one_dir  


nbox, overlap, xyzmin, xyzmax = int(nbox), float(overlap), float(xyzmin), float(xyzmax)
#outputfile = filename+'_merged'
outputfile =  filename+'.nbox'+str(nbox)+'_overlap%.1f'%overlap+'_xyz%.1f'%xyzmin+'to%.1f'%xyzmax+'.iboxALL_rockstar_halo.ascii.xyzvxvyvz_mvir_vmax' 

if not in_one_dir:
    files = [dirname+str(ibox)+'/'+filename for ibox in range(1,nbox**3+1)]
else:
    files = [dirname+'/'+filename+'.nbox'+str(nbox)+'_overlap%.1f'%overlap+'_xyz%.1f'%xyzmin+'to%.1f'%xyzmax+'.ibox'+str(ibox)+'_rockstar_halo.ascii.xyzvxvyvz_mvir_vmax' for ibox in range(1,nbox**3+1)]
    outputfile = dirname + '/' + outputfile

print 'Write results to ', outputfile
dxyz = (xyzmax - xyzmin) / float(nbox)


zranges = [[xyzmin + dxyz*(i%nbox), xyzmin+dxyz*((i)%nbox+1)] for i in range(0,nbox**3)]
yranges = [[xyzmin + dxyz*((i//nbox)%nbox), xyzmin+dxyz*((i//nbox)%nbox+1)] for i in range(0,nbox**3)]
xranges = [[xyzmin + dxyz*((i//(nbox**2))%nbox), xyzmin+dxyz*((i//(nbox**2))%nbox+1)] for i in range(0,nbox**3)]

print 'Processing files with enforced ranges:'
for ibox in range(nbox*nbox*nbox):
    print '\tbox', ibox, ':', xranges[ibox], yranges[ibox], zranges[ibox]

print '###########################################'
fout = open(outputfile, 'w')
for ibox in range(nbox**3):
    if shell_output: 
        if xranges[ibox][1] < 0 or yranges[ibox][1] < 0 or zranges[ibox][1] < 0:
           continue
    print '\tbox', ibox, ':', xranges[ibox], yranges[ibox], zranges[ibox]
    print '\t\tRead in ', files[ibox], '...'
    nowf = open(files[ibox], 'r')
    i1, i2 = 0, 0
    while True:
        nowstr = nowf.readline()
        if nowstr == '': break
        nowstrs = nowstr.split()
        x, y, z = float(nowstrs[0]), float(nowstrs[1]), float(nowstrs[2])
        if shell_output:
            if x<0 or y<0 or z<0: continue
        if xranges[ibox][0] < x and x < xranges[ibox][1] and yranges[ibox][0] < y and y < yranges[ibox][1] and zranges[ibox][0] < z and z < zranges[ibox][1]:
            fout.write(nowstr)
            i2 += 1
        i1 += 1
    nowf.close()
    print '\t\t\t\t * Lines   read-in / write    =    %10i'%i1, '%10i'%i2
fout.close()
